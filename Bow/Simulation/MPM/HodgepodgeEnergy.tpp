#include "HodgepodgeEnergy.h"

#include <Bow/Physics/FixedCorotated.h>

namespace Bow::MPM {

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::precompute(const Field<Vector<T, dim>>& x)
{
    Field<TM> m_gradXp(m_X.size(), TM::Zero());
    tbb::parallel_for(size_t(0), m_X.size(), [&](size_t i) {
        const Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        Matrix<T, dim, dim> gradXp = Matrix<T, dim, dim>::Zero();
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            gradXp.noalias() += x[g.idx] * dw.transpose();
        });
        m_gradXp[i] = gradXp;
    });
    for (auto& model : elasticity_models)
        model->trial_strain(m_gradXp);
}

template <class T, int dim, class StorageIndex>
T HodgepodgeEnergy<T, dim, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
{
    T total_energy = 0;
    grid.iterateGridSerial([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        T m = g.v_and_m(dim);
        TV x_n = node.template cast<T>() * dx;
        TV v_n = g.v_and_m.template topLeftCorner<dim, 1>();
        TV a_n = g.a;
        TV x_tilde = x_n + v_n * dt + tsParam[tsMethod][0] * (1 - 2 * tsParam[tsMethod][1]) * a_n * dt * dt;
        total_energy += (T)0.5 * m * (x[g.idx] - x_tilde).squaredNorm();
        total_energy -= 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * m * gravity.dot(x[g.idx] - x_n);
    });
    Field<T> t_energy(m_X.size(), (T)0);
    for (auto& model : elasticity_models) {
        model->trial_energy(t_energy);
    }
    for (size_t i = 0; i < m_X.size(); ++i) {
        total_energy += 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * t_energy[i];
    }
    return this->energy_scale * total_energy;
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
{
    grad.assign(grid.num_nodes, Vector<T, dim>::Zero());
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        if (g.idx < 0) return;
        auto idx = g.idx;
        T m = g.v_and_m(dim);
        TV x_n = node.template cast<T>() * dx;
        TV v_n = g.v_and_m.template topLeftCorner<dim, 1>();
        TV a_n = g.a;
        TV x_tilde = x_n + v_n * dt + tsParam[tsMethod][0] * (1 - 2 * tsParam[tsMethod][1]) * a_n * dt * dt;
        grad[idx] += this->energy_scale * m * (x[idx] - x_tilde);
        grad[idx] -= this->energy_scale * 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * m * gravity;
    });
    Field<TM> t_gradient(m_X.size(), TM::Zero());
    for (auto& model : elasticity_models) {
        model->trial_gradient(t_gradient);
    }
    grid.colored_for([&](int i) {
        Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        Matrix<T, dim, dim> stress = t_gradient[i];
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            auto idx = g.idx;
            grad[idx] += this->energy_scale * 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * stress * dw;
        });
    });
    for (int i = 0; i < grid.num_nodes; ++i) {
        grad[i] = BC_basis[i].transpose() * grad[i];
        // overwrite rhs 0 for empty row and col in I - S^TS
        for (int d = 0; d < BC_order[i]; ++d)
            grad[i](d) = 0;
    }
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::multiply(const Field<Vector<T, dim>>& _x, Field<Vector<T, dim>>& Ax, bool project_pd)
{
    Field<Vector<T, dim>> x = _x;
    // right multiply
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        int idx = g.idx;
        x[idx] = BC_basis[idx] * x[idx];
    });
    Ax.assign(grid.num_nodes, Vector<T, dim>::Zero());
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        auto idx = g.idx;
        T m = g.v_and_m(dim);
        Ax[idx] += m * x[idx];
    });
    Field<TM> d_F(m_X.size(), TM::Zero());
    Field<TM> t_differential(m_X.size(), TM::Zero());
    grid.parallel_for([&](int i) {
        const Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        TM gradVp = TM::Zero();
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            gradVp.noalias() += x[g.idx] * dw.transpose();
        });
        d_F[i] = gradVp;
    });
    for (auto& model : elasticity_models) {
        model->trial_differential(d_F, t_differential, project_pd);
    }
    grid.colored_for([&](int i) {
        const Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            auto idx = g.idx;
            Ax[idx] += this->energy_scale * 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * t_differential[i] * dw;
        });
    });
    // left multiply
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        int idx = g.idx;
        Ax[idx] = BC_basis[idx].transpose() * Ax[idx];
    });
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::precondition(Field<Vector<T, dim>>& diagonal) {
    diagonal.assign(grid.num_nodes, Vector<T, dim>::Zero());
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        int idx = g.idx;
        diagonal[idx] += g.v_and_m(dim) * TV::Ones();
    });
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::project(Field<Vector<T, dim>>& b) {
    for (int i = 0; i < grid.num_nodes; ++i)
        for (int d = 0; d < BC_order[i]; ++d)
            b[i](d) = 0;
}

// TODO: optimize neasted loop
// https://github.com/penn-graphics-research/ziran/blob/LBFGSAMG/Projects/multigrid/ImplicitSolver.h
// https://github.com/penn-graphics-research/ziran/blob/LBFGSAMG/Projects/multigrid/ImplicitSolver_prev.h
template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd)
{
    std::vector<int> entryRow(grid.num_nodes * kernelSize(), 0);
    std::vector<int> entryCol(grid.num_nodes * kernelSize(), 0);
    std::vector<TM> entryVal(grid.num_nodes * kernelSize(), TM::Zero());
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        T m = g.v_and_m(dim);
        entryRow[g.idx * kernelSize() + kernelOffset(Vector<int, dim>::Zero())] = g.idx;
        entryCol[g.idx * kernelSize() + kernelOffset(Vector<int, dim>::Zero())] = g.idx;
        entryVal[g.idx * kernelSize() + kernelOffset(Vector<int, dim>::Zero())] += m * TM::Identity();
    });
    Field<Matrix<T, dim * dim, dim * dim>> t_hessian(m_X.size(), Matrix<T, dim * dim, dim * dim>::Zero());
    for (auto& model : elasticity_models) {
        model->trial_hessian(t_hessian, project_pd);
    }
    grid.colored_for([&](int i) {
        Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        Matrix<T, dim * dim, dim * dim> deformed_dPdF = t_hessian[i];
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            grid.iterateKernel(spline, [&](const Vector<int, dim>& _node, T _w, Vector<T, dim> _dw, GridState<T, dim>& _g) {
                if (_g.idx < 0) return;
                TM dFdX = TM::Zero();
                for (int u = 0; u < dim; ++u)
                    for (int p = 0; p < dim; ++p)
                        for (int x = 0; x < dim; ++x)
                            for (int y = 0; y < dim; ++y)
                                dFdX(u, p) += deformed_dPdF(u + x * dim, p + y * dim) * dw(x) * _dw(y);
                entryRow[g.idx * kernelSize() + kernelOffset(node - _node)] = g.idx;
                entryCol[g.idx * kernelSize() + kernelOffset(node - _node)] = _g.idx;
                entryVal[g.idx * kernelSize() + kernelOffset(node - _node)] += 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt * dFdX;
            });
        });
    });
    using IJK = Eigen::Triplet<T>;
    std::vector<IJK> coeffs;
    for (int i = 0; i < grid.num_nodes * kernelSize(); ++i)
        entryVal[i] = BC_basis[entryRow[i]].transpose() * entryVal[i] * BC_basis[entryCol[i]];
    for (int i = 0; i < grid.num_nodes * kernelSize(); ++i)
        for (int u = 0; u < dim; ++u)
            for (int v = 0; v < dim; ++v) {
                // Eliminate I - S^TS
                if (u < BC_order[entryRow[i]] || v < BC_order[entryCol[i]])
                    continue;
                coeffs.push_back(IJK(entryRow[i] * dim + u, entryCol[i] * dim + v, this->energy_scale * entryVal[i](u, v)));
            }
    // add diagonal 1 for empty row and col in I - S^TS
    for (int i = 0; i < grid.num_nodes; ++i)
        for (int d = 0; d < BC_order[i]; ++d)
            coeffs.push_back(IJK(i * dim + d, i * dim + d, 1));
    hess.derived().resize(grid.num_nodes * dim, grid.num_nodes * dim);
    hess.derived().setZero();
    hess.derived().setFromTriplets(coeffs.begin(), coeffs.end());
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::postprocess(Field<Vector<T, dim>>& direction)
{
    for (int i = 0; i < grid.num_nodes; ++i)
        direction[i] = BC_basis[i] * direction[i];
}

template <class T, int dim, class StorageIndex>
T HodgepodgeEnergy<T, dim, StorageIndex>::stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& direction)
{
    T alpha = 1.0;
    Field<TM> m_gradDVp(m_X.size(), TM::Zero());
    tbb::parallel_for(size_t(0), m_X.size(), [&](size_t i) {
        const Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        Matrix<T, dim, dim> gradDVp = Matrix<T, dim, dim>::Zero();
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            TV new_v = direction[g.idx] / dt;
            gradDVp.noalias() += new_v * dw.transpose();
        });
        m_gradDVp[i] = gradDVp;
    });
    for (auto& model : elasticity_models) {
        alpha = std::min(alpha, model->stepsize_upperbound(m_gradDVp, dt));
    }
    return alpha;
}

template <class T, int dim, class StorageIndex>
void HodgepodgeEnergy<T, dim, StorageIndex>::internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force)
{
    force.assign(grid.num_nodes, Vector<T, dim>::Zero());
    grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
        if (g.idx < 0) return;
        auto idx = g.idx;
        T m = g.v_and_m(dim);
        force[idx] -= m * gravity;
    });
    precompute(xn);
    Field<TM> t_gradient(m_X.size(), TM::Zero());
    for (auto& model : elasticity_models) {
        model->trial_gradient(t_gradient);
    }
    grid.colored_for([&](int i) {
        Vector<T, dim>& Xp = m_X[i];
        BSplineWeights<T, dim> spline(Xp, dx);
        Matrix<T, dim, dim> stress = t_gradient[i];
        grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            auto idx = g.idx;
            force[idx] += stress * dw;
        });
    });
}

}