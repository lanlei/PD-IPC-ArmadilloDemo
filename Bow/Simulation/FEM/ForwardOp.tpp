#include "ForwardOp.h"
#include "Utils.h"
#include <Bow/Optimization/Newton.h>
#include <iostream>
#include <Bow/Utils/Timer.h>
#include <oneapi/tbb.h>
#include <fstream>
#include <Bow/Math/Utils.h>

namespace Bow {
namespace FEM {

/* InitializeOp */
template <class T, int dim>
inline void InitializeOp<T, dim>::operator()()
{
    m_vol.resize(m_elem.size());
    m_IB.resize(m_elem.size());
    tbb::parallel_for(size_t(0), m_elem.size(), [&](size_t e) {
        const auto& indices = m_elem[e];
        const auto& X0 = m_X[indices[0]];
        const auto& X1 = m_X[indices[1]];
        const auto& X2 = m_X[indices[2]];
        Matrix<T, dim, dim> B;
        B.col(0) = X1 - X0;
        B.col(1) = X2 - X0;
        if constexpr (dim == 3) {
            const auto& X3 = m_X[indices[3]];
            B.col(2) = X3 - X0;
            m_vol[e] = B.determinant() / T(6);
        }
        else {
            m_vol[e] = B.determinant() / T(2);
        }
        m_IB[e] = B.inverse();
    });

    // compute node mass
    m_mass.resize(m_X.size());
    tbb::parallel_for(size_t(0), m_X.size(), [&](size_t i) {
        m_mass[i] = 0;
    });
    for (size_t e = 0; e < m_elem.size(); ++e) {
        const auto& indices = m_elem[e];
        T mass_contrib = m_density[e] * m_vol[e] / T(dim + 1);
        for (int i = 0; i < dim + 1; ++i) {
            m_mass[indices[i]] += mass_contrib;
        }
    }
}

/* BackwardEulerUpdateOp */

template <class T, int dim, class _StorageIndex>
BackwardEulerUpdateOp<T, dim, _StorageIndex>::BackwardEulerUpdateOp(const Field<Matrix<T, dim, dim>>& BC_basis, const std::vector<int>& BC_order, const Field<Vector<T, dim>>& BC_target, Field<T>& mass, Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& v, Field<Vector<T, dim>>& a, Field<Vector<T, dim>>& x_tilde)
    : BC_basis(BC_basis)
    , BC_order(BC_order)
    , BC_target(BC_target)
    , m_x(x)
    , m_v(v)
    , m_a(a)
    , m_x_tilde(x_tilde)
    , m_mass(mass)
{
    update_transformation_matrix();
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::update_transformation_matrix()
{
    // construct transform matrix
    std::vector<StorageIndex> ptr(BC_basis.size() * dim + 1);
    ptr[0] = 0;
    std::vector<StorageIndex> row(BC_basis.size() * dim * dim);
    for (size_t i = 0; i < BC_basis.size(); ++i) {
        for (int d1 = 0; d1 < dim; ++d1) {
            ptr[i * dim + d1 + 1] = ptr[i * dim + d1] + dim;
            for (int d2 = 0; d2 < dim; ++d2)
                row[i * dim * dim + d1 * dim + d2] = i * dim + d2;
        }
    }
    std::vector<T> val(dim * dim * BC_basis.size());
    memcpy(val.data(), reinterpret_cast<const T*>(BC_basis.data()), sizeof(T) * dim * dim * BC_basis.size());
    m_transform_matrix.setZero();
   // Math::sparse_from_csr(ptr, row, val, BC_basis.size() * dim, BC_basis.size() * dim, m_transform_matrix);
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::initialize_acceleration()
{
    Field<Vector<T, dim>> force(m_x.size(), Vector<T, dim>::Zero());
    for (auto energy : this->m_energy_terms) {
        Field<Vector<T, dim>> sub_force;
        energy->internal_force(m_x, m_v, sub_force);
        force -= sub_force;
    }
    tbb::parallel_for((size_t)0, m_x.size(), [&](size_t i) {
        if (BC_order[i] > 0) {
            force[i] = BC_basis[i].transpose() * force[i];
            for (int d = 0; d < BC_order[i]; ++d)
                force[i](d) = 0;
            force[i] = BC_basis[i] * force[i];
        }
        m_a[i] = force[i] / m_mass[i];
    });
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::update_predictive_pos()
{
    m_x_tilde.resize(m_x.size());
    tbb::parallel_for((size_t)0, m_x_tilde.size(), [&](size_t i) {
        m_x_tilde[i] = m_x[i] + m_v[i] * dt + tsParam[tsMethod][0] * (1 - 2 * tsParam[tsMethod][1]) * m_a[i] * dt * dt;
        if (BC_order[i] > 0) {
            m_x[i] = BC_basis[i].transpose() * m_x[i];
            for (int d = 0; d < BC_order[i]; ++d)
                m_x[i](d) = BC_target[i](d);
            m_x[i] = BC_basis[i] * m_x[i];
        }
    });
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::set_ts_weights()
{
    for (auto energy : this->m_energy_terms)
        if (dynamic_cast<InertialEnergyOp<T, dim>*>(energy)) {
            energy->energy_scale = 1;
        }
        else {
            energy->energy_scale = 2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt;
        }
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::gradient(const Vec& x, Vec& grad)
{
    Base::gradient(x, grad);
    grad = m_transform_matrix.transpose() * grad;
    if (project_dirichlet) {
        tbb::parallel_for((size_t)0, m_x.size(), [&](size_t i) {
            if (BC_order[i] > 0) {
                for (int d = 0; d < BC_order[i]; ++d)
                    grad(i * dim + d) = 0;
            }
        });
    }
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::hessian(const Vec& x, Mat& hess, const bool project_pd)
{
    Base::hessian(x, hess, project_pd);
    hess = m_transform_matrix.transpose() * hess * m_transform_matrix;
    if (project_dirichlet) {
        tbb::parallel_for(0, (int)m_x.size(), [&](int i) {
            for (int d = 0; d < dim; ++d) {
                bool clear_col = d < BC_order[i];
                for (typename Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>::InnerIterator it(hess, dim * i + d); it; ++it) {
                    bool clear_row = (it.row() % dim) < BC_order[it.row() / dim];
                    if (clear_row || clear_col) {
                        if (it.col() == it.row())
                            it.valueRef() = 1;
                        else
                            it.valueRef() = 0;
                    }
                }
            }
        });
    }
}

template <class T, int dim, class StorageIndex>
T BackwardEulerUpdateOp<T, dim, StorageIndex>::residual(const Vec& x, const Vec& grad, const Vec& direction)
{
    return direction.cwiseAbs().maxCoeff() / this->dt;
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::postprocess(Vec& direction)
{
    direction = m_transform_matrix * direction;
}

template <class T, int dim, class StorageIndex>
void BackwardEulerUpdateOp<T, dim, StorageIndex>::advance(const Vec& new_x)
{
    tbb::parallel_for(0, (int)m_x.size(), [&](int i) {
        m_x[i] = new_x.template segment<dim>(i * dim);
    });

    tbb::parallel_for(0, (int)m_x.size(), [&](int i) {
        Bow::Vector<T, dim> new_ai = (m_x[i] - m_x_tilde[i]) / (2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt);
        m_v[i] += dt * ((1 - tsParam[tsMethod][2]) * m_a[i] + tsParam[tsMethod][2] * new_ai);
        // m_a[i] = new_ai;
    });

}

template <class T, int dim, class _StorageIndex>
inline void BackwardEulerUpdateOp<T, dim, _StorageIndex>::operator()()
{
    update_transformation_matrix();
    initialize_acceleration();
    update_predictive_pos();
    set_ts_weights();
    Vector<T, Eigen::Dynamic> new_x = Eigen::Map<Vector<T, Eigen::Dynamic>>(reinterpret_cast<T*>(m_x.data()), dim * m_x.size());
    this->optimize(new_x);
    advance(new_x);
}
}
} // namespace Bow::FEM