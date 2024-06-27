#pragma once

#include "MPMTransfer.h"
#include <Bow/Simulation/MPM/MPMGrid.h>
#include <Bow/Math/LinearSolver/SparseQR.h>
#include <Bow/Math/LinearSolver/SparseCholesky.h>
#include <Bow/Physics/FixedCorotated.h>
#include <Bow/Utils/Timer.h>
#include <Bow/Utils/Logging.h>
#include <Bow/Utils/FiniteDiff.h>
#include <Bow/Optimization/Newton.h>
#include <tbb/tbb.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Page_Map.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

using namespace SPGrid;

namespace Bow {
namespace MPM {

template <class T, int dim, class StorageIndex = int, bool MatrixFree = false>
class BackwardEulerUpdateOp : public Optimization::Newton<T, dim, StorageIndex, MatrixFree> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    using Vec = Bow::Vector<T, Eigen::Dynamic>;
    using Mat = Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>;

    MPMGrid<T, dim>& grid;
    Geometry::BoundaryConditionManager<T, dim>& BC;

    T dx;
    T dt;

    enum TSMethod {
        BE,
        NM
    };
    TSMethod tsMethod = BE;
    T tsParam[2][3] = {
        { 1, 0.5, 1 },
        { 0.5, 0.25, 0.5 }
    };

    BackwardEulerUpdateOp(MPMGrid<T, dim>& grid, Geometry::BoundaryConditionManager<T, dim>& BC, T dx, T dt)
        : grid(grid), BC(BC), dx(dx), dt(dt) {}

    // https://www.overleaf.com/project/601dcb598a7b7415a29351c7
    // BC_basis is V in doc, which is [v_n, v_m, v_l]
    Field<Matrix<T, dim, dim>> BC_basis;
    std::vector<int> BC_order;

    void generate_BC()
    {
        BC_basis.assign(grid.num_nodes, TM::Identity());
        BC_order.assign(grid.num_nodes, 0);
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            BC.mpm_implicit_update(node.template cast<T>() * dx, BC_basis[g.idx], BC_order[g.idx]);
        });
    }

    void applyNewtonResult(const Field<Vector<T, dim>>& x_hat)
    {
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            g.x = x_hat[g.idx];
            TV x_n = node.template cast<T>() * dx;
            TV v_n = g.v_and_m.template topLeftCorner<dim, 1>();
            TV a_n = g.a;
            TV x_tilde = x_n + v_n * dt + tsParam[tsMethod][0] * (1 - 2 * tsParam[tsMethod][1]) * a_n * dt * dt;
            Bow::Vector<T, dim> new_ai = (g.x - x_tilde) / (2 * tsParam[tsMethod][0] * tsParam[tsMethod][1] * dt * dt);
            g.v_and_m.template topLeftCorner<dim, 1>() += dt * ((1 - tsParam[tsMethod][2]) * a_n + tsParam[tsMethod][2] * new_ai);
        });
    }

    T residual(const Vec& x, const Vec& grad, const Vec& direction)
    {
        return direction.cwiseAbs().maxCoeff() / this->dt;
    }

    void diff_test_with_matrix(const Vec& x)
    {
        const auto f = [&](const Eigen::VectorXd& x_vec) -> double {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            this->m_energy_terms[0]->precompute(x);
            return this->m_energy_terms[0]->energy(x);
        };
        const auto g = [&](const Eigen::VectorXd& x_vec, Eigen::VectorXd& grad_vec) {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            Bow::Field<Bow::Vector<T, dim>> grad;
            this->m_energy_terms[0]->precompute(x);
            this->m_energy_terms[0]->gradient(x, grad);
            grad_vec = Bow::to_vec(grad);
        };
        const auto h = [&](const Eigen::VectorXd& x_vec, Eigen::SparseMatrix<T>& hess) {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            this->m_energy_terms[0]->precompute(x);
            this->m_energy_terms[0]->hessian(x, hess, false);
        };
        const auto project = [&](Eigen::VectorXd& step_vec) {
            Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
            grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
                TV x_n = node.template cast<T>() * dx;
                for (int d = 0; d < BC_order[g.idx]; ++d) {
                    TV n = BC_basis[g.idx].col(d);
                    step[g.idx] -= step[g.idx].dot(n) * n;
                }
            });
            step_vec = to_vec(step);
        };
        const auto transform = [&](Eigen::VectorXd& step_vec) {
            Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
            grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
                step[g.idx] = BC_basis[g.idx].transpose() * step[g.idx];
            });
            step_vec = to_vec(step);
        };
        FiniteDiff::ziran_check_false(x, f, g, h, project, transform);
    }

    void diff_test_matrix_free(const Vec& x)
    {
        const auto f = [&](const Eigen::VectorXd& x_vec) -> double {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            this->m_energy_terms[0]->precompute(x);
            return this->m_energy_terms[0]->energy(x);
        };
        const auto g = [&](const Eigen::VectorXd& x_vec, Eigen::VectorXd& grad_vec) {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            Bow::Field<Bow::Vector<T, dim>> grad;
            this->m_energy_terms[0]->precompute(x);
            this->m_energy_terms[0]->gradient(x, grad);
            grad_vec = Bow::to_vec(grad);
        };
        const auto h = [&](const Eigen::VectorXd& x_vec, const Eigen::VectorXd& dx_vec, Eigen::VectorXd& Adx_vec) {
            Bow::Field<Bow::Vector<T, dim>> x = Bow::to_field<dim>(x_vec);
            Bow::Field<Bow::Vector<T, dim>> dx = Bow::to_field<dim>(dx_vec);
            Bow::Field<Bow::Vector<T, dim>> Adx;
            this->m_energy_terms[0]->precompute(x);
            this->m_energy_terms[0]->multiply(dx, Adx, false);
            this->m_energy_terms[0]->project(Adx);
            Adx_vec = to_vec(Adx);
        };
        const auto project = [&](Eigen::VectorXd& step_vec) {
            Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
            grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
                TV x_n = node.template cast<T>() * dx;
                for (int d = 0; d < BC_order[g.idx]; ++d) {
                    TV n = BC_basis[g.idx].col(d);
                    step[g.idx] -= step[g.idx].dot(n) * n;
                }
            });
            step_vec = to_vec(step);
        };
        const auto transform = [&](Eigen::VectorXd& step_vec) {
            Bow::Field<Bow::Vector<T, dim>> step = Bow::to_field<dim>(step_vec);
            grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
                step[g.idx] = BC_basis[g.idx].transpose() * step[g.idx];
            });
            step_vec = to_vec(step);
        };
        FiniteDiff::ziran_check_true(x, f, g, h, project, transform);
    }

    void initialize_acceleration()
    {
        Field<Vector<T, dim>> xn(grid.num_nodes, Vector<T, dim>::Zero());
        Field<Vector<T, dim>> vn(grid.num_nodes, Vector<T, dim>::Zero());
        Field<Vector<T, dim>> force(grid.num_nodes, Vector<T, dim>::Zero());
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            xn[g.idx] = node.template cast<T>() * dx;
            vn[g.idx] = g.v_and_m.template segment<dim>(0);
        });

        for (auto energy : this->m_energy_terms) {
            Field<Vector<T, dim>> sub_force;
            energy->internal_force(xn, vn, sub_force);
            force -= sub_force;
        }

        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            if (g.idx < 0) return;
            if (BC_order[g.idx] > 0) {
                force[g.idx] = BC_basis[g.idx].transpose() * force[g.idx];
                for (int d = 0; d < BC_order[g.idx]; ++d)
                    force[g.idx](d) = 0;
                force[g.idx] = BC_basis[g.idx] * force[g.idx];
            }
            g.a = force[g.idx] / g.v_and_m(dim);
        });
    }

    void operator()()
    {
        BOW_TIMER_FLAG("newton one step");
        Field<Vector<T, dim>> x_hat(grid.num_nodes, Vector<T, dim>::Zero());
        generate_BC();
        initialize_acceleration();
        // modify initial v_n to satisfy BC
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            TV x_n = node.template cast<T>() * dx;
            x_hat[g.idx] = x_n + g.v_and_m.template head<dim>() * dt;
            for (int d = 0; d < BC_order[g.idx]; ++d) {
                TV v = (x_hat[g.idx] - x_n) / dt;
                TV n = BC_basis[g.idx].col(d);
                v -= v.dot(n) * n;
                x_hat[g.idx] = x_n + v * dt;
            }
        });
        Vec x = to_vec(x_hat);
        this->optimize(x);
        x_hat = to_field<dim>(x);
        applyNewtonResult(x_hat);
    }
};
}
} // namespace Bow::MPM