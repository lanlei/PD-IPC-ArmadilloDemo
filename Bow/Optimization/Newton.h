#ifndef OPTIMIZER_NEWTON_H
#define OPTIMIZER_NEWTON_H

#include <Bow/Macros.h>
#include <Eigen/Eigen>
#include <Bow/Utils/Timer.h>
#include <Eigen/SparseCholesky>
#include <Bow/Utils/Logging.h>
#include <Bow/Math/LinearSolver/SparseCholesky.h>
#include <Bow/Math/LinearSolver/ConjugateGradient.h>
#include <Bow/Simulation/Prototypes.h>
#include <functional>
#include <Bow/Types.h>

namespace Bow {
namespace Optimization {
template <class Scalar, int dim, class StorageIndex, bool MatrixFree = false>
class Newton {
public:
    using Vec = Bow::Vector<Scalar, Eigen::Dynamic>;
    using Mat = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex>;

    bool line_search = true;
    int max_iter = 1000;
    Scalar tol = 1e-3;
    std::vector<EnergyOp<Scalar, dim>*> m_energy_terms;

    virtual Scalar energy(const Vec& x_vec)
    {
        Scalar total_energy = 0.0;
        auto x = to_field<dim>(x_vec);
		for (auto e : m_energy_terms)
		{
			total_energy += e->energy(x);
		}
        return total_energy;
    }

    virtual void gradient(const Vec& x_vec, Vec& grad_vec)
    {
        grad_vec.resize(x_vec.size());
        grad_vec.setZero();
        auto x = to_field<dim>(x_vec);
        for (auto e : m_energy_terms) {
            Field<Vector<Scalar, dim>> sub_grad;
            e->gradient(x, sub_grad);
            grad_vec += to_vec(sub_grad);
        }
    }

    virtual void hessian(const Vec& x_vec, Mat& hess, const bool project_pd)
    {
        hess.derived().resize(x_vec.size(), x_vec.size());
        hess.derived().setZero();
        auto x = to_field<dim>(x_vec);
        for (auto e : m_energy_terms) {
            Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex> sub_hess;
            e->hessian(x, sub_hess, project_pd);
            hess += sub_hess;
        }
    }

    virtual void multiply(const Vec& x_vec, Vec& Ax_vec) const
    {
        Ax_vec.resize(x_vec.size());
        Ax_vec.setZero();
        auto x = to_field<dim>(x_vec);
        for (auto e : m_energy_terms) {
            Field<Vector<Scalar, dim>> sub_Ax;
            e->multiply(x, sub_Ax);
            Ax_vec += to_vec(sub_Ax);
        }
    }
    virtual void precondition(const Vec& in_vec, Vec& out_vec) const
    {
        Vec diagonal_vec = in_vec;
        diagonal_vec.setZero();
        for (auto e : m_energy_terms) {
            Field<Vector<Scalar, dim>> diagonal;
            e->precondition(diagonal);
            diagonal_vec += to_vec(diagonal);
        }
        out_vec = in_vec.cwiseQuotient(diagonal_vec);
    }
    virtual void project(Vec& b_vec) const
    {
        auto b = to_field<dim>(b_vec);
        for (auto e : m_energy_terms)
            e->project(b);
        b_vec = to_vec(b);
    }

    virtual Scalar residual(const Vec& x, const Vec& grad, const Vec& direction) = 0;

    virtual void postprocess(Vec& direction_vec)
    {
        auto direction = to_field<dim>(direction_vec);
        for (auto e : m_energy_terms)
            e->postprocess(direction);
        direction_vec = to_vec(direction);
    }

    virtual void callback(const Vec& x_vec)
    {
        auto x = to_field<dim>(x_vec);
        for (auto e : m_energy_terms)
            e->callback(x);
    }

    virtual Scalar initial_stepsize(const Vec& x_vec, const Vec& dx_vec)
    {
        auto x = to_field<dim>(x_vec);
        auto dx = to_field<dim>(dx_vec);
        Scalar upper_bound = 1.0;
        for (auto e : m_energy_terms)
            upper_bound = std::min(upper_bound, e->stepsize_upperbound(x, dx));
        return upper_bound;
    }

    virtual void precompute(const Vec& x_vec)
    {
        auto x = to_field<dim>(x_vec);
        for (auto e : m_energy_terms)
            e->precompute(x);
    }

    Newton()
    {
    }
    /**
     * x is modified in place.
     */
    template <class DerivedX>
    int optimize(Eigen::MatrixBase<DerivedX>& x)
    {
        BOW_TIMER_FLAG("Newton");
        int newton_iter = 0;
        Bow::Vector<Scalar, Eigen::Dynamic> xn = x;
        precompute(xn);
        for (; newton_iter < max_iter; ++newton_iter) {
            callback(xn);

            Mat hess;
            if constexpr (!MatrixFree) {
                BOW_TIMER_FLAG("Compute Hessian");
                hessian(xn, hess, true);
            }

            Vec grad;
            gradient(xn, grad);

            Vec direction;
            if constexpr (MatrixFree) {
                BOW_TIMER_FLAG("Matrix Free Solve");
                Scalar linear_solve_tolerance_scale = 1;
                Scalar residual_norm = grad.norm();
                Scalar tolerance = 1e-5;
                Scalar linear_solve_relative_tolerance = std::min((Scalar)0.5, linear_solve_tolerance_scale * std::sqrt(std::max(residual_norm, tolerance)));
                ConjugateGradient<Scalar, Newton<Scalar, dim, StorageIndex, MatrixFree>, Vec> cg(10000);
                cg.setTolerance(1);
                cg.setRelativeTolerance(linear_solve_relative_tolerance);
                direction = x;
                direction.setZero();
                cg.solve(*this, direction, -grad);
            }
            else {
#ifdef BOW_SUITESPARSE
                BOW_TIMER_FLAG("Build Matrix Cholmod Solve");
                Math::LinearSolver::CholmodLLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex>> solver;
#else
                BOW_TIMER_FLAG("Build Matrix Eigen Solve");
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex>> solver;
#endif
                solver.compute(hess);
                direction = -solver.solve(grad);
            }
            postprocess(direction);

            Scalar res = residual(xn, grad, direction);
            if (newton_iter > 0 && res < tol) {
                static int total_newton_iter = 0;
                total_newton_iter += newton_iter;
                Bow::Logging::info("Converged Newton iter: ", newton_iter, ",\tResidual: ", res, ",\tTotal iter: ", total_newton_iter);
                break;
            }

            {
                BOW_TIMER_FLAG("Line Search");
                Scalar E0 = energy(xn);
                Scalar alpha = initial_stepsize(xn, direction);
                Scalar alpha0 = alpha;
                Vec new_x = xn + alpha * direction;
                precompute(new_x);
                if (line_search) {
                    Scalar E = energy(new_x);
                    while (E > E0) {
                        alpha *= 0.5;
                        new_x = xn + alpha * direction;
                        precompute(new_x);
                        E = energy(new_x);
                    }
                }
                if (alpha < 1e-10) {
                    new_x = xn + direction;
                }
                xn = new_x;
                Bow::Logging::info("Newton iter: ", newton_iter, ",\tResidual: ", res, ",\tLine search: ", alpha, ",\tInitial search: ", alpha0);
            }
        }
        x = xn;
        return newton_iter;
    }
};
}
} // namespace Bow::Optimization

#endif