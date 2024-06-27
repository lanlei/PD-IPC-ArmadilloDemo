#pragma once

#include <Bow/Simulation/MPM/MPMGrid.h>
#include <Bow/IO/ply.h>
#include <Bow/Simulation/MPM/ParticlesLevelSet.h>
#include <Bow/Math/LinearSolver/ConjugateGradient.h>
#include <Bow/Physics/FixedCorotated.h>
#include <Bow/Utils/Timer.h>
#include <Bow/Utils/Logging.h>
#include <Bow/Geometry/AnalyticalLevelSet.h>
#include <tbb/tbb.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Page_Map.h>

using namespace SPGrid;

namespace Bow {
namespace MPM {

class AbstractOp {
};

template <class T, int dim, bool symplectic = true>
class ParticlesToGridOp : public AbstractOp {
public:
    using SparseMask = typename MPMGrid<T, dim>::SparseMask;
    Field<Vector<T, dim>>& m_X;
    Field<Vector<T, dim>>& m_V;
    std::vector<T>& m_mass;
    Field<Matrix<T, dim, dim>>& m_C;
    Field<Matrix<T, dim, dim>>& stress;

    MPMGrid<T, dim>& grid;
    T dx;
    T dt;

    void operator()()
    {
        BOW_TIMER_FLAG("P2G");
        grid.colored_for([&](int i) {
            const Vector<T, dim> pos = m_X[i];
            const Vector<T, dim> v = m_V[i];
            const T mass = m_mass[i];
            const Matrix<T, dim, dim> C = m_C[i] * mass;
            const Vector<T, dim> momentum = mass * v;
            const Matrix<T, dim, dim> delta_t_tmp_force = -dt * stress[i];
            BSplineWeights<T, dim> spline(pos, dx);
            grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, const Vector<T, dim>& dw, GridState<T, dim>& g) {
                Vector<T, dim> xi_minus_xp = node.template cast<T>() * dx - pos;
                Vector<T, dim + 1> velocity_term = Vector<T, dim + 1>::Zero();
                velocity_term.template topLeftCorner<dim, 1>() = momentum + C * xi_minus_xp;
                velocity_term(dim) = mass;
                Vector<T, dim + 1> stress_term_dw = Vector<T, dim + 1>::Zero();
                if constexpr (symplectic)
                    stress_term_dw.template topLeftCorner<dim, 1>() = delta_t_tmp_force * dw;
                Vector<T, dim + 1> delta = w * velocity_term + stress_term_dw;
                g.v_and_m += delta;
            });
        });
        grid.countNumNodes();
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            T mass = g.v_and_m(dim);
            Vector<T, dim + 1> alpha;
            alpha.template topLeftCorner<dim, 1>() = Vector<T, dim>::Ones() * ((T)1 / mass);
            alpha(dim) = T(1);
            g.v_and_m = g.v_and_m.cwiseProduct(alpha);
            g.old_v = g.v_and_m.template segment<dim>(0);
        });
    }
};

template <class T, int dim>
class BoundaryConditionUpdateOp : public AbstractOp {
public:
    using SparseMask = typename MPMGrid<T, dim>::SparseMask;
    MPMGrid<T, dim>& grid;
    Vector<T, dim>& gravity;
    Geometry::BoundaryConditionManager<T, dim>& BC;
    T dx;
    T dt;

    void operator()()
    {
        BOW_TIMER_FLAG("grid update");
        grid.iterateGrid([&](const Vector<int, dim>& node, GridState<T, dim>& g) {
            g.v_and_m.template head<dim>() += gravity * dt;
            Vector<T, dim> new_v = g.v_and_m.template head<dim>();
            BC.mpm_explicit_update(node.template cast<T>() * dx, new_v);
            g.v_and_m.template head<dim>() = new_v;
        });
    }
};

template <class T, int dim>
class GridToParticlesOp : public AbstractOp {
public:
    Field<Vector<T, dim>>& m_X;
    Field<Vector<T, dim>>& m_V;
    Field<Matrix<T, dim, dim>>& m_C;

    MPMGrid<T, dim>& grid;
    T dx;
    T dt;
    T flip_pic_ratio = 0.98;

    Field<Matrix<T, dim, dim>> m_gradXp;

    template <bool apic>
    void grid_to_particle()
    {
        BOW_TIMER_FLAG("G2P");
        T D_inverse = (T)4 / (dx * dx);
        m_gradXp.assign(m_X.size(), Matrix<T, dim, dim>());
        grid.parallel_for([&](int i) {
            Vector<T, dim>& Xp = m_X[i];
            Vector<T, dim> picV = Vector<T, dim>::Zero();
            BSplineWeights<T, dim> spline(Xp, dx);

            Matrix<T, dim, dim> Bp = Matrix<T, dim, dim>::Zero();
            Matrix<T, dim, dim> gradXp = Matrix<T, dim, dim>::Zero();
            Vector<T, dim> oldV = Vector<T, dim>::Zero();
            Vector<T, dim> newX = Vector<T, dim>::Zero();
            grid.iterateKernel(spline, [&](const Vector<int, dim>& node, T w, Vector<T, dim> dw, GridState<T, dim>& g) {
                Vector<T, dim> new_v = g.v_and_m.template topLeftCorner<dim, 1>();
                picV += w * new_v;
                oldV += w * g.old_v;
                newX += w * g.x;
                if constexpr (apic) {
                    Vector<T, dim> xi_minus_xp = node.template cast<T>() * dx - Xp;
                    Bp.noalias() += w * new_v * xi_minus_xp.transpose();
                }
                gradXp.noalias() += g.x * dw.transpose();
            });
            if constexpr (apic) {
                m_C[i] = Bp * D_inverse;
                m_V[i] = picV;
                m_X[i] += picV * dt;
            }
            else {
                m_C[i].setZero();
                m_V[i] *= flip_pic_ratio;
                m_V[i] += picV - flip_pic_ratio * oldV;
                m_X[i] = newX;
            }
            m_gradXp[i] = gradXp;
        });
    }

    void operator()(bool apic = true)
    {
        if (apic)
            grid_to_particle<true>();
        else
            grid_to_particle<false>();
    }
};
}
} // namespace Bow::MPM