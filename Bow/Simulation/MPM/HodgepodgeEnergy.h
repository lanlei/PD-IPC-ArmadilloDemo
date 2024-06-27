#ifndef HODGEPODGE_ENERGY_H
#define HODGEPODGE_ENERGY_H

#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <tbb/tbb.h>
#include <Eigen/Sparse>
#include <Bow/Simulation/Prototypes.h>
#include <Bow/Simulation/MPM/ElasticityOp.h>

namespace Bow::MPM {

template <class T, int dim, class StorageIndex = int>
class HodgepodgeEnergy : public EnergyOp<T, dim, StorageIndex> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    MPMGrid<T, dim>& grid;
    Vector<T, dim> gravity;
    T dx;
    T dt;

    Field<Matrix<T, dim, dim>>& BC_basis;
    std::vector<int>& BC_order;

    // particles
    Field<Vector<T, dim>>& m_X;
    std::vector<std::shared_ptr<ElasticityOp<T, dim>>>& elasticity_models;

    enum TSMethod {
        BE,
        NM
    };
    TSMethod tsMethod = BE;
    T tsParam[2][3] = {
        { 1, 0.5, 1 },
        { 0.5, 0.25, 0.5 }
    };

    HodgepodgeEnergy(MPMGrid<T, dim>& grid, Vector<T, dim> gravity, T dx, T dt,
        Field<Matrix<T, dim, dim>>& BC_basis, std::vector<int>& BC_order, Field<Vector<T, dim>>& m_X,
        std::vector<std::shared_ptr<ElasticityOp<T, dim>>>& elasticity_models)
        : grid(grid), gravity(gravity), dx(dx), dt(dt), BC_basis(BC_basis), BC_order(BC_order), m_X(m_X), elasticity_models(elasticity_models) {}

    static inline int kernelSize()
    {
        if constexpr (dim == 2) {
            return 25;
        }
        else {
            return 125;
        }
    }

    static inline int kernelOffset(const Vector<int, dim>& dnode)
    {
        if constexpr (dim == 2) {
            return (dnode(0) + 2) * 5 + (dnode(1) + 2);
        }
        else {
            return (dnode(0) + 2) * 25 + (dnode(1) + 2) * 5 + (dnode(2) + 2);
        }
    }

    void precompute(const Field<Vector<T, dim>>& x) override;
    T energy(const Field<Vector<T, dim>>& x) override;
    void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
    void multiply(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& Ax, bool project_pd) override;
    void precondition(Field<Vector<T, dim>>& diagonal) override;
    void project(Field<Vector<T, dim>>& b) override;
    void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd) override;
    void postprocess(Field<Vector<T, dim>>& direction) override;
    void internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force) override;
    T stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& dx) override;
};

} // namespace Bow::MPM

#include "HodgepodgeEnergy.tpp"

#endif
