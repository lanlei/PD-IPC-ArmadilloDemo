#ifndef FEM_INERTIAL_ENERGY_H
#define FEM_INERTIAL_ENERGY_H

#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <oneapi/tbb.h>
#include <Eigen/Sparse>
#include "../Prototypes.h"

namespace Bow::FEM {
template <class T, int dim, class StorageIndex = int>
class InertialEnergyOp : public EnergyOp<T, dim, StorageIndex> {
public:
    const Field<T>& m_mass;
    const Field<Vector<T, dim>>& m_x_tilde;
    T energy_scale;
    InertialEnergyOp(const Field<T>& mass, const Field<Vector<T, dim>>& x_tilde, const T energy_scale = 1.0);
    T energy(const Field<Vector<T, dim>>& x) override;
    void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
    void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd = true) override;
};
} // namespace Bow::FEM

#include "InertialEnergy.tpp"

#endif
