#ifndef FEM_NODAL_FORCE_ENERGY
#define FEM_NODAL_FORCE_ENERGY
#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <oneapi/tbb.h>
#include <Eigen/Sparse>
#include "../Prototypes.h"

namespace Bow::FEM {

template <class T, int dim, class StorageIndex = int>
class StaticForceEnergy : public EnergyOp<T, dim, StorageIndex> {
public:
    const Field<Vector<T, dim>>& m_f;
    const Field<T>& m_mass;
    const Vector<T, dim>& gravity;
    Field<Vector<T, dim>> zero_potential_ref;
    StaticForceEnergy(const Field<Vector<T, dim>>& f, const Field<T>& mass, const Vector<T, dim>& gravity, const T energy_scale = 1.0);
    T energy(const Field<Vector<T, dim>>& x) override;
    void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
    void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd = true) override;
    void internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force)
    {
        gradient(xn, force);
        force *= T(1) / this->energy_scale;
    }
    void callback(const Field<Vector<T, dim>>& xn) override {};
};
} // namespace Bow::FEM

#include "ExternalForceEnergy.tpp"
#endif