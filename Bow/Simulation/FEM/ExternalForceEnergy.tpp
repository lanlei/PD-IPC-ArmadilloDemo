#include "ExternalForceEnergy.h"

namespace Bow::FEM {
template <class T, int dim, class StorageIndex>
StaticForceEnergy<T, dim, StorageIndex>::StaticForceEnergy(const Field<Vector<T, dim>>& f, const Field<T>& mass, const Vector<T, dim>& gravity, const T energy_scale)
    : m_f(f), m_mass(mass), gravity(gravity)
{
    this->energy_scale = energy_scale;
    zero_potential_ref.resize(m_f.size(), Vector<T, dim>::Zero());
}

template <class T, int dim, class StorageIndex>
T StaticForceEnergy<T, dim, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
{
    T total_energy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        total_energy -= m_mass[i] * gravity.dot(x[i] - zero_potential_ref[i]);
    }
    total_energy -= dotProduct(x, m_f);
    return this->energy_scale * total_energy;
}
template <class T, int dim, class StorageIndex>
void StaticForceEnergy<T, dim, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
{
    grad.resize(x.size());
    tbb::parallel_for(size_t(0), x.size(), [&](size_t i) {
        grad[i] = -this->energy_scale * m_mass[i] * gravity;
    });
    grad -= this->energy_scale * m_f;
}
template <class T, int dim, class StorageIndex>
void StaticForceEnergy<T, dim, StorageIndex>::hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd)
{
    hess.derived().resize(x.size() * dim, x.size() * dim);
    hess.derived().setZero();
}
} // namespace Bow::FEM