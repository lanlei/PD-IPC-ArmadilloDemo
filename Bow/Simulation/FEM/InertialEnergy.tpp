#include "InertialEnergy.h"
#include <oneapi/tbb.h>

namespace Bow::FEM {
template <class T, int dim, class StorageIndex>
InertialEnergyOp<T, dim, StorageIndex>::InertialEnergyOp(const Field<T>& mass, const Field<Vector<T, dim>>& x_tilde, const T energy_scale)
    : m_mass(mass), m_x_tilde(x_tilde), energy_scale(energy_scale)
{
}

template <class T, int dim, class StorageIndex>
T InertialEnergyOp<T, dim, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
{
    T total_energy = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        total_energy += 0.5 * m_mass[i] * (x[i] - m_x_tilde[i]).squaredNorm();
    }
    return energy_scale * total_energy;
}

template <class T, int dim, class StorageIndex>
void InertialEnergyOp<T, dim, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
{
    grad.resize(x.size());
    std::fill(grad.begin(), grad.end(), Vector<T, dim>::Zero());
    tbb::parallel_for(size_t(0), x.size(), [&](size_t i) {
        grad[i] = energy_scale * m_mass[i] * (x[i] - m_x_tilde[i]);
    });
}

template <class T, int dim, class StorageIndex>
void InertialEnergyOp<T, dim, StorageIndex>::hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd)
{
    int nrows = x.size() * dim;
    hess.setZero();
    hess.derived().resize(nrows, nrows);
    hess.derived().reserve(nrows);
    tbb::parallel_for(0, nrows, [&](int i) {
        hess.derived().outerIndexPtr()[i + 1] = i + 1;
        hess.derived().innerIndexPtr()[i] = i;
        hess.derived().valuePtr()[i] = energy_scale * m_mass[i / dim];
    });
    hess.derived().finalize();
}

} // namespace Bow::FEM