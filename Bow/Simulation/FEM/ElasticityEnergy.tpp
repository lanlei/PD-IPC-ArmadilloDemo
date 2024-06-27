#include "ElasticityEnergy.h"
#include "Utils.h"
#include<numeric>

namespace Bow::FEM {
template <class T, int dim, class StorageIndex>
ElasticityEnergyOp<T, dim, StorageIndex>::ElasticityEnergyOp(const Field<Vector<int, dim + 1>>& elem, const Field<T>& vol, const Field<T>& mu, const Field<T>& lam, const Field<Matrix<T, dim, dim>>& IB, std::vector<std::pair<int, int>>& offsets, T energy_scale)
    : m_elem(elem)
    , m_vol(vol)
    , m_mu(mu)
    , m_lam(lam)
    , m_IB(IB)
    , m_offsets(offsets)
{
    this->energy_scale = energy_scale;
}

template <class T, int dim, class StorageIndex>
T ElasticityEnergyOp<T, dim, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
{
    T total_energy = 0.0;
    for (auto range : m_offsets) {
        Field<T> elem_energy(range.second - range.first);
        tbb::parallel_for(range.first, range.second, [&](int e) {
            Matrix<T, dim, dim> F;
            const auto& vertices = m_elem[e];
            if constexpr (dim == 2)
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], m_IB[e], F);
            else
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], x[vertices[3]], m_IB[e], F);
            elem_energy[e - range.first] = m_vol[e] * this->psi(F, m_mu[e], m_lam[e]);
        });
        total_energy += this->energy_scale * std::accumulate(elem_energy.begin(), elem_energy.end(), T(0));
    }
    return total_energy;
}

template <class T, int dim, class StorageIndex>
void ElasticityEnergyOp<T, dim, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
{
    grad.resize(x.size());
    std::fill(grad.begin(), grad.end(), Vector<T, dim>::Zero());

    for (auto range : m_offsets) {
        // compute local gradients in parallel
        Field<Vector<T, dim*(dim + 1)>> local_grad(range.second - range.first);
        tbb::parallel_for(range.first, range.second, [&](int e) {
            // deformation gradient
            Matrix<T, dim, dim> F;
            const auto& vertices = m_elem[e];
            if constexpr (dim == 2)
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], m_IB[e], F);
            else
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], x[vertices[3]], m_IB[e], F);
            // first piola
            Matrix<T, dim, dim> P;
            this->first_piola(F, m_mu[e], m_lam[e], P);
            // backpropagate
            backpropagate_element_gradient(m_IB[e], P, local_grad[e - range.first]);
            local_grad[e - range.first] *= this->energy_scale * m_vol[e];
        });

        // assemble global gradient
        for (int e = range.first; e < range.second; ++e) {
            const auto& vertices = m_elem[e];
            for (int local_index = 0; local_index < dim + 1; ++local_index) {
                grad[vertices(local_index)] += local_grad[e - range.first].template segment<dim>(local_index * dim);
            }
        }
    }
}

template <class T, int dim, class StorageIndex>
template <bool project_pd>
void ElasticityEnergyOp<T, dim, StorageIndex>::hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess)
{
    using IJK = Eigen::Triplet<T, StorageIndex>;
    static const int nrow = (dim + 1) * dim;
    static const int nelem = nrow * nrow;
    std::vector<IJK> coeffs;
    for (auto range : m_offsets) {
        const int coeff_ind0 = coeffs.size();
        coeffs.resize(coeff_ind0 + (range.second - range.first) * nelem);
        tbb::parallel_for(range.first, range.second, [&](int e) {
            // deformation gradient
            Matrix<T, dim, dim> F;
            const auto& vertices = m_elem[e];
            if constexpr (dim == 2)
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], m_IB[e], F);
            else
                deformation_gradient(x[vertices[0]], x[vertices[1]], x[vertices[2]], x[vertices[3]], m_IB[e], F);
            // first piola derivative
            Matrix<T, dim * dim, dim * dim> dPdF;
            this->first_piola_derivative<project_pd>(F, m_mu[e], m_lam[e], dPdF);
            Matrix<T, nrow, nrow> local_hessian;
            backpropagate_element_hessian(m_IB[e], dPdF, local_hessian);
            local_hessian *= this->energy_scale * m_vol[e];
            int indMap[nrow];
            for (int i = 0; i < dim + 1; ++i)
                for (int d = 0; d < dim; ++d)
                    indMap[i * dim + d] = vertices[i] * dim + d;
            for (int row = 0; row < nrow; ++row)
                for (int col = 0; col < nrow; ++col) {
                    coeffs[coeff_ind0 + (e - range.first) * nelem + row * nrow + col] = std::move(IJK(indMap[row], indMap[col], local_hessian(row, col)));
                }
        });
    }
    hess.derived().resize(x.size() * dim, x.size() * dim);
    hess.derived().setZero();
    hess.derived().setFromTriplets(coeffs.begin(), coeffs.end());
}

} // namespace Bow::FEM