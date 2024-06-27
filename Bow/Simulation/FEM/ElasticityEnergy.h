#ifndef FEM_ELASTICITY_ENERGY_H
#define FEM_ELASTICITY_ENERGY_H

#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <Eigen/Sparse>
#include "../Prototypes.h"
#include <Bow/Physics/ConstitutiveModel.h>
#include <Bow/Physics/FixedCorotated.h>
#include <Bow/Physics/NeoHookean.h>
#include <Bow/Physics/LinearElasticity.h>
#include <Bow/Physics/AsRigidAsPossible.h>
#include <Bow/Math/Utils.h>

namespace Bow::FEM {
template <class T, int dim, class StorageIndex = int>
class ElasticityEnergyOp : public EnergyOp<T, dim, StorageIndex> {
public:
    const Field<Vector<int, dim + 1>>& m_elem;
    const Field<T>&m_vol, m_mu, m_lam;
    const Field<Matrix<T, dim, dim>>& m_IB;
    std::vector<std::pair<int, int>>& m_offsets;
    ElasticityEnergyOp(const Field<Vector<int, dim + 1>>& elem, const Field<T>& vol, const Field<T>& mu, const Field<T>& lam, const Field<Matrix<T, dim, dim>>& IB, std::vector<std::pair<int, int>>& offsets, T energy_scale = 1.0);
    T energy(const Field<Vector<T, dim>>& x) override;
    void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
    template <bool project_pd = true>
    void hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess);
    void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, const bool project_pd = true) override
    {
        if (project_pd)
            hessian_impl<true>(x, hess);
        else
            hessian_impl<false>(x, hess);
    }
    void internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force)
    {
        force.resize(xn.size());
        std::fill(force.begin(), force.end(), Vector<T, dim>::Zero());
        this->energy_scale = 1.0;
        gradient(xn, force);
        force *= T(1) / this->energy_scale;
    }

protected:
    virtual T psi(const Matrix<T, dim, dim>& F, const T mu, const T lam) { return 0; }
    virtual void first_piola(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim, dim>& P) {}
    template <bool project_pd = true>
    void first_piola_derivative(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP)
    {
        if constexpr (project_pd)
            first_piola_derivative_pd(F, mu, lam, dP);
        else
            first_piola_derivative_npd(F, mu, lam, dP);
    }
    virtual void first_piola_derivative_pd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP)
    {
    }
    virtual void first_piola_derivative_npd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP)
    {
    }
};

template <class T, int dim, class StorageIndex = int>
class FixedCorotatedEnergyOp : public ElasticityEnergyOp<T, dim, StorageIndex> {
public:
    using ElasticityEnergyOp<T, dim, StorageIndex>::ElasticityEnergyOp;
    T psi(const Matrix<T, dim, dim>& F, const T mu, const T lam) override
    {
        return ConstitutiveModel::FixedCorotated::psi(F, mu, lam);
    }
    void first_piola(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim, dim>& P) override
    {
        ConstitutiveModel::FixedCorotated::first_piola(F, mu, lam, P);
    }
    void first_piola_derivative_pd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::FixedCorotated::first_piola_derivative<true>(F, mu, lam, dP);
    }
    void first_piola_derivative_npd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::FixedCorotated::first_piola_derivative<false>(F, mu, lam, dP);
    }
};

template <class T, int dim, class StorageIndex = int>
class NeoHookeanEnergyOp : public ElasticityEnergyOp<T, dim, StorageIndex> {
public:
    using ElasticityEnergyOp<T, dim, StorageIndex>::ElasticityEnergyOp;
    using ElasticityEnergyOp<T, dim, StorageIndex>::m_elem;
    using ElasticityEnergyOp<T, dim, StorageIndex>::m_offsets;
    T psi(const Matrix<T, dim, dim>& F, const T mu, const T lam) override
    {
        return ConstitutiveModel::NeoHookean::psi(F, mu, lam);
    }
    void first_piola(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim, dim>& P) override
    {
        ConstitutiveModel::NeoHookean::first_piola(F, mu, lam, P);
    }
    void first_piola_derivative_pd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::NeoHookean::first_piola_derivative<true>(F, mu, lam, dP);
    }
    void first_piola_derivative_npd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::NeoHookean::first_piola_derivative<false>(F, mu, lam, dP);
    }
    T stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& dx) override
    {
        // https://www.dropbox.com/scl/fi/pdlkk5uzvlrpp9bd85ygo/Non-invertable-stepsize.paper?dl=0&rlkey=pm8gy8wqxset4rnt22rrfw1my
        Vector<T, Eigen::Dynamic> alphas(m_elem.size());
        alphas.setOnes();
        for (auto range : m_offsets)
            for (int e = range.first; e < range.second; ++e) {
                const auto& vert = m_elem[e];
                Matrix<T, dim, dim> basis, dbasis;
                basis.row(0) = x[vert[1]] - x[vert[0]];
                dbasis.row(0) = dx[vert[1]] - dx[vert[0]];
                basis.row(1) = x[vert[2]] - x[vert[0]];
                dbasis.row(1) = dx[vert[2]] - dx[vert[0]];
                if constexpr (dim == 3) {
                    basis.row(2) = x[vert[3]] - x[vert[0]];
                    dbasis.row(2) = dx[vert[3]] - dx[vert[0]];
                }
                if (dbasis.norm() == 0) continue;
                Matrix<T, dim, dim> A = basis.partialPivLu().solve(dbasis);

                T a, b, c, d;
                if constexpr (dim == 2) {
                    a = 0;
                    b = A.determinant();
                }
                else {
                    a = A.determinant();
                    b = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0) + A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) + A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
                }
                c = A.diagonal().sum();
                d = 0.9;

                T t = Math::get_smallest_positive_real_cubic_root(a, b, c, d);
                if (t < 0 || t > 1) t = 1;
                alphas(e) = t;
            }
        return alphas.minCoeff();
    }
};

template <class T, int dim, class StorageIndex = int>
class LinearElasticityEnergyOp : public ElasticityEnergyOp<T, dim, StorageIndex> {
public:
    using ElasticityEnergyOp<T, dim, StorageIndex>::ElasticityEnergyOp;
    T psi(const Matrix<T, dim, dim>& F, const T mu, const T lam) override
    {
        return ConstitutiveModel::LinearElasticity::psi(F, mu, lam);
    }
    void first_piola(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim, dim>& P) override
    {
        ConstitutiveModel::LinearElasticity::first_piola(F, mu, lam, P);
    }
    void first_piola_derivative_pd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::LinearElasticity::first_piola_derivative<true>(F, mu, lam, dP);
    }
    void first_piola_derivative_npd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
    {
        ConstitutiveModel::LinearElasticity::first_piola_derivative<false>(F, mu, lam, dP);
    }
};

template <class T, int dim, class StorageIndex = int>
class AsRigidAsPossibleEnergyOp : public ElasticityEnergyOp<T, dim, StorageIndex> {
public:
	using ElasticityEnergyOp<T, dim, StorageIndex>::ElasticityEnergyOp;
	T psi(const Matrix<T, dim, dim>& F, const T mu, const T lam) override
	{
		return ConstitutiveModel::AsRigidAsPossible::psi(F, mu, T(0.0));
	}
	void first_piola(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim, dim>& P) override
	{
		ConstitutiveModel::AsRigidAsPossible::first_piola(F, mu, T(0.0), P);
	}
	void first_piola_derivative_pd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
	{
		ConstitutiveModel::AsRigidAsPossible::first_piola_derivative<true>(F, mu, T(0.0), dP);
	}
	void first_piola_derivative_npd(const Matrix<T, dim, dim>& F, const T mu, const T lam, Matrix<T, dim * dim, dim * dim>& dP) override
	{
		ConstitutiveModel::AsRigidAsPossible::first_piola_derivative<false>(F, mu, T(0.0), dP);
	}
};



} // namespace Bow::FEM

#include "ElasticityEnergy.tpp"

#endif
