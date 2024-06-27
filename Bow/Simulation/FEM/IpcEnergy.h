#ifndef IPC_OP_H
#define IPC_OP_H
#include <Bow/Types.h>
#include <Bow/Simulation/FEM/FrictionUtils.h>
#include "../Prototypes.h"

namespace Bow::FEM::IPC {
	template <class T, class StorageIndex = int>
	class IpcEnergyOp2D : public EnergyOp<T, 2, StorageIndex> {
	public:
		static const int dim = 2;
		using TV = Vector<T, dim>;

		// inputs
		const Field<int>& m_boundary_points;
		const Field<Vector<int, 2>>& m_boundary_edges;
		const Field<T>& m_mass;
		const std::map<int, T>& m_snode_area;

		// parameters
		T dHat = 1e-3;
		T kappa = 1e4;

		// intermediate variables;
		Field<Vector<int, 2>> PP;
		Field<Vector<int, 3>> PE;

		// friction related
		void initialize_friction(T mu_input, T epsv_input, T dt);
		virtual void update_weight_and_xhat(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, const Field<Vector<T, dim>>& an, const T dt, const T tsParam[]);
		T mu = 0, epsvh = 0, dt = 0;
		bool update_basis = true;
		Field<Vector<T, 2>> x_hat;
		Field<Vector<int, 2>> PP_friction;
		Field<Vector<int, 3>> PE_friction;
		Field<T> PP_normalForce;
		Field<T> PE_normalForce;
		Field<TV> PP_tanBasis;
		Field<TV> PE_tanBasis;
		Field<T> PE_yita;
		T x_weight = 1.0; // 2.0 for Newmark

		IpcEnergyOp2D(const Field<int>& boundary_points, const Field<Vector<int, 2>>& boundary_edges, const Field<T>& mass, const std::map<int, T>& snode_area, T energy_scale = 1.0);
		/* find constraint set*/
		void precompute(const Field<Vector<T, dim>>& x) override;
		/* adaptive kappa */
		void callback(const Field<Vector<T, dim>>& x) override;
		T stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& dx) override;
		T energy(const Field<Vector<T, dim>>& x) override;
		void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
		template <bool project_pd = true>
		void hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess) const;
		void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, const bool project_pd = true) override
		{
			if (project_pd)
				hessian_impl<true>(x, hess);
			else
				hessian_impl<false>(x, hess);
		}
		void internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force) override;
	};

	template <class T, class StorageIndex = int>
	class IpcEnergyOp3D : public EnergyOp<T, 3, StorageIndex> {
	public:
		static const int dim = 3;
		using TV = Vector<T, dim>;

		// inputs
		const Field<int>& m_boundary_points;
		const Field<Vector<int, 2>>& m_boundary_edges;
		const Field<Vector<int, 3>>& m_boundary_faces;
		const Field<Vector<T, dim>>& m_X;
		const Field<T>& m_mass;
		const std::map<int, T>& m_snode_area;
		std::map<std::pair<int, int>, T> m_sedge_area;

		// parameters
		T dHat = 1e-3;
		T kappa = 1e4;

		// intermediate variables;
		Field<Vector<T, 3>> x_hat;
		Field<Vector<int, 2>> PP; // PP from PT only
		Field<Vector<int, 3>> PE; // PE from PT only
		Field<Vector<int, 4>> PT;
		Field<Vector<int, 4>> PPM; // PP from EE only， the 0st and 2nd digits are actual pair
		Field<Vector<int, 4>> PEM; // PE from EE only, the 0st, 2nd, 3rd digits are actual pair
		Field<Vector<int, 4>> EEM;

		// friction related
		void initialize_friction(T mu_input, T epsv_input, T dt);
		virtual void update_weight_and_xhat(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, const Field<Vector<T, dim>>& an, const T dt, const T tsParam[]);

		void friction_precompute(const Field<Vector<T, dim>>& x);
		T mu = 0, epsvh = 0, epsv = 0, dt = 0;
		bool update_basis = true;
		Field<Vector<T, 3>> xhat;
		T x_weight = 1.0;
		Field<Vector<int, 4>> PP_friction;
		Field<Vector<int, 5>> PE_friction;
		Field<Vector<int, 4>> PT_friction;
		Field<T> PP_normalForce;
		Field<T> PE_normalForce;
		Field<T> PT_normalForce;
		Field<Matrix<T, dim, dim - 1>> PP_tanBasis;
		Field<Matrix<T, dim, dim - 1>> PE_tanBasis;
		Field<Matrix<T, dim, dim - 1>> PT_tanBasis;
		Field<T> PE_yita;
		Field<Vector<T, 2>> PT_yita;

		IpcEnergyOp3D(const Field<int>& boundary_points, const Field<Vector<int, 2>>& boundary_edges, const Field<Vector<int, 3>>& boundary_faces, const Field<Vector<T, dim>>& X, const Field<T>& mass, const std::map<int, T>& snode_area, T energy_scale = 1.0);
		/* find constraint set*/
		void precompute(const Field<Vector<T, dim>>& x) override;
		/* adaptive kappa */
		void callback(const Field<Vector<T, dim>>& x) override;
		T stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& dx) override;
		T energy(const Field<Vector<T, dim>>& x) override;
		void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad) override;
		template <bool project_pd = true>
		void hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess) const;
		void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, const bool project_pd = true) override
		{
			if (project_pd)
				hessian_impl<true>(x, hess);
			else
				hessian_impl<false>(x, hess);
		}

	protected:
		void point_triangle_constraints(const Field<Vector<T, 3>>& x);
		void edge_edge_constraints(const Field<Vector<T, 3>>& x);
	};

	template <class T, int dim, class StorageIndex = int>
	using IpcEnergyOp = typename std::conditional<dim == 2, IpcEnergyOp2D<T, StorageIndex>, IpcEnergyOp3D<T, StorageIndex>>::type;
} // namespace Bow::FEM::IPC

#include "IpcEnergy.tpp"

#endif
