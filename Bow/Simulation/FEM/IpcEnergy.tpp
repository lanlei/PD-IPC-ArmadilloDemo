#include "IpcEnergy.h"
#include <Bow/Geometry/IpcTookit/CCD.h>
#include <Bow/Geometry/IpcTookit/PointPointDistance.h>
#include <Bow/Geometry/IpcTookit/PointEdgeDistance.h>
#include <Bow/Geometry/IpcTookit/PointTriangleDistance.h>
#include <Bow/Geometry/IpcTookit/EdgeEdgeDistance.h>
#include <Bow/Geometry/IpcTookit/DistanceType.h>
#include <Bow/Math/Barrier.h>
#include <Bow/Math/Utils.h>
#include <Bow/Utils/Logging.h>
#include <oneapi/tbb.h>

namespace Bow::FEM::IPC {

	template <class T, class StorageIndex>
	IpcEnergyOp2D<T, StorageIndex>::IpcEnergyOp2D(const Field<int>& boundary_points, const Field<Vector<int, 2>>& boundary_edges, const Field<T>& mass, const std::map<int, T>& snode_area, T energy_scale)
		: m_boundary_points(boundary_points)
		, m_boundary_edges(boundary_edges)
		, m_mass(mass)
		, m_snode_area(snode_area)
	{
		this->energy_scale = energy_scale;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::precompute(const Field<Vector<T, dim>>& x)
	{
		PP.resize(0);
		PE.resize(0);
		T dHat2 = dHat * dHat;
		for (const auto p : m_boundary_points) {
			for (const auto& edge : m_boundary_edges) {
				if (p == edge[0] || p == edge[1])
					continue;
				const auto& e0 = x[edge[0]];
				const auto& e1 = x[edge[1]];
				if (!Geometry::IPC::point_edge_cd_broadphase(x[p], e0, e1, dHat))
					continue;
				switch (Geometry::IPC::point_edge_distance_type(x[p], e0, e1)) {
				case 0:
					if (Geometry::IPC::point_point_distance(x[p], e0) < dHat2) {
						Vector<int, 2> pair(p, edge[0]);
						PP.emplace_back(pair);
					}
					break;
				case 1:
					if (Geometry::IPC::point_point_distance(x[p], e1) < dHat2) {
						Vector<int, 2> pair(p, edge[1]);
						PP.emplace_back(pair);
					}
					break;
				case 2:
					if (Geometry::IPC::point_edge_distance(x[p], e0, e1) < dHat2) {
						PE.emplace_back(Vector<int, 3>(p, edge[0], edge[1]));
					}
					break;
				}
			}
		}
		//    Logging::info("# PP constraint: ", PP.size());
		//    Logging::info("# PE constraint: ", PE.size());
		if (mu > 0 && update_basis) {
			update_basis = false;
			PP_friction = PP;
			PE_friction = PE;
			PP_normalForce.clear();
			PE_normalForce.clear();
			PP_tanBasis.clear();
			PE_tanBasis.clear();
			PE_yita.clear();
			for (const auto& pp_pair : PP_friction) {
				TV p0 = x[pp_pair[0]], p1 = x[pp_pair[1]];
				T dist2 = Geometry::IPC::point_point_distance(p0, p1);
				T bGrad = Math::barrier_gradient(dist2, dHat2, kappa);
				const auto area = m_snode_area.find(pp_pair[0]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				PP_normalForce.push_back(-bGrad * 2 * area->second * dHat * std::sqrt(dist2));

				TV m = (p1 - p0).normalized();
				PP_tanBasis.push_back(TV(m(1), -m(0)));
			}
			for (const auto& pe_pair : PE_friction) {
				TV p = x[pe_pair[0]], e0 = x[pe_pair[1]], e1 = x[pe_pair[2]];
				T dist2 = Geometry::IPC::point_edge_distance(p, e0, e1);
				T bGrad = Math::barrier_gradient(dist2, dHat2, kappa);
				const auto area = m_snode_area.find(pe_pair[0]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				PE_normalForce.push_back(-bGrad * 2 * area->second * dHat * std::sqrt(dist2));

				TV m = (e1 - e0).normalized();
				PE_tanBasis.push_back(m);
				PE_yita.push_back(m.dot(p - e0) / (e1 - e0).norm());
			}
		}
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::initialize_friction(T mu_input, T epsv_input, T dt_input)
	{
		mu = mu_input;
		epsvh = epsv_input * dt_input;
		dt = dt_input;
		update_basis = true;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::internal_force(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force)
	{
		bool update_basis_bk = update_basis;
		update_basis = true;
		precompute(x);
		update_basis = update_basis_bk;
		T dHat2 = dHat * dHat;
		force.resize(x.size());
		std::fill(force.begin(), force.end(), Vector<T, dim>::Zero());
		for (const auto& pp_pair : PP) {
			Vector<T, 2 * 2> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[1]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PP_grad *= area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			force[pp_pair[0]] += PP_grad.template segment<2>(0);
			force[pp_pair[1]] += PP_grad.template segment<2>(2);
		}

		for (const auto& pe_pair : PE) {
			Vector<T, 2 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PE_grad *= area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			force[pe_pair[0]] += PE_grad.template segment<2>(0);
			force[pe_pair[1]] += PE_grad.template segment<2>(2);
			force[pe_pair[2]] += PE_grad.template segment<2>(4);
		}

		if (mu > 0) {
			Field<Vector<T, dim>> friction(x.size(), Vector<T, dim>::Zero());
			for (size_t i = 0; i < PP_friction.size(); ++i) {
				TV p0 = vn[PP_friction[i][0]] * dt, p1 = vn[PP_friction[i][1]] * dt;
				TV relDX2D;
				Point_Point_RelDX_2D(p0, p1, relDX2D);
				T relDX = PP_tanBasis[i].transpose() * relDX2D;
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				relDX *= f1_div_relDXNorm * mu * PP_normalForce[i];
				Eigen::Matrix<T, 4, 1> TTTDX;
				Point_Point_RelDXTan_To_Mesh_2D(relDX, PP_tanBasis[i], TTTDX);
				force[PP_friction[i][0]] += TTTDX.template segment<2>(0);
				force[PP_friction[i][1]] += TTTDX.template segment<2>(2);
			}

			for (size_t i = 0; i < PE_friction.size(); ++i) {
				TV p = vn[PE_friction[i][0]] * dt, e0 = vn[PE_friction[i][1]] * dt, e1 = vn[PE_friction[i][2]] * dt;
				TV relDX2D;
				Point_Edge_RelDX(p, e0, e1, PE_yita[i], relDX2D);
				T relDX = PE_tanBasis[i].transpose() * relDX2D;
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				relDX *= f1_div_relDXNorm * mu * PE_normalForce[i];
				Eigen::Matrix<T, 6, 1> TTTDX;
				Point_Edge_RelDXTan_To_Mesh_2D(relDX, PE_tanBasis[i], PE_yita[i], TTTDX);
				force[PE_friction[i][0]] += TTTDX.template segment<2>(0);
				force[PE_friction[i][1]] += TTTDX.template segment<2>(2);
				force[PE_friction[i][2]] += TTTDX.template segment<2>(4);
			}
		}
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::callback(const Field<Vector<T, dim>>& x)
	{
		//NOTE: kappa lower bound (adaptive kappa) is not necessary for the new constitutive IPC,
		// instead kappa can simply be set like a Young's modulus in elasticity (~10x that of density should be great)
		// T dHat2 = dHat * dHat;
		// T H_b = Math::barrier_hessian(1.0e-16, dHat2, 1.0);
		// T total_mass = Eigen::Map<const Vector<T, Eigen::Dynamic>>(m_mass.data(), m_mass.size()).sum();
		// // kappa = 1.0e16 * total_mass / T(m_mass.size()) / (4.0e-16 * H_b) / 0.516693 * 400;

		//TODO: implement kappa adjustment when distances become too tiny for numerical fail-safe
	}

	template <class T, class StorageIndex>
	T IpcEnergyOp2D<T, StorageIndex>::stepsize_upperbound(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& dx)
	{
		T alpha = 1.0;
		for (const int p : m_boundary_points) {
			for (const auto& edge : m_boundary_edges) {
				if (p == edge[0] || p == edge[1])
					continue;
				if (Geometry::IPC::point_edge_ccd_broadphase(xn[p], xn[edge[0]], xn[edge[1]], dx[p], dx[edge[0]], dx[edge[1]], dHat))
					alpha = std::min(alpha, Geometry::IPC::point_edge_ccd(xn[p], xn[edge[0]], xn[edge[1]], dx[p], dx[edge[0]], dx[edge[1]], 0.1));
			}
		}
		return alpha;
	}

	template <class T, class StorageIndex>
	T IpcEnergyOp2D<T, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
	{
		T energy = 0;
		T dHat2 = dHat * dHat;
		for (const auto& pp_pair : PP) {
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			energy += this->energy_scale * area->second * dHat * Math::barrier(dist2, dHat2, kappa);
		}
		for (const auto& pe_pair : PE) {
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			energy += this->energy_scale * area->second * dHat * Math::barrier(dist2, dHat2, kappa);
		}

		if (mu > 0) {
			for (size_t i = 0; i < PP_friction.size(); ++i) {
				TV p0 = x_weight * x[PP_friction[i][0]] - x_hat[PP_friction[i][0]], p1 = x_weight * x[PP_friction[i][1]] - x_hat[PP_friction[i][1]];
				TV relDX2D;
				Point_Point_RelDX_2D(p0, p1, relDX2D);
				T relDX = PP_tanBasis[i].transpose() * relDX2D;
				energy += this->energy_scale * f0_SF(relDX * relDX, epsvh) * mu * PP_normalForce[i] / x_weight;
			}
			for (size_t i = 0; i < PE_friction.size(); ++i) {
				TV p = x_weight * x[PE_friction[i][0]] - x_hat[PE_friction[i][0]], e0 = x_weight * x[PE_friction[i][1]] - x_hat[PE_friction[i][1]], e1 = x_weight * x[PE_friction[i][2]] - x_hat[PE_friction[i][2]];
				TV relDX2D;
				Point_Edge_RelDX(p, e0, e1, PE_yita[i], relDX2D);
				T relDX = PE_tanBasis[i].transpose() * relDX2D;
				energy += this->energy_scale * f0_SF(relDX * relDX, epsvh) * mu * PE_normalForce[i] / x_weight;
			}
		}
		return energy;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
	{
		T dHat2 = dHat * dHat;
		grad.resize(x.size());
		std::fill(grad.begin(), grad.end(), Vector<T, dim>::Zero());
		for (const auto& pp_pair : PP) {
			Vector<T, 2 * 2> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[1]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PP_grad *= this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			grad[pp_pair[0]] += PP_grad.template segment<2>(0);
			grad[pp_pair[1]] += PP_grad.template segment<2>(2);
		}
		for (const auto& pe_pair : PE) {
			Vector<T, 2 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PE_grad *= this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			grad[pe_pair[0]] += PE_grad.template segment<2>(0);
			grad[pe_pair[1]] += PE_grad.template segment<2>(2);
			grad[pe_pair[2]] += PE_grad.template segment<2>(4);
		}
		if (mu > 0) {
			for (size_t i = 0; i < PP_friction.size(); ++i) {
				TV p0 = x_weight * x[PP_friction[i][0]] - x_hat[PP_friction[i][0]], p1 = x_weight * x[PP_friction[i][1]] - x_hat[PP_friction[i][1]];
				TV relDX2D;
				Point_Point_RelDX_2D(p0, p1, relDX2D);
				T relDX = PP_tanBasis[i].transpose() * relDX2D;
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				relDX *= f1_div_relDXNorm * mu * PP_normalForce[i];
				Eigen::Matrix<T, 4, 1> TTTDX;
				Point_Point_RelDXTan_To_Mesh_2D(relDX, PP_tanBasis[i], TTTDX);
				TTTDX *= this->energy_scale;
				grad[PP_friction[i][0]] += TTTDX.template segment<2>(0);
				grad[PP_friction[i][1]] += TTTDX.template segment<2>(2);
			}
			for (size_t i = 0; i < PE_friction.size(); ++i) {
				TV p = x_weight * x[PE_friction[i][0]] - x_hat[PE_friction[i][0]], e0 = x_weight * x[PE_friction[i][1]] - x_hat[PE_friction[i][1]], e1 = x_weight * x[PE_friction[i][2]] - x_hat[PE_friction[i][2]];
				TV relDX2D;
				Point_Edge_RelDX(p, e0, e1, PE_yita[i], relDX2D);
				T relDX = PE_tanBasis[i].transpose() * relDX2D;
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				relDX *= f1_div_relDXNorm * mu * PE_normalForce[i];
				Eigen::Matrix<T, 6, 1> TTTDX;
				Point_Edge_RelDXTan_To_Mesh_2D(relDX, PE_tanBasis[i], PE_yita[i], TTTDX);
				TTTDX *= this->energy_scale;
				grad[PE_friction[i][0]] += TTTDX.template segment<2>(0);
				grad[PE_friction[i][1]] += TTTDX.template segment<2>(2);
				grad[PE_friction[i][2]] += TTTDX.template segment<2>(4);
			}
		}
	}

	template <class T, class StorageIndex>
	template <bool project_pd>
	void IpcEnergyOp2D<T, StorageIndex>::hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess) const
	{
		T dHat2 = dHat * dHat;
		using IJK = Eigen::Triplet<T, StorageIndex>;
		std::vector<IJK> coeffs;
		for (const auto& pp_pair : PP) {
			Matrix<T, 2 * 2, 2 * 2> PP_hess;
			Geometry::IPC::point_point_distance_hessian(x[pp_pair[0]], x[pp_pair[1]], PP_hess);
			Vector<T, 2 * 2> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[1]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PP_hess = this->energy_scale * area->second * dHat * Math::barrier_hessian(dist2, dHat2, kappa) * PP_grad * PP_grad.transpose() + this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa) * PP_hess;
			if constexpr (project_pd)
				Math::make_pd(PP_hess);
			int indMap[] = { 2 * pp_pair[0], 2 * pp_pair[0] + 1, 2 * pp_pair[1], 2 * pp_pair[1] + 1 };
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PP_hess(i, j)));
		}
		for (const auto& pe_pair : PE) {
			Matrix<T, 2 * 3, 2 * 3> PE_hess;
			Geometry::IPC::point_edge_distance_hessian(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_hess);
			Vector<T, 2 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PE_hess = this->energy_scale * area->second * dHat * Math::barrier_hessian(dist2, dHat2, kappa) * PE_grad * PE_grad.transpose() + this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa) * PE_hess;
			if constexpr (project_pd)
				Math::make_pd(PE_hess);
			int indMap[] = { 2 * pe_pair[0], 2 * pe_pair[0] + 1, 2 * pe_pair[1], 2 * pe_pair[1] + 1, 2 * pe_pair[2], 2 * pe_pair[2] + 1 };
			for (int i = 0; i < 6; ++i)
				for (int j = 0; j < 6; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PE_hess(i, j)));
		}

		if (mu > 0) {
			for (size_t i = 0; i < PP_friction.size(); ++i) {
				TV p0 = x_weight * x[PP_friction[i][0]] - x_hat[PP_friction[i][0]], p1 = x_weight * x[PP_friction[i][1]] - x_hat[PP_friction[i][1]];
				TV relDX2D;
				Point_Point_RelDX_2D(p0, p1, relDX2D);
				T relDX = PP_tanBasis[i].transpose() * relDX2D;

				Eigen::Matrix<T, 1, 4> TT;
				Point_Point_TT_2D(PP_tanBasis[i], TT);
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				T f2_term = f2_SF_Term(relDX * relDX, epsvh);
				Eigen::Matrix<T, 4, 4> HessianI = Eigen::Matrix<T, 4, 4>::Zero();
				if (relDX * relDX >= epsvh * epsvh) {
					// zero
				}
				else {
					if (relDX == 0) {
						// no SPD projection needed
						HessianI = ((mu * PP_normalForce[i] * f1_div_relDXNorm) * TT.transpose()) * TT;
					}
					else {
						// only need to project the inner 2x2 matrix to SPD
						T innerMtr = ((f2_term / std::sqrt(relDX * relDX)) * relDX) * relDX;
						innerMtr += f1_div_relDXNorm;
						innerMtr *= mu * PP_normalForce[i];
						if constexpr (project_pd)
							innerMtr = std::max(innerMtr, (T)0);
						// tensor product:
						HessianI = TT.transpose() * innerMtr * TT;
					}
				}
				HessianI *= x_weight * this->energy_scale;
				int cIVInd[2] = { PP_friction[i][0], PP_friction[i][1] };
				for (int i = 0; i < 2; ++i)
					for (int j = 0; j < 2; ++j)
						for (int idI = 0; idI < 2; ++idI)
							for (int jdI = 0; jdI < 2; ++jdI)
								coeffs.emplace_back(
									cIVInd[i] * 2 + idI,
									cIVInd[j] * 2 + jdI,
									HessianI(i * 2 + idI, j * 2 + jdI));
			}

			for (size_t i = 0; i < PE_friction.size(); ++i) {
				TV p = x_weight * x[PE_friction[i][0]] - x_hat[PE_friction[i][0]], e0 = x_weight * x[PE_friction[i][1]] - x_hat[PE_friction[i][1]], e1 = x_weight * x[PE_friction[i][2]] - x_hat[PE_friction[i][2]];
				TV relDX2D;
				Point_Edge_RelDX(p, e0, e1, PE_yita[i], relDX2D);
				T relDX = PE_tanBasis[i].transpose() * relDX2D;

				Eigen::Matrix<T, 1, 6> TT;
				Point_Edge_TT_2D(PE_tanBasis[i], PE_yita[i], TT);
				T f1_div_relDXNorm = f1_SF_Div_RelDXNorm(relDX * relDX, epsvh);
				T f2_term = f2_SF_Term(relDX * relDX, epsvh);
				Eigen::Matrix<T, 6, 6> HessianI = Eigen::Matrix<T, 6, 6>::Zero();
				if (relDX * relDX >= epsvh * epsvh) {
					// zero
				}
				else {
					if (relDX == 0) {
						// no SPD projection needed
						HessianI = ((mu * PE_normalForce[i] * f1_div_relDXNorm) * TT.transpose()) * TT;
					}
					else {
						// only need to project the inner 2x2 matrix to SPD
						T innerMtr = ((f2_term / std::sqrt(relDX * relDX)) * relDX) * relDX;
						innerMtr += f1_div_relDXNorm;
						innerMtr *= mu * PE_normalForce[i];
						if constexpr (project_pd)
							innerMtr = std::max(innerMtr, (T)0);
						// tensor product:
						HessianI = TT.transpose() * innerMtr * TT;
					}
				}
				HessianI *= x_weight * this->energy_scale;
				int cIVInd[3] = { PE_friction[i][0], PE_friction[i][1], PE_friction[i][2] };
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						for (int idI = 0; idI < 2; ++idI)
							for (int jdI = 0; jdI < 2; ++jdI)
								coeffs.emplace_back(
									cIVInd[i] * 2 + idI,
									cIVInd[j] * 2 + jdI,
									HessianI(i * 2 + idI, j * 2 + jdI));
			}
		}

		hess.resize(x.size() * dim, x.size() * dim);
		hess.setZero();
		hess.setFromTriplets(coeffs.begin(), coeffs.end());
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp2D<T, StorageIndex>::update_weight_and_xhat(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, const Field<Vector<T, dim>>& an, const T dt, const T tsParam[])
	{
		x_weight = tsParam[2] / (2 * tsParam[0] * tsParam[1]);
		x_hat.resize(xn.size());
		tbb::parallel_for(size_t(0), x_hat.size(), [&](size_t i) {
			TV x_tilde = xn[i] + vn[i] * dt + tsParam[0] * (1 - 2 * tsParam[1]) * an[i] * dt * dt;
			x_hat[i] = x_weight * x_tilde - vn[i] * dt - (1 - tsParam[2]) * an[i] * dt * dt;
		});
	}

	/** ############### IPC 3D ################# */

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::initialize_friction(T mu_input, T epsv_input, T dt_input)
	{
		mu = mu_input;
		epsvh = epsv_input * dt_input;
		epsv = epsv_input;
		update_basis = true;
		dt = dt_input;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::friction_precompute(const Field<Vector<T, dim>>& x)
	{
		// update_basis = false;
		// PP_friction = PP;
		// PE_friction = PE;
		// PT_friction = PT;
		// PP_normalForce.clear();
		// PE_normalForce.clear();
		// PT_normalForce.clear();
		// PP_tanBasis.clear();
		// PE_tanBasis.clear();
		// PT_tanBasis.clear();
		// PE_yita.clear();
		// T dHat2 = dHat * dHat;
		// for (const auto& pp_pair : PP_friction) {
		//     const auto& p0 = x[pp_pair[0]], p1 = x[pp_pair[1]];
		//     T dist2 = Geometry::IPC::point_point_distance(p0, p1);
		//     T bGrad = Math::barrier_gradient(dist2, dHat2, kappa);
		//     const auto area = m_snode_area.find(pp_pair[0]);
		//     BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
		//     PP_normalForce.push_back(-bGrad * 2 * area->second * dHat * std::sqrt(dist2));
		//     Matrix<T, dim, dim - 1> basis;
		//     Point_Point_Tangent_Basis(p0, p1, basis);
		//     PP_tanBasis.push_back(basis);
		// }
		// for (const auto& pe_pair : PE_friction) {
		//     const auto& p = x[pe_pair[0]], e0 = x[pe_pair[1]], e1 = x[pe_pair[2]];
		//     T dist2 = Geometry::IPC::point_edge_distance(p, e0, e1);
		//     T bGrad = Math::barrier_gradient(dist2, dHat2, kappa);
		//     const auto area = m_snode_area.find(pe_pair[0]);
		//     BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
		//     PE_normalForce.push_back(-bGrad * 2 * area->second * dHat * std::sqrt(dist2));
		//     Matrix<T, dim, dim - 1> basis;
		//     Point_Edge_Tangent_Basis(p, e0, e1, basis);
		//     PE_tanBasis.push_back(basis);
		//     PE_yita.push_back(basis.col(0).dot(p - e0) / (e1 - e0).norm());
		// }
	}

	template <class T, class StorageIndex>
	IpcEnergyOp3D<T, StorageIndex>::IpcEnergyOp3D(const Field<int>& boundary_points, const Field<Vector<int, 2>>& boundary_edges, const Field<Vector<int, 3>>& boundary_faces, const Field<Vector<T, dim>>& X, const Field<T>& mass, const std::map<int, T>& snode_area, T energy_scale)
		: m_boundary_points(boundary_points)
		, m_boundary_edges(boundary_edges)
		, m_boundary_faces(boundary_faces)
		, m_X(X)
		, m_mass(mass)
		, m_snode_area(snode_area)
	{
		this->energy_scale = energy_scale;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::precompute(const Field<Vector<T, dim>>& x)
	{
		PP.resize(0);
		PE.resize(0);
		PT.resize(0);
		PPM.resize(0);
		PEM.resize(0);
		EEM.resize(0);
		point_triangle_constraints(x);
		edge_edge_constraints(x);

		// Logging::info("# PP constraint: ", PP.size());
		// Logging::info("# PE constraint: ", PE.size());
		// Logging::info("# PT constraint: ", PT.size());
		// Logging::info("# PPM constraint: ", PPM.size());
		// Logging::info("# PEM constraint: ", PEM.size());
		// Logging::info("# EEM constraint: ", EEM.size());
		if (mu > 0 && update_basis) {
			// friction_precompute(x);
		}
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::callback(const Field<Vector<T, dim>>& x)
	{
		//NOTE: kappa lower bound (adaptive kappa) is not necessary for the new constitutive IPC,
		// instead kappa can simply be set like a Young's modulus in elasticity (~10x that of density should be great)
		// T dHat2 = dHat * dHat;
		// T H_b = Math::barrier_hessian(1.0e-16, dHat2, 1.0);
		// T total_mass = Eigen::Map<const Vector<T, Eigen::Dynamic>>(m_mass.data(), m_mass.size()).sum();
		// // kappa = 1.0e16 * total_mass / T(m_mass.size()) / (4.0e-16 * H_b) / 0.516693 * 400;

		//TODO: implement kappa adjustment when distances become too tiny for numerical fail-safe
	}

	template <class T, class StorageIndex>
	T IpcEnergyOp3D<T, StorageIndex>::stepsize_upperbound(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& dx)
	{
		T alpha = 1.0;
		for (const int p : m_boundary_points) {
			for (const auto& face : m_boundary_faces) {
				if (p == face[0] || p == face[1] || p == face[2])
					continue;
				if (Geometry::IPC::point_triangle_ccd_broadphase(xn[p], xn[face[0]], xn[face[1]], xn[face[2]], dx[p], dx[face[0]], dx[face[1]], dx[face[2]], dHat))
					alpha = std::min(alpha, Geometry::IPC::point_triangle_ccd(xn[p], xn[face[0]], xn[face[1]], xn[face[2]], dx[p], dx[face[0]], dx[face[1]], dx[face[2]], T(0.1), T(0.0)));
			}
		}
		for (size_t i = 0; i < m_boundary_edges.size() - 1; ++i) {
			auto edge0 = m_boundary_edges[i];
			for (size_t j = i + 1; j < m_boundary_edges.size(); ++j) {
				auto edge1 = m_boundary_edges[j];
				if (edge0[0] == edge1[0] || edge1[1] == edge0[1] || edge0[0] == edge1[1] || edge1[0] == edge0[1])
					continue;
				if (Geometry::IPC::edge_edge_ccd_broadphase(xn[edge0[0]], xn[edge0[1]], xn[edge1[0]], xn[edge1[1]], dx[edge0[0]], dx[edge0[1]], dx[edge1[0]], dx[edge1[1]], dHat))
					alpha = std::min(alpha, Geometry::IPC::edge_edge_ccd(xn[edge0[0]], xn[edge0[1]], xn[edge1[0]], xn[edge1[1]], dx[edge0[0]], dx[edge0[1]], dx[edge1[0]], dx[edge1[1]], T(0.1), T(0.0)));
			}
		}
		return alpha;
	}

	template <class T, class StorageIndex>
	T IpcEnergyOp3D<T, StorageIndex>::energy(const Field<Vector<T, dim>>& x)
	{
		T energy = 0;
		T dHat2 = dHat * dHat;

		// PP
		for (const auto& pp_pair : PP) {
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			energy += this->energy_scale * area->second * dHat * Math::barrier(dist2, dHat2, kappa);
		}
		for (const auto& pe_pair : PE) {
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			energy += this->energy_scale * area->second * dHat * Math::barrier(dist2, dHat2, kappa);
		}
		for (const auto& pt_pair : PT) {
			T dist2 = Geometry::IPC::point_triangle_distance(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]]);
			const auto area = m_snode_area.find(pt_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			energy += this->energy_scale * area->second * dHat * Math::barrier(dist2, dHat2, kappa);
		}

		// EE
		auto mollifier_info = [&](const Vector<int, 4>& pair_info) {
			const Vector<T, dim>& ea0_rest = m_X[pair_info[0]];
			const Vector<T, dim>& ea1_rest = m_X[pair_info[1]];
			const Vector<T, dim>& eb0_rest = m_X[pair_info[2]];
			const Vector<T, dim>& eb1_rest = m_X[pair_info[3]];
			T eps_x = Geometry::IPC::edge_edge_mollifier_threshold(ea0_rest, ea1_rest, eb0_rest, eb1_rest);
			return Geometry::IPC::edge_edge_mollifier(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x);
		};
		for (const auto& pp_pair : PPM) {
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[2]]);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pp_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			energy += this->energy_scale * ee_area * dHat * mollifier_info(pp_pair) * Math::barrier(dist2, dHat2, kappa);
		}
		for (const auto& pe_pair : PEM) {
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]]);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pe_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			energy += this->energy_scale * ee_area * dHat * mollifier_info(pe_pair) * Math::barrier(dist2, dHat2, kappa);
		}
		for (const auto& ee_pair : EEM) {
			T dist2 = Geometry::IPC::edge_edge_distance(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]]);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(ee_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			energy += this->energy_scale * ee_area * dHat * mollifier_info(ee_pair) * Math::barrier(dist2, dHat2, kappa);
		}
		return energy;
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad)
	{
		T dHat2 = dHat * dHat;
		grad.resize(x.size());
		std::fill(grad.begin(), grad.end(), Vector<T, dim>::Zero());

		// PP
		for (const auto& pp_pair : PP) {
			Vector<T, 2 * 3> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[1]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PP_grad *= this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			grad[pp_pair[0]] += PP_grad.template segment<3>(0);
			grad[pp_pair[1]] += PP_grad.template segment<3>(3);
		}
		for (const auto& pe_pair : PE) {
			Vector<T, 3 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PE_grad *= this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			grad[pe_pair[0]] += PE_grad.template segment<3>(0);
			grad[pe_pair[1]] += PE_grad.template segment<3>(3);
			grad[pe_pair[2]] += PE_grad.template segment<3>(6);
		}
		for (const auto& pt_pair : PT) {
			Vector<T, 4 * 3> PT_grad;
			Geometry::IPC::point_triangle_distance_gradient(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]], PT_grad);
			T dist2 = Geometry::IPC::point_triangle_distance(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]]);
			const auto area = m_snode_area.find(pt_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PT_grad *= this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa);
			grad[pt_pair[0]] += PT_grad.template segment<3>(0);
			grad[pt_pair[1]] += PT_grad.template segment<3>(3);
			grad[pt_pair[2]] += PT_grad.template segment<3>(6);
			grad[pt_pair[3]] += PT_grad.template segment<3>(9);
		}

		// EE
		auto mollifier_info = [&](const Vector<int, 4>& pair_info, T& m, Vector<T, dim * 4>& gm) {
			const Vector<T, dim>& ea0_rest = m_X[pair_info[0]];
			const Vector<T, dim>& ea1_rest = m_X[pair_info[1]];
			const Vector<T, dim>& eb0_rest = m_X[pair_info[2]];
			const Vector<T, dim>& eb1_rest = m_X[pair_info[3]];
			T eps_x = Geometry::IPC::edge_edge_mollifier_threshold(ea0_rest, ea1_rest, eb0_rest, eb1_rest);
			m = Geometry::IPC::edge_edge_mollifier(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x);
			Geometry::IPC::edge_edge_mollifier_gradient(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x, gm);
		};
		for (const auto& pp_pair : PPM) {
			Vector<T, 2 * 3> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[2]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[2]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			mollifier_info(pp_pair, mollifier, mollifier_grad);
			PP_grad *= Math::barrier_gradient(dist2, dHat2, kappa);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pp_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			T scale = this->energy_scale * ee_area * dHat;
			grad[pp_pair[0]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(0) + mollifier * PP_grad.template segment<3>(0));
			grad[pp_pair[1]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(3));
			grad[pp_pair[2]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(6) + mollifier * PP_grad.template segment<3>(3));
			grad[pp_pair[3]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(9));
		}
		for (const auto& pe_pair : PEM) {
			Vector<T, 3 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			mollifier_info(pe_pair, mollifier, mollifier_grad);
			PE_grad *= Math::barrier_gradient(dist2, dHat2, kappa);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pe_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			T scale = this->energy_scale * ee_area * dHat;
			grad[pe_pair[0]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(0) + mollifier * PE_grad.template segment<3>(0));
			grad[pe_pair[1]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(3));
			grad[pe_pair[2]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(6) + mollifier * PE_grad.template segment<3>(3));
			grad[pe_pair[3]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(9) + mollifier * PE_grad.template segment<3>(6));
		}
		for (const auto& ee_pair : EEM) {
			Vector<T, 4 * 3> EE_grad;
			Geometry::IPC::edge_edge_distance_gradient(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]], EE_grad);
			T dist2 = Geometry::IPC::edge_edge_distance(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			mollifier_info(ee_pair, mollifier, mollifier_grad);
			EE_grad *= Math::barrier_gradient(dist2, dHat2, kappa);
			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(ee_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			T scale = this->energy_scale * ee_area * dHat;
			grad[ee_pair[0]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(0) + mollifier * EE_grad.template segment<3>(0));
			grad[ee_pair[1]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(3) + mollifier * EE_grad.template segment<3>(3));
			grad[ee_pair[2]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(6) + mollifier * EE_grad.template segment<3>(6));
			grad[ee_pair[3]] += scale * (barrier_dist2 * mollifier_grad.template segment<3>(9) + mollifier * EE_grad.template segment<3>(9));
		}
	}

	template <class T, class StorageIndex>
	template <bool project_pd>
	void IpcEnergyOp3D<T, StorageIndex>::hessian_impl(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess) const
	{
		T dHat2 = dHat * dHat;
		using IJK = Eigen::Triplet<T, StorageIndex>;
		std::vector<IJK> coeffs;

		// PT
		for (const auto& pp_pair : PP) {
			Matrix<T, 2 * 3, 2 * 3> PP_hess;
			Geometry::IPC::point_point_distance_hessian(x[pp_pair[0]], x[pp_pair[1]], PP_hess);
			Vector<T, 2 * 3> PP_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[1]], PP_grad);
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[1]]);
			const auto area = m_snode_area.find(pp_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PP_hess = this->energy_scale * area->second * dHat * Math::barrier_hessian(dist2, dHat2, kappa) * PP_grad * PP_grad.transpose() + this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa) * PP_hess;
			if constexpr (project_pd)
				Math::make_pd(PP_hess);
			int indMap[] = { 3 * pp_pair[0], 3 * pp_pair[0] + 1, 3 * pp_pair[0] + 2,
				3 * pp_pair[1], 3 * pp_pair[1] + 1, 3 * pp_pair[1] + 2 };
			for (int i = 0; i < 2 * 3; ++i)
				for (int j = 0; j < 2 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PP_hess(i, j)));
		}
		for (const auto& pe_pair : PE) {
			Matrix<T, 3 * 3, 3 * 3> PE_hess;
			Geometry::IPC::point_edge_distance_hessian(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_hess);
			Vector<T, 3 * 3> PE_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]], PE_grad);
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[1]], x[pe_pair[2]]);
			const auto area = m_snode_area.find(pe_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PE_hess = this->energy_scale * area->second * dHat * Math::barrier_hessian(dist2, dHat2, kappa) * PE_grad * PE_grad.transpose() + this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa) * PE_hess;
			if constexpr (project_pd)
				Math::make_pd(PE_hess);
			int indMap[] = { 3 * pe_pair[0], 3 * pe_pair[0] + 1, 3 * pe_pair[0] + 2,
				3 * pe_pair[1], 3 * pe_pair[1] + 1, 3 * pe_pair[1] + 2,
				3 * pe_pair[2], 3 * pe_pair[2] + 1, 3 * pe_pair[2] + 2 };
			for (int i = 0; i < 3 * 3; ++i)
				for (int j = 0; j < 3 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PE_hess(i, j)));
		}
		for (const auto& pt_pair : PT) {
			Matrix<T, 4 * 3, 4 * 3> PT_hess;
			Geometry::IPC::point_triangle_distance_hessian(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]], PT_hess);
			Vector<T, 4 * 3> PT_grad;
			Geometry::IPC::point_triangle_distance_gradient(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]], PT_grad);
			T dist2 = Geometry::IPC::point_triangle_distance(x[pt_pair[0]], x[pt_pair[1]], x[pt_pair[2]], x[pt_pair[3]]);
			const auto area = m_snode_area.find(pt_pair[0]);
			BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
			PT_hess = this->energy_scale * area->second * dHat * Math::barrier_hessian(dist2, dHat2, kappa) * PT_grad * PT_grad.transpose() + this->energy_scale * area->second * dHat * Math::barrier_gradient(dist2, dHat2, kappa) * PT_hess;
			if constexpr (project_pd)
				Math::make_pd(PT_hess);
			int indMap[] = { 3 * pt_pair[0], 3 * pt_pair[0] + 1, 3 * pt_pair[0] + 2,
				3 * pt_pair[1], 3 * pt_pair[1] + 1, 3 * pt_pair[1] + 2,
				3 * pt_pair[2], 3 * pt_pair[2] + 1, 3 * pt_pair[2] + 2,
				3 * pt_pair[3], 3 * pt_pair[3] + 1, 3 * pt_pair[3] + 2 };
			for (int i = 0; i < 4 * 3; ++i)
				for (int j = 0; j < 4 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PT_hess(i, j)));
		}

		// EE
		auto mollifier_info = [&](const Vector<int, 4>& pair_info, T& m, Vector<T, dim * 4>& gm, Matrix<T, dim * 4, dim * 4>& hm) {
			const Vector<T, dim>& ea0_rest = m_X[pair_info[0]];
			const Vector<T, dim>& ea1_rest = m_X[pair_info[1]];
			const Vector<T, dim>& eb0_rest = m_X[pair_info[2]];
			const Vector<T, dim>& eb1_rest = m_X[pair_info[3]];
			T eps_x = Geometry::IPC::edge_edge_mollifier_threshold(ea0_rest, ea1_rest, eb0_rest, eb1_rest);
			m = Geometry::IPC::edge_edge_mollifier(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x);
			Geometry::IPC::edge_edge_mollifier_gradient(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x, gm);
			Geometry::IPC::edge_edge_mollifier_hessian(x[pair_info[0]], x[pair_info[1]], x[pair_info[2]], x[pair_info[3]], eps_x, hm);
		};
		for (const auto& pp_pair : PPM) {
			T dist2 = Geometry::IPC::point_point_distance(x[pp_pair[0]], x[pp_pair[2]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			Vector<T, 2 * 3> dist2_grad;
			Geometry::IPC::point_point_distance_gradient(x[pp_pair[0]], x[pp_pair[2]], dist2_grad);
			Vector<T, 2 * 3> barrier_dist2_grad;
			barrier_dist2_grad = dist2_grad * Math::barrier_gradient(dist2, dHat2, kappa);
			Vector<T, 4 * 3> barrier_dist2_grad_extended;
			barrier_dist2_grad_extended.setZero();
			barrier_dist2_grad_extended.template segment<dim>(0) = barrier_dist2_grad.template segment<dim>(0);
			barrier_dist2_grad_extended.template segment<dim>(6) = barrier_dist2_grad.template segment<dim>(3);
			Matrix<T, 2 * 3, 2 * 3> barrier_dist2_hess;
			Geometry::IPC::point_point_distance_hessian(x[pp_pair[0]], x[pp_pair[2]], barrier_dist2_hess);
			barrier_dist2_hess = Math::barrier_hessian(dist2, dHat2, kappa) * dist2_grad * dist2_grad.transpose() + Math::barrier_gradient(dist2, dHat2, kappa) * barrier_dist2_hess;

			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			Matrix<T, 4 * 3, 4 * 3> mollifier_hess;
			mollifier_info(pp_pair, mollifier, mollifier_grad, mollifier_hess);

			Matrix<T, 4 * 3, 4 * 3> PP_hess = barrier_dist2 * mollifier_hess
				+ mollifier_grad * barrier_dist2_grad_extended.transpose()
				+ barrier_dist2_grad_extended * mollifier_grad.transpose();
			PP_hess.template block<dim, dim>(0, 0) += mollifier * barrier_dist2_hess.template block<dim, dim>(0, 0);
			PP_hess.template block<dim, dim>(0, 6) += mollifier * barrier_dist2_hess.template block<dim, dim>(0, 3);
			PP_hess.template block<dim, dim>(6, 0) += mollifier * barrier_dist2_hess.template block<dim, dim>(3, 0);
			PP_hess.template block<dim, dim>(6, 6) += mollifier * barrier_dist2_hess.template block<dim, dim>(3, 3);

			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pp_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			PP_hess *= this->energy_scale * ee_area * dHat;

			if constexpr (project_pd)
				Math::make_pd(PP_hess);

			int indMap[] = { 3 * pp_pair[0], 3 * pp_pair[0] + 1, 3 * pp_pair[0] + 2,
				3 * pp_pair[1], 3 * pp_pair[1] + 1, 3 * pp_pair[1] + 2,
				3 * pp_pair[2], 3 * pp_pair[2] + 1, 3 * pp_pair[2] + 2,
				3 * pp_pair[3], 3 * pp_pair[3] + 1, 3 * pp_pair[3] + 2 };
			for (int i = 0; i < 4 * 3; ++i)
				for (int j = 0; j < 4 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PP_hess(i, j)));
		}
		for (const auto& pe_pair : PEM) {
			T dist2 = Geometry::IPC::point_edge_distance(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			Vector<T, 3 * 3> dist2_grad;
			Geometry::IPC::point_edge_distance_gradient(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]], dist2_grad);
			Vector<T, 3 * 3> barrier_dist2_grad;
			barrier_dist2_grad = dist2_grad * Math::barrier_gradient(dist2, dHat2, kappa);
			Vector<T, 4 * 3> barrier_dist2_grad_extended;
			barrier_dist2_grad_extended.setZero();
			barrier_dist2_grad_extended.template segment<dim>(0) = barrier_dist2_grad.template segment<dim>(0);
			barrier_dist2_grad_extended.template segment<2 * dim>(6) = barrier_dist2_grad.template segment<2 * dim>(3);
			Matrix<T, 3 * 3, 3 * 3> barrier_dist2_hess;
			Geometry::IPC::point_edge_distance_hessian(x[pe_pair[0]], x[pe_pair[2]], x[pe_pair[3]], barrier_dist2_hess);
			barrier_dist2_hess = Math::barrier_hessian(dist2, dHat2, kappa) * dist2_grad * dist2_grad.transpose() + Math::barrier_gradient(dist2, dHat2, kappa) * barrier_dist2_hess;

			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			Matrix<T, 4 * 3, 4 * 3> mollifier_hess;
			mollifier_info(pe_pair, mollifier, mollifier_grad, mollifier_hess);
			Matrix<T, 4 * 3, 4 * 3> PE_hess = barrier_dist2 * mollifier_hess
				+ mollifier_grad * barrier_dist2_grad_extended.transpose()
				+ barrier_dist2_grad_extended * mollifier_grad.transpose();
			PE_hess.template block<dim, dim>(0, 0) += mollifier * barrier_dist2_hess.template block<dim, dim>(0, 0);
			PE_hess.template block<dim, 2 * dim>(0, 6) += mollifier * barrier_dist2_hess.template block<dim, 2 * dim>(0, 3);
			PE_hess.template block<2 * dim, dim>(6, 0) += mollifier * barrier_dist2_hess.template block<2 * dim, dim>(3, 0);
			PE_hess.template block<2 * dim, 2 * dim>(6, 6) += mollifier * barrier_dist2_hess.template block<2 * dim, 2 * dim>(3, 3);

			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(pe_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}

			PE_hess *= this->energy_scale * ee_area * dHat;
			if constexpr (project_pd)
				Math::make_pd(PE_hess);

			int indMap[] = { 3 * pe_pair[0], 3 * pe_pair[0] + 1, 3 * pe_pair[0] + 2,
				3 * pe_pair[1], 3 * pe_pair[1] + 1, 3 * pe_pair[1] + 2,
				3 * pe_pair[2], 3 * pe_pair[2] + 1, 3 * pe_pair[2] + 2,
				3 * pe_pair[3], 3 * pe_pair[3] + 1, 3 * pe_pair[3] + 2 };
			for (int i = 0; i < 4 * 3; ++i)
				for (int j = 0; j < 4 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], PE_hess(i, j)));
		}
		for (const auto& ee_pair : EEM) {
			T dist2 = Geometry::IPC::edge_edge_distance(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]]);
			T barrier_dist2 = Math::barrier(dist2, dHat2, kappa);
			Vector<T, 4 * 3> dist2_grad;
			Geometry::IPC::edge_edge_distance_gradient(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]], dist2_grad);
			Vector<T, 4 * 3> barrier_dist2_grad;
			barrier_dist2_grad = dist2_grad * Math::barrier_gradient(dist2, dHat2, kappa);
			Matrix<T, 4 * 3, 4 * 3> barrier_dist2_hess;
			Geometry::IPC::edge_edge_distance_hessian(x[ee_pair[0]], x[ee_pair[1]], x[ee_pair[2]], x[ee_pair[3]], barrier_dist2_hess);
			barrier_dist2_hess = Math::barrier_hessian(dist2, dHat2, kappa) * dist2_grad * dist2_grad.transpose() + Math::barrier_gradient(dist2, dHat2, kappa) * barrier_dist2_hess;

			T mollifier;
			Vector<T, 4 * 3> mollifier_grad;
			Matrix<T, 4 * 3, 4 * 3> mollifier_hess;
			mollifier_info(ee_pair, mollifier, mollifier_grad, mollifier_hess);

			Matrix<T, 4 * 3, 4 * 3> EE_hess = barrier_dist2 * mollifier_hess + mollifier * barrier_dist2_hess
				+ mollifier_grad * barrier_dist2_grad.transpose()
				+ barrier_dist2_grad * mollifier_grad.transpose();

			T ee_area = 0;
			for (int ii = 0; ii < 4; ++ii) {
				const auto area = m_snode_area.find(ee_pair[ii]);
				BOW_ASSERT_INFO(area != m_snode_area.end(), "cannot find surface node area for IPC weighting");
				ee_area += area->second / 4.0;
			}
			EE_hess *= this->energy_scale * ee_area * dHat;

			if constexpr (project_pd)
				Math::make_pd(EE_hess);

			int indMap[] = { 3 * ee_pair[0], 3 * ee_pair[0] + 1, 3 * ee_pair[0] + 2,
				3 * ee_pair[1], 3 * ee_pair[1] + 1, 3 * ee_pair[1] + 2,
				3 * ee_pair[2], 3 * ee_pair[2] + 1, 3 * ee_pair[2] + 2,
				3 * ee_pair[3], 3 * ee_pair[3] + 1, 3 * ee_pair[3] + 2 };
			for (int i = 0; i < 4 * 3; ++i)
				for (int j = 0; j < 4 * 3; ++j)
					coeffs.push_back(IJK(indMap[i], indMap[j], EE_hess(i, j)));
		}

		hess.resize(x.size() * dim, x.size() * dim);
		hess.setZero();
		hess.setFromTriplets(coeffs.begin(), coeffs.end());
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::point_triangle_constraints(const Field<Vector<T, 3>>& x)
	{
		using namespace Geometry::IPC;
		T dHat2 = dHat * dHat;
		for (const auto p : m_boundary_points) {
			Vector<T, 3> p0 = x[p];
			for (const auto& face : m_boundary_faces) {
				if (p == face[0] || p == face[1] || p == face[2])
					continue;
				const Vector<T, 3>& t0 = x[face[0]];
				const Vector<T, 3>& t1 = x[face[1]];
				const Vector<T, 3>& t2 = x[face[2]];
				if (!Geometry::IPC::point_triangle_cd_broadphase(p0, t0, t1, t2, dHat))
					continue;
				switch (Geometry::IPC::point_triangle_distance_type(p0, t0, t1, t2)) {
				case 0:
					if (Geometry::IPC::point_point_distance(p0, t0) < dHat2)
						PP.push_back(Eigen::Vector2i(p, face[0]));
					break;
				case 1:
					if (Geometry::IPC::point_point_distance(p0, t1) < dHat2)
						PP.push_back(Eigen::Vector2i(p, face[1]));
					break;
				case 2:
					if (Geometry::IPC::point_point_distance(p0, t2) < dHat2)
						PP.push_back(Eigen::Vector2i(p, face[2]));
					break;
				case 3:
					if (Geometry::IPC::point_edge_distance(p0, t0, t1) < dHat2)
						PE.push_back(Eigen::Vector3i(p, face[0], face[1]));
					break;
				case 4:
					if (Geometry::IPC::point_edge_distance(p0, t1, t2) < dHat2)
						PE.push_back(Eigen::Vector3i(p, face[1], face[2]));
					break;
				case 5:
					if (Geometry::IPC::point_edge_distance(p0, t2, t0) < dHat2)
						PE.push_back(Eigen::Vector3i(p, face[2], face[0]));
					break;
				case 6:
					if (Geometry::IPC::point_triangle_distance(p0, t0, t1, t2) < dHat2)
						PT.push_back(Eigen::Vector4i(p, face[0], face[1], face[2]));
					break;
				default:
					break;
				}
			}
		}
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::edge_edge_constraints(const Field<Vector<T, 3>>& x)
	{
		using namespace Geometry::IPC;
		T dHat2 = dHat * dHat;
		for (size_t i = 0; i < m_boundary_edges.size() - 1; ++i) {
			auto edge0 = m_boundary_edges[i];
			const Vector<T, dim>& ea0 = x[edge0[0]];
			const Vector<T, dim>& ea1 = x[edge0[1]];
			for (size_t j = i + 1; j < m_boundary_edges.size(); ++j) {
				auto edge1 = m_boundary_edges[j];
				if (edge0[0] == edge1[0] || edge1[1] == edge0[1] || edge0[0] == edge1[1] || edge1[0] == edge0[1])
					continue;
				const Vector<T, dim>& eb0 = x[edge1[0]];
				const Vector<T, dim>& eb1 = x[edge1[1]];
				if (!edge_edge_cd_broadphase(ea0, ea1, eb0, eb1, dHat))
					continue;
				switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
				case 0:
					if (point_point_distance(ea0, eb0) < dHat2)
						PPM.push_back(Vector<int, 4>(edge0[0], edge0[1], edge1[0], edge1[1]));
					break;
				case 1:
					if (point_point_distance(ea0, eb1) < dHat2)
						PPM.push_back(Vector<int, 4>(edge0[0], edge0[1], edge1[1], edge1[0]));
					break;
				case 2:
					if (point_edge_distance(ea0, eb0, eb1) < dHat2)
						PEM.push_back(Vector<int, 4>(edge0[0], edge0[1], edge1[0], edge1[1]));
					break;
				case 3:
					if (point_point_distance(ea1, eb0) < dHat2)
						PPM.push_back(Vector<int, 4>(edge0[1], edge0[0], edge1[0], edge1[1]));
					break;
				case 4:
					if (point_point_distance(ea1, eb1) < dHat2)
						PPM.push_back(Vector<int, 4>(edge0[1], edge0[0], edge1[1], edge1[0]));
					break;
				case 5:
					if (point_edge_distance(ea1, eb0, eb1) < dHat2)
						PEM.push_back(Vector<int, 4>(edge0[1], edge0[0], edge1[0], edge1[1]));
					break;
				case 6:
					if (point_edge_distance(eb0, ea0, ea1) < dHat2)
						PEM.push_back(Vector<int, 4>(edge1[0], edge1[1], edge0[0], edge0[1]));
					break;
				case 7:
					if (point_edge_distance(eb1, ea0, ea1) < dHat2)
						PEM.push_back(Vector<int, 4>(edge1[1], edge1[0], edge0[0], edge0[1]));
					break;
				case 8:
					if (edge_edge_distance(ea0, ea1, eb0, eb1) < dHat2)
						EEM.push_back(Vector<int, 4>(edge0[0], edge0[1], edge1[0], edge1[1]));
					break;
				default:
					break;
				}
			}
		}
	}

	template <class T, class StorageIndex>
	void IpcEnergyOp3D<T, StorageIndex>::update_weight_and_xhat(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, const Field<Vector<T, dim>>& an, const T dt, const T tsParam[])
	{
		x_weight = tsParam[2] / (2 * tsParam[0] * tsParam[1]);
		x_hat.resize(xn.size());
		tbb::parallel_for(size_t(0), x_hat.size(), [&](size_t i) {
			TV x_tilde = xn[i] + vn[i] * dt + tsParam[0] * (1 - 2 * tsParam[1]) * an[i] * dt * dt;
			x_hat[i] = x_weight * x_tilde - vn[i] * dt - (1 - tsParam[2]) * an[i] * dt * dt;
		});
	}

} // namespace Bow::FEM::IPC