#pragma once

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow::FEM {

	// C1 clamping
	template <class T>
	inline T f0_SF(T x2, T epsvh)
	{
		if (x2 >= epsvh * epsvh) {
			return std::sqrt(x2);
		}
		else {
			return x2 * (-std::sqrt(x2) / 3.0 + epsvh) / (epsvh * epsvh) + epsvh / 3.0;
		}
	}

	template <class T>
	inline T f1_SF_Div_RelDXNorm(T x2, T epsvh)
	{
		if (x2 >= epsvh * epsvh) {
			return 1 / std::sqrt(x2);
		}
		else {
			return (-std::sqrt(x2) + 2.0 * epsvh) / (epsvh * epsvh);
		}
	}

	template <class T>
	inline T f2_SF_Term(T x2, T epsvh)
	{
		return -1 / (epsvh * epsvh);
		// same for x2 >= epsvh * epsvh for C1 clamped friction
	}

	template <class T>
	void compute_PE_yita_3D(const Vector<T, 3>& v0, const Vector<T, 3>& v1, const Vector<T, 3>& v2, T& yita)
	{
		Eigen::Matrix<T, 3, 1> e12 = v2 - v1;
		yita = (v0 - v1).dot(e12) / e12.squaredNorm();
	}

	template <class T>
	void compute_PT_yita_3D(const Vector<T, 3>& v0, const Vector<T, 3>& v1, const Vector<T, 3>& v2, const Vector<T, 3>& v3, Vector<T, 2>& yita)
	{
		Eigen::Matrix<T, 2, 3> basis;
		basis.row(0) = (v2 - v1).transpose();
		basis.row(1) = (v3 - v1).transpose();
		yita = (basis * basis.transpose()).ldlt().solve(basis * (v0 - v1));
	}


	template <class T>
	void compute_EE_yita_3D(const Vector<T, 3>& v0, const Vector<T, 3>& v1, const Vector<T, 3>& v2, const Vector<T, 3>& v3, Vector<T, 2>& yita)
	{
		Eigen::Matrix<T, 1, 3> e20 = (v0 - v2).transpose();
		Eigen::Matrix<T, 1, 3> e01 = (v1 - v0).transpose();
		Eigen::Matrix<T, 1, 3> e23 = (v3 - v2).transpose();

		Eigen::Matrix<T, 2, 2> coefMtr;
		coefMtr(0, 0) = e01.squaredNorm();
		coefMtr(0, 1) = coefMtr(1, 0) = -e23.dot(e01);
		coefMtr(1, 1) = e23.squaredNorm();

		Eigen::Matrix<T, 2, 1> rhs;
		rhs[0] = -e20.dot(e01);
		rhs[1] = e20.dot(e23);

		yita = coefMtr.ldlt().solve(rhs);
	}


	// Point - Point
	template <class T, int dim>
	inline void Point_Point_Tangent_Basis(
		const Vector<T, dim>& p0,
		const Vector<T, dim>& p1,
		Matrix<T, dim, dim - 1>& basis)
	{
		if constexpr (dim == 2) {
			auto normal = (p1 - p0).normalized();
			basis(0) = normal(1);
			basis(1) = -normal(0);
		}
		else if constexpr (dim == 3) {
			auto v01 = (p1 - p0);
			auto xCross = Vector<T, dim>::UnitX().cross(v01);
			auto yCross = Vector<T, dim>::UnitY().cross(v01);
			if (xCross.squaredNorm() > yCross.squaredNorm()) {
				basis.col(0) = xCross.normalized();
				basis.col(1) = v01.cross(xCross).normalized();
			}
			else {
				basis.col(0) = yCross.normalized();
				basis.col(1) = v01.cross(yCross).normalized();
			}
		}
	}

	template <class T>
	inline void Point_Point_RelDX_2D(
		const Eigen::Matrix<T, 2, 1>& dx0,
		const Eigen::Matrix<T, 2, 1>& dx1,
		Eigen::Matrix<T, 2, 1>& relDX)
	{
		relDX = dx0 - dx1;
	}

	template <class T>
	inline void Point_Point_RelDX_3D(
		const Eigen::Matrix<T, 3, 1>& dx0,
		const Eigen::Matrix<T, 3, 1>& dx1,
		Eigen::Matrix<T, 3, 1>& relDX)
	{
		relDX = dx0 - dx1;
	}

	template <class T>
	inline void Point_Point_RelDXTan_To_Mesh_2D(
		const T& relDXTan,
		const Eigen::Matrix<T, 2, 1>& basis,
		Eigen::Matrix<T, 4, 1>& TTTDX)
	{
		TTTDX.template segment<2>(0) = basis * relDXTan;
		TTTDX.template segment<2>(2) = -TTTDX.template segment<2>(0);
	}

	template <class T>
	inline void Point_Point_RelDXTan_To_Mesh_3D(
		const Eigen::Matrix<T, 2, 1>& relDXTan,
		const Eigen::Matrix<T, 3, 2>& basis,
		Eigen::Matrix<T, 6, 1>& TTTDX)
	{
		TTTDX.template segment<3>(0) = basis * relDXTan;
		TTTDX.template segment<3>(3) = -TTTDX.template segment<3>(0);
	}

	template <class T>
	inline void Point_Point_TT_2D(
		const Eigen::Matrix<T, 2, 1>& basis,
		Eigen::Matrix<T, 1, 4>& TT)
	{
		TT.template block<1, 2>(0, 0) = basis.transpose();
		TT.template block<1, 2>(0, 2) = -basis.transpose();
	}

	template <class T>
	inline void Point_Point_TT_3D(
		const Eigen::Matrix<T, 3, 2>& basis,
		Eigen::Matrix<T, 2, 6>& TT)
	{
		TT.template block<2, 3>(0, 0) = basis.transpose();
		TT.template block<2, 3>(0, 3) = -basis.transpose();
	}

	//


	// Point - Edge
	template <class T, int dim>
	inline void Point_Edge_Tangent_Basis(
		const Vector<T, dim>& v0,
		const Vector<T, dim>& v1,
		const Vector<T, dim>& v2,
		Matrix<T, dim, dim - 1>& basis)
	{
		if constexpr (dim == 2) {
			basis = (v1 - v0).normalized();
		}
		else if constexpr (dim == 3) {
			Matrix<T, 3, 1> v12 = v2 - v1;
			basis.col(0) = v12.normalized();
			basis.col(1) = v12.cross(v0 - v1).normalized();
		}
	}

	template <class T, int dim>
	inline void Point_Edge_RelDX(
		const Eigen::Matrix<T, dim, 1>& dx0,
		const Eigen::Matrix<T, dim, 1>& dx1,
		const Eigen::Matrix<T, dim, 1>& dx2,
		T yita,
		Eigen::Matrix<T, dim, 1>& relDX)
	{
		relDX = dx0 - (dx1 + yita * (dx2 - dx1));
	}

	template <class T>
	inline void Point_Edge_RelDXTan_To_Mesh_2D(
		T relDXTan,
		const Eigen::Matrix<T, 2, 1>& basis,
		T yita,
		Eigen::Matrix<T, 6, 1>& TTTDX)
	{
		TTTDX.template segment<2>(0) = basis * relDXTan;
		TTTDX.template segment<2>(2) = (yita - 1.0) * TTTDX.template segment<2>(0);
		TTTDX.template segment<2>(4) = -yita * TTTDX.template segment<2>(0);
	}

	template <class T>
	inline void Point_Edge_RelDXTan_To_Mesh_3D(
		const Eigen::Matrix<T, 2, 1>& relDXTan,
		const Eigen::Matrix<T, 3, 2>& basis,
		T yita,
		Eigen::Matrix<T, 9, 1>& TTTDX)
	{
		TTTDX.template segment<3>(0) = basis * relDXTan;
		TTTDX.template segment<3>(3) = (yita - 1.0) * TTTDX.template segment<3>(0);
		TTTDX.template segment<3>(6) = -yita * TTTDX.template segment<3>(0);
	}

	template <class T>
	inline void Point_Edge_TT_2D(
		const Eigen::Matrix<T, 2, 1>& basis,
		T yita,
		Eigen::Matrix<T, 1, 6>& TT)
	{
		TT.template block<1, 2>(0, 0) = basis.transpose();
		TT.template block<1, 2>(0, 2) = (yita - 1.0) * basis.transpose();
		TT.template block<1, 2>(0, 4) = -yita * basis.transpose();
	}

	template <class T>
	inline void Point_Edge_TT_3D(
		const Eigen::Matrix<T, 3, 2>& basis,
		T yita,
		Eigen::Matrix<T, 2, 9>& TT)
	{
		TT.template block<2, 3>(0, 0) = basis.transpose();
		TT.template block<2, 3>(0, 3) = (yita - 1.0) * basis.transpose();
		TT.template block<2, 3>(0, 6) = -yita * basis.transpose();
	}

	// Point - Triangle
	template <class T>
	inline void Point_Triangle_Tangent_Basis(
		const Matrix<T, 3, 1>& v0,
		const Matrix<T, 3, 1>& v1,
		const Matrix<T, 3, 1>& v2,
		const Matrix<T, 3, 1>& v3,
		Matrix<T, 3, 2>& basis)
	{
		auto v12 = v2 - v1;
		basis.col(0) = v12.normalized();
		basis.col(1) = v12.cross(v3 - v1).cross(v12).normalized();
	}

	template <class T>
	inline void Point_Triangle_RelDX(
		const Eigen::Matrix<T, 3, 1>& dx0,
		const Eigen::Matrix<T, 3, 1>& dx1,
		const Eigen::Matrix<T, 3, 1>& dx2,
		const Eigen::Matrix<T, 3, 1>& dx3,
		T beta1, T beta2,
		Eigen::Matrix<T, 3, 1>& relDX)
	{
		relDX = dx0 - (dx1 + beta1 * (dx2 - dx1) + beta2 * (dx3 - dx1));
	}

	template <class T>
	inline void Point_Triangle_RelDXTan_To_Mesh(
		const Eigen::Matrix<T, 2, 1>& relDXTan,
		const Eigen::Matrix<T, 3, 2>& basis,
		T beta1, T beta2,
		Eigen::Matrix<T, 12, 1>& TTTDX)
	{
		TTTDX.template segment<3>(0) = basis * relDXTan;
		TTTDX.template segment<3>(3) = (-1 + beta1 + beta2) * TTTDX.template segment<3>(0);
		TTTDX.template segment<3>(6) = -beta1 * TTTDX.template segment<3>(0);
		TTTDX.template segment<3>(9) = -beta2 * TTTDX.template segment<3>(0);
	}

	template <class T>
	inline void Point_Triangle_TT(
		const Eigen::Matrix<T, 3, 2>& basis,
		T beta1, T beta2,
		Eigen::Matrix<T, 2, 12>& TT)
	{
		TT.template block<2, 3>(0, 0) = basis.transpose();
		TT.template block<2, 3>(0, 3) = (-1 + beta1 + beta2) * basis.transpose();
		TT.template block<2, 3>(0, 6) = -beta1 * basis.transpose();
		TT.template block<2, 3>(0, 9) = -beta2 * basis.transpose();
	}

	 //Edge - Edge
	template <class T>
	inline void Edge_Edge_Tangent_Basis(
		const Matrix<T, 3, 1>& v0,
		const Matrix<T, 3, 1>& v1,
		const Matrix<T, 3, 1>& v2,
		const Matrix<T, 3, 1>& v3,
		Matrix<T, 3, 2>& basis)
	{
		Eigen::Matrix<T, 3, 1> v01 = v1 - v0;
		basis.col(0) = v01.normalized();
		basis.col(1) = v01.cross(v3 - v2).cross(v01).normalized();
	}

	template <class T>
	inline void Edge_Edge_RelDX(
		const Eigen::Matrix<T, 3, 1>& dx0,
		const Eigen::Matrix<T, 3, 1>& dx1,
		const Eigen::Matrix<T, 3, 1>& dx2,
		const Eigen::Matrix<T, 3, 1>& dx3,
		T gamma1, T gamma2,
		Eigen::Matrix<T, 3, 1>& relDX)
	{
		relDX = dx0 + gamma1 * (dx1 - dx0) - (dx2 + gamma2 * (dx3 - dx2));
	}

	template <class T>
	inline void Edge_Edge_RelDXTan_To_Mesh(
		const Eigen::Matrix<T, 2, 1>& relDXTan,
		const Eigen::Matrix<T, 3, 2>& basis,
		T gamma1, T gamma2,
		Eigen::Matrix<T, 12, 1>& TTTDX)
	{
		Eigen::Matrix<T, 3, 1> relDXTan3D = basis * relDXTan;
		TTTDX.template segment<3>(0) = (1.0 - gamma1) * relDXTan3D;
		TTTDX.template segment<3>(3) = gamma1 * relDXTan3D;
		TTTDX.template segment<3>(6) = (gamma2 - 1.0) * relDXTan3D;
		TTTDX.template segment<3>(9) = -gamma2 * relDXTan3D;
	}

	template <class T>
	inline void Edge_Edge_TT(
		const Eigen::Matrix<T, 3, 2>& basis,
		T gamma1, T gamma2,
		Eigen::Matrix<T, 2, 12>& TT)
	{
		TT.template block<2, 3>(0, 0) = (1.0 - gamma1) * basis.transpose();
		TT.template block<2, 3>(0, 3) = gamma1 * basis.transpose();
		TT.template block<2, 3>(0, 6) = (gamma2 - 1.0) * basis.transpose();
		TT.template block<2, 3>(0, 9) = -gamma2 * basis.transpose();
	}


} // namespace Bow::FEM
