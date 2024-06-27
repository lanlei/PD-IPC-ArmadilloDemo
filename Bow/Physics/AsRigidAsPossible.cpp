#include "AsRigidAsPossible.h"
#include <Bow/Math/SVD.h>
#include <Bow/Math/PolarDecomposition.h>
#include <Bow/Physics/ConstitutiveModel.h>
#include <Bow/Math/Utils.h>
#include <Bow/Math/MathTools.h>

namespace Bow {
	namespace ConstitutiveModel {
		namespace AsRigidAsPossible {
			namespace internal {
				template <class DerivedV, class DerivedW, class Scalar>
				inline void dpsi_dsigma(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedW>& de_dsigma)
				{
					EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedV);
					using T = Scalar;
					static const int dim = DerivedV::SizeAtCompileTime;
					de_dsigma = T(2) * mu * (sigma - Bow::Vector<T, dim>::Ones());
				}

				template <class DerivedV, class DerivedM, class Scalar>
				inline void d2psi_dsigma2(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedM>& d2e_dsigma2)
				{
					EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedV);
					using T = Scalar;
					static const int dim = DerivedV::SizeAtCompileTime;

					T _2mu = mu * T(2);
					T zero = T(0);
					if constexpr (dim == 2) {
						d2e_dsigma2(0, 0) = _2mu;
						d2e_dsigma2(1, 1) = _2mu;
						d2e_dsigma2(0, 1) = d2e_dsigma2(1, 0) = zero;
					}
					else {
						Bow::Vector<T, dim> sigma_prod_noI;
						sigma_prod_noI[0] = sigma[1] * sigma[2];
						sigma_prod_noI[1] = sigma[0] * sigma[2];
						sigma_prod_noI[2] = sigma[0] * sigma[1];
						d2e_dsigma2(0, 0) = _2mu;
						d2e_dsigma2(1, 1) = _2mu;
						d2e_dsigma2(2, 2) = _2mu;
						d2e_dsigma2(0, 1) = d2e_dsigma2(1, 0) = zero;
						d2e_dsigma2(0, 2) = d2e_dsigma2(2, 0) = zero;
						d2e_dsigma2(1, 2) = d2e_dsigma2(2, 1) = zero;
					}
				}

				template <class DerivedV, class DerivedW, class Scalar>
				inline void B_left_coeff(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedW>& left_coeff)
				{
					// https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
					// Eq 77
					using T = typename DerivedV::Scalar;
					const int dim = DerivedV::SizeAtCompileTime;

					if constexpr (dim == 2)
						left_coeff[0] = mu;
					else {
						left_coeff[0] = mu;
						left_coeff[1] = mu;
						left_coeff[2] = mu;
					}
				}
			} // namespace internal

			template <class DerivedM, class Scalar>
			BOW_INLINE Scalar psi(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam)
			{
				using T = Scalar;
				static const int dim = DerivedM::RowsAtCompileTime;

				Bow::Matrix<T, dim, dim> U, V;
				Bow::Vector<T, dim> sigma;
				Math::svd(F, U, sigma, V);
				return mu * (sigma - Bow::Vector<T, dim>::Ones()).squaredNorm();
			}

			template <class DerivedM, class DerivedN, class Scalar>
			BOW_INLINE void first_piola(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedN>& P)
			{
				using T = Scalar;
				static const int dim = DerivedM::RowsAtCompileTime;

				Bow::Matrix<T, dim, dim> R, S;
				Math::polar_decomposition(F, R, S);
				P = T(2) * mu * (F - R);
			}

			template <bool projectPD, class DerivedM, class DerivedH, class Scalar>
			BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedM>& F, Scalar mu, Scalar lam, Eigen::MatrixBase<DerivedH>& dPdF)
			{
				using T = Scalar;
				static const int dim = DerivedM::RowsAtCompileTime;

				dPdF.setZero();
				Bow::Matrix<T, dim, dim> U, V;
				Bow::Vector<T, dim> sigma;
				Math::svd(F, U, sigma, V);
				Bow::Vector<T, dim> de_dsigma;
				Bow::Matrix<T, dim, dim> d2e_dsigma2;
				internal::dpsi_dsigma(sigma, mu, lam, de_dsigma);
				internal::d2psi_dsigma2(sigma, mu, lam, d2e_dsigma2);
				Bow::Vector<T, 2 * dim - 3> left_coeff;
				internal::B_left_coeff(sigma, mu, lam, left_coeff);
				Bow::ConstitutiveModel::first_piola_derivative<projectPD>(U, sigma, V, de_dsigma, left_coeff, d2e_dsigma2, dPdF);
			}

			// TODO: not redo SVD every time
			template <bool projectPD, class DerivedM, class Scalar>
			BOW_INLINE void first_piola_differential(const Eigen::MatrixBase<DerivedM>& F, const Eigen::MatrixBase<DerivedM>& dF, Scalar mu, Scalar lam, Eigen::MatrixBase<DerivedM>& dP)
			{
				using T = Scalar;
				static const int dim = DerivedM::RowsAtCompileTime;

				dP.setZero();
				Bow::Matrix<T, dim, dim> U, V;
				Bow::Vector<T, dim> sigma;
				Math::svd(F, U, sigma, V);
				/*	DerivedM JFinvT;
					Math::cofactor(F, JFinvT);
					T J = F.determinant();*/
				if constexpr (!projectPD) {
					DerivedM R = U * V.transpose();
					DerivedM S = V * sigma.asDiagonal() * V.transpose();

					dP.noalias() = DerivedM::Zero();
					dP += 2 * mu * dF;
					MATH_TOOLS::addScaledRotationalDifferential(R, S, dF, -2 * mu, dP);
					MATH_TOOLS::addScaledCofactorMatrixDifferential(F, dF, T(0), dP);
				}
				else {
					Bow::Matrix<T, dim, dim> D = U.transpose() * dF * V;
					Bow::Matrix<T, dim, dim> K;
					auto clamp_small_magnitude = [&](T x, T eps) {
						if (x < -eps)
							return x;
						else if (x < 0)
							return -eps;
						else if (x < eps)
							return eps;
						else
							return x;
					};
					if constexpr (dim == 2) {
						Bow::Matrix<T, dim, dim> Aij;
						Bow::Matrix<T, 2, 2> B01;
						T _2mu = mu * 2;
						T Sprod[2] = { sigma(1), sigma(0) };
						T psi0 = _2mu * (sigma(0) - 1);
						T psi1 = _2mu * (sigma(1) - 1);
						T psi00 = _2mu;
						T psi11 = _2mu;
						T psi01 = T(0);
						T m01 = _2mu;
						T p01 = (psi0 + psi1) / clamp_small_magnitude(sigma(0) + sigma(1), 1e-6);

						Aij(0, 0) = psi00;
						Aij(0, 1) = Aij(1, 0) = psi01;
						Aij(1, 1) = psi11;
						B01(0, 0) = B01(1, 1) = (m01 + p01) * 0.5;
						B01(0, 1) = B01(1, 0) = (m01 - p01) * 0.5;

						Math::make_pd(Aij);
						Math::make_pd(B01);

						K(0, 0) = Aij(0, 0) * D(0, 0) + Aij(0, 1) * D(1, 1);
						K(1, 1) = Aij(1, 0) * D(0, 0) + Aij(1, 1) * D(1, 1);
						K(0, 1) = B01(0, 0) * D(0, 1) + B01(0, 1) * D(1, 0);
						K(1, 0) = B01(1, 0) * D(0, 1) + B01(1, 1) * D(1, 0);
					}
					else {
						Bow::Matrix<T, dim, dim> Aij;
						Bow::Matrix<T, 2, 2> B01, B12, B20;
						T _2mu = mu * 2;
						T Sprod[3] = { sigma(1) * sigma(2), sigma(0) * sigma(2), sigma(0) * sigma(1) };
						T psi0 = _2mu * (sigma(0) - 1);
						T psi1 = _2mu * (sigma(1) - 1);
						T psi2 = _2mu * (sigma(2) - 1);
						T psi00 = _2mu;
						T psi11 = _2mu;
						T psi22 = _2mu;
						T psi01 = T(0);
						T psi02 = T(0);
						T psi12 = T(0);
						T m01 = _2mu; // i = 0
						T m02 = _2mu; // i = 2
						T m12 = _2mu; // i = 1
						T p01 = (psi0 + psi1) / clamp_small_magnitude(sigma(0) + sigma(1), 1e-6);
						T p02 = (psi0 + psi2) / clamp_small_magnitude(sigma(0) + sigma(2), 1e-6);
						T p12 = (psi1 + psi2) / clamp_small_magnitude(sigma(1) + sigma(2), 1e-6);

						Aij(0, 0) = psi00;
						Aij(1, 1) = psi11;
						Aij(2, 2) = psi22;
						Aij(0, 1) = Aij(1, 0) = psi01;
						Aij(0, 2) = Aij(2, 0) = psi02;
						Aij(1, 2) = Aij(2, 1) = psi12;
						B01(0, 0) = B01(1, 1) = (m01 + p01) * 0.5;
						B01(0, 1) = B01(1, 0) = (m01 - p01) * 0.5;
						B12(0, 0) = B12(1, 1) = (m12 + p12) * 0.5;
						B12(0, 1) = B12(1, 0) = (m12 - p12) * 0.5;
						B20(0, 0) = B20(1, 1) = (m02 + p02) * 0.5;
						B20(0, 1) = B20(1, 0) = (m02 - p02) * 0.5;

						Math::make_pd(Aij);
						Math::make_pd(B01);
						Math::make_pd(B12);
						Math::make_pd(B20);

						K(0, 0) = Aij(0, 0) * D(0, 0) + Aij(0, 1) * D(1, 1) + Aij(0, 2) * D(2, 2);
						K(1, 1) = Aij(1, 0) * D(0, 0) + Aij(1, 1) * D(1, 1) + Aij(1, 2) * D(2, 2);
						K(2, 2) = Aij(2, 0) * D(0, 0) + Aij(2, 1) * D(1, 1) + Aij(2, 2) * D(2, 2);
						K(0, 1) = B01(0, 0) * D(0, 1) + B01(0, 1) * D(1, 0);
						K(1, 0) = B01(1, 0) * D(0, 1) + B01(1, 1) * D(1, 0);
						K(0, 2) = B20(0, 0) * D(0, 2) + B20(0, 1) * D(2, 0);
						K(2, 0) = B20(1, 0) * D(0, 2) + B20(1, 1) * D(2, 0);
						K(1, 2) = B12(0, 0) * D(1, 2) + B12(0, 1) * D(2, 1);
						K(2, 1) = B12(1, 0) * D(1, 2) + B12(1, 1) * D(2, 1);
					}

					dP = U * K * V.transpose();
				}
			}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
			template float psi(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam);
			template void first_piola(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& P);
			template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
			template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
			template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const Eigen::MatrixBase<Matrix<float, 2, 2>>& dF, float mu, float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& dP);
			template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const Eigen::MatrixBase<Matrix<float, 2, 2>>& dF, float mu, float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& dP);
#endif
#ifdef BOW_COMPILE_3D
			template float psi(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam);
			template void first_piola(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& P);
			template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
			template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
			template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const Eigen::MatrixBase<Matrix<float, 3, 3>>& dF, float mu, float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& dP);
			template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const Eigen::MatrixBase<Matrix<float, 3, 3>>& dF, float mu, float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& dP);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
			template double psi(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam);
			template void first_piola(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& P);
			template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
			template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
			template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const Eigen::MatrixBase<Matrix<double, 2, 2>>& dF, double mu, double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& dP);
			template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const Eigen::MatrixBase<Matrix<double, 2, 2>>& dF, double mu, double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& dP);
#endif
#ifdef BOW_COMPILE_3D
			template double psi(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam);
			template void first_piola(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& P);
			template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
			template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
			template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const Eigen::MatrixBase<Matrix<double, 3, 3>>& dF, double mu, double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& dP);
			template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const Eigen::MatrixBase<Matrix<double, 3, 3>>& dF, double mu, double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& dP);
#endif
#endif
#endif

		}
	}
} // namespace Bow::ConstitutiveModel::FixedCorotated