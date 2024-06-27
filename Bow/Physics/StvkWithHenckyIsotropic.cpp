#include "StvkWithHenckyIsotropic.h"
#include <Bow/Math/SVD.h>
#include <Bow/Math/PolarDecomposition.h>
#include "ConstitutiveModel.h"
#include <Bow/Math/Utils.h>
#include <Bow/Math/MathTools.h>
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace ConstitutiveModel {
namespace StvkWithHenckyIsotropic {

namespace internal {
template <class DerivedV, class DerivedW, class Scalar>
inline void dpsi_dsigma(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedW>& de_dsigma)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedV);
    using T = Scalar;
    static const int dim = DerivedV::SizeAtCompileTime;
    Vector<T, dim> log_sigma = sigma.array().abs().log();
    T sum_log_sigma = log_sigma.sum();
    const T inv0 = 1.0 / sigma[0];
    de_dsigma[0] = (2 * mu * log_sigma(0) + lam * sum_log_sigma) * inv0;
    const T inv1 = 1.0 / sigma[1];
    de_dsigma[1] = (2 * mu * log_sigma(1) + lam * sum_log_sigma) * inv1;
    if constexpr (dim == 3) {
        const T inv2 = 1.0 / sigma[2];
        de_dsigma[2] = (2 * mu * log_sigma(2) + lam * sum_log_sigma) * inv2;
    }
}

template <class DerivedV, class DerivedM, class Scalar>
inline void d2psi_dsigma2(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedM>& d2e_dsigma2)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedV);
    using T = Scalar;
    static const int dim = DerivedV::SizeAtCompileTime;
    Vector<T, dim> log_sigma = sigma.array().abs().log();
    T sum_log_sigma = log_sigma.sum();
    const double inv2_0 = T(1) / (sigma[0] * sigma[0]);
    d2e_dsigma2(0, 0) = (2 * mu * (1 - log_sigma(0)) + lam * (1 - sum_log_sigma)) * inv2_0;
    const double inv2_1 = T(1) / (sigma[1] * sigma[1]);
    d2e_dsigma2(1, 1) = (2 * mu * (1 - log_sigma(1)) + lam * (1 - sum_log_sigma)) * inv2_1;
    d2e_dsigma2(0, 1) = d2e_dsigma2(1, 0) = lam / sigma[0] / sigma[1];
    if constexpr (dim == 3) {
        const double inv2_2 = T(1) / sigma[2] / sigma[2];
        d2e_dsigma2(2, 2) = (2 * mu * (1 - log_sigma(2)) + lam * (1 - sum_log_sigma)) * inv2_2;
        d2e_dsigma2(1, 2) = d2e_dsigma2(2, 1) = lam / sigma[1] / sigma[2];
        d2e_dsigma2(2, 0) = d2e_dsigma2(0, 2) = lam / sigma[2] / sigma[0];
    }
}

template <class DerivedV, class DerivedW, class Scalar>
inline void B_left_coeff(const Eigen::MatrixBase<DerivedV>& sigma, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedW>& left_coeff)
{
    // https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
    // Eq 77
    // (psiA-psiB)/(sigmaA-sigmaB)
    using T = typename DerivedV::Scalar;
    const int dim = DerivedV::SizeAtCompileTime;
    Vector<T, dim> log_sigma = sigma.array().abs().log();
    T eps = 1e-6;
    if constexpr (dim == 2) {
        T q = std::max(sigma(0) / sigma(1) - 1, -1 + eps);
        T h = (std::abs(q) < eps) ? 1 : (std::log1p(q) / q);
        T t = h / sigma(1);
        T z = log_sigma(1) - t * sigma(1);
        left_coeff[0] = -(lam * (log_sigma(0) + log_sigma(1)) + 2 * mu * z) / sigma.prod() / T(2);
    }
    else {
        T sum_log_sigma = log_sigma.sum();
        left_coeff[0] = -(lam * sum_log_sigma + 2 * mu * MATH_TOOLS::diff_interlock_log_over_diff(sigma(0), std::abs(sigma(1)), log_sigma(1), eps)) / (sigma[0] * sigma[1]) / T(2);
        left_coeff[1] = -(lam * sum_log_sigma + 2 * mu * MATH_TOOLS::diff_interlock_log_over_diff(sigma(1), std::abs(sigma(2)), log_sigma(2), eps)) / (sigma[1] * sigma[2]) / T(2);
        left_coeff[2] = -(lam * sum_log_sigma + 2 * mu * MATH_TOOLS::diff_interlock_log_over_diff(sigma(0), std::abs(sigma(2)), log_sigma(2), eps)) / (sigma[0] * sigma[2]) / T(2);
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
    Vector<T, dim> log_sigma_squared = sigma.array().abs().log().square();
    T trace_log_sigma = sigma.array().abs().log().sum();
    return mu * log_sigma_squared.sum() + (T).5 * lam * trace_log_sigma * trace_log_sigma;
}

template <class DerivedM, class DerivedN, class Scalar>
BOW_INLINE void first_piola(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedN>& P)
{
    using T = Scalar;
    static const int dim = DerivedM::RowsAtCompileTime;
    Bow::Matrix<T, dim, dim> U, V;
    Bow::Vector<T, dim> sigma;
    Math::svd(F, U, sigma, V);
    Bow::Vector<T, dim> de_dsigma;
    internal::dpsi_dsigma(sigma, mu, lam, de_dsigma);
    P = U * de_dsigma.asDiagonal() * V.transpose();
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

template <bool projectPD, class DerivedM, class Scalar>
BOW_INLINE void first_piola_differential(const Eigen::MatrixBase<DerivedM>& F, const Eigen::MatrixBase<DerivedM>& dF, Scalar mu, Scalar lam, Eigen::MatrixBase<DerivedM>& dP)
{
    BOW_NOT_IMPLEMENTED
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template float psi(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam);
template void first_piola(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& P);
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const Eigen::MatrixBase<Matrix<float, 2, 2>>& dF, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& dP);
template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& F, const Eigen::MatrixBase<Matrix<float, 2, 2>>& dF, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 2, 2>>& dP);
#endif
#ifdef BOW_COMPILE_3D
template float psi(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam);
template void first_piola(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& P);
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const Eigen::MatrixBase<Matrix<float, 3, 3>>& dF, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& dP);
template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& F, const Eigen::MatrixBase<Matrix<float, 3, 3>>& dF, const float mu, const float lam, Eigen::MatrixBase<Matrix<float, 3, 3>>& dP);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template double psi(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam);
template void first_piola(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& P);
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const Eigen::MatrixBase<Matrix<double, 2, 2>>& dF, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& dP);
template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& F, const Eigen::MatrixBase<Matrix<double, 2, 2>>& dF, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 2, 2>>& dP);
#endif
#ifdef BOW_COMPILE_3D
template double psi(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam);
template void first_piola(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& P);
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
template void first_piola_differential<true>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const Eigen::MatrixBase<Matrix<double, 3, 3>>& dF, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& dP);
template void first_piola_differential<false>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& F, const Eigen::MatrixBase<Matrix<double, 3, 3>>& dF, const double mu, const double lam, Eigen::MatrixBase<Matrix<double, 3, 3>>& dP);
#endif
#endif
#endif
}
}
} // namespace Bow::ConstitutiveModel::StvkWithHenckyIsotropic
