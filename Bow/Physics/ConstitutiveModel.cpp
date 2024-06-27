#include "ConstitutiveModel.h"
#include <Bow/Math/SVD.h>
#include <Bow/Types.h>
#include <Bow/Math/Utils.h>

namespace Bow {
namespace ConstitutiveModel {
template <class T>
BOW_INLINE std::pair<T, T> lame_paramters(T E, T nu)
{
    T mu = 0.5 * E / (1 + nu);
    T lam = E * nu / ((1 + nu) * (1 - 2 * nu));
    return std::make_pair(mu, lam);
}

template <class T>
BOW_INLINE std::pair<T, T> E_nu(T mu, T lam)
{
    T mu_plus_lam = mu + lam;
    T nu = lam / (2 * mu_plus_lam);
    T E = mu * (2 + lam / mu_plus_lam);
    return std::make_pair(E, nu);
}

template <bool projectPD, class DerivedU, class DerivedS, class DerivedV, class DerivedG, class DerivedB, class DerivedHS, class DerivedHF>
BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedU>& U, const Eigen::MatrixBase<DerivedS>& sigma, const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedG>& de_dsigma, const Eigen::MatrixBase<DerivedB>& B_left_coeff, const Eigen::MatrixBase<DerivedHS>& _d2e_dsigma2, Eigen::MatrixBase<DerivedHF>& dPdF)
{
    using T = typename DerivedS::Scalar;
    const int dim = DerivedS::SizeAtCompileTime;
    Bow::Matrix<T, dim, dim> d2e_dsigma2 = _d2e_dsigma2;
    if constexpr (projectPD)
        Math::make_pd(d2e_dsigma2);
    if constexpr (dim == 2) {
        T left_coef = B_left_coeff[0];
        T right_coef = de_dsigma[0] + de_dsigma[1];
        T sum_sigma = std::max(sigma[0] + sigma[1], T(0.000001));
        right_coef /= (T(2) * sum_sigma);
        Bow::Matrix<T, 2, 2> B;
        B << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        if constexpr (projectPD) {
            Math::make_pd(B);
        }
        Bow::Matrix<T, dim * dim, dim * dim> M;
        M.setZero();
        M(0, 0) = d2e_dsigma2(0, 0);
        M(0, 3) = d2e_dsigma2(0, 1);
        M(1, 1) = B(0, 0);
        M(1, 2) = B(0, 1);
        M(2, 1) = B(1, 0);
        M(2, 2) = B(1, 1);
        M(3, 0) = d2e_dsigma2(1, 0);
        M(3, 3) = d2e_dsigma2(1, 1);
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                for (int s = 0; s < 2; ++s)
                    for (int r = 0; r < 2; ++r) {
                        int ij = j * 2 + i;
                        int rs = s * 2 + r;
                        dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0)
                            + M(0, 3) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1)
                            + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1)
                            + M(1, 2) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0)
                            + M(2, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1)
                            + M(2, 2) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0)
                            + M(3, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0)
                            + M(3, 3) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1);
                    }
    }
    else {
        T left_coef = B_left_coeff[0];
        T right_coef = de_dsigma[0] + de_dsigma[1];
        T sum_sigma = std::max(sigma[0] + sigma[1], T(0.000001));
        right_coef /= (T(2) * sum_sigma);
        Bow::Matrix<T, 2, 2> B0;
        B0 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        if constexpr (projectPD)
            Math::make_pd(B0);

        left_coef = B_left_coeff[1];
        right_coef = de_dsigma[1] + de_dsigma[2];
        sum_sigma = std::max(sigma[1] + sigma[2], T(0.000001));
        right_coef /= (T(2) * sum_sigma);
        Bow::Matrix<T, 2, 2> B1;
        B1 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        if constexpr (projectPD)
            Math::make_pd(B1);

        left_coef = B_left_coeff[2];
        right_coef = de_dsigma[2] + de_dsigma[0];
        sum_sigma = std::max(sigma[2] + sigma[0], T(0.000001));
        right_coef /= (T(2) * sum_sigma);
        Bow::Matrix<T, 2, 2> B2;
        B2 << left_coef + right_coef, left_coef - right_coef, left_coef - right_coef, left_coef + right_coef;
        if constexpr (projectPD)
            Math::make_pd(B2);

        Bow::Matrix<T, dim * dim, dim * dim> M;
        M.setZero();
        M(0, 0) = d2e_dsigma2(0, 0);
        M(0, 4) = d2e_dsigma2(0, 1);
        M(0, 8) = d2e_dsigma2(0, 2);
        M(4, 0) = d2e_dsigma2(1, 0);
        M(4, 4) = d2e_dsigma2(1, 1);
        M(4, 8) = d2e_dsigma2(1, 2);
        M(8, 0) = d2e_dsigma2(2, 0);
        M(8, 4) = d2e_dsigma2(2, 1);
        M(8, 8) = d2e_dsigma2(2, 2);
        M(1, 1) = B0(0, 0);
        M(1, 3) = B0(0, 1);
        M(3, 1) = B0(1, 0);
        M(3, 3) = B0(1, 1);
        M(5, 5) = B1(0, 0);
        M(5, 7) = B1(0, 1);
        M(7, 5) = B1(1, 0);
        M(7, 7) = B1(1, 1);
        M(2, 2) = B2(1, 1);
        M(2, 6) = B2(1, 0);
        M(6, 2) = B2(0, 1);
        M(6, 6) = B2(0, 0);

        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                for (int s = 0; s < 3; ++s)
                    for (int r = 0; r < 3; ++r) {
                        int ij = j * 3 + i;
                        int rs = s * 3 + r;
                        dPdF(ij, rs) = M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0)
                            + M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1)
                            + M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2)
                            + M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0)
                            + M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1)
                            + M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2)
                            + M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0)
                            + M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1)
                            + M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2)
                            + M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1)
                            + M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0)
                            + M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1)
                            + M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0)
                            + M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2)
                            + M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1)
                            + M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2)
                            + M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1)
                            + M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2)
                            + M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0)
                            + M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2)
                            + M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);
                    }
    }
}

template <bool projectPD, class DerivedU, class DerivedS, class DerivedV, class DerivedG, class DerivedB, class DerivedHS, class DerivedHF>
BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedU>& U, const Eigen::MatrixBase<DerivedS>& sigma, const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedG>& de_dsigma, const Eigen::MatrixBase<DerivedB>& B_left_coeff, const Eigen::MatrixBase<DerivedHS>& d2e_dsigma2, Eigen::MatrixBase<DerivedHF>& dPdF);

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template std::pair<float, float> lame_paramters(float E, float nu);
template std::pair<float, float> E_nu(float E, float nu);
#ifdef BOW_COMPILE_2D
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& U, const Eigen::MatrixBase<Vector<float, 2>>& sigma, const Eigen::MatrixBase<Matrix<float, 2, 2>>& V, const Eigen::MatrixBase<Vector<float, 2>>& de_dsigma, const Eigen::MatrixBase<Vector<float, 1>>& B_left_coeff, const Eigen::MatrixBase<Matrix<float, 2, 2>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 2, 2>>& U, const Eigen::MatrixBase<Vector<float, 2>>& sigma, const Eigen::MatrixBase<Matrix<float, 2, 2>>& V, const Eigen::MatrixBase<Vector<float, 2>>& de_dsigma, const Eigen::MatrixBase<Vector<float, 1>>& B_left_coeff, const Eigen::MatrixBase<Matrix<float, 2, 2>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<float, 4, 4>>& dPdF);
#endif
#ifdef BOW_COMPILE_3D
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& U, const Eigen::MatrixBase<Vector<float, 3>>& sigma, const Eigen::MatrixBase<Matrix<float, 3, 3>>& V, const Eigen::MatrixBase<Vector<float, 3>>& de_dsigma, const Eigen::MatrixBase<Vector<float, 3>>& B_left_coeff, const Eigen::MatrixBase<Matrix<float, 3, 3>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<float, 3, 3>>& U, const Eigen::MatrixBase<Vector<float, 3>>& sigma, const Eigen::MatrixBase<Matrix<float, 3, 3>>& V, const Eigen::MatrixBase<Vector<float, 3>>& de_dsigma, const Eigen::MatrixBase<Vector<float, 3>>& B_left_coeff, const Eigen::MatrixBase<Matrix<float, 3, 3>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<float, 9, 9>>& dPdF);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
template std::pair<double, double> lame_paramters(double E, double nu);
template std::pair<double, double> E_nu(double E, double nu);
#ifdef BOW_COMPILE_2D
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& U, const Eigen::MatrixBase<Vector<double, 2>>& sigma, const Eigen::MatrixBase<Matrix<double, 2, 2>>& V, const Eigen::MatrixBase<Vector<double, 2>>& de_dsigma, const Eigen::MatrixBase<Vector<double, 1>>& B_left_coeff, const Eigen::MatrixBase<Matrix<double, 2, 2>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 2, 2>>& U, const Eigen::MatrixBase<Vector<double, 2>>& sigma, const Eigen::MatrixBase<Matrix<double, 2, 2>>& V, const Eigen::MatrixBase<Vector<double, 2>>& de_dsigma, const Eigen::MatrixBase<Vector<double, 1>>& B_left_coeff, const Eigen::MatrixBase<Matrix<double, 2, 2>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<double, 4, 4>>& dPdF);
#endif
#ifdef BOW_COMPILE_3D
template void first_piola_derivative<true>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& U, const Eigen::MatrixBase<Vector<double, 3>>& sigma, const Eigen::MatrixBase<Matrix<double, 3, 3>>& V, const Eigen::MatrixBase<Vector<double, 3>>& de_dsigma, const Eigen::MatrixBase<Vector<double, 3>>& B_left_coeff, const Eigen::MatrixBase<Matrix<double, 3, 3>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
template void first_piola_derivative<false>(const Eigen::MatrixBase<Matrix<double, 3, 3>>& U, const Eigen::MatrixBase<Vector<double, 3>>& sigma, const Eigen::MatrixBase<Matrix<double, 3, 3>>& V, const Eigen::MatrixBase<Vector<double, 3>>& de_dsigma, const Eigen::MatrixBase<Vector<double, 3>>& B_left_coeff, const Eigen::MatrixBase<Matrix<double, 3, 3>>& d2e_dsigma2, Eigen::MatrixBase<Matrix<double, 9, 9>>& dPdF);
#endif
#endif
#endif
}
} // namespace Bow::ConstitutiveModel