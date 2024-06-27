#include "LinearElasticity.h"
#include <Bow/Math/SVD.h>
#include <Bow/Math/PolarDecomposition.h>
#include "ConstitutiveModel.h"
#include <Bow/Math/Utils.h>
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace ConstitutiveModel {
namespace LinearElasticity {

template <class DerivedM, class Scalar>
BOW_INLINE Scalar psi(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam)
{
    using T = Scalar;
    static const int dim = DerivedM::RowsAtCompileTime;

    Bow::Matrix<T, dim, dim> R = Bow::Matrix<T, dim, dim>::Identity(); // can be updated per time step for lagged corotational
    Bow::Matrix<T, dim, dim> RtF = R.transpose() * F;
    Bow::Matrix<T, dim, dim> smallStrain = 0.5 * (RtF + RtF.transpose()) - Bow::Matrix<T, dim, dim>::Identity();
    T tr_smallStrain = smallStrain.diagonal().sum();

    return mu * smallStrain.squaredNorm() + lam * 0.5 * tr_smallStrain * tr_smallStrain;
}

template <class DerivedM, class DerivedN, class Scalar>
BOW_INLINE void first_piola(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedN>& P)
{
    using T = Scalar;
    static const int dim = DerivedM::RowsAtCompileTime;

    Bow::Matrix<T, dim, dim> R = Bow::Matrix<T, dim, dim>::Identity(); // can be updated per time step for lagged corotational
    Bow::Matrix<T, dim, dim> RtF = R.transpose() * F;
    Bow::Matrix<T, dim, dim> smallStrain = 0.5 * (RtF + RtF.transpose()) - Bow::Matrix<T, dim, dim>::Identity();
    T tr_smallStrain = smallStrain.diagonal().sum();

    P.noalias() = 2 * mu * R * smallStrain + lam * tr_smallStrain * R;
}

template <bool projectPD, class DerivedM, class DerivedH, class Scalar>
BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedM>& F, Scalar mu, Scalar lam, Eigen::MatrixBase<DerivedH>& dPdF)
{
    using T = Scalar;
    static const int dim = DerivedM::RowsAtCompileTime;

    Bow::Matrix<T, dim, dim> R = Bow::Matrix<T, dim, dim>::Identity(); // can be updated per time step for lagged corotational

    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) {
            int row_idx = i + j * dim;
            for (int a = 0; a < dim; ++a)
                for (int b = 0; b < dim; ++b) {
                    int col_idx = a + b * dim;
                    int ia = (i == a);
                    int jb = (j == b);
                    dPdF(row_idx, col_idx) = mu * (ia * jb + R(i, b) * R(a, j)) + lam * R(i, j) * R(a, b);
                }
        }
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
} // namespace Bow::ConstitutiveModel::LinearElasticity
