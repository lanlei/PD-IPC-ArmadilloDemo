#include "PolarDecomposition.h"
#include "GivensRotation.h"
#include "SVD.h"
#include <Bow/Types.h>

namespace Bow {
namespace Math {
template <class DerivedA, class DerivedR, class DerivedS>
BOW_INLINE void polar_decomposition(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedR>& R,
    Eigen::MatrixBase<DerivedS>& S)
{
    using T = typename DerivedA::Scalar;
    static const int dim = DerivedA::RowsAtCompileTime;
    if constexpr (dim == 1) {
        S = A;
        R(0, 0) = 1.;
    }
    else if constexpr (dim == 2) {
        using Bow::Math::GivensRotation;
        GivensRotation<T> r(0, 1);
        Bow::Vector<T, 2> x(A(0, 0) + A(1, 1), A(1, 0) - A(0, 1));
        T denominator = x.norm();
        r.m_c = (T)1;
        r.m_s = (T)0;
        if (denominator != 0) {
            /*
            No need to use a tolerance here because x(0) and x(1) always have
            smaller magnitude then denominator, therefore overflow never happens.
            */
            r.m_c = x(0) / denominator;
            r.m_s = -x(1) / denominator;
        }
        S = A;
        r.row_rotation(S);
        r.fill(R);
    }
    else if constexpr (dim == 3) {
        Bow::Matrix<T, 3, 3> U, V;
        Bow::Vector<T, 3> sigma;
        svd(A, U, sigma, V);
        R = U * V.transpose();
        S = V * sigma.asDiagonal() * V.transpose();
    }
    else {
        static_assert(dim <= 3, "polar_decomposition is only implemented for 1x1, 2x2 and 3x3 matrices.");
    }
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template void polar_decomposition(const Eigen::MatrixBase<Matrix<float, 2, 2>>& A, Eigen::MatrixBase<Matrix<float, 2, 2>>& R, Eigen::MatrixBase<Matrix<float, 2, 2>>& S);
#endif
#ifdef BOW_COMPILE_3D
template void polar_decomposition(const Eigen::MatrixBase<Matrix<float, 3, 3>>& A, Eigen::MatrixBase<Matrix<float, 3, 3>>& R, Eigen::MatrixBase<Matrix<float, 3, 3>>& S);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template void polar_decomposition(const Eigen::MatrixBase<Matrix<double, 2, 2>>& A, Eigen::MatrixBase<Matrix<double, 2, 2>>& R, Eigen::MatrixBase<Matrix<double, 2, 2>>& S);
#endif
#ifdef BOW_COMPILE_3D
template void polar_decomposition(const Eigen::MatrixBase<Matrix<double, 3, 3>>& A, Eigen::MatrixBase<Matrix<double, 3, 3>>& R, Eigen::MatrixBase<Matrix<double, 3, 3>>& S);
#endif
#endif
#endif

}} // namespace Bow::Math