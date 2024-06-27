#include "MatrixDerivative.h"

namespace Bow {
namespace Math {

template <class DerivedA, class DerivedX, class DerivedJ>
void dAX(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, Eigen::MatrixBase<DerivedJ>& J)
{
    if constexpr (DerivedJ::RowsAtCompileTime == Eigen::Dynamic)
        J.derived().resize(A.rows() * X.cols(), X.size());
    J.setZero();
    for (int i = 0; i < A.rows(); i++)
        for (int j = 0; j < X.cols(); ++j)
            for (int m = 0; m < X.rows(); ++m)
                J(j * A.rows() + i, j * X.rows() + m) = A(i, m);
}

template <class DerivedA, class DerivedX, class DerivedJ>
void dXA(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, Eigen::MatrixBase<DerivedJ>& J)
{
    if constexpr (DerivedJ::RowsAtCompileTime == Eigen::Dynamic)
        J.derived().resize(X.rows() * A.cols(), X.size());
    J.setZero();
    for (int i = 0; i < X.rows(); i++)
        for (int j = 0; j < A.cols(); ++j)
            for (int n = 0; n < A.rows(); ++n)
                J(j * X.rows() + i, n * X.rows() + i) = A(n, j);
}

template <class DerivedA, class DerivedJ>
void dXinv(const Eigen::MatrixBase<DerivedA>& X, Eigen::MatrixBase<DerivedJ>& J)
{
    if constexpr (DerivedJ::RowsAtCompileTime == Eigen::Dynamic)
        J.derived().resize(X.size(), X.size());
    J.setZero();
    Bow::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic> Xinv = X.inverse();
    int rows = X.rows();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < rows; ++j)
            for (int m = 0; m < rows; ++m)
                for (int n = 0; n < rows; ++n)
                    J(n * rows + m, j * rows + i) = -Xinv(m, i) * Xinv(j, n);
}

template <class DerivedA, class DerivedX, class DerivedP, class DerivedJ>
void df_AXinv(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, const Eigen::MatrixBase<DerivedP>& df, Eigen::MatrixBase<DerivedJ>& J)
{
    Bow::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic> Xinv = X.inverse();
    Bow::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic> F = A * Xinv;
    J = -F.transpose() * df * Xinv.transpose();
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void dAX(
    const Eigen::MatrixBase<Matrix<float, 2, 2>>&,
    const Eigen::MatrixBase<Matrix<float, 2, 2>>&,
    Eigen::MatrixBase<Matrix<float, 4, 4>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<float, 2, 2>>&,
    const Eigen::MatrixBase<Matrix<float, 2, 2>>&,
    Eigen::MatrixBase<Matrix<float, 4, 4>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<float, 2, 2>>&,
    Eigen::MatrixBase<Matrix<float, 4, 4>>&);

template void dAX(
    const Eigen::MatrixBase<Matrix<float, 3, 3>>&,
    const Eigen::MatrixBase<Matrix<float, 3, 3>>&,
    Eigen::MatrixBase<Matrix<float, 9, 9>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<float, 3, 3>>&,
    const Eigen::MatrixBase<Matrix<float, 3, 3>>&,
    Eigen::MatrixBase<Matrix<float, 9, 9>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<float, 3, 3>>&,
    Eigen::MatrixBase<Matrix<float, 9, 9>>&);

template void dAX(
    const Eigen::MatrixBase<Matrix<float, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<float, -1, -1>>&,
    Eigen::MatrixBase<Matrix<float, -1, -1>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<float, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<float, -1, -1>>&,
    Eigen::MatrixBase<Matrix<float, -1, -1>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<float, -1, -1>>&,
    Eigen::MatrixBase<Matrix<float, -1, -1>>&);
#endif
#ifdef BOW_COMPILE_DOUBLE
template void dAX(
    const Eigen::MatrixBase<Matrix<double, 2, 2>>&,
    const Eigen::MatrixBase<Matrix<double, 2, 2>>&,
    Eigen::MatrixBase<Matrix<double, 4, 4>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<double, 2, 2>>&,
    const Eigen::MatrixBase<Matrix<double, 2, 2>>&,
    Eigen::MatrixBase<Matrix<double, 4, 4>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<double, 2, 2>>&,
    Eigen::MatrixBase<Matrix<double, 4, 4>>&);

template void dAX(
    const Eigen::MatrixBase<Matrix<double, 3, 3>>&,
    const Eigen::MatrixBase<Matrix<double, 3, 3>>&,
    Eigen::MatrixBase<Matrix<double, 9, 9>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<double, 3, 3>>&,
    const Eigen::MatrixBase<Matrix<double, 3, 3>>&,
    Eigen::MatrixBase<Matrix<double, 9, 9>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<double, 3, 3>>&,
    Eigen::MatrixBase<Matrix<double, 9, 9>>&);

template void dAX(
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    Eigen::MatrixBase<Matrix<double, -1, -1>>&);
template void dXA(
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    Eigen::MatrixBase<Matrix<double, -1, -1>>&);
template void dXinv(
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    Eigen::MatrixBase<Matrix<double, -1, -1>>&);
template void df_AXinv(const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    const Eigen::MatrixBase<Matrix<double, -1, -1>>&,
    Eigen::MatrixBase<Matrix<double, -1, -1>>&);
#endif
#endif
}
} // namespace Bow::Math