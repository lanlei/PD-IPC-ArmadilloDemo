#ifndef MATRIX_DERIVATIVE_H
#define MATRIX_DERIVATIVE_H

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow {
namespace Math {

/**
 * \brief dAX/dX. AX and X are flattened in column-major order.
 * \param[in] A a matrix of size (m,n)
 * \param[in] X a matrix of size (n,l)
 * \param[out] dAX/dX as a matrix of size (m*l, n*l)
 */
template <class DerivedA, class DerivedX, class DerivedJ>
void dAX(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, Eigen::MatrixBase<DerivedJ>& J);

/**
 * \brief dAX/dX. XA and X are flattened in column-major order.
 * \param[in] A a matrix of size (n,l)
 * \param[in] X a matrix of size (m,n)
 * \param[out] dAX/dX as a matrix of size (m*l, m*n)
 */
template <class DerivedA, class DerivedX, class DerivedJ>
void dXA(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, Eigen::MatrixBase<DerivedJ>& J);

/**
 * \brief dX^{-1}/dX. X^{-1} and X are flattened in column-major order.
 * \param[in] X a matrix of size (n,n)
 * \param[out] dX^{-1}/dX as a matrix of size (n*n, n*n)
 */
template <class DerivedA, class DerivedJ>
void dXinv(const Eigen::MatrixBase<DerivedA>& X, Eigen::MatrixBase<DerivedJ>& J);

/**
 * \brief df(AX^{-1})/dX.
 * \param[in] A a matrix of size (m,n)
 * \param[in] X a matrix of size (n,n)
 * \param[in] df df(F)/dF as a matrix of size (m, n)
 * \param[out] df(AX^{-1})/dX as a matrix of size (n, n)
 */
template <class DerivedA, class DerivedX, class DerivedP, class DerivedJ>
void df_AXinv(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedX>& X, const Eigen::MatrixBase<DerivedP>& df, Eigen::MatrixBase<DerivedJ>& J);
}
} // namespace Bow::Math

#ifndef BOW_STATIC_LIBRARY
#include "MatrixDerivative.cpp"
#endif

#endif