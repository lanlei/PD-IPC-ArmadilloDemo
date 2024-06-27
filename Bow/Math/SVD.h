#ifndef SVD_H
#define SVD_H

#include <Eigen/Dense>
#include <Bow/Macros.h>

namespace Bow {
namespace Math {
/**
   \brief SVD (singular value decomposition) A=USV'
   \param[in] A Input square matrix.
   \param[out] U Robustly a rotation matrix.
   \param[out] Sigma Vector of singular values sorted with decreasing magnitude. The second one can be negative.
   \param[out] V Robustly a rotation matrix.
*/
template <class DerivedA, class DerivedU, class DerivedSigma, class DerivedV>
BOW_INLINE void svd(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedU>& U,
    Eigen::MatrixBase<DerivedSigma>& Sigma,
    Eigen::MatrixBase<DerivedV>& V);
}
} // namespace Bow::Math

#ifndef BOW_STATIC_LIBRARY
#include "SVD.cpp"
#endif

#endif