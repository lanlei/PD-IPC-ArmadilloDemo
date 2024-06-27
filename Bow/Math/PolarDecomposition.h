#ifndef POLAR_DECOMPOSITION_H
#define POLAR_DECOMPOSITION_H

#include <Eigen/Eigen>
#include <Bow/Macros.h>

namespace Bow {
namespace Math {

/**
   \brief polar decomposition. Only 1x1, 2 x 2 or 3 x 3 matrices are supported.
   \param[in] A matrix.
   \param[out] R Robustly a rotation matrix.
   \param[out] S Symmetric.

   Polar guarantees negative sign is on the small magnitude singular value.
   S is guaranteed to be the closest one to identity.
   R is guaranteed to be the closest rotation to A.
*/
template <class DerivedA, class DerivedR, class DerivedS>
void polar_decomposition(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedR>& R,
    Eigen::MatrixBase<DerivedS>& S);
}
} // namespace Bow::Math

#ifndef BOW_STATIC_LIBRARY
#include "PolarDecomposition.cpp"
#endif

#endif