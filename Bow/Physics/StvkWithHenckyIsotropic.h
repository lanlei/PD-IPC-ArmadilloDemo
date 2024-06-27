#ifndef STVK_HENCKY_ISOTROPIC_H
#define STVK_HENCKY_ISOTROPIC_H
#include <Bow/Macros.h>
#include <Bow/Types.h>
namespace Bow::ConstitutiveModel::StvkWithHenckyIsotropic {
template <class DerivedM, class Scalar>
BOW_INLINE Scalar psi(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam);

template <class DerivedM, class DerivedN, class Scalar>
BOW_INLINE void first_piola(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedN>& P);

template <bool projectPD = true, class DerivedM, class DerivedH, class Scalar>
BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedM>& F, const Scalar mu, const Scalar lam, Eigen::MatrixBase<DerivedH>& dPdF);

template <bool projectPD, class DerivedM, class Scalar>
BOW_INLINE void first_piola_differential(const Eigen::MatrixBase<DerivedM>& F, const Eigen::MatrixBase<DerivedM>& dF, Scalar mu, Scalar lam, Eigen::MatrixBase<DerivedM>& dP);

} // namespace Bow::ConstitutiveModel::StvkWithHenckyIsotropic

#ifndef BOW_STATIC_LIBRARY
#include "StvkWithHenckyIsotropic.cpp"
#endif

#endif