#ifndef CONSTITUTIVE_MODEL_H
#define CONSTITUTIVE_MODEL_H

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow {
namespace ConstitutiveModel {
enum Type {
    FIXED_COROTATED,
    NEO_HOOKEAN,
    LINEAR_ELASTICITY,
	AS_RIGID_AS_POSSIBLE
};

template <class T>
BOW_INLINE std::pair<T, T> lame_paramters(T E, T nu);

template <class T>
BOW_INLINE std::pair<T, T> E_nu(T mu, T lam);

template <bool projectPD = true, class DerivedU, class DerivedS, class DerivedV, class DerivedG, class DerivedB, class DerivedHS, class DerivedHF>
BOW_INLINE void first_piola_derivative(const Eigen::MatrixBase<DerivedU>& U, const Eigen::MatrixBase<DerivedS>& sigma, const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedG>& de_dsigma, const Eigen::MatrixBase<DerivedB>& B_left_coeff, const Eigen::MatrixBase<DerivedHS>& d2e_dsigma2, Eigen::MatrixBase<DerivedHF>& dPdF);
}
} // namespace Bow::ConstitutiveModel

#ifndef BOW_STATIC_LIBRARY
#include "ConstitutiveModel.cpp"
#endif

#endif