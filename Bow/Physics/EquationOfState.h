#pragma once

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow::ConstitutiveModel::EquationOfState {

template <class Scalar>
BOW_INLINE Scalar psi(const Scalar& J, const Scalar bulk, const Scalar gamma);

template <class Scalar>
BOW_INLINE void first_piola(const Scalar& J, const Scalar bulk, const Scalar gamma, Scalar& P);

template <class Scalar>
BOW_INLINE void first_piola_derivative(const Scalar& J, const Scalar bulk, const Scalar gamma, Scalar& dPdF);

} // namespace Bow::ConstitutiveModel::EquationOfState

#ifndef BOW_STATIC_LIBRARY
#include "EquationOfState.cpp"
#endif
