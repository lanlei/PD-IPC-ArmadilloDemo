#include "EquationOfState.h"
#include <Bow/Utils/Logging.h>

namespace Bow::ConstitutiveModel::EquationOfState {

template <class Scalar>
BOW_INLINE Scalar psi(const Scalar& J, const Scalar bulk, const Scalar gamma)
{
    using T = Scalar;
    T J2 = J * J;
    T J6 = J2 * J2 * J2;
    BOW_ASSERT(gamma == 7);
    return -bulk * (1. / J6 / (-6.) - J);
}

template <class Scalar>
BOW_INLINE void first_piola(const Scalar& J, const Scalar bulk, const Scalar gamma, Scalar& P)
{
    using T = Scalar;
    T J2 = J * J;
    T J4 = J2 * J2;
    T J7 = J4 * J2 * J;
    BOW_ASSERT(gamma == 7);
    P = -bulk * (1. / J7 - 1.);
}

template <class Scalar>
BOW_INLINE void first_piola_derivative(const Scalar& J, const Scalar bulk, const Scalar gamma, Scalar& dPdJ)
{
    using T = Scalar;
    T J2 = J * J;
    T J4 = J2 * J2;
    T J8 = J4 * J4;
    BOW_ASSERT(gamma == 7);
    dPdJ = bulk * (1. / J8 * 7.);
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template float psi(const float& J, const float bulk, const float gamma);
template void first_piola(const float& J, const float bulk, const float gamma, float& P);
template void first_piola_derivative(const float& J, const float bulk, const float gamma, float& dPdJ);
#endif
#ifdef BOW_COMPILE_DOUBLE
template double psi(const double& J, const double bulk, const double gamma);
template void first_piola(const double& J, const double bulk, const double gamma, double& P);
template void first_piola_derivative(const double& J, const double bulk, const double gamma, double& dPdJ);
#endif
#endif

} // namespace Bow::ConstitutiveModel::EquationOfState