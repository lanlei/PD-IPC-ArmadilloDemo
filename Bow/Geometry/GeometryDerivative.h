#ifndef GEOMETRY_DERIVATIVE_H
#define GEOMETRY_DERIVATIVE_H
#include <Bow/Types.h>
#include <Bow/Macros.h>
namespace Bow {
namespace Geometry {
template <class T>
T simplex_volume(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2);

template <class T, class Derived>
void simplex_volume_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, Eigen::MatrixBase<Derived>& grad);

template <class T>
T simplex_volume(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3);

template <class T, class Derived>
void simplex_volume_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3, Eigen::MatrixBase<Derived>& grad);

template <class T, int dim>
void center_of_mass(const Field<Vector<T, dim>>& x, const Field<T>& mass, Vector<T, dim>& cm);

/**
 * \brief compute df(cm)/dx, given df/dcm
 */
template <class T, class DerivedX, class DerivedY>
void center_of_mass_gradient(const Field<T>& mass, const Eigen::MatrixBase<DerivedX>& dLdcm, Eigen::MatrixBase<DerivedY>& grad);
}
} // namespace Bow::Geometry

#ifndef BOW_STATIC_LIBRARY
#include "GeometryDerivative.cpp"
#endif

#endif