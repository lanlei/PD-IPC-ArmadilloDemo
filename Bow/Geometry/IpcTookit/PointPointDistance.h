#ifndef PP_DISTANCE_H
#define PP_DISTANCE_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace Geometry {
namespace IPC {

template <class T, int dim>
BOW_INLINE T point_point_distance(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <class T, int dim>
BOW_INLINE void point_point_distance_gradient(const Vector<T, dim>& a, const Vector<T, dim>& b,
    Vector<T, dim * 2>& grad);

template <class T, int dim>
BOW_INLINE void point_point_distance_hessian(const Vector<T, dim>& a, const Vector<T, dim>& b,
    Matrix<T, dim * 2, dim * 2>& hess);

}}} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "PointPointDistance.cpp"
#endif

#endif