#ifndef PT_DISTANCE_H
#define PT_DISTANCE_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace Geometry {
namespace IPC {

template <class T>
BOW_INLINE T point_triangle_distance(const Vector<T, 3>& p, const Vector<T, 3>& t0, const Vector<T, 3>& t1, const Vector<T, 3>& t2);

template <class T>
BOW_INLINE void point_triangle_distance_gradient(const Vector<T, 3>& p, const Vector<T, 3>& t0, const Vector<T, 3>& t1, const Vector<T, 3>& t2,
    Vector<T, 12>& grad);

template <class T>
BOW_INLINE void point_triangle_distance_hessian(const Vector<T, 3>& p, const Vector<T, 3>& t0, const Vector<T, 3>& t1, const Vector<T, 3>& t2,
    Matrix<T, 12, 12>& hess);

}}} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "PointTriangleDistance.cpp"
#endif

#endif