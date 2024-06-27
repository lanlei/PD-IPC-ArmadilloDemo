#ifndef PE_DISTANCE_H
#define PE_DISTANCE_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow {
namespace Geometry {
namespace IPC {
template <class T, int dim>
BOW_INLINE T point_edge_distance(const Vector<T, dim>& p, const Vector<T, dim>& e0, const Vector<T, dim>& e1);

template <class T, int dim>
BOW_INLINE void point_edge_distance_gradient(const Vector<T, dim>& p,
    const Vector<T, dim>& e0,
    const Vector<T, dim>& e1,
    Vector<T, 3 * dim>& grad);

template <class T, int dim>
BOW_INLINE void point_edge_distance_hessian(const Vector<T, dim>& p,
    const Vector<T, dim>& e0,
    const Vector<T, dim>& e1,
    Matrix<T, dim * 3, dim * 3>& hess);

}}} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "PointEdgeDistance.cpp"
#endif

#endif