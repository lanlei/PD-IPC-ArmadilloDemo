#ifndef Edge_Edge_Distance_H
#define Edge_Edge_Distance_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow::Geometry::IPC {

template <class T>
BOW_INLINE T edge_edge_cross_norm2(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1);

template <class T>
BOW_INLINE void edge_edge_cross_norm2_gradient(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, Vector<T, 12>& grad);

template <class T>
BOW_INLINE void edge_edge_cross_norm2_hessian(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, Matrix<T, 12, 12>& hessian);

template <class T>
BOW_INLINE T edge_edge_mollifier_threshold(const Vector<T, 3>& ea0_rest, const Vector<T, 3>& ea1_rest, const Vector<T, 3>& eb0_rest, const Vector<T, 3>& eb1_rest);

template <class T>
BOW_INLINE T edge_edge_mollifier(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, const T eps_x);

template <class T>
BOW_INLINE void edge_edge_mollifier_gradient(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, const T eps_x, Vector<T, 12>& grad);

template <class T>
BOW_INLINE void edge_edge_mollifier_hessian(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, const T eps_x, Matrix<T, 12, 12>& hessian);

template <class T>
BOW_INLINE T edge_edge_distance(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1);

template <class T>
BOW_INLINE void edge_edge_distance_gradient(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, Vector<T, 12>& grad);

template <class T>
BOW_INLINE void edge_edge_distance_hessian(const Vector<T, 3>& ea0, const Vector<T, 3>& ea1, const Vector<T, 3>& eb0, const Vector<T, 3>& eb1, Matrix<T, 12, 12>& hessian);

} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "EdgeEdgeDistance.cpp"
#endif

#endif