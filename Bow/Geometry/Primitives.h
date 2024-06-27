#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow {
namespace Geometry {
/**
     * Triangulate a 2D/3D cube with identical triangles.
     */
template <class T, int dim>
BOW_INLINE void cube(const Vector<int, dim> resolution, const T dx, Field<Vector<T, dim>>& points, Field<Vector<int, dim + 1>>& elements, const Vector<T, dim> center = Vector<T, dim>::Zero());
}
} // namespace Bow::Geometry

#ifndef BOW_STATIC_LIBRARY
#include "Primitives.cpp"
#endif

#endif