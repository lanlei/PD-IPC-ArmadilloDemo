#ifndef EDGE_EDGE_CCD_H
#define EDGE_EDGE_CCD_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow::Geometry::IPC {
template <class T>
BOW_INLINE bool edge_edge_cd_broadphase(
    const Vector<T, 3>& ea0,
    const Vector<T, 3>& ea1,
    const Vector<T, 3>& eb0,
    const Vector<T, 3>& eb1,
    T dist);

template <class T>
BOW_INLINE bool edge_edge_ccd_broadphase(
    const Vector<T, 3>& ea0,
    const Vector<T, 3>& ea1,
    const Vector<T, 3>& eb0,
    const Vector<T, 3>& eb1,
    const Vector<T, 3>& dea0,
    const Vector<T, 3>& dea1,
    const Vector<T, 3>& deb0,
    const Vector<T, 3>& deb1,
    T dist);

template <class T>
BOW_INLINE T edge_edge_ccd(
    const Vector<T, 3>& ea0,
    const Vector<T, 3>& ea1,
    const Vector<T, 3>& eb0,
    const Vector<T, 3>& eb1,
    const Vector<T, 3>& dea0,
    const Vector<T, 3>& dea1,
    const Vector<T, 3>& deb0,
    const Vector<T, 3>& deb1,
    T eta, T thickness);
} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "EdgeEdgeCcd.cpp"
#endif

#endif