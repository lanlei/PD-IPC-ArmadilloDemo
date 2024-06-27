#ifndef POINT_EDGE_CCD_H
#define POINT_EDGE_CCD_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow::Geometry::IPC {

template <class T, int dim>
bool point_edge_cd_broadphase(
    const Vector<T, dim>& x0,
    const Vector<T, dim>& x1,
    const Vector<T, dim>& x2,
    T dist);

template <class T, int dim>
bool point_edge_ccd_broadphase(
    const Vector<T, dim>& p,
    const Vector<T, dim>& e0,
    const Vector<T, dim>& e1,
    const Vector<T, dim>& dp,
    const Vector<T, dim>& de0,
    const Vector<T, dim>& de1,
    const T dist);

template <class T>
T point_edge_ccd(
    const Vector<T, 2>& x0,
    const Vector<T, 2>& x1,
    const Vector<T, 2>& x2,
    const Vector<T, 2>& d0,
    const Vector<T, 2>& d1,
    const Vector<T, 2>& d2,
    T eta);
} // namespace Bow::Geometry::IPC

#ifndef BOW_STATIC_LIBRARY
#include "PointEdgeCcd.cpp"
#endif

#endif