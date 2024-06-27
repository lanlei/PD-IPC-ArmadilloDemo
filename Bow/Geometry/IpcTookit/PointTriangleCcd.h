#ifndef POINT_TRIANGLE_CCD_H
#define POINT_TRIANGLE_CCD_H

#include <Bow/Types.h>
#include <Bow/Macros.h>

namespace Bow::Geometry::IPC {
template <class T>
bool point_triangle_cd_broadphase(
    const Vector<T, 3>& p,
    const Vector<T, 3>& t0,
    const Vector<T, 3>& t1,
    const Vector<T, 3>& t2,
    T dist);

template <class T>
bool point_triangle_ccd_broadphase(
    const Vector<T, 3>& p,
    const Vector<T, 3>& t0,
    const Vector<T, 3>& t1,
    const Vector<T, 3>& t2,
    const Vector<T, 3>& dp,
    const Vector<T, 3>& dt0,
    const Vector<T, 3>& dt1,
    const Vector<T, 3>& dt2,
    T dist);

template <class T>
T point_triangle_ccd(
    const Vector<T, 3>& p,
    const Vector<T, 3>& t0,
    const Vector<T, 3>& t1,
    const Vector<T, 3>& t2,
    const Vector<T, 3>& dp,
    const Vector<T, 3>& dt0,
    const Vector<T, 3>& dt1,
    const Vector<T, 3>& dt2,
    T eta, T thickness);

template <class T>
T point_triangle_ccd2(
	int tid,
	const Vector<T, 3>& p,
	const Vector<T, 3>& t0,
	const Vector<T, 3>& t1,
	const Vector<T, 3>& t2,
	const Vector<T, 3>& dp,
	const Vector<T, 3>& dt0,
	const Vector<T, 3>& dt1,
	const Vector<T, 3>& dt2,
	T eta, T thickness);
} // namespace Bow::Geometry::IPC






#ifndef BOW_STATIC_LIBRARY
#include "PointTriangleCcd.cpp"
#endif

#endif