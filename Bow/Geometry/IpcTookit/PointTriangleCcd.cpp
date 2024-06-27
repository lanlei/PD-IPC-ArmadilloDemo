#include "PointTriangleCcd.h"
#include "PointPointDistance.h"
#include "PointEdgeDistance.h"
#include "PointTriangleDistance.h"
#include "DistanceType.h"

namespace Bow::Geometry::IPC {

namespace internal {
template <class T>
T point_triangle_distance_unclassified(
    const Vector<T, 3>& p,
    const Vector<T, 3>& t0,
    const Vector<T, 3>& t1,
    const Vector<T, 3>& t2)
{
    switch (point_triangle_distance_type(p, t0, t1, t2)) {
    case 0:
        return point_point_distance(p, t0);
    case 1:
        return point_point_distance(p, t1);
    case 2:
        return point_point_distance(p, t2);
    case 3:
        return point_edge_distance(p, t0, t1);
    case 4:
        return point_edge_distance(p, t1, t2);
    case 5:
        return point_edge_distance(p, t2, t0);
    case 6:
        return point_triangle_distance(p, t0, t1, t2);
    default:
        return std::numeric_limits<T>::max();
    }
}
} // namespace internal

template <class T>
bool point_triangle_cd_broadphase(
    const Vector<T, 3>& p,
    const Vector<T, 3>& t0,
    const Vector<T, 3>& t1,
    const Vector<T, 3>& t2,
    T dist)
{
    auto max_tri = t0.array().max(t1.array()).max(t2.array());
    auto min_tri = t0.array().min(t1.array()).min(t2.array());
    if ((p.array() - max_tri > dist).any() || (min_tri - p.array() > dist).any())
        return false;
    else
        return true;
}

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
    T dist)
{
    auto max_p = p.array().max((p + dp).array());
    auto min_p = p.array().min((p + dp).array());
    auto max_tri = t0.array().max(t1.array()).max(t2.array()).max((t0 + dt0).array()).max((t1 + dt1).array()).max((t2 + dt2).array());
    auto min_tri = t0.array().min(t1.array()).min(t2.array()).min((t0 + dt0).array()).min((t1 + dt1).array()).min((t2 + dt2).array());
    if ((min_p - max_tri > dist).any() || (min_tri - max_p > dist).any())
        return false;
    else
        return true;
}

template <class T>
T point_triangle_ccd(
    const Vector<T, 3>& _p,
    const Vector<T, 3>& _t0,
    const Vector<T, 3>& _t1,
    const Vector<T, 3>& _t2,
    const Vector<T, 3>& _dp,
    const Vector<T, 3>& _dt0,
    const Vector<T, 3>& _dt1,
    const Vector<T, 3>& _dt2,
    T eta, T thickness)
{
    Vector<T, 3> p = _p, t0 = _t0, t1 = _t1, t2 = _t2, dp = _dp, dt0 = _dt0, dt1 = _dt1, dt2 = _dt2;
    Vector<T, 3> mov = (dt0 + dt1 + dt2 + dp) / 4;
    dt0 -= mov;
    dt1 -= mov;
    dt2 -= mov;
    dp -= mov;
    std::vector<T> disp_mag2_vec{ dt0.squaredNorm(), dt1.squaredNorm(), dt2.squaredNorm() };

	T sd = *std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end());
    T max_disp_mag = dp.norm() + std::sqrt(*std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end()));
    if (max_disp_mag == 0)
        return 1.0;

    T dist2_cur = internal::point_triangle_distance_unclassified(p, t0, t1, t2);

    T dist_cur = std::sqrt(dist2_cur);
    T gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);
    T toc = 0.0;
    while (true) {
        T toc_lower_bound = (1 - eta) * (dist2_cur - thickness * thickness) / ((dist_cur + thickness) * max_disp_mag);
		//printf("cpu: %f, %f, %f, %f, %f, %f\n", toc_lower_bound, dist2_cur, eta, thickness, dist_cur, max_disp_mag);
        p += toc_lower_bound * dp;
        t0 += toc_lower_bound * dt0;
        t1 += toc_lower_bound * dt1;
        t2 += toc_lower_bound * dt2;
        dist2_cur = internal::point_triangle_distance_unclassified(p, t0, t1, t2);
        dist_cur = std::sqrt(dist2_cur);
        if (toc && ((dist2_cur - thickness * thickness) / (dist_cur + thickness) < gap)) {
            break;
        }

        toc += toc_lower_bound;
        if (toc > 1.0) {
            return 1.0;
        }
    }
    return toc;
}

template <class T>
T point_triangle_ccd2(
	int tid,
	const Vector<T, 3>& _p,
	const Vector<T, 3>& _t0,
	const Vector<T, 3>& _t1,
	const Vector<T, 3>& _t2,
	const Vector<T, 3>& _dp,
	const Vector<T, 3>& _dt0,
	const Vector<T, 3>& _dt1,
	const Vector<T, 3>& _dt2,
	T eta, T thickness)
{
	Vector<T, 3> p = _p, t0 = _t0, t1 = _t1, t2 = _t2, dp = _dp, dt0 = _dt0, dt1 = _dt1, dt2 = _dt2;
	Vector<T, 3> mov = (dt0 + dt1 + dt2 + dp) / 4;
	dt0 -= mov;
	dt1 -= mov;
	dt2 -= mov;
	dp -= mov;
	std::vector<T> disp_mag2_vec{ dt0.squaredNorm(), dt1.squaredNorm(), dt2.squaredNorm() };

	T sd = *std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end());
	T max_disp_mag = dp.norm() + std::sqrt(*std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end()));



	if (max_disp_mag == 0)
		return 1.0;

	T dist2_cur = internal::point_triangle_distance_unclassified(p, t0, t1, t2);

	T dist_cur = std::sqrt(dist2_cur);
	T gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);
	T toc = 0.0;

	//if (tid == 0)
	//{
	//	printf("cpu %f, %f, %f, %f\n", max_disp_mag, dist2_cur, dist_cur, gap);
	//}

	while (true) {
		T toc_lower_bound = (1 - eta) * (dist2_cur - thickness * thickness) / ((dist_cur + thickness) * max_disp_mag);
		//printf("cpu: %f, %f, %f, %f, %f, %f\n", toc_lower_bound, dist2_cur, eta, thickness, dist_cur, max_disp_mag);
		p += toc_lower_bound * dp;
		t0 += toc_lower_bound * dt0;
		t1 += toc_lower_bound * dt1;
		t2 += toc_lower_bound * dt2;
		dist2_cur = internal::point_triangle_distance_unclassified(p, t0, t1, t2);
		dist_cur = std::sqrt(dist2_cur);
		if (toc && ((dist2_cur - thickness * thickness) / (dist_cur + thickness) < gap)) {
			break;
		}

		toc += toc_lower_bound;
		if (toc > 1.0) {
			return 1.0;
		}
	}
	return toc;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_3D
template bool point_triangle_cd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float);
template bool point_triangle_ccd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float);
template float point_triangle_ccd(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float, float);
template float point_triangle_ccd2(int tid, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float, float);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_3D
template bool point_triangle_cd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double);
template bool point_triangle_ccd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double);
template double point_triangle_ccd(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double, double);
template double point_triangle_ccd2(int tid, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double, double);
#endif
#endif
#endif

} // namespace Bow::Geometry::IPC