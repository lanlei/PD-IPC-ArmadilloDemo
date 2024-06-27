#include "EdgeEdgeCcd.h"
#include "EdgeEdgeDistance.h"
#include "PointPointDistance.h"
#include "PointEdgeDistance.h"
#include "DistanceType.h"

namespace Bow::Geometry::IPC {

namespace internal {
template <class T>
T edge_edge_distance_unclassified(
    const Vector<T, 3>& ea0,
    const Vector<T, 3>& ea1,
    const Vector<T, 3>& eb0,
    const Vector<T, 3>& eb1)
{
    switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
    case 0:
        return point_point_distance(ea0, eb0);
    case 1:
        return point_point_distance(ea0, eb1);
    case 2:
        return point_edge_distance(ea0, eb0, eb1);
    case 3:
        return point_point_distance(ea1, eb0);
    case 4:
        return point_point_distance(ea1, eb1);
    case 5:
        return point_edge_distance(ea1, eb0, eb1);
    case 6:
        return point_edge_distance(eb0, ea0, ea1);
    case 7:
        return point_edge_distance(eb1, ea0, ea1);
    case 8:
        return edge_edge_distance(ea0, ea1, eb0, eb1);
    default:
        return std::numeric_limits<T>::max();
    }
}
} // namespace internal

template <class T>
BOW_INLINE bool edge_edge_cd_broadphase(
    const Vector<T, 3>& ea0,
    const Vector<T, 3>& ea1,
    const Vector<T, 3>& eb0,
    const Vector<T, 3>& eb1,
    T dist)
{
    auto max_a = ea0.array().max(ea1.array());
    auto min_a = ea0.array().min(ea1.array());
    auto max_b = eb0.array().max(eb1.array());
    auto min_b = eb0.array().min(eb1.array());
    if ((min_a - max_b > dist).any() || (min_b - max_a > dist).any())
        return false;
    else
        return true;
}

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
    T dist)
{
    auto max_a = ea0.array().max(ea1.array()).max((ea0 + dea0).array()).max((ea1 + dea1).array());
    auto min_a = ea0.array().min(ea1.array()).min((ea0 + dea0).array()).min((ea1 + dea1).array());
    auto max_b = eb0.array().max(eb1.array()).max((eb0 + deb0).array()).max((eb1 + deb1).array());
    auto min_b = eb0.array().min(eb1.array()).min((eb0 + deb0).array()).min((eb1 + deb1).array());
    if ((min_a - max_b > dist).any() || (min_b - max_a > dist).any())
        return false;
    else
        return true;
}

template <class T>
BOW_INLINE T edge_edge_ccd(
    const Vector<T, 3>& _ea0,
    const Vector<T, 3>& _ea1,
    const Vector<T, 3>& _eb0,
    const Vector<T, 3>& _eb1,
    const Vector<T, 3>& _dea0,
    const Vector<T, 3>& _dea1,
    const Vector<T, 3>& _deb0,
    const Vector<T, 3>& _deb1,
    T eta, T thickness)
{
    Vector<T, 3> ea0 = _ea0, ea1 = _ea1, eb0 = _eb0, eb1 = _eb1, dea0 = _dea0, dea1 = _dea1, deb0 = _deb0, deb1 = _deb1;
    Vector<T, 3> mov = (dea0 + dea1 + deb0 + deb1) / 4;
    dea0 -= mov;
    dea1 -= mov;
    deb0 -= mov;
    deb1 -= mov;
    T max_disp_mag = std::sqrt(std::max(dea0.squaredNorm(), dea1.squaredNorm())) + std::sqrt(std::max(deb0.squaredNorm(), deb1.squaredNorm()));
    if (max_disp_mag == 0)
        return 1.0;

    T dist2_cur = internal::edge_edge_distance_unclassified(ea0, ea1, eb0, eb1);

    T dFunc = dist2_cur - thickness * thickness;
    if (dFunc <= 0) {
        // since we ensured other place that all dist smaller than dHat are positive,
        // this must be some far away nearly parallel edges
        std::vector<T> dists{ (ea0 - eb0).squaredNorm(), (ea0 - eb1).squaredNorm(), (ea1 - eb0).squaredNorm(), (ea1 - eb1).squaredNorm() };
        dist2_cur = *std::min_element(dists.begin(), dists.end());
        dFunc = dist2_cur - thickness * thickness;
    }
    T dist_cur = std::sqrt(dist2_cur);
    T gap = eta * dFunc / (dist_cur + thickness);
    T toc = 0.0;

	int i = 0;
    while (true) {
        T toc_lower_bound = (1 - eta) * dFunc / ((dist_cur + thickness) * max_disp_mag);
        ea0 += toc_lower_bound * dea0;
        ea1 += toc_lower_bound * dea1;
        eb0 += toc_lower_bound * deb0;
        eb1 += toc_lower_bound * deb1;
        dist2_cur = internal::edge_edge_distance_unclassified(ea0, ea1, eb0, eb1);
        dFunc = dist2_cur - thickness * thickness;
        if (dFunc <= 0) {
            // since we ensured other place that all dist smaller than dHat are positive,
            // this must be some far away nearly parallel edges
            std::vector<T> dists{ (ea0 - eb0).squaredNorm(), (ea0 - eb1).squaredNorm(), (ea1 - eb0).squaredNorm(), (ea1 - eb1).squaredNorm() };
            dist2_cur = *std::min_element(dists.begin(), dists.end());
            dFunc = dist2_cur - thickness * thickness;
        }
        dist_cur = std::sqrt(dist2_cur);
        if (toc && (dFunc / (dist_cur + thickness) < gap)) {
            break;
        }
        toc += toc_lower_bound;
        if (toc > 1.0)
            return 1.0;
    }
    return toc;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_3D
template bool edge_edge_cd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float);
template bool edge_edge_ccd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float);
template float edge_edge_ccd(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float, float);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_3D
template bool edge_edge_cd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double);
template bool edge_edge_ccd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double);
template double edge_edge_ccd(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double, double);
#endif
#endif
#endif

} // namespace Bow::Geometry::IPC