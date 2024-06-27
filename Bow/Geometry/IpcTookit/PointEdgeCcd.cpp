#include "PointEdgeCcd.h"

namespace Bow::Geometry::IPC {

namespace internal {
template <class T, int dim>
bool check_overlap(
    const Vector<T, dim>& x0,
    const Vector<T, dim>& x1,
    const Vector<T, dim>& x2,
    const Vector<T, dim>& d0,
    const Vector<T, dim>& d1,
    const Vector<T, dim>& d2,
    T root)
{
    Vector<T, dim> p0 = x0 + d0 * root;
    Vector<T, dim> e0 = x1 + d1 * root;
    Vector<T, dim> e1 = x2 + d2 * root;
    Vector<T, dim> e = e1 - e0;
    T ratio = e.dot(p0 - e0) / e.squaredNorm();
    return 0 <= ratio && ratio <= 1;
}
} // namespace internal

template <class T, int dim>
bool point_edge_cd_broadphase(
    const Vector<T, dim>& x0,
    const Vector<T, dim>& x1,
    const Vector<T, dim>& x2,
    T dist)
{
    const Array<T, dim, 1> max_e = x1.array().max(x2.array());
    const Array<T, dim, 1> min_e = x1.array().min(x2.array());
    if ((x0.array() - max_e > dist).any() || (min_e - x0.array() > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T, int dim>
bool point_edge_ccd_broadphase(
    const Vector<T, dim>& p,
    const Vector<T, dim>& e0,
    const Vector<T, dim>& e1,
    const Vector<T, dim>& dp,
    const Vector<T, dim>& de0,
    const Vector<T, dim>& de1,
    const T dist)
{
    const Array<T, dim, 1> max_p = p.array().max((p + dp).array());
    const Array<T, dim, 1> min_p = p.array().min((p + dp).array());
    const Array<T, dim, 1> max_e = e0.array().max(e1.array()).max((e0 + de0).array()).max((e1 + de1).array());
    const Array<T, dim, 1> min_e = e0.array().min(e1.array()).min((e0 + de0).array()).min((e1 + de1).array());
    // check overlap of bounding box
    if ((min_p - max_e > dist).any() || (min_e - max_p > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
T point_edge_ccd(
    const Vector<T, 2>& x0,
    const Vector<T, 2>& x1,
    const Vector<T, 2>& x2,
    const Vector<T, 2>& d0,
    const Vector<T, 2>& d1,
    const Vector<T, 2>& d2,
    T eta)
{
    T toc = 1;
    T a = d0[0] * (d2[1] - d1[1]) + d0[1] * (d1[0] - d2[0]) + d2[0] * d1[1] - d2[1] * d1[0];
    T b = x0[0] * (d2[1] - d1[1]) + d0[0] * (x2[1] - x1[1]) + d0[1] * (x1[0] - x2[0]) + x0[1] * (d1[0] - d2[0]) + d1[1] * x2[0] + d2[0] * x1[1] - d1[0] * x2[1] - d2[1] * x1[0];
    T c = x0[0] * (x2[1] - x1[1]) + x0[1] * (x1[0] - x2[0]) + x2[0] * x1[1] - x2[1] * x1[0];

    if (a == 0) {
        if (b == 0) {
            // parallel motion, only need to handle colinear case
            if (c == 0) {
                // colinear PP CCD
                if ((x0 - x1).dot(d0 - d1) < 0) {
                    T root = std::sqrt((x0 - x1).squaredNorm() / (d0 - d1).squaredNorm());
                    if (root > 0 && root <= 1)
                        toc = std::min(toc, root * (1 - eta));
                }
                if ((x0 - x2).dot(d0 - d2) < 0) {
                    T root = std::sqrt((x0 - x2).squaredNorm() / (d0 - d2).squaredNorm());
                    if (root > 0 && root <= 1)
                        toc = std::min(toc, root * (1 - eta));
                }
            }
        }
        else {
            T root = -c / b;
            if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                toc = std::min(toc, root * (1 - eta));
        }
    }
    else {
        T delta = b * b - 4 * a * c;
        if (delta == 0) {
            T root = -b / (2 * a);
            if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                toc = std::min(toc, root * (1 - eta));
        }
        else if (delta > 0) {
            // accurate expression differs in b's sign
            if (b > 0) {
                T root = (-b - std::sqrt(delta)) / (2 * a);
                if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                    toc = std::min(toc, root * (1 - eta));
                root = 2 * c / (-b - std::sqrt(delta));
                if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                    toc = std::min(toc, root * (1 - eta));
            }
            else {
                T root = 2 * c / (-b + std::sqrt(delta));
                if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                    toc = std::min(toc, root * (1 - eta));
                root = (-b + std::sqrt(delta)) / (2 * a);
                if (root > 0 && root <= 1 && internal::check_overlap(x0, x1, x2, d0, d1, d2, root))
                    toc = std::min(toc, root * (1 - eta));
            }
        }
    }
    return toc;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template bool point_edge_cd_broadphase(const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, float dist);
template bool point_edge_ccd_broadphase(const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const float);
template float point_edge_ccd(const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, const Vector<float, 2>&, float);
#endif
#ifdef BOW_COMPILE_3D
template bool point_edge_cd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, float dist);
bool point_edge_ccd_broadphase(const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const Vector<float, 3>&, const float);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template bool point_edge_cd_broadphase(const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, double dist);
template bool point_edge_ccd_broadphase(const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const double);
template double point_edge_ccd(const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, const Vector<double, 2>&, double);
#endif
#ifdef BOW_COMPILE_3D
template bool point_edge_cd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, double dist);
template bool point_edge_ccd_broadphase(const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const Vector<double, 3>&, const double);
#endif
#endif
#endif

} // namespace Bow::Geometry::IPC