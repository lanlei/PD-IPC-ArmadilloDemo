#include "PointPointDistance.h"

namespace Bow::Geometry::IPC {
template <class T, int dim>
BOW_INLINE T point_point_distance(const Vector<T, dim>& a, const Vector<T, dim>& b)
{
    return (a - b).squaredNorm();
}

template <class T, int dim>
BOW_INLINE void point_point_distance_gradient(const Vector<T, dim>& a, const Vector<T, dim>& b,
    Vector<T, dim * 2>& grad)
{
    grad.template segment<dim>(0) = 2.0 * (a - b);
    grad.template segment<dim>(dim) = -grad.template segment<dim>(0);
}

template <class T, int dim>
BOW_INLINE void point_point_distance_hessian(const Vector<T, dim>& a, const Vector<T, dim>& b,
    Matrix<T, dim * 2, dim * 2>& hess)
{
    hess.setZero();
    hess.diagonal().setConstant(2.0);
    if constexpr (dim == 2) {
        hess(0, 2) = hess(1, 3) = hess(2, 0) = hess(3, 1) = -2.0;
    }
    else {
        hess(0, 3) = hess(1, 4) = hess(2, 5) = hess(3, 0) = hess(4, 1) = hess(5, 2) = -2.0;
    }
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template float point_point_distance(const Vector<float, 2>&, const Vector<float, 2>&);
template void point_point_distance_gradient(const Vector<float, 2>&, const Vector<float, 2>&, Vector<float, 4>&);
template void point_point_distance_hessian(const Vector<float, 2>&, const Vector<float, 2>&, Matrix<float, 4, 4>&);
#endif
#ifdef BOW_COMPILE_3D
template float point_point_distance(const Vector<float, 3>&, const Vector<float, 3>&);
template void point_point_distance_gradient(const Vector<float, 3>&, const Vector<float, 3>&, Vector<float, 6>&);
template void point_point_distance_hessian(const Vector<float, 3>&, const Vector<float, 3>&, Matrix<float, 6, 6>&);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template double point_point_distance(const Vector<double, 2>&, const Vector<double, 2>&);
template void point_point_distance_gradient(const Vector<double, 2>&, const Vector<double, 2>&, Vector<double, 4>&);
template void point_point_distance_hessian(const Vector<double, 2>&, const Vector<double, 2>&, Matrix<double, 4, 4>&);
#endif
#ifdef BOW_COMPILE_3D
template double point_point_distance(const Vector<double, 3>&, const Vector<double, 3>&);
template void point_point_distance_gradient(const Vector<double, 3>&, const Vector<double, 3>&, Vector<double, 6>&);
template void point_point_distance_hessian(const Vector<double, 3>&, const Vector<double, 3>&, Matrix<double, 6, 6>&);
#endif
#endif
#endif

} // namespace Bow::Geometry::IPC