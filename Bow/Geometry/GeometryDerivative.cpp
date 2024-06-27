#include "GeometryDerivative.h"

namespace Bow {
namespace Geometry {
namespace internal {
template <class T>
inline Vector<T, 2> area_weighted_normal(const Vector<T, 2>& x0, const Vector<T, 2>& x1)
{
    Vector<T, 2> e = x1 - x0;
    return Vector<T, 2>(-e[1], e[0]);
}
template <class T>
inline Vector<T, 3> area_weighted_normal(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2)
{
    Vector<T, 3> e1 = x1 - x0;
    Vector<T, 3> e2 = x2 - x0;
    return 0.5 * e1.cross(e2);
}
} // namespace internal

template <class T>
T simplex_volume(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2)
{
    Matrix<T, 2, 2> local_basis;
    local_basis.col(0) = x1 - x0;
    local_basis.col(1) = x2 - x0;
    return local_basis.determinant() / T(2);
}

template <class T, class Derived>
void simplex_volume_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, Eigen::MatrixBase<Derived>& grad)
{
    EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
    grad.setZero();
    grad.template segment<2>(0) = internal::area_weighted_normal(x1, x2);
    grad.template segment<2>(2) = internal::area_weighted_normal(x2, x0);
    grad.template segment<2>(4) = internal::area_weighted_normal(x0, x1);
    grad *= T(0.5);
}

template <class T>
T simplex_volume(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3)
{
    Matrix<T, 3, 3> local_basis;
    local_basis.col(0) = x1 - x0;
    local_basis.col(1) = x2 - x0;
    local_basis.col(2) = x3 - x0;
    return local_basis.determinant() / T(6);
}

template <class T, class Derived>
void simplex_volume_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3, Eigen::MatrixBase<Derived>& grad)
{
    EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
    grad.setZero();
    grad.template segment<3>(0) = internal::area_weighted_normal(x2, x1, x3);
    grad.template segment<3>(3) = internal::area_weighted_normal(x0, x2, x3);
    grad.template segment<3>(6) = internal::area_weighted_normal(x0, x3, x1);
    grad.template segment<3>(9) = internal::area_weighted_normal(x0, x1, x2);
    grad /= T(3);
}

template <class T, int dim>
void center_of_mass(const Field<Vector<T, dim>>& x, const Field<T>& mass, Vector<T, dim>& cm)
{
    Eigen::Map<const Vector<T, Eigen::Dynamic>> mass_vec(mass.data(), mass.size());
    Eigen::Map<const Matrix<T, dim, Eigen::Dynamic>> x_stack(reinterpret_cast<const T*>(x.data()), dim, x.size());
    cm = x_stack * mass_vec / mass_vec.sum();
}

template <class T, class DerivedX, class DerivedY>
void center_of_mass_gradient(const Field<T>& mass, const Eigen::MatrixBase<DerivedX>& dLdcm, Eigen::MatrixBase<DerivedY>& grad)
{
    Eigen::Map<const Vector<T, Eigen::Dynamic>> mass_vec(mass.data(), mass.size());
    const int dim = DerivedX::RowsAtCompileTime;
    Matrix<T, dim, Eigen::Dynamic> grad_mat = dLdcm * mass_vec.transpose();
    grad = Eigen::Map<Vector<T, Eigen::Dynamic>>(reinterpret_cast<T*>(grad_mat.data()), dim * mass.size()) / mass_vec.sum();
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template float simplex_volume(const Vector<float, 2>& x0, const Vector<float, 2>& x1, const Vector<float, 2>& x2);
template float simplex_volume(const Vector<float, 3>& x0, const Vector<float, 3>& x1, const Vector<float, 3>& x2, const Vector<float, 3>& x3);
template void simplex_volume_gradient(const Vector<float, 2>& x0, const Vector<float, 2>& x1, const Vector<float, 2>& x2, Eigen::MatrixBase<Vector<float, 6>>&);
template void simplex_volume_gradient(const Vector<float, 3>& x0, const Vector<float, 3>& x1, const Vector<float, 3>& x2, const Vector<float, 3>& x3, Eigen::MatrixBase<Vector<float, 12>>&);
template void center_of_mass(const Field<Vector<float, 2>>& x, const Field<float>& mass, Vector<float, 2>& cm);
template void center_of_mass(const Field<Vector<float, 3>>& x, const Field<float>& mass, Vector<float, 3>& cm);
template void center_of_mass_gradient(const Field<float>& mass, const Eigen::MatrixBase<Eigen::Vector2f>& dLdcm, Eigen::MatrixBase<Eigen::VectorXf>& grad);
template void center_of_mass_gradient(const Field<float>& mass, const Eigen::MatrixBase<Eigen::Vector3f>& dLdcm, Eigen::MatrixBase<Eigen::VectorXf>& grad);

#endif
#ifdef BOW_COMPILE_DOUBLE
template double simplex_volume(const Vector<double, 2>& x0, const Vector<double, 2>& x1, const Vector<double, 2>& x2);
template double simplex_volume(const Vector<double, 3>& x0, const Vector<double, 3>& x1, const Vector<double, 3>& x2, const Vector<double, 3>& x3);
template void simplex_volume_gradient(const Vector<double, 2>& x0, const Vector<double, 2>& x1, const Vector<double, 2>& x2, Eigen::MatrixBase<Vector<double, 6>>&);
template void simplex_volume_gradient(const Vector<double, 3>& x0, const Vector<double, 3>& x1, const Vector<double, 3>& x2, const Vector<double, 3>& x3, Eigen::MatrixBase<Vector<double, 12>>&);
template void center_of_mass(const Field<Vector<double, 2>>& x, const Field<double>& mass, Vector<double, 2>& cm);
template void center_of_mass(const Field<Vector<double, 3>>& x, const Field<double>& mass, Vector<double, 3>& cm);
template void center_of_mass_gradient(const Field<double>& mass, const Eigen::MatrixBase<Eigen::Vector2d>& dLdcm, Eigen::MatrixBase<Eigen::VectorXd>& grad);
template void center_of_mass_gradient(const Field<double>& mass, const Eigen::MatrixBase<Eigen::Vector3d>& dLdcm, Eigen::MatrixBase<Eigen::VectorXd>& grad);
#endif
#endif

}} // namespace Bow::Geometry