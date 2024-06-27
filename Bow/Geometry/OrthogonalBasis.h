#ifndef ORTHONORMAL_BASIS_H
#include <Bow/Types.h>

namespace Bow::Geometry {

template <class T>
Matrix<T, 3, 3> extend_basis(const Vector<T, 3>& v1, const Vector<T, 3>& v2)
{
    Matrix<T, 3, 3> basis;
    basis.col(0) = v1;
    basis.col(1) = v2;
    basis.col(2) = v1.cross(v2);
    return basis;
}

template <class T, int dim>
Matrix<T, dim, dim> extend_basis(const Vector<T, dim>& v1)
{
    Matrix<T, dim, dim> basis;
    basis.col(0) = v1;
    if constexpr (dim == 2) {
        basis(0, 1) = -basis(1, 0);
        basis(1, 1) = basis(0, 0);
    }
    else {
        if (basis.col(0).dot(Vector<T, dim>(1, 0, 0)) > 0.5) {
            basis.col(1) = basis.col(0).cross(Vector<T, dim>(0, 1, 0));
        }
        else {
            basis.col(1) = basis.col(0).cross(Vector<T, dim>(1, 0, 0));
        }
        basis.col(1).normalize();
        basis.col(2) = basis.col(0).cross(basis.col(1));
    }
    return basis;
}
} // namespace Bow::Geometry

#endif