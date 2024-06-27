#include "GivensRotation.h"
#include <Bow/Types.h>

namespace Bow {
namespace Math {
template <class T>
BOW_INLINE void GivensRotation<T>::compute(T a, T b)
{
    using std::sqrt;
    T d = a * a + b * b;
    m_c = 1;
    m_s = 0;
    T sqrtd = sqrt(d);
    if (sqrtd) {
        T t = 1 / sqrtd;
        m_c = a * t;
        m_s = -b * t;
    }
}

template <class T>
BOW_INLINE void GivensRotation<T>::compute_unconventional(T a, T b)
{
    using std::sqrt;
    T d = a * a + b * b;
    m_c = 0;
    m_s = 1;
    T sqrtd = sqrt(d);
    if (sqrtd) {
        T t = T(1) / sqrtd;
        m_s = a * t;
        m_c = b * t;
    }
}

template <class T>
template <class DerivedR>
BOW_INLINE void GivensRotation<T>::fill(Eigen::MatrixBase<DerivedR>& R) const
{
    R.setIdentity();
    R(m_rowi, m_rowi) = m_c;
    R(m_rowk, m_rowi) = -m_s;
    R(m_rowi, m_rowk) = m_s;
    R(m_rowk, m_rowk) = m_c;
}

template <class T>
template <class DerivedA>
BOW_INLINE void GivensRotation<T>::row_rotation(Eigen::MatrixBase<DerivedA>& A) const
{
    for (int j = 0; j < DerivedA::RowsAtCompileTime; j++) {
        T tau1 = A(m_rowi, j);
        T tau2 = A(m_rowk, j);
        A(m_rowi, j) = m_c * tau1 - m_s * tau2;
        A(m_rowk, j) = m_s * tau1 + m_c * tau2;
    }
    //not type safe :/
}

template <class T>
template <class DerivedA>
BOW_INLINE void GivensRotation<T>::column_rotation(Eigen::MatrixBase<DerivedA>& A) const
{
    for (int j = 0; j < DerivedA::RowsAtCompileTime; j++) {
        T tau1 = A(j, m_rowi);
        T tau2 = A(j, m_rowk);
        A(j, m_rowi) = m_c * tau1 - m_s * tau2;
        A(j, m_rowk) = m_s * tau1 + m_c * tau2;
    }
    //not type safe :/
}

template <class T>
BOW_INLINE void GivensRotation<T>::operator*=(const GivensRotation<T>& A)
{
    T new_c = m_c * A.m_c - m_s * A.m_s;
    T new_s = m_s * A.m_c + m_c * A.m_s;
    m_c = new_c;
    m_s = new_s;
}

template <class T>
BOW_INLINE GivensRotation<T> GivensRotation<T>::operator*(const GivensRotation<T>& A) const
{
    GivensRotation<T> r(*this);
    r *= A;
    return r;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template class GivensRotation<float>;
#ifdef BOW_COMPILE_2D
template void GivensRotation<float>::fill(Eigen::MatrixBase<Matrix<float, 2, 2>>& R) const;
template void GivensRotation<float>::row_rotation(Eigen::MatrixBase<Matrix<float, 2, 2>>& A) const;
template void GivensRotation<float>::column_rotation(Eigen::MatrixBase<Matrix<float, 2, 2>>& A) const;
#endif
#ifdef BOW_COMPILE_3D
template void GivensRotation<float>::fill(Eigen::MatrixBase<Matrix<float, 3, 3>>& R) const;
template void GivensRotation<float>::row_rotation(Eigen::MatrixBase<Matrix<float, 3, 3>>& A) const;
template void GivensRotation<float>::column_rotation(Eigen::MatrixBase<Matrix<float, 3, 3>>& A) const;
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
template class GivensRotation<double>;
#ifdef BOW_COMPILE_2D
template void GivensRotation<double>::fill(Eigen::MatrixBase<Matrix<double, 2, 2>>& R) const;
template void GivensRotation<double>::row_rotation(Eigen::MatrixBase<Matrix<double, 2, 2>>& A) const;
template void GivensRotation<double>::column_rotation(Eigen::MatrixBase<Matrix<double, 2, 2>>& A) const;
#endif
#ifdef BOW_COMPILE_3D
template void GivensRotation<double>::fill(Eigen::MatrixBase<Matrix<double, 3, 3>>& R) const;
template void GivensRotation<double>::row_rotation(Eigen::MatrixBase<Matrix<double, 3, 3>>& A) const;
template void GivensRotation<double>::column_rotation(Eigen::MatrixBase<Matrix<double, 3, 3>>& A) const;
#endif
#endif
#endif
}
} // namespace Bow::Math
