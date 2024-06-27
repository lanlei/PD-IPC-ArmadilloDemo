#ifndef GIVENS_ROTATION_H
#define GIVENS_ROTATION_H

#include <Eigen/Dense>
#include <Bow/Macros.h>

namespace Bow {
namespace Math {
template <class T>
class GivensRotation {
public:
    int m_rowi;
    int m_rowk;
    T m_c;
    T m_s;

    GivensRotation(int rowi, int rowk)
        : m_rowi(rowi), m_rowk(rowk), m_c(1), m_s(0) {}
    GivensRotation(T a, T b, int rowi, int rowk)
        : m_rowi(rowi), m_rowk(rowk) { compute(a, b); }
    ~GivensRotation() {}
    void reset()
    {
        m_c = 1;
        m_s = 0;
    }
    void transpose() { m_s = -m_s; }
    /**
        Compute c and s from a and b so that
        ( c -s ) ( a )  =  ( * )
        ( s  c ) ( b )     ( 0 )
    */
    BOW_INLINE void compute(T a, T b);
    /**
        Compute c and s from a and b so that
        ( c -s ) ( a )  =  ( 0 )
        ( s  c ) ( b )     ( * )
    */
    BOW_INLINE void compute_unconventional(T a, T b);
    /**
        Fill the R with the entries of this rotation
    */
    template <class DerivedR>
    BOW_INLINE void fill(Eigen::MatrixBase<DerivedR>& R) const;
    /**
        This function does something like Q^T A -> A
        [ c -s  0 ]
        [ s  c  0 ] A -> A
        [ 0  0  1 ]
        It only affects row i and row k of A.
    */
    template <class DerivedA>
    BOW_INLINE void row_rotation(Eigen::MatrixBase<DerivedA>& A) const;
    /**
        This function does something like A Q -> A
           [ c  s  0 ]
        A  [-s  c  0 ]  -> A
           [ 0  0  1 ]
        It only affects column i and column k of A.
    */
    template <class DerivedA>
    BOW_INLINE void column_rotation(Eigen::MatrixBase<DerivedA>& A) const;
    /**
        Multiply givens must be for same row and column
    **/
    BOW_INLINE void operator*=(const GivensRotation<T>& A);
    /**
        Multiply givens must be for same row and column
    **/
    BOW_INLINE GivensRotation<T> operator*(const GivensRotation<T>& A) const;
};
}
} // namespace Bow::Math

#ifndef BOW_STATIC_LIBRARY
#include "GivensRotation.cpp"
#endif

#endif