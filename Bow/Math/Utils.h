#ifndef MATH_UTILS
#define MATH_UTILS

#include <Bow/Macros.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Bow {
namespace Math {
template <class DerivedA, class DerivedCoA>
BOW_INLINE void cofactor(const Eigen::MatrixBase<DerivedA>& A, Eigen::MatrixBase<DerivedCoA>& coA);
template <class DerivedA>
BOW_INLINE void make_pd(Eigen::MatrixBase<DerivedA>& symA);
/**
 * \brief CSR arraies to Eigen sparse matrix
 */
template <class Derived>
BOW_INLINE void sparse_from_csr(const std::vector<typename Derived::StorageIndex>& outer_indices,
    const std::vector<typename Derived::StorageIndex>& inner_indices,
    const std::vector<typename Derived::Scalar>& val,
    const typename Derived::StorageIndex nrows,
    const typename Derived::StorageIndex ncols,
    Eigen::SparseMatrixBase<Derived>& A);
/** 
 * \brief Select cols (for ColMajor)
*/
template <class T, class StorageIndex>
BOW_INLINE void select_cols(const Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& A,
    const std::vector<StorageIndex>& selection,
    Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& A_cols);

template <class T, class StorageIndex>
BOW_INLINE void shift(const Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& out);

inline double get_smallest_positive_real_quad_root(double a, double b, double c, double tol = 1e-6)
{
    // return negative value if no positive real root is found
    using std::abs;
    using std::sqrt;
    double t;
    if (abs(a) <= tol) {
        if (abs(b) <= tol) // f(x) = c > 0 for all x
            t = -1;
        else
            t = -c / b;
    }
    else {
        double desc = b * b - 4 * a * c;
        if (desc > 0) {
            t = (-b - sqrt(desc)) / (2 * a);
            if (t < 0)
                t = (-b + sqrt(desc)) / (2 * a);
        }
        else // desv<0 ==> imag
            t = -1;
    }
    return t;
}

inline double get_smallest_positive_real_cubic_root(double a, double b, double c, double d, double tol = 1e-6)
{
    // return negative value if no positive real root is found
    using std::abs;
    using std::complex;
    using std::pow;
    using std::sqrt;
    double t = -1;
    if (abs(a) <= tol)
        t = get_smallest_positive_real_quad_root(b, c, d, tol);
    else {
        complex<double> i(0, 1);
        complex<double> delta0(b * b - 3 * a * c, 0);
        complex<double> delta1(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 0);
        complex<double> C = pow((delta1 + sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
        if (abs(C) < tol)
            C = pow((delta1 - sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
        complex<double> u2 = (-1.0 + sqrt(3.0) * i) / 2.0;
        complex<double> u3 = (-1.0 - sqrt(3.0) * i) / 2.0;
        complex<double> t1 = (b + C + delta0 / C) / (-3.0 * a);
        complex<double> t2 = (b + u2 * C + delta0 / (u2 * C)) / (-3.0 * a);
        complex<double> t3 = (b + u3 * C + delta0 / (u3 * C)) / (-3.0 * a);
        if ((abs(imag(t1)) < tol) && (real(t1) > 0))
            t = real(t1);
        if ((abs(imag(t2)) < tol) && (real(t2) > 0) && ((real(t2) < t) || (t < 0)))
            t = real(t2);
        if ((abs(imag(t3)) < tol) && (real(t3) > 0) && ((real(t3) < t) || (t < 0)))
            t = real(t3);
    }
    return t;
}
}
} // namespace Bow::Math

#ifndef BOW_STATIC_LIBRARY
#include "Utils.cpp"
#endif

#endif