#include <Bow/Types.h>
#include "SVD.h"
#include "GivensRotation.h"
#include "Utils.h"

namespace Bow {
namespace Math {
namespace internal {
template <class DerivedM, class DerivedV>
inline void eigen_values(const Eigen::MatrixBase<DerivedM>& A_Sym, Eigen::MatrixBase<DerivedV>& lambda)
{
    static_assert(std::is_same<typename DerivedM::Scalar, double>::value, "eigen_values only supports double precision");
    static_assert(DerivedV::SizeAtCompileTime == 3, "eigen_values only supports 3 x 3 matrices");
    using T = double;
    using std::max;
    using std::sqrt;
    using std::swap;
    T m = ((T)1 / 3) * (A_Sym(0, 0) + A_Sym(1, 1) + A_Sym(2, 2));
    T a00 = A_Sym(0, 0) - m;
    T a11 = A_Sym(1, 1) - m;
    T a22 = A_Sym(2, 2) - m;
    T a12_sqr = A_Sym(0, 1) * A_Sym(0, 1);
    T a13_sqr = A_Sym(0, 2) * A_Sym(0, 2);
    T a23_sqr = A_Sym(1, 2) * A_Sym(1, 2);
    T p = ((T)1 / 6) * (a00 * a00 + a11 * a11 + a22 * a22 + 2 * (a12_sqr + a13_sqr + a23_sqr));
    T q = (T).5 * (a00 * (a11 * a22 - a23_sqr) - a11 * a13_sqr - a22 * a12_sqr) + A_Sym(0, 1) * A_Sym(0, 2) * A_Sym(1, 2);
    T sqrt_p = sqrt(p);
    T disc = p * p * p - q * q;
    T phi = ((T)1 / 3) * atan2(sqrt(max((T)0, disc)), q);
    T c = cos(phi), s = sin(phi);
    T sqrt_p_cos = sqrt_p * c;
    T root_three_sqrt_p_sin = sqrt((T)3) * sqrt_p * s;

    lambda(0) = m + 2 * sqrt_p_cos;
    lambda(1) = m - sqrt_p_cos - root_three_sqrt_p_sin;
    lambda(2) = m - sqrt_p_cos + root_three_sqrt_p_sin;

    if (lambda(0) < lambda(1))
        swap(lambda(0), lambda(1));
    if (lambda(1) < lambda(2))
        swap(lambda(1), lambda(2));
    if (lambda(0) < lambda(1))
        swap(lambda(0), lambda(1));
}

template <class DerivedM, class DerivedV>
inline void eigen_vectors(const Eigen::MatrixBase<DerivedM>& A_Sym, const Eigen::MatrixBase<DerivedV>& lambda, Eigen::MatrixBase<DerivedM>& V)
{
    static_assert(std::is_same<typename DerivedM::Scalar, double>::value, "eigen_vectors only supports double precision");
    static_assert(DerivedV::SizeAtCompileTime == 3, "eigen_vectors only supports 3 x 3 matrices");
    using T = double;
    using std::sqrt;
    using std::swap;

    bool flipped = false;
    Vector<T, 3> lambda_flip(lambda);
    if (lambda(0) - lambda(1) < lambda(1) - lambda(2)) {
        swap(lambda_flip(0), lambda_flip(2));
        flipped = true;
    }

    // get first eigenvector
    Matrix<T, 3, 3> C1;
    Matrix<T, 3, 3> temp1 = -lambda_flip(0) * DerivedM::Identity() + A_Sym;
    cofactor(temp1, C1);

    typename Vector<T, 3>::Index i;
    T norm2 = C1.colwise().squaredNorm().maxCoeff(&i);

    Vector<T, 3> v1;
    if (norm2 != 0) {
        T one_over_sqrt = (T)1 / sqrt(norm2);
        v1 = C1.col(i) * one_over_sqrt;
    }
    else
        v1 << 1, 0, 0;
    // form basis for orthogonal complement to v1, and reduce A to this space
    Vector<T, 3> v1_orthogonal = v1.unitOrthogonal(); // 6m+2a+1d+1s (tweak: 5m+1a+1d+1s)
    Eigen::Matrix<T, 3, 2, 0, 3, 2> other_v;
    other_v.col(0) = v1_orthogonal;
    other_v.col(1) = v1.cross(v1_orthogonal); // 6m+3a (tweak: 4m+1a)
    Matrix<T, 2, 2> A_reduced = other_v.transpose() * A_Sym * other_v; // 21m+12a (tweak: 18m+9a)

    // find third eigenvector from A_reduced, and fill in second via cross product
    Matrix<T, 2, 2> C3;
    Matrix<T, 2, 2> temp3 = -lambda_flip(2) * Matrix<T, 2, 2>::Identity() + A_reduced;
    cofactor(temp3, C3);
    Vector<T, 2>::Index j;
    norm2 = C3.colwise().squaredNorm().maxCoeff(&j); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s
    Vector<T, 3> v3;
    if (norm2 != 0) {
        T one_over_sqrt = (T)1 / sqrt(norm2);
        v3 = other_v * C3.col(j) * one_over_sqrt;
    }
    else
        v3 = other_v.col(0);

    Vector<T, 3> v2 = v3.cross(v1); // 6m+3a

    // finish
    if (flipped) {
        V.col(0) = v3;
        V.col(1) = v2;
        V.col(2) = -v1;
    }
    else {
        V.col(0) = v1;
        V.col(1) = v2;
        V.col(2) = v3;
    }
}

template <class DerivedA, class DerivedU, class DerivedSigma, class DerivedV>
inline void _svd(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedU>& U,
    Eigen::MatrixBase<DerivedSigma>& sigma,
    Eigen::MatrixBase<DerivedV>& V)
{
    static_assert(std::is_same<typename DerivedA::Scalar, double>::value, "_svd only supports double precision");
    static_assert(DerivedA::RowsAtCompileTime == 3, "_svd only supports 3 x 3 matrices");
    using T = double;
    static const int dim = 3;
    Vector<T, dim> lambda;
    Matrix<T, dim, dim> A_sym = A.transpose() * A;
    eigen_values(A_sym, lambda);
    eigen_vectors(A_sym, lambda, V);
    // compute singular values
    if (lambda(2) < 0)
        lambda = (lambda.array() >= (T)0).select(lambda, (T)0);
    sigma = lambda.array().sqrt();
    if (A.determinant() < 0)
        sigma(2) = -sigma(2);

    // compute singular vectors
    U.col(0) = A * V.col(0);
    T norm = U.col(0).norm();
    if (norm != 0) {
        T one_over_norm = (T)1 / norm;
        U.col(0) = U.col(0) * one_over_norm;
    }
    else
        U.col(0) << 1, 0, 0;
    Vector<T, 3> v1_orthogonal = U.col(0).unitOrthogonal();
    Eigen::Matrix<T, 3, 2, 0, 3, 2> other_v;
    other_v.col(0) = v1_orthogonal;
    other_v.col(1) = U.col(0).cross(v1_orthogonal);
    Vector<T, 2> w = other_v.transpose() * A * V.col(1);
    norm = w.norm();
    if (norm != 0) {
        T one_over_norm = (T)1 / norm;
        w = w * one_over_norm;
    }
    else
        w << 1, 0;
    U.col(1) = other_v * w;
    U.col(2) = U.col(0).cross(U.col(1));
}
} // namespace internal

template <class DerivedA, class DerivedU, class DerivedSigma, class DerivedV>
BOW_INLINE void svd(
    const Eigen::MatrixBase<DerivedA>& A,
    Eigen::MatrixBase<DerivedU>& U,
    Eigen::MatrixBase<DerivedSigma>& sigma,
    Eigen::MatrixBase<DerivedV>& V)
{
    using T = typename DerivedA::Scalar;
    const int dim = DerivedA::RowsAtCompileTime;
    if constexpr (dim == 1) {
        sigma = A;
        U.array() = 1;
        V.array() = 1;
    }
    else if constexpr (dim == 2) {
        using std::sqrt;
        GivensRotation<T> gv(0, 1);
        GivensRotation<T> gu(0, 1);
        // 2d polar decomposition
        Vector<T, dim> xx(A(0, 0) + A(1, 1), A(1, 0) - A(0, 1));
        T denominator = xx.norm();
        gu.m_c = (T)1;
        gu.m_s = (T)0;
        if (denominator != 0) {
            gu.m_c = xx(0) / denominator;
            gu.m_s = -xx(1) / denominator;
        }
        Matrix<T, dim, dim> S(A);
        gu.row_rotation(S);

        T cosine, sine;
        T x = S(0, 0);
        T y = S(0, 1);
        T z = S(1, 1);
        T y2 = y * y;
        if (y2 == 0) {
            // S is already diagonal
            cosine = 1;
            sine = 0;
            sigma(0) = x;
            sigma(1) = z;
        }
        else {
            T tau = T(0.5) * (x - z);
            T w = sqrt(tau * tau + y2);
            // w > y > 0
            T t;
            if (tau > 0) {
                // tau + w > w > y > 0 ==> division is safe
                t = y / (tau + w);
            }
            else {
                // tau - w < -w < -y < 0 ==> division is safe
                t = y / (tau - w);
            }
            cosine = T(1) / sqrt(t * t + T(1));
            sine = -t * cosine;
            /*
            V = [cosine -sine; sine cosine]
            Sigma = V'SV. Only compute the diagonals for efficiency.
            Also utilize symmetry of S and don't form V yet.
            */
            T c2 = cosine * cosine;
            T csy = 2 * cosine * sine * y;
            T s2 = sine * sine;
            sigma(0) = c2 * x - csy + s2 * z;
            sigma(1) = s2 * x + csy + c2 * z;
        }
        // Sorting
        // Polar already guarantees negative sign is on the small magnitude singular value.
        if (sigma(0) < sigma(1)) {
            std::swap(sigma(0), sigma(1));
            gv.m_c = -sine;
            gv.m_s = cosine;
        }
        else {
            gv.m_c = cosine;
            gv.m_s = sine;
        }
        gu *= gv;
        gu.fill(U);
        gv.fill(V);
    }
    else if constexpr (dim == 3) {
        if constexpr (std::is_same<typename DerivedA::Scalar, double>::value) {
            internal::_svd(A, U, sigma, V);
        }
        else {
            Matrix<double, 3, 3> Ud = U.template cast<double>();
            Vector<double, 3> sd = sigma.template cast<double>();
            Matrix<double, 3, 3> Vd = V.template cast<double>();
            Matrix<double, 3, 3> Ad = A.template cast<double>();
            internal::_svd(Ad, Ud, sd, Vd);
            V = Vd.template cast<T>();
            U = Ud.template cast<T>();
            sigma = sd.template cast<T>();
        }
    }
    else {
        // Not implemented!!
        static_assert(dim <= 3, "svd is only implemented for 1x1, 2x2 and 3x3 matrices.");
    }
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template void svd(
    const Eigen::MatrixBase<Matrix<float, 2, 2>>& A,
    Eigen::MatrixBase<Matrix<float, 2, 2>>& U,
    Eigen::MatrixBase<Vector<float, 2>>& Sigma,
    Eigen::MatrixBase<Matrix<float, 2, 2>>& V);
#endif
#ifdef BOW_COMPILE_3D
template void svd(
    const Eigen::MatrixBase<Matrix<float, 3, 3>>& A,
    Eigen::MatrixBase<Matrix<float, 3, 3>>& U,
    Eigen::MatrixBase<Vector<float, 3>>& Sigma,
    Eigen::MatrixBase<Matrix<float, 3, 3>>& V);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template void svd(
    const Eigen::MatrixBase<Matrix<double, 2, 2>>& A,
    Eigen::MatrixBase<Matrix<double, 2, 2>>& U,
    Eigen::MatrixBase<Vector<double, 2>>& Sigma,
    Eigen::MatrixBase<Matrix<double, 2, 2>>& V);
#endif
#ifdef BOW_COMPILE_3D
template void svd(
    const Eigen::MatrixBase<Matrix<double, 3, 3>>& A,
    Eigen::MatrixBase<Matrix<double, 3, 3>>& U,
    Eigen::MatrixBase<Vector<double, 3>>& Sigma,
    Eigen::MatrixBase<Matrix<double, 3, 3>>& V);
#endif
#endif
#endif
}
} // namespace Bow::Math
