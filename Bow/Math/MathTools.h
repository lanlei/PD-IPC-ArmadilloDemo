#pragma once

#include <cmath>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <Bow/Types.h>
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace MATH_TOOLS {

template <class T>
void addScaledRotationalDifferential(
    const Eigen::MatrixBase<Matrix<T, 1, 1>>& R,
    const Eigen::MatrixBase<Matrix<T, 1, 1>>& S,
    const Eigen::MatrixBase<Matrix<T, 1, 1>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 1, 1>>& M)
{
}
template <class T>
void addScaledRotationalDifferential(
    const Eigen::MatrixBase<Matrix<T, 2, 2>>& R,
    const Eigen::MatrixBase<Matrix<T, 2, 2>>& S,
    const Eigen::MatrixBase<Matrix<T, 2, 2>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 2, 2>>& M)
{
    T trace_s = S.trace();
    BOW_ASSERT_INFO(trace_s != 0, "dR computation encountered division by zero");
    //Majorly exploited but basically derived from "rotationalDifferential'
    T omega = (R(0, 0) * dF(0, 1) + R(1, 0) * dF(1, 1) - R(0, 1) * dF(0, 0) - R(1, 1) * dF(1, 0)) / trace_s * scale;
    M(0, 0) -= omega * R(0, 1);
    M(1, 0) -= omega * R(1, 1);
    M(0, 1) += omega * R(0, 0);
    M(1, 1) += omega * R(1, 0);
}
template <class T>
void addScaledRotationalDifferential(
    const Eigen::MatrixBase<Matrix<T, 3, 3>>& R,
    const Eigen::MatrixBase<Matrix<T, 3, 3>>& S,
    const Eigen::MatrixBase<Matrix<T, 3, 3>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 3, 3>>& M)
{
    using TM = Matrix<T, 3, 3>;

    TM S_hat = -S;
    S_hat.diagonal().array() += S.trace();
    T b = S_hat.determinant();
    BOW_ASSERT_INFO(b != 0, "dRdF computation encountered division by zero");
    TM A = R.transpose() * dF;
    M.noalias() += (scale / b) * R * S_hat * (A - A.transpose()) * S_hat;
}

template <class T>
void addScaledCofactorMatrixDifferential(
    const Eigen::MatrixBase<Matrix<T, 1, 1>>& F,
    const Eigen::MatrixBase<Matrix<T, 1, 1>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 1, 1>>& M)
{
    M += scale * dF;
}
template <class T>
void addScaledCofactorMatrixDifferential(
    const Eigen::MatrixBase<Matrix<T, 2, 2>>& F,
    const Eigen::MatrixBase<Matrix<T, 2, 2>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 2, 2>>& M)
{
    M(0, 0) += scale * dF(1, 1);
    M(1, 0) -= scale * dF(0, 1);
    M(0, 1) -= scale * dF(1, 0);
    M(1, 1) += scale * dF(0, 0);
}
template <class T>
void addScaledCofactorMatrixDifferential(
    const Eigen::MatrixBase<Matrix<T, 3, 3>>& F,
    const Eigen::MatrixBase<Matrix<T, 3, 3>>& dF,
    T scale,
    Eigen::MatrixBase<Matrix<T, 3, 3>>& M)
{
    M(0, 0) += scale * (dF(1, 1) * F(2, 2) + F(1, 1) * dF(2, 2) - dF(2, 1) * F(1, 2) - F(2, 1) * dF(1, 2));
    M(1, 0) += scale * (dF(2, 1) * F(0, 2) + F(2, 1) * dF(0, 2) - dF(0, 1) * F(2, 2) - F(0, 1) * dF(2, 2));
    M(2, 0) += scale * (dF(0, 1) * F(1, 2) + F(0, 1) * dF(1, 2) - dF(1, 1) * F(0, 2) - F(1, 1) * dF(0, 2));
    M(0, 1) += scale * (dF(2, 0) * F(1, 2) + F(2, 0) * dF(1, 2) - dF(1, 0) * F(2, 2) - F(1, 0) * dF(2, 2));
    M(1, 1) += scale * (dF(0, 0) * F(2, 2) + F(0, 0) * dF(2, 2) - dF(2, 0) * F(0, 2) - F(2, 0) * dF(0, 2));
    M(2, 1) += scale * (dF(1, 0) * F(0, 2) + F(1, 0) * dF(0, 2) - dF(0, 0) * F(1, 2) - F(0, 0) * dF(1, 2));
    M(0, 2) += scale * (dF(1, 0) * F(2, 1) + F(1, 0) * dF(2, 1) - dF(2, 0) * F(1, 1) - F(2, 0) * dF(1, 1));
    M(1, 2) += scale * (dF(2, 0) * F(0, 1) + F(2, 0) * dF(0, 1) - dF(0, 0) * F(2, 1) - F(0, 0) * dF(2, 1));
    M(2, 2) += scale * (dF(0, 0) * F(1, 1) + F(0, 0) * dF(1, 1) - dF(1, 0) * F(0, 1) - F(1, 0) * dF(0, 1));
}

/**
   int floor much faster than (int)floor(x)
 */

inline static int int_floor(float x)
{
    int i = (int)x; /* truncate */
    return i - (i > x); /* convert trunc to floor */
}

inline static int int_floor(double x)
{
    int i = (int)x; /* truncate */
    return i - (i > x); /* convert trunc to floor */
}

template <class T>
inline static Vector<int, 2> int_floor(const Vector<T, 2>& v)
{
    return Vector<int, 2>(
        int_floor(v(0)),
        int_floor(v(1)));
}

template <class T>
inline static Vector<int, 3> int_floor(const Vector<T, 3>& v)
{
    return Vector<int, 3>(
        int_floor(v(0)),
        int_floor(v(1)),
        int_floor(v(2)));
}

// Generate pseudorand float between 0 and 1
// only 22 significant bits
inline static float pseudorand_f(int seed)
{
    //shuffle bits
    seed ^= seed >> 16;
    seed *= 0x85ebca6b;
    seed ^= seed >> 13;
    seed *= 0xc2b2ae35;
    seed ^= seed >> 16;
    // keep 22 bits
    float f = seed >> 10;
    return scalbnf(f, -22);
}

template <class Derived>
inline static Vector<float, 2> pseudorand_2f(const Eigen::MatrixBase<Derived>& seed)
{
    // multiply by invertible matrix mod 2^32
    // to ensure each pseudo random seed depends on
    // all 3 coordinates
    return Vector<float, 2>(
        pseudorand_f(12 * seed(0) + 89 * seed(1)),
        pseudorand_f(32 * seed(0) + 35 * seed(1)));
}

template <class Derived>
inline static Vector<float, 3> pseudorand_3f(const Eigen::MatrixBase<Derived>& seed)
{
    // multiply by invertible matrix mod 2^32
    // to ensure each pseudo random seed depends on
    // all 3 coordinates
    return Vector<float, 3>(
        pseudorand_f(12 * seed(0) + 89 * seed(1) + 63 * seed(2)),
        pseudorand_f(32 * seed(0) + 35 * seed(1) + 14 * seed(2)),
        pseudorand_f(51 * seed(0) + 76 * seed(1) + 45 * seed(2)));
}

/**
    \brief Permutation Tensor of rank 3
    \param i first index.
    \param j second index.
    \param k third index.

    return the permutation tensor of index (i,j,k)
*/

inline static int permutationTensor(int i, int j, int k)
{
    int r = i + 3 * j + 9 * k;
    return ((0x200880 >> r) & 1) - ((0x088020 >> r) & 1);
}

/**
  \brief Approximate inverse square root

  A fast approximation to the inverse sqrt
  The relative error is less than  1.5*2^-12
*/
inline float approx_rsqrt(float a)
{
    return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(a)));
}

/**
  \brief Inverse square root
  computed from approx_rsqrt and one newton step
*/
inline float rsqrt(float a)
{
    return (float)1 / std::sqrt(a);

    // float b = approx_rsqrt(a);
    // // Newton step with f(x) = a - 1/x^2
    // b = 0.5f * b * (3.0f - a * (b * b));
    // return b;
}

/**
  \brief Inverse square root
  computed from 1/std::sqrt
*/
inline double rsqrt(double a)
{
    using std::sqrt;
    return 1 / sqrt(a);
}

/**
  \brief Integer power
  constexpr (can be used in template parameters)
*/
constexpr int power(int base, unsigned int exponent, int result = 1)
{
    return (exponent < 1) ? result : power(base * base, exponent / 2, (exponent % 2) ? result * base : result);
}

/**
 square
*/
template <class T>
inline T sqr(T a)
{
    return a * a;
}

template <class T>
inline T cube(T a)
{
    return a * a * a;
}

template <class T>
inline T to_the_fourth(T a)
{
    return a * a * a * a;
}

template <class T>
inline T clamp_small_magnitude(const T x, const T eps)
{
    assert(eps >= 0);
    if (x < -eps)
        return x;
    else if (x < 0)
        return -eps;
    else if (x < eps)
        return eps;
    else
        return x;
}

template <class T>
inline void clamp(T& x, const T min, const T max)
{
    assert(max >= min);
    if (x < min)
        x = min;
    else if (x > max)
        x = max;
}

template <class T, int dim>
inline void clamp(Vector<T, dim>& x, const T min, const T max)
{
    assert(max >= min);
    for (int d = 0; d < dim; ++d) clamp(x(d), min, max);
}

/**
   Robustly computing log(x+1)/x
 */
template <class T>
inline T log_1px_over_x(const T x, const T eps)
{
    assert(eps > 0);
    if (std::fabs(x) < eps)
        return (T)1;
    else
        return std::log1p(x) / x;
}

/**
   Robustly computing (logx-logy)/(x-y)
 */
template <class T>
inline T diff_log_over_diff(const T x, const T y, const T eps)
{
    assert(eps > 0);
    T p = x / y - 1;
    return log_1px_over_x(p, eps) / y;
}

/**
   Robustly computing (expx-1)/x
 */
template <class T>
inline T exp_m1x_over_x(const T x, const T eps)
{
    assert(eps > 0);
    if (std::fabs(x) < eps)
        return (T)1;
    else
        return std::expm1(x) / x;
}

/**
   Robustly computing (expx-expy)/(x-y)
 */
template <class T>
inline T diff_exp_over_diff(const T x, const T y, const T eps)
{
    assert(eps > 0);
    T p = x - y;
    return exp_m1x_over_x(p, eps) * std::exp(y);
}

/**
   Robustly computing sin(x)/x
 */
template <class T>
inline T sinx_over_x(const T x, const T eps)
{
    assert(eps > 0);
    if (std::fabs(x) < eps)
        return 1;
    else
        return std::sin(x) / x;
}

/**
   Robustly computing (x logy- y logx)/(x-y)
 */
template <class T>
inline T diff_interlock_log_over_diff(const T x, const T y, const T logy, const T eps)
{
    assert(eps > 0);
    return logy - y * diff_log_over_diff(x, y, eps);
}

template <class T>
inline T ziran_cross_product(const Vector<T, 2>& a, const Vector<T, 2>& b)
{
    return a(0) * b(1) - a(1) * b(0);
}

template <class T>
inline T ziran_cross_product(const Vector<T, 3>& a, const Vector<T, 3>& b)
{
    return a.cross(b);
}

/* Compute macaulay brackets (+) */
template <class T>
inline T macaulay(const T x)
{
    return (x + std::abs(x)) / 2;
}
}
}
