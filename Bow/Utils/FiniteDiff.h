#ifndef FINITE_DIFF_H
#define FINITE_DIFF_H

#include <Bow/Macros.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Bow {
namespace FiniteDiff {

/**
 * Adopted from https://github.com/zfergus/finite-diff
 */

/**
 * @brief Enumeration of available orders of accuracy for finite differences.
 *
 * The corresponding integer values are used internally and should be ignored.
 */
enum AccuracyOrder {
    SECOND = 0, ///< @brief Second order accuracy.
    FOURTH = 1, ///< @brief Fourth order accuracy.
    SIXTH = 2, ///< @brief Sixth order accuracy.
    EIGHTH = 3 ///< @brief Eighth order accuracy.
};

// MatrixFree = false
BOW_INLINE void ziran_check_false(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& energy,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& gradient,
    const std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&)>& hessian,
    const std::function<void(Eigen::VectorXd&)>& project = nullptr,
    const std::function<void(Eigen::VectorXd&)>& transform = nullptr,
    double diff_test_perturbation_scale = 1.);

// MatrixFree = true
BOW_INLINE void ziran_check_true(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& energy,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& gradient,
    const std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::VectorXd&)>& hessian,
    const std::function<void(Eigen::VectorXd&)>& project = nullptr,
    const std::function<void(Eigen::VectorXd&)>& transform = nullptr,
    double diff_test_perturbation_scale = 1.);

/**
 * @brief Compute the gradient of a function using finite differences.
 *
 * @param[in]  x         Point at which to compute the gradient.
 * @param[in]  f         Compute the gradient of this function.
 * @param[out] grad      Computed gradient.
 * @param[in]  eps       Value of the finite difference step.
 */
BOW_INLINE void finite_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    Eigen::VectorXd& grad,
    const AccuracyOrder accuracy = SECOND,
    const double eps = 1.0e-8);

/**
 * @brief Compute the jacobian of a function using finite differences.
 *
 * @param[in]  x         Point at which to compute the jacobian.
 * @param[in]  f         Compute the jacobian of this function.
 * @param[out] jac       Computed jacobian.
 * @param[in]  eps       Value of the finite difference step.
 */
BOW_INLINE void finite_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    Eigen::MatrixXd& jac,
    const AccuracyOrder accuracy = SECOND,
    const double eps = 1.0e-8);

/**
 * https://www.cs.ucr.edu/~craigs/papers/2019-derivatives/course.pdf
 */
BOW_INLINE bool check_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const double eps = 1e-4, const double pass_ratio = 1e-3);

/**
 * https://www.cs.ucr.edu/~craigs/papers/2019-derivatives/course.pdf
 */
template <class Jacobian>
BOW_INLINE bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Jacobian&)>& g,
    const double eps = 1e-4, const double pass_ratio = 1e-3);

BOW_INLINE bool check_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio = 1e-3);

template <class Jacobian>
BOW_INLINE bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Jacobian&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio = 1e-3);

/**
 * https://www.cs.ucr.edu/~craigs/papers/2019-derivatives/course.pdf
 */
BOW_INLINE bool check_jacobian_matrix_free(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const double eps = 1e-4, const double pass_ratio = 1e-3);
}
} // namespace Bow::FiniteDiff

#ifndef BOW_STATIC_LIBRARY
#include "FiniteDiff.cpp"
#endif

#endif