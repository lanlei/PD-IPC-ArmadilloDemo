#include "FiniteDiff.h"
#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <Eigen/Sparse>
#include <Bow/Utils/Logging.h>

namespace Bow {
namespace FiniteDiff {

template <bool MatrixFree, class Jacobian>
BOW_INLINE void ziran_check_impl(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& energy,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& gradient,
    const Jacobian& hessian,
    const std::function<void(Eigen::VectorXd&)>& project,
    const std::function<void(Eigen::VectorXd&)>& transform,
    double diff_test_perturbation_scale)
{
    Logging::info("Running diff test with perturbation scale: ", diff_test_perturbation_scale);
    using Vec = Eigen::VectorXd;
    // Eigen::VectorXd step, f0, df0, f1, df1;
    Vec step = Vec::Random(x.size());
    if (project != nullptr) {
        project(step);
    }
    step.normalize();
    Vec raw_step = step;
    if (transform != nullptr) {
        transform(raw_step);
    }
    double e0 = energy(x);
    Vec f0;
    gradient(x, f0);
    f0 = -f0;
    Vec df0;
    if constexpr (!MatrixFree) {
        Eigen::SparseMatrix<double> jacobian0;
        hessian(x, jacobian0);
        df0 = jacobian0 * raw_step;
    }
    else {
        hessian(x, raw_step, df0);
    }

    const int DIFF_SIZE = 20;
    std::vector<double> energy_difference(DIFF_SIZE);
    std::vector<double> energy_differential(DIFF_SIZE);
    std::vector<double> energy_err(DIFF_SIZE);
    std::vector<double> energy_log_err(DIFF_SIZE);
    std::vector<double> force_difference_norm(DIFF_SIZE);
    std::vector<double> force_differential_norm(DIFF_SIZE);
    std::vector<double> force_err(DIFF_SIZE);
    std::vector<double> force_log_err(DIFF_SIZE);

    std::setprecision(12);
    std::cout << "e0\t=" << std::setw(20) << e0 << std::endl;

    for (int i = 1; i <= DIFF_SIZE; ++i) {
        double h = diff_test_perturbation_scale * std::pow((double)(2), -i);
        Vec x_new = x + h * step;
        double e1 = energy(x_new);
        std::cout << "e1\t=" << std::setw(20) << e1 << "\th = " << h << std::endl;
        Vec f1;
        gradient(x_new, f1);
        f1 = -f1;
        double difference = (e0 - e1) / h;
        double differential = (f0 + f1).dot(raw_step) / 2;

        double err = (difference - differential);
        double log_err = std::log(std::abs(err));
        energy_difference[i - 1] = difference;
        energy_differential[i - 1] = differential;
        energy_err[i - 1] = err;
        energy_log_err[i - 1] = log_err;

        Vec df1;
        if constexpr (!MatrixFree) {
            Eigen::SparseMatrix<double> jacobian1;
            hessian(x_new, jacobian1);
            df1 = jacobian1 * raw_step;
        }
        else {
            hessian(x_new, raw_step, df1);
        }
        Vec force_difference = (1 / h) * (f0 - f1);
        Vec force_differential = 0.5 * (df0 + df1);
        double err_force = (force_difference - force_differential).norm();
        double log_err_force = std::log(std::abs(err_force));
        force_difference_norm[i - 1] = force_difference.norm();
        force_differential_norm[i - 1] = force_differential.norm();
        force_err[i - 1] = err_force;
        force_log_err[i - 1] = log_err_force;
    }

    std::cout << std::setprecision(12) << "energy["
              << "i"
              << "] = " << std::setw(20) << "difference"
              << std::setw(20) << "differential"
              << std::setw(20) << "err"
              << std::setw(20) << "log_err"
              << std::endl;
    for (int i = 0; i < DIFF_SIZE; ++i) {
        std::cout << std::setprecision(12) << "energy[" << i << "] = " << std::setw(20) << energy_difference[i]
                  << std::setw(20) << energy_differential[i]
                  << std::setw(20) << energy_err[i]
                  << std::setw(20) << energy_log_err[i]
                  << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::setprecision(12) << "force["
              << "i"
              << "] = " << std::setw(20) << "difference_norm"
              << std::setw(20) << "differential_norm"
              << std::setw(20) << "err"
              << std::setw(20) << "log_err"
              << std::endl;
    for (int i = 0; i < DIFF_SIZE; ++i) {
        std::cout << std::setprecision(12) << "force[" << i << "] = " << std::setw(20) << force_difference_norm[i]
                  << std::setw(20) << force_differential_norm[i]
                  << std::setw(20) << force_err[i]
                  << std::setw(20) << force_log_err[i]
                  << std::endl;
    }
    std::cin.get();
}

// fewer DoF preferred
// transform is from original space to reduced space
BOW_INLINE void ziran_check_false(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& energy,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& gradient,
    const std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&)>& hessian,
    const std::function<void(Eigen::VectorXd&)>& project,
    const std::function<void(Eigen::VectorXd&)>& transform,
    double diff_test_perturbation_scale)
{
    ziran_check_impl<false>(x, energy, gradient, hessian, project, transform, diff_test_perturbation_scale);
}
BOW_INLINE void ziran_check_true(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& energy,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& gradient,
    const std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::VectorXd&)>& hessian,
    const std::function<void(Eigen::VectorXd&)>& project,
    const std::function<void(Eigen::VectorXd&)>& transform,
    double diff_test_perturbation_scale)
{
    ziran_check_impl<true>(x, energy, gradient, hessian, project, transform, diff_test_perturbation_scale);
}

BOW_INLINE void finite_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    Eigen::VectorXd& grad,
    const AccuracyOrder accuracy,
    const double eps)
{
    // Create an array of the coefficients for finite differences.
    // See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    // clang-format off
    // The external coefficients, c1, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    // The internal coefficients, c2, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    // clang-format on
    // The denominators of the finite difference.
    static const std::array<double, 4> dd = { { 2, 12, 60, 840 } };

    grad.resize(x.size());

    const size_t innerSteps = 2 * (accuracy + 1);
    const double ddVal = dd[accuracy] * eps;

    Eigen::VectorXd xx = x;
    for (long d = 0; d < x.rows(); d++) {
        grad[d] = 0;
        for (size_t s = 0; s < innerSteps; ++s) {
            xx[d] += coeff2[accuracy][s] * eps;
            grad[d] += coeff[accuracy][s] * f(xx);
            xx[d] = x[d];
        }
        grad[d] /= ddVal;
    }
}

BOW_INLINE void finite_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    Eigen::MatrixXd& jac,
    const AccuracyOrder accuracy,
    const double eps)
{
    // Create an array of the coefficients for finite differences.
    // See: https://en.wikipedia.org/wiki/Finite_difference_coefficient
    // clang-format off
    // The external coefficients, c1, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    // The internal coefficients, c2, in c1 * f(x + c2).
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    // clang-format on
    // The denominators of the finite difference.
    static const std::array<double, 4> dd = { { 2, 12, 60, 840 } };

    Eigen::VectorXd fv;
    f(x, fv);
    jac.resize(fv.size(), x.size());

    const size_t innerSteps = 2 * (accuracy + 1);
    const double ddVal = dd[accuracy] * eps;

    Eigen::VectorXd xx = x;
    for (long d = 0; d < x.rows(); d++) {
        jac.col(d).setZero();
        for (size_t s = 0; s < innerSteps; ++s) {
            xx[d] += coeff2[accuracy][s] * eps;
            Eigen::VectorXd fv;
            f(xx, fv);
            jac.col(d) += coeff[accuracy][s] * fv;
            xx[d] = x[d];
        }
        jac.col(d) /= ddVal;
    }
}

BOW_INLINE bool check_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const double eps, const double pass_ratio)
{
    std::random_device rd;
    Eigen::VectorXd direction = Eigen::VectorXd::Random(x.size());
    direction /= direction.cwiseAbs().maxCoeff();
    Eigen::VectorXd dx = eps * direction;
    return check_gradient(x, f, g, dx, pass_ratio);
}

template <class Jacobian>
BOW_INLINE bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Jacobian&)>& g,
    const double eps, const double pass_ratio)
{
    std::random_device rd;
    Eigen::VectorXd direction = Eigen::VectorXd::Random(x.size());
    direction /= direction.norm();
    Eigen::VectorXd dx = eps * direction;
    return check_jacobian(x, f, g, dx, pass_ratio);
}

BOW_INLINE bool check_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio)
{
    double eps = dx.norm();
    Eigen::VectorXd x0 = x - dx;
    Eigen::VectorXd x1 = x + dx;
    double f0 = f(x0);
    double f1 = f(x1);
    Eigen::VectorXd g0, g1;
    g(x0, g0);
    g(x1, g1);
    double true_value = std::abs(f1 - f0 - (g1 + g0).transpose() * dx) / eps;
    double fake_value = std::abs(f1 - f0 - 2 * (g1 + g0).transpose() * dx) / eps;
    Bow::Logging::info("[Check Gradient] real_value: ", true_value, ",\tfake_value: ", fake_value, ",\tratio: ", true_value / fake_value);
    return true_value / fake_value < pass_ratio;
}

template <class Jacobian>
BOW_INLINE bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Jacobian&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio)
{
    double eps = dx.cwiseAbs().maxCoeff();
    Eigen::VectorXd x0 = x - dx;
    Eigen::VectorXd x1 = x + dx;
    Eigen::VectorXd f0, f1;
    f(x0, f0);
    f(x1, f1);
    Jacobian g0, g1;
    g(x0, g0);
    g(x1, g1);
    double true_value = (f1 - f0 - (g1 + g0) * dx).cwiseAbs().maxCoeff() / eps;
    double fake_value = (f1 - f0 - 2 * (g1 + g0) * dx).cwiseAbs().maxCoeff() / eps;
    Bow::Logging::info("[Check Jacobian] real_value: ", true_value, ",\tfake_value: ", fake_value, ",\tratio: ", true_value / fake_value);
    return true_value / fake_value < pass_ratio;
}

BOW_INLINE bool check_jacobian_matrix_free(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, const Eigen::VectorXd&, Eigen::VectorXd&)>& g,
    const double eps, const double pass_ratio)
{
    std::random_device rd;
    Eigen::VectorXd dx = eps * Eigen::VectorXd::Random(x.size());
    Eigen::VectorXd x0 = x - dx;
    Eigen::VectorXd x1 = x + dx;
    Eigen::VectorXd f0, f1;
    f(x0, f0);
    f(x1, f1);
    Eigen::VectorXd g0, g1;
    g(x0, dx, g0);
    g(x1, dx, g1);
    double true_value = (f1 - f0 - (g1 + g0)).cwiseAbs().maxCoeff() / eps;
    double fake_value = (f1 - f0 - 2 * (g1 + g0)).cwiseAbs().maxCoeff() / eps;
    Bow::Logging::info("[Check Jacobian] real_value: ", true_value, ",\tfake_value: ", fake_value, ",\tratio: ", true_value / fake_value);
    return true_value / fake_value < pass_ratio;
}

#ifdef BOW_STATIC_LIBRARY
template bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&)>& g,
    const double eps, const double pass_ratio);

template bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::MatrixXd&)>& g,
    const double eps, const double pass_ratio);

template bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::SparseMatrix<double>&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio);

template bool check_jacobian(
    const Eigen::VectorXd& x,
    const std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& f,
    const std::function<void(const Eigen::VectorXd&, Eigen::MatrixXd&)>& g,
    const Eigen::VectorXd& dx, const double pass_ratio);

#endif
}
} // namespace Bow::FiniteDiff
