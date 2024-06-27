#ifdef BOW_AMGCL
#include "Amgcl.h"
#include <Bow/Utils/Timer.h>
#include <Eigen/SparseCore>
#include <Bow/Utils/Logging.h>
#include <boost/foreach.hpp>

namespace Bow {
namespace Math {
namespace LinearSolver {

namespace internal {
/* https://stackoverflow.com/questions/15904896/range-based-for-loop-on-a-dynamic-array */
template <typename T>
struct WrappedArray {
    WrappedArray(const T* first, const T* last)
        : begin_{ first }, end_{ last } {}
    WrappedArray(const T* first, std::ptrdiff_t size)
        : WrappedArray{ first, first + size } {}

    const T* begin() const noexcept { return begin_; }
    const T* end() const noexcept { return end_; }
    const T& operator[](const size_t i) const { return begin_[i]; }

    const T* begin_;
    const T* end_;
};
} // namespace internal

template <class Derived, bool cuda_backend>
AmgclSolver<Derived, cuda_backend>::AmgclSolver()
{
    // prm.put("precond.class", "amg");
    // prm.put("precond.coarsening.type", "smoothed_aggregation");
    // prm.put("solver.type", "lgmres");
    // prm.put("solver.M", 100);
    // prm.put("precond.coarsening.aggr.eps_strong", 0.0);
    // prm.put("precond.relax.type", "gauss_seidel");
    // prm.put("solver.tol", T(1e-8));
    // prm.put("solver.maxiter", 100);
    prm.put("solver.tol", 1e-7); // relative
    prm.put("solver.maxiter", 1000);
    prm.put("precond.class", "amg");
    prm.put("precond.relax.type", "chebyshev");
    prm.put("precond.relax.degree", 16);
    prm.put("precond.relax.power_iters", 100);
    prm.put("precond.relax.higher", 2.0f);
    prm.put("precond.relax.lower", 1.0f / 120.0f);
    prm.put("precond.relax.scale", true);
    prm.put("precond.max_levels", 6);
    prm.put("precond.direct_coarse", false);
    prm.put("precond.ncycle", 2);
    prm.put("precond.coarsening.type", "smoothed_aggregation");
    prm.put("precond.coarsening.estimate_spectral_radius", true);
    prm.put("precond.coarsening.relax", 1.0f);
    prm.put("precond.coarsening.aggr.eps_strong", 0.0);
    // prm.put("precond.coarsening.aggr.block_size", block_size);
    prm.put("solver.type", "lgmres");
    prm.put("solver.M", 100);
}
template <class Derived, bool cuda_backend>
AmgclSolver<Derived, cuda_backend>::AmgclSolver(boost::property_tree::ptree prm_in)
    : AmgclSolver()
{
    BOOST_FOREACH (boost::property_tree::ptree::value_type& v, prm_in) {
        prm.put(v.first, v.second.data());
    }
}

template <class Derived, bool cuda_backend>
AmgclSolver<Derived, cuda_backend>::~AmgclSolver()
{
}

template <class Derived, bool cuda_backend>
bool AmgclSolver<Derived, cuda_backend>::compute(const Eigen::SparseMatrixBase<Derived>& mat)
{
    n_rows = mat.rows();
    n_cols = mat.cols();
    internal::WrappedArray<StorageIndex> ptr(mat.derived().outerIndexPtr(), n_rows + 1);
    internal::WrappedArray<StorageIndex> col(mat.derived().innerIndexPtr(), mat.derived().nonZeros());
    internal::WrappedArray<typename Derived::Scalar> val(mat.derived().valuePtr(), mat.derived().nonZeros());
    auto A = std::make_tuple(n_rows, ptr, col, val);
#ifdef ENABLE_AMGCL_CUDA
    if constexpr (cuda_backend) {
        vex::Context ctx(vex::Filter::Env);
        Logging::info("Computation Context: ", ctx);
        bprm.q = ctx;
        solver = std::make_unique<Solver>(A, prm, bprm);
    }
    else
        solver = std::make_unique<Solver>(A, prm);
#else
    solver = std::make_unique<Solver>(A, prm);
#endif
    Logging::info(solver->precond());
    return true;
}

template <class Derived, bool cuda_backend>
template <class DerivedB>
Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> AmgclSolver<Derived, cuda_backend>::solve(const Eigen::MatrixBase<DerivedB>& rhs) const
{
    std::vector<T> F(rhs.derived().data(), rhs.derived().data() + rhs.size());
    std::vector<T> X(n_cols, T(0.0));
    auto f_b = Backend::copy_vector(F, bprm);
    auto x_b = Backend::copy_vector(X, bprm);

    size_t iters;
    double resid;
    std::tie(iters, resid) = (*solver)(*f_b, *x_b);
    Logging::info("Iterations: ", iters);
    Logging::info("Error:      ", resid);
#ifdef ENABLE_AMGCL_CUDA
    if constexpr (cuda_backend)
        vex::copy(*x_b, X);
    else
        std::copy(&(*x_b)[0], &(*x_b)[0] + X.size(), X.data());
#else
    std::copy(&(*x_b)[0], &(*x_b)[0] + X.size(), X.data());
#endif
    Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> result = Eigen::Map<Bow::Vector<typename Derived::Scalar, Eigen::Dynamic>>(X.data(), n_cols);
    return result;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template class AmgclSolver<Eigen::SparseMatrix<float>>;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template class AmgclSolver<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
#endif
#ifdef BOW_COMPILE_DOUBLE
template class AmgclSolver<Eigen::SparseMatrix<double>>;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template class AmgclSolver<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
#endif

#if ENABLE_AMGCL_CUDA
#ifdef BOW_COMPILE_FLOAT
template class AmgclSolver<Eigen::SparseMatrix<float>, true>;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float>, true>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float>, true>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXf>>& rhs) const;
template class AmgclSolver<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>, true>;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>, true>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template Eigen::VectorXf AmgclSolver<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>, true>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXf>>& rhs) const;
#endif
#ifdef BOW_COMPILE_DOUBLE
template class AmgclSolver<Eigen::SparseMatrix<double>, true>;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double>, true>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double>, true>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXd>>& rhs) const;
template class AmgclSolver<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>, true>;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>, true>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template Eigen::VectorXd AmgclSolver<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>, true>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXd>>& rhs) const;
#endif
#endif

#endif
}
}} // namespace Bow::Math::LinearSolver
#endif