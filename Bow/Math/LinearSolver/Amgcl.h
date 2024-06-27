#ifdef BOW_AMGCL
#ifndef BOW_AMGCL_H
#define BOW_AMGCL_H
#include <Bow/Macros.h>
#include <Eigen/Core>
#include <Bow/Types.h>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/plain_aggregates.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/io/mm.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#ifdef ENABLE_AMGCL_CUDA
#include <amgcl/backend/vexcl.hpp>
#endif
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/reorder.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/profiler.hpp>
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <type_traits>

/**
 * CAUTION: The matrix is assumed to be in row-major format, since AMGCL assumes 
 *          that the outer index is for row. If the matrix is symmetric, you 
 *          are fine, because CSR and CSC are the same. If the matrix is not 
 *          symmetric and you pass in a column-major matrix, the solver will 
 *          actually solve A^T x = b.
 */

namespace Bow {
namespace Math {
namespace LinearSolver {
template <class Derived, bool cuda_backend = false>
class AmgclSolver {
    using T = typename Derived::Scalar;
#ifdef ENABLE_AMGCL_CUDA
    typedef typename std::conditional<cuda_backend, amgcl::backend::vexcl<T>, amgcl::backend::builtin<T>>::type Backend;
#else
    typedef amgcl::backend::builtin<T> Backend;
#endif
    using Solver = amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>>;
    typename Backend::params bprm;
    std::unique_ptr<Solver> solver;
    int n_rows;
    int n_cols;

public:
    boost::property_tree::ptree prm;
    using StorageIndex = typename Derived::StorageIndex;
    AmgclSolver();
    AmgclSolver(boost::property_tree::ptree prm);
    ~AmgclSolver();
    bool compute(const Eigen::SparseMatrixBase<Derived>& mat);
    template <class DerivedB>
    Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> solve(const Eigen::MatrixBase<DerivedB>& rhs) const;
};
}
}
} // namespace Bow::Math::LinearSolver

#ifndef BOW_STATIC_LIBRARY
#include "Amgcl.cpp"
#endif

#endif
#endif