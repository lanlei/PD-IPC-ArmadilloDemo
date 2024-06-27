#ifdef BOW_SUITESPARSE
#include "SparseCholesky.h"
#include <Bow/Utils/Timer.h>
#include <Eigen/SparseCore>

namespace Bow {
namespace Math {
namespace LinearSolver {
template <class Derived>
CholmodLLT<Derived>::CholmodLLT(const int supernodal)
{
    cholmod_start(&cm);
    if constexpr (std::is_same<typename Derived::Scalar, double>::value)
        cm.dtype = CHOLMOD_DOUBLE;
    else
        cm.dtype = CHOLMOD_SINGLE;
    cm.supernodal = supernodal;
    if constexpr (std::is_same<typename Derived::StorageIndex, long>::value)
        cm.itype = CHOLMOD_LONG;
    else
        cm.itype = CHOLMOD_INT;
}

template <class Derived>
CholmodLLT<Derived>::~CholmodLLT()
{
    if (A)
        cholmod_free_sparse(&A, &cm);
    if (L)
        cholmod_free_factor(&L, &cm);
    cholmod_finish(&cm);
}

template <class Derived>
bool CholmodLLT<Derived>::compute(const Eigen::SparseMatrixBase<Derived>& mat)
{
    if (A)
        cholmod_free_sparse(&A, &cm);
    A = cholmod_allocate_sparse(mat.rows(), mat.cols(), mat.derived().nonZeros(), true, true, -1, CHOLMOD_REAL, &cm);
    memcpy(A->p, mat.derived().outerIndexPtr(), (mat.rows() + 1) * sizeof(typename Derived::StorageIndex));
    memcpy(A->i, mat.derived().innerIndexPtr(), mat.derived().nonZeros() * sizeof(typename Derived::StorageIndex));
    memcpy(A->x, mat.derived().valuePtr(), mat.derived().nonZeros() * sizeof(typename Derived::Scalar));
    if (L)
        cholmod_free_factor(&L, &cm);
    L = cholmod_analyze(A, &cm);
    cholmod_factorize(A, L, &cm);
    return cm.status != CHOLMOD_NOT_POSDEF;
}

template <class Derived>
template <class DerivedB>
Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> CholmodLLT<Derived>::solve(const Eigen::MatrixBase<DerivedB>& rhs) const
{
    cholmod_common _cm = cm;
    cholmod_dense* b = cholmod_allocate_dense(rhs.size(), 1, rhs.size(), CHOLMOD_REAL, &_cm);
    memcpy(b->x, rhs.derived().data(), rhs.size() * sizeof(typename Derived::Scalar));
    cholmod_dense* x = cholmod_solve(CHOLMOD_A, L, b, &_cm);
    Bow::Vector<typename DerivedB::Scalar, Eigen::Dynamic> result(x->nrow);
    memcpy(result.data(), x->x, result.size() * sizeof(typename Derived::Scalar));
    cholmod_free_dense(&b, &_cm);
    cholmod_free_dense(&x, &_cm);
    return result;
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template class CholmodLLT<Eigen::SparseMatrix<float>>;
template Eigen::VectorXf CholmodLLT<Eigen::SparseMatrix<float>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template Eigen::VectorXf CholmodLLT<Eigen::SparseMatrix<float>>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXf>>& rhs) const;
template class CholmodLLT<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>;
template Eigen::VectorXf CholmodLLT<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template Eigen::VectorXf CholmodLLT<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXf>>& rhs) const;
#endif
#ifdef BOW_COMPILE_DOUBLE
template class CholmodLLT<Eigen::SparseMatrix<double>>;
template Eigen::VectorXd CholmodLLT<Eigen::SparseMatrix<double>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template Eigen::VectorXd CholmodLLT<Eigen::SparseMatrix<double>>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXd>>& rhs) const;
template class CholmodLLT<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>;
template Eigen::VectorXd CholmodLLT<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template Eigen::VectorXd CholmodLLT<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::Map<const Eigen::VectorXd>>& rhs) const;
#endif
#endif
}
}} // namespace Bow::Math::LinearSolver
#endif