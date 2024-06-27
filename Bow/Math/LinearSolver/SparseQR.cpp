#ifdef BOW_SUITESPARSE
#include "SparseQR.h"
#include <Bow/Utils/Timer.h>
#include <Eigen/SparseCore>
#include <oneapi/tbb.h>

namespace Bow {
namespace Math {
namespace LinearSolver {
template <class Derived>
SparseQR<Derived>::SparseQR()
{
    // if constexpr (std::is_same<typename Derived::Scalar, double>::value)
    cholmod_l_start(&cm);
    cm.dtype = CHOLMOD_DOUBLE;
    // else
    // cm.dtype = CHOLMOD_SINGLE;
}

template <class Derived>
SparseQR<Derived>::~SparseQR()
{
    if (A)
        cholmod_l_free_sparse(&A, &cm);
    if (QR)
        SuiteSparseQR_free(&QR, &cm);
    cholmod_l_finish(&cm);
}

template <class Derived>
bool SparseQR<Derived>::compute(const Eigen::SparseMatrixBase<Derived>& mat)
{
    // the Sparse QR routines are only implemented for "long" indices and double precision.
    // https://forum.kde.org/viewtopic.php?f=74&t=96087
    if (A)
        cholmod_l_free_sparse(&A, &cm);
    A = cholmod_l_allocate_sparse(mat.rows(), mat.cols(), mat.derived().nonZeros(), true, true, 0, CHOLMOD_REAL, &cm);

    tbb::parallel_for(0, (int)mat.rows() + 1, [&](int i) {
        static_cast<long*>(A->p)[i] = (long)mat.derived().outerIndexPtr()[i];
    });
    tbb::parallel_for(0, (int)mat.derived().nonZeros(), [&](int i) {
        static_cast<long*>(A->i)[i] = (long)mat.derived().innerIndexPtr()[i];
    });
    tbb::parallel_for(0, (int)mat.derived().nonZeros(), [&](int i) {
        static_cast<double*>(A->x)[i] = (double)mat.derived().valuePtr()[i];
    });
    if (QR)
        SuiteSparseQR_free<double>(&QR, &cm);
    QR = SuiteSparseQR_factorize<double>(SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, A, &cm);
    return true;
}

template <class Derived>
template <class DerivedB>
Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> SparseQR<Derived>::solve(const Eigen::MatrixBase<DerivedB>& _rhs) const
{
    Eigen::VectorXd rhs = _rhs.template cast<double>();
    cholmod_common _cm = cm;
    cholmod_dense* b = cholmod_l_allocate_dense(rhs.size(), 1, rhs.size(), CHOLMOD_REAL, &_cm);
    memcpy(b->x, rhs.derived().data(), rhs.size() * sizeof(double));
    cholmod_dense* y = SuiteSparseQR_qmult<double>(SPQR_QTX, QR, b, &_cm);
    cholmod_dense* x = SuiteSparseQR_solve<double>(SPQR_RETX_EQUALS_B, QR, y, &_cm);
    Eigen::VectorXd result(x->nrow);
    memcpy(result.data(), x->x, result.size() * sizeof(typename Derived::Scalar));
    cholmod_l_free_dense(&b, &_cm);
    cholmod_l_free_dense(&x, &_cm);
    cholmod_l_free_dense(&y, &_cm);
    return result.template cast<typename Derived::Scalar>();
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template class SparseQR<Eigen::SparseMatrix<float>>;
template Eigen::VectorXf SparseQR<Eigen::SparseMatrix<float>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
template class SparseQR<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>;
template Eigen::VectorXf SparseQR<Eigen::SparseMatrix<float, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXf>& rhs) const;
#endif
#ifdef BOW_COMPILE_DOUBLE
template class SparseQR<Eigen::SparseMatrix<double>>;
template Eigen::VectorXd SparseQR<Eigen::SparseMatrix<double>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
template class SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>;
template Eigen::VectorXd SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor, long int>>::solve(const Eigen::MatrixBase<Eigen::VectorXd>& rhs) const;
#endif
#endif
}
}} // namespace Bow::Math::LinearSolver
#endif