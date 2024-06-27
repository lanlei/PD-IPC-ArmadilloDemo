#include "Utils.h"
#include <Bow/Types.h>
#include <oneapi/tbb.h>

namespace Bow {
namespace Math {
template <class DerivedA, class DerivedCoA>
BOW_INLINE void cofactor(const Eigen::MatrixBase<DerivedA>& A, Eigen::MatrixBase<DerivedCoA>& coA)
{
    static const int dim = DerivedA::RowsAtCompileTime;
    if constexpr (dim == 2) {
        coA(0, 0) = A(1, 1);
        coA(0, 1) = -A(1, 0);
        coA(1, 0) = -A(0, 1);
        coA(1, 1) = A(0, 0);
    }
    else if constexpr (dim == 3) {
        coA(0, 0) = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
        coA(0, 1) = A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2);
        coA(0, 2) = A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0);
        coA(1, 0) = A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2);
        coA(1, 1) = A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0);
        coA(1, 2) = A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1);
        coA(2, 0) = A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1);
        coA(2, 1) = A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2);
        coA(2, 2) = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    }
    else {
        static_assert(dim == 2 || dim == 3, "cofactor is only implemented for 2x2 and 3x3 matrices!!");
    }
}

template <class DerivedA>
BOW_INLINE void make_pd(Eigen::MatrixBase<DerivedA>& sym_A)
{
    using T = typename DerivedA::Scalar;
    static const int dim = DerivedA::RowsAtCompileTime;
    Eigen::SelfAdjointEigenSolver<DerivedA> eigenSolver(sym_A);
    if (eigenSolver.eigenvalues()[0] >= 0) {
        return;
    }
    Vector<T, dim> D = eigenSolver.eigenvalues();
    for (int i = 0; i < dim; ++i) {
        if (D[i] < 0) {
            D[i] = 0;
        }
        else {
            break;
        }
    }
    sym_A = eigenSolver.eigenvectors() * D.asDiagonal() * eigenSolver.eigenvectors().transpose();
}

template <class Derived>
BOW_INLINE void sparse_from_csr(const std::vector<typename Derived::StorageIndex>& outer_indices,
    const std::vector<typename Derived::StorageIndex>& inner_indices,
    const std::vector<typename Derived::Scalar>& val,
    const typename Derived::StorageIndex nrows,
    const typename Derived::StorageIndex ncols,
    Eigen::SparseMatrixBase<Derived>& A)
{
    A.derived().setZero();
    A.derived().resize(nrows, ncols);
    A.derived().reserve(val.size());
    memcpy(A.derived().valuePtr(), val.data(),
        val.size() * sizeof(typename Derived::Scalar));
    memcpy(A.derived().innerIndexPtr(), inner_indices.data(),
        inner_indices.size() * sizeof(typename Derived::StorageIndex));
    memcpy(A.derived().outerIndexPtr(), outer_indices.data(),
        outer_indices.size() * sizeof(typename Derived::StorageIndex));
    A.derived().finalize();
}

template <class T, class StorageIndex>
BOW_INLINE void select_cols(const Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& A,
    const std::vector<StorageIndex>& selection,
    Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& A_cols)
{
    using Scalar = T;
    std::vector<StorageIndex> outer_indices;
    std::vector<StorageIndex> inner_indices;
    std::vector<Scalar> val;
    const auto* A_outer_index_ptr = A.derived().outerIndexPtr();
    const auto* A_inner_index_ptr = A.derived().innerIndexPtr();
    const auto* A_val_ptr = A.derived().valuePtr();
    outer_indices.clear();
    inner_indices.clear();
    val.clear();
    outer_indices.push_back(0);
    for (auto index : selection) {
        int start = A_outer_index_ptr[index];
        int end = A_outer_index_ptr[index + 1];
        int inner_length = end - start;
        inner_indices.resize(inner_indices.size() + inner_length);
        memcpy(&inner_indices[outer_indices.back()], &A_inner_index_ptr[start], inner_length * sizeof(StorageIndex));
        val.resize(val.size() + inner_length);
        memcpy(&val[outer_indices.back()], &A_val_ptr[start], inner_length * sizeof(Scalar));
        outer_indices.push_back(inner_length + outer_indices.back());
    }
    //sparse_from_csr(outer_indices, inner_indices, val, A.rows(), selection.size(), A_cols);
}

template <class T, class StorageIndex>
BOW_INLINE void shift(const Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& out)
{
    out.resize(new_rows, new_cols);
    out.reserve(in.nonZeros());
    memcpy(out.outerIndexPtr() + dcol, in.outerIndexPtr(), sizeof(StorageIndex) * (in.cols() + 1));
    memcpy(out.valuePtr(), in.valuePtr(), sizeof(T) * in.nonZeros());
    tbb::parallel_for(long(0), long(in.nonZeros()), [&](long i) {
        out.innerIndexPtr()[i] = in.innerIndexPtr()[i] + drow;
    });
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
template void make_pd(Eigen::MatrixBase<Matrix<float, 2, 2>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<float, 3, 3>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<float, 4, 4>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<float, 6, 6>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<float, 9, 9>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<float, 12, 12>>& symA);
template void shift(const Eigen::SparseMatrix<float, Eigen::ColMajor, int>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<float, Eigen::ColMajor, int>& out);
template void shift(const Eigen::SparseMatrix<float, Eigen::ColMajor, long>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<float, Eigen::ColMajor, long>& out);
#ifdef BOW_COMPILE_2D
template void cofactor(const Eigen::MatrixBase<Matrix<float, 2, 2>>& A, Eigen::MatrixBase<Matrix<float, 2, 2>>& coA);
#endif
#ifdef BOW_COMPILE_3D
template void cofactor(const Eigen::MatrixBase<Matrix<float, 3, 3>>& A, Eigen::MatrixBase<Matrix<float, 3, 3>>& coA);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
template void make_pd(Eigen::MatrixBase<Matrix<double, 2, 2>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<double, 3, 3>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<double, 4, 4>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<double, 6, 6>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<double, 9, 9>>& symA);
template void make_pd(Eigen::MatrixBase<Matrix<double, 12, 12>>& symA);
template void shift(const Eigen::SparseMatrix<double, Eigen::ColMajor, int>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<double, Eigen::ColMajor, int>& out);
template void shift(const Eigen::SparseMatrix<double, Eigen::ColMajor, long>& in, const int drow, const int dcol, const int new_rows, const int new_cols, Eigen::SparseMatrix<double, Eigen::ColMajor, long>& out);
//template void sparse_from_csr(const std::vector<int>&, const std::vector<int>&, const std::vector<double>&, const int, const int, Eigen::SparseMatrixBase<Eigen::SparseMatrix<double>>&);
//template void sparse_from_csr(const std::vector<int>&, const std::vector<int>&, const std::vector<double>&, const int, const int, Eigen::SparseMatrixBase<Eigen::SparseMatrix<double, Eigen::RowMajor, int>>&);
template void select_cols(const Eigen::SparseMatrix<double>&, const std::vector<int>&, Eigen::SparseMatrix<double>&);

#ifdef BOW_COMPILE_2D
template void cofactor(const Eigen::MatrixBase<Matrix<double, 2, 2>>& A, Eigen::MatrixBase<Matrix<double, 2, 2>>& coA);
#endif
#ifdef BOW_COMPILE_3D
template void cofactor(const Eigen::MatrixBase<Matrix<double, 3, 3>>& A, Eigen::MatrixBase<Matrix<double, 3, 3>>& coA);
#endif
#endif
#endif

}} // namespace Bow::Math
