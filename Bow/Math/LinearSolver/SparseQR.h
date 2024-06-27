#ifdef BOW_SUITESPARSE
#ifndef BOW_SPARSE_QR_H
#define BOW_SPARSE_QR_H
#include <Bow/Macros.h>
#include <Eigen/Core>
#include <Bow/Types.h>
#include <cholmod.h>
#include <SuiteSparseQR.hpp>

namespace Bow {
namespace Math {
namespace LinearSolver {
template <class Derived>
class SparseQR {
protected:
    cholmod_common cm;
    cholmod_sparse* A = NULL;
    SuiteSparseQR_factorization<double>* QR = NULL;

public:
    using StorageIndex = typename Derived::StorageIndex;
    SparseQR();
    ~SparseQR();
    bool compute(const Eigen::SparseMatrixBase<Derived>& mat);
    template <class DerivedB>
    Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> solve(const Eigen::MatrixBase<DerivedB>& rhs) const;
};
}
}
} // namespace Bow::Math::LinearSolver

#ifndef BOW_STATIC_LIBRARY
#include "SparseQR.cpp"
#endif

#endif
#endif