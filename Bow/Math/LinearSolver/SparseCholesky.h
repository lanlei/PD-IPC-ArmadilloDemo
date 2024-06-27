#ifdef BOW_SUITESPARSE
#ifndef BOW_SPARSE_CHOLESKY_H
#define BOW_SPARSE_CHOLESKY_H
#include <Bow/Macros.h>
#include <Eigen/Core>
#include <Bow/Types.h>
#include <cholmod.h>

namespace Bow {
namespace Math {
namespace LinearSolver {
template <class Derived>
class CholmodLLT {
protected:
    cholmod_common cm;
    cholmod_sparse* A = NULL;
    cholmod_factor* L = NULL;

public:
    using StorageIndex = typename Derived::StorageIndex;
    CholmodLLT(const int supernodal = CHOLMOD_SUPERNODAL); //CHOLMOD_SUPERNODAL or CHOLMOD_SIMPLICIAL or CHOLMOD_AUTO
    ~CholmodLLT();
    bool compute(const Eigen::SparseMatrixBase<Derived>& mat);
    template <class DerivedB>
    Bow::Vector<typename Derived::Scalar, Eigen::Dynamic> solve(const Eigen::MatrixBase<DerivedB>& rhs) const;
};
}
}
} // namespace Bow::Math::LinearSolver

#ifndef BOW_STATIC_LIBRARY
#include "SparseCholesky.cpp"
#endif

#endif
#endif