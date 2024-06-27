#ifndef FEM_DEFORMATION_GRADIENT
#define FEM_DEFORMATION_GRADIENT

#include <Bow/Macros.h>
#include <Bow/Types.h>

namespace Bow {
namespace FEM {

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, const Matrix<T, 2, 2>& IB, Matrix<T, 2, 2>& F);

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Matrix<T, 3, 3>& IB, Matrix<T, 3, 3>& F);


template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3, const Matrix<T, 3, 3>& IB, Matrix<T, 3, 3>& F);

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, const Vector<T, 2>& x3, const Matrix<T, 2, 2>& IB, Matrix<T, 2, 2>& F);

template <class DerivedB, class DerivedP, class DerivedG>
BOW_INLINE void backpropagate_element_gradient(const Eigen::MatrixBase<DerivedB>& IB, const Eigen::MatrixBase<DerivedP>& de_dF, Eigen::MatrixBase<DerivedG>& de_dX);

template <class DerivedB, class DerivedP, class DerivedH>
BOW_INLINE void backpropagate_element_hessian(const Eigen::MatrixBase<DerivedB>& IB, const Eigen::MatrixBase<DerivedP>& d2e_dF2, Eigen::MatrixBase<DerivedH>& d2e_dX2);
}
} // namespace Bow::FEM

#ifndef BOW_STATIC_LIBRARY
#include "Utils.cpp"
#endif

#endif