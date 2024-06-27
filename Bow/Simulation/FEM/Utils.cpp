#include "Utils.h"
#include <oneapi/tbb.h>

namespace Bow {
namespace FEM {

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, const Matrix<T, 2, 2>& IB, Matrix<T, 2, 2>& F)
{
    F.col(0) = x1 - x0;
    F.col(1) = x2 - x0;
    F *= IB;
}

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Matrix<T, 3, 3>& IB, Matrix<T, 3, 3>& F)
{

}

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 3>& x0, const Vector<T, 3>& x1, const Vector<T, 3>& x2, const Vector<T, 3>& x3, const Matrix<T, 3, 3>& IB, Matrix<T, 3, 3>& F)
{
    F.col(0) = x1 - x0;
    F.col(1) = x2 - x0;
    F.col(2) = x3 - x0;
    F *= IB;
}

template <class T>
BOW_INLINE void deformation_gradient(const Vector<T, 2>& x0, const Vector<T, 2>& x1, const Vector<T, 2>& x2, const Vector<T, 2>& x3, const Matrix<T, 2, 2>& IB, Matrix<T, 2, 2>& F)
{

}

template <class DerivedB, class DerivedP, class DerivedG>
BOW_INLINE void backpropagate_element_gradient(const Eigen::MatrixBase<DerivedB>& IB, const Eigen::MatrixBase<DerivedP>& de_dF, Eigen::MatrixBase<DerivedG>& de_dX)
{
    using T = typename DerivedB::Scalar;
    const int dim = DerivedB::RowsAtCompileTime;
    de_dX.setZero();
    if constexpr (dim == 2) {
        de_dX(1 * 2 + 0) = de_dF(0, 0) * IB(0, 0) + de_dF(0, 1) * IB(0, 1);
        de_dX(1 * 2 + 1) = de_dF(1, 0) * IB(0, 0) + de_dF(1, 1) * IB(0, 1);
        de_dX(2 * 2 + 0) = de_dF(0, 0) * IB(1, 0) + de_dF(0, 1) * IB(1, 1);
        de_dX(2 * 2 + 1) = de_dF(1, 0) * IB(1, 0) + de_dF(1, 1) * IB(1, 1);
        de_dX(0 * 2 + 0) = -de_dF(0, 0) * IB(0, 0) - de_dF(0, 1) * IB(0, 1) - de_dF(0, 0) * IB(1, 0) - de_dF(0, 1) * IB(1, 1);
        de_dX(0 * 2 + 1) = -de_dF(1, 0) * IB(0, 0) - de_dF(1, 1) * IB(0, 1) - de_dF(1, 0) * IB(1, 0) - de_dF(1, 1) * IB(1, 1);
    }
    else {
        T R10 = IB(0, 0) * de_dF(0, 0) + IB(0, 1) * de_dF(0, 1) + IB(0, 2) * de_dF(0, 2);
        T R11 = IB(0, 0) * de_dF(1, 0) + IB(0, 1) * de_dF(1, 1) + IB(0, 2) * de_dF(1, 2);
        T R12 = IB(0, 0) * de_dF(2, 0) + IB(0, 1) * de_dF(2, 1) + IB(0, 2) * de_dF(2, 2);
        T R20 = IB(1, 0) * de_dF(0, 0) + IB(1, 1) * de_dF(0, 1) + IB(1, 2) * de_dF(0, 2);
        T R21 = IB(1, 0) * de_dF(1, 0) + IB(1, 1) * de_dF(1, 1) + IB(1, 2) * de_dF(1, 2);
        T R22 = IB(1, 0) * de_dF(2, 0) + IB(1, 1) * de_dF(2, 1) + IB(1, 2) * de_dF(2, 2);
        T R30 = IB(2, 0) * de_dF(0, 0) + IB(2, 1) * de_dF(0, 1) + IB(2, 2) * de_dF(0, 2);
        T R31 = IB(2, 0) * de_dF(1, 0) + IB(2, 1) * de_dF(1, 1) + IB(2, 2) * de_dF(1, 2);
        T R32 = IB(2, 0) * de_dF(2, 0) + IB(2, 1) * de_dF(2, 1) + IB(2, 2) * de_dF(2, 2);
        de_dX(1 * 3 + 0) = R10;
        de_dX(1 * 3 + 1) = R11;
        de_dX(1 * 3 + 2) = R12;
        de_dX(2 * 3 + 0) = R20;
        de_dX(2 * 3 + 1) = R21;
        de_dX(2 * 3 + 2) = R22;
        de_dX(3 * 3 + 0) = R30;
        de_dX(3 * 3 + 1) = R31;
        de_dX(3 * 3 + 2) = R32;
        de_dX(0 * 3 + 0) = -R10 - R20 - R30;
        de_dX(0 * 3 + 1) = -R11 - R21 - R31;
        de_dX(0 * 3 + 2) = -R12 - R22 - R32;
    }
}

template <class DerivedB, class DerivedP, class DerivedH>
BOW_INLINE void backpropagate_element_hessian(const Eigen::MatrixBase<DerivedB>& IB, const Eigen::MatrixBase<DerivedP>& d2e_dF2, Eigen::MatrixBase<DerivedH>& d2e_dX2)
{
    using T = typename DerivedB::Scalar;
    const int dim = DerivedB::RowsAtCompileTime;
    d2e_dX2.setZero();
    Matrix<T, dim*(dim + 1), dim * dim> intermediate;
    intermediate.setZero();
    if constexpr (dim == 2) {
        for (int colI = 0; colI < 4; ++colI) {
            T _000 = d2e_dF2(0, colI) * IB(0, 0);
            T _010 = d2e_dF2(0, colI) * IB(1, 0);
            T _101 = d2e_dF2(2, colI) * IB(0, 1);
            T _111 = d2e_dF2(2, colI) * IB(1, 1);
            T _200 = d2e_dF2(1, colI) * IB(0, 0);
            T _210 = d2e_dF2(1, colI) * IB(1, 0);
            T _301 = d2e_dF2(3, colI) * IB(0, 1);
            T _311 = d2e_dF2(3, colI) * IB(1, 1);
            intermediate(2, colI) = _000 + _101;
            intermediate(3, colI) = _200 + _301;
            intermediate(4, colI) = _010 + _111;
            intermediate(5, colI) = _210 + _311;
            intermediate(0, colI) = -intermediate(2, colI) - intermediate(4, colI);
            intermediate(1, colI) = -intermediate(3, colI) - intermediate(5, colI);
        }
        for (int colI = 0; colI < 6; ++colI) {
            T _000 = intermediate(colI, 0) * IB(0, 0);
            T _010 = intermediate(colI, 0) * IB(1, 0);
            T _101 = intermediate(colI, 2) * IB(0, 1);
            T _111 = intermediate(colI, 2) * IB(1, 1);
            T _200 = intermediate(colI, 1) * IB(0, 0);
            T _210 = intermediate(colI, 1) * IB(1, 0);
            T _301 = intermediate(colI, 3) * IB(0, 1);
            T _311 = intermediate(colI, 3) * IB(1, 1);
            d2e_dX2(2, colI) = _000 + _101;
            d2e_dX2(3, colI) = _200 + _301;
            d2e_dX2(4, colI) = _010 + _111;
            d2e_dX2(5, colI) = _210 + _311;
            d2e_dX2(0, colI) = -_000 - _101 - _010 - _111;
            d2e_dX2(1, colI) = -_200 - _301 - _210 - _311;
        }
    }
    else {
        for (int colI = 0; colI < 9; ++colI) {
            intermediate(3, colI) = IB(0, 0) * d2e_dF2(0, colI) + IB(0, 1) * d2e_dF2(3, colI) + IB(0, 2) * d2e_dF2(6, colI);
            intermediate(4, colI) = IB(0, 0) * d2e_dF2(1, colI) + IB(0, 1) * d2e_dF2(4, colI) + IB(0, 2) * d2e_dF2(7, colI);
            intermediate(5, colI) = IB(0, 0) * d2e_dF2(2, colI) + IB(0, 1) * d2e_dF2(5, colI) + IB(0, 2) * d2e_dF2(8, colI);
            intermediate(6, colI) = IB(1, 0) * d2e_dF2(0, colI) + IB(1, 1) * d2e_dF2(3, colI) + IB(1, 2) * d2e_dF2(6, colI);
            intermediate(7, colI) = IB(1, 0) * d2e_dF2(1, colI) + IB(1, 1) * d2e_dF2(4, colI) + IB(1, 2) * d2e_dF2(7, colI);
            intermediate(8, colI) = IB(1, 0) * d2e_dF2(2, colI) + IB(1, 1) * d2e_dF2(5, colI) + IB(1, 2) * d2e_dF2(8, colI);
            intermediate(9, colI) = IB(2, 0) * d2e_dF2(0, colI) + IB(2, 1) * d2e_dF2(3, colI) + IB(2, 2) * d2e_dF2(6, colI);
            intermediate(10, colI) = IB(2, 0) * d2e_dF2(1, colI) + IB(2, 1) * d2e_dF2(4, colI) + IB(2, 2) * d2e_dF2(7, colI);
            intermediate(11, colI) = IB(2, 0) * d2e_dF2(2, colI) + IB(2, 1) * d2e_dF2(5, colI) + IB(2, 2) * d2e_dF2(8, colI);
            intermediate(0, colI) = -intermediate(3, colI) - intermediate(6, colI) - intermediate(9, colI);
            intermediate(1, colI) = -intermediate(4, colI) - intermediate(7, colI) - intermediate(10, colI);
            intermediate(2, colI) = -intermediate(5, colI) - intermediate(8, colI) - intermediate(11, colI);
        }
        for (int rowI = 0; rowI < 12; ++rowI) {
            T _000 = IB(0, 0) * intermediate(rowI, 0);
            T _013 = IB(0, 1) * intermediate(rowI, 3);
            T _026 = IB(0, 2) * intermediate(rowI, 6);
            T _001 = IB(0, 0) * intermediate(rowI, 1);
            T _014 = IB(0, 1) * intermediate(rowI, 4);
            T _027 = IB(0, 2) * intermediate(rowI, 7);
            T _002 = IB(0, 0) * intermediate(rowI, 2);
            T _015 = IB(0, 1) * intermediate(rowI, 5);
            T _028 = IB(0, 2) * intermediate(rowI, 8);
            T _100 = IB(1, 0) * intermediate(rowI, 0);
            T _113 = IB(1, 1) * intermediate(rowI, 3);
            T _126 = IB(1, 2) * intermediate(rowI, 6);
            T _101 = IB(1, 0) * intermediate(rowI, 1);
            T _114 = IB(1, 1) * intermediate(rowI, 4);
            T _127 = IB(1, 2) * intermediate(rowI, 7);
            T _102 = IB(1, 0) * intermediate(rowI, 2);
            T _115 = IB(1, 1) * intermediate(rowI, 5);
            T _128 = IB(1, 2) * intermediate(rowI, 8);
            T _200 = IB(2, 0) * intermediate(rowI, 0);
            T _213 = IB(2, 1) * intermediate(rowI, 3);
            T _226 = IB(2, 2) * intermediate(rowI, 6);
            T _201 = IB(2, 0) * intermediate(rowI, 1);
            T _214 = IB(2, 1) * intermediate(rowI, 4);
            T _227 = IB(2, 2) * intermediate(rowI, 7);
            T _202 = IB(2, 0) * intermediate(rowI, 2);
            T _215 = IB(2, 1) * intermediate(rowI, 5);
            T _228 = IB(2, 2) * intermediate(rowI, 8);
            d2e_dX2(rowI, 3) = _000 + _013 + _026; // 1
            d2e_dX2(rowI, 4) = _001 + _014 + _027; // 2
            d2e_dX2(rowI, 5) = _002 + _015 + _028; // 3
            d2e_dX2(rowI, 6) = _100 + _113 + _126; // 1
            d2e_dX2(rowI, 7) = _101 + _114 + _127; // 2
            d2e_dX2(rowI, 8) = _102 + _115 + _128; // 3
            d2e_dX2(rowI, 9) = _200 + _213 + _226; // 1
            d2e_dX2(rowI, 10) = _201 + _214 + _227; // 2
            d2e_dX2(rowI, 11) = _202 + _215 + _228; // 3
            d2e_dX2(rowI, 0) = -_200 - _213 - _226 - _100 - _113 - _126 - _000 - _013 - _026;
            d2e_dX2(rowI, 1) = -_001 - _014 - _027 - _101 - _114 - _127 - _201 - _214 - _227;
            d2e_dX2(rowI, 2) = -_002 - _015 - _028 - _102 - _115 - _128 - _202 - _215 - _228;
        }
    }
}

#ifdef BOW_STATIC_LIBRARY
#ifdef BOW_COMPILE_FLOAT
#ifdef BOW_COMPILE_2D
template void deformation_gradient(const Vector<float, 2>& x0, const Vector<float, 2>& x1, const Vector<float, 2>& x2, const Matrix<float, 2, 2>& IB, Matrix<float, 2, 2>& F);
template void backpropagate_element_gradient(const Eigen::MatrixBase<Matrix<float, 2, 2>>& IB, const Eigen::MatrixBase<Matrix<float, 2, 2>>& de_dF, Eigen::MatrixBase<Vector<float, 6>>& de_dX);
template void backpropagate_element_hessian(const Eigen::MatrixBase<Matrix<float, 2, 2>>& IB, const Eigen::MatrixBase<Matrix<float, 4, 4>>& d2e_dF2, Eigen::MatrixBase<Matrix<float, 6, 6>>& d2e_dX2);
#endif
#ifdef BOW_COMPILE_3D
template void deformation_gradient(const Vector<float, 3>& x0, const Vector<float, 3>& x1, const Vector<float, 3>& x2, const Vector<float, 3>& x3, const Matrix<float, 3, 3>& IB, Matrix<float, 3, 3>& F);
template void backpropagate_element_gradient(const Eigen::MatrixBase<Matrix<float, 3, 3>>& IB, const Eigen::MatrixBase<Matrix<float, 3, 3>>& de_dF, Eigen::MatrixBase<Vector<float, 12>>& de_dX);
template void backpropagate_element_hessian(const Eigen::MatrixBase<Matrix<float, 3, 3>>& IB, const Eigen::MatrixBase<Matrix<float, 9, 9>>& d2e_dF2, Eigen::MatrixBase<Matrix<float, 12, 12>>& d2e_dX2);
#endif
#endif
#ifdef BOW_COMPILE_DOUBLE
#ifdef BOW_COMPILE_2D
template void deformation_gradient(const Vector<double, 2>& x0, const Vector<double, 2>& x1, const Vector<double, 2>& x2, const Matrix<double, 2, 2>& IB, Matrix<double, 2, 2>& F);
template void backpropagate_element_gradient(const Eigen::MatrixBase<Matrix<double, 2, 2>>& IB, const Eigen::MatrixBase<Matrix<double, 2, 2>>& de_dF, Eigen::MatrixBase<Vector<double, 6>>& de_dX);
template void backpropagate_element_hessian(const Eigen::MatrixBase<Matrix<double, 2, 2>>& IB, const Eigen::MatrixBase<Matrix<double, 4, 4>>& d2e_dF2, Eigen::MatrixBase<Matrix<double, 6, 6>>& d2e_dX2);
#endif
#ifdef BOW_COMPILE_3D
template void deformation_gradient(const Vector<double, 3>& x0, const Vector<double, 3>& x1, const Vector<double, 3>& x2, const Vector<double, 3>& x3, const Matrix<double, 3, 3>& IB, Matrix<double, 3, 3>& F);
template void backpropagate_element_gradient(const Eigen::MatrixBase<Matrix<double, 3, 3>>& IB, const Eigen::MatrixBase<Matrix<double, 3, 3>>& de_dF, Eigen::MatrixBase<Vector<double, 12>>& de_dX);
template void backpropagate_element_hessian(const Eigen::MatrixBase<Matrix<double, 3, 3>>& IB, const Eigen::MatrixBase<Matrix<double, 9, 9>>& d2e_dF2, Eigen::MatrixBase<Matrix<double, 12, 12>>& d2e_dX2);
#endif
#endif
#endif
}
} // namespace Bow::FEM
