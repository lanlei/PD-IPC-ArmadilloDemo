#pragma once
#ifndef MATRIX_CORE_H
#define MATRIX_CORE_H
#ifndef  EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
//#define EIGEN_VECTORIZE_SSE4_2
#endif // ! EIGEN_USE_MKL_ALL

#include <Eigen/Eigen>
#include <Eigen/PardisoSupport>
#include "DataCore.h"

typedef Eigen::Matrix<qeal, 2, 1, Eigen::ColMajor> Vector2;
typedef Eigen::Matrix<qeal, 1, 2, Eigen::RowMajor> VectorR2;
typedef Eigen::Matrix<qeal, 3, 1, Eigen::ColMajor> Vector3;
typedef Eigen::Matrix<qeal, 1, 3, Eigen::RowMajor> VectorR3;
typedef Eigen::Matrix<qeal, 4, 1, Eigen::ColMajor> Vector4;
typedef Eigen::Matrix<qeal, 1, 4, Eigen::RowMajor> VectorR4;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 1, Eigen::ColMajor> VectorX;
typedef Eigen::Matrix<qeal, 1, Eigen::Dynamic, Eigen::RowMajor> VectorXR;

typedef Eigen::Matrix<qeal, 2, 2, Eigen::ColMajor> Matrix2;
typedef Eigen::Matrix<qeal, 3, 3, Eigen::ColMajor> Matrix3;
typedef Eigen::Matrix<qeal, 4, 4, Eigen::ColMajor> Matrix4;
typedef Eigen::Matrix<qeal, 2, Eigen::Dynamic, Eigen::ColMajor> Matrix2X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 2, Eigen::RowMajor> MatrixXR2;
typedef Eigen::Matrix<qeal, 3, Eigen::Dynamic, Eigen::ColMajor> Matrix3X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixXR3;
typedef Eigen::Matrix<qeal, 4, Eigen::Dynamic, Eigen::ColMajor> Matrix4X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 4, Eigen::ColMajor> MatrixXR4;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixX;

typedef Eigen::Matrix<int, 2, 1, Eigen::ColMajor> Vector2i;
typedef Eigen::Matrix<int, 1, 2, Eigen::RowMajor> VectorR2i;
typedef Eigen::Matrix<int, 3, 1, Eigen::ColMajor> Vector3i;
typedef Eigen::Matrix<int, 1, 3, Eigen::RowMajor> VectorR3i;
typedef Eigen::Matrix<int, 4, 1, Eigen::ColMajor> Vector4i;
typedef Eigen::Matrix<int, 1, 4, Eigen::RowMajor> VectorR4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1, Eigen::ColMajor> VectorXi;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic, Eigen::RowMajor> VectorXRi;

typedef Eigen::Matrix<int, 2, 2, Eigen::ColMajor> Matrix2i;
typedef Eigen::Matrix<int, 3, 3, Eigen::ColMajor> Matrix3i;
typedef Eigen::Matrix<int, 4, 4, Eigen::ColMajor> Matrix4i;
typedef Eigen::Matrix<int, 2, Eigen::Dynamic, Eigen::ColMajor> Matrix2Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor> MatrixXR2i;
typedef Eigen::Matrix<int, 3, Eigen::Dynamic, Eigen::ColMajor> Matrix3Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixXR3i;
typedef Eigen::Matrix<int, 4, Eigen::Dynamic, Eigen::ColMajor> Matrix4Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> MatrixXR4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXi;

typedef Eigen::SparseMatrix<qeal/*, Eigen::RowMajor*/, Eigen::ColMajor> SparseMatrix;
typedef Eigen::Triplet<qeal> TripletX;
typedef Eigen::Triplet<int> TripletXi;



#endif