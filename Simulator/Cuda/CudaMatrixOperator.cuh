#ifndef CUDA_MATRIX_OPERATOR_CUH
#define CUDA_MATRIX_OPERATOR_CUH

#include "CudaHeader.cuh"

__device__ __forceinline__
void setVectorValue(qeal* v, qeal val, int dim);

__device__ __forceinline__
void scaleVector(qeal* v, qeal s, int dim = 3);

__device__ __forceinline__
qeal getVectorNorm(qeal* v, int dim = 3);

__device__ __forceinline__
qeal getVectorSqNorm(qeal* v, int dim = 3);

__device__ __forceinline__
qeal getVectorDot(qeal* v1, qeal* v2, int dim = 3);

__device__ __forceinline__
void getVectorNormalize(qeal* v, int dim = 3);

__device__ __forceinline__
void getVectorNormalized(qeal* v, qeal* result, int dim = 3);

__device__ __forceinline__
void getVector3Cross(qeal* v1, qeal* v2, qeal* result);

__device__ __forceinline__
void getVectorAdd(qeal* v1, qeal* v2, qeal* result, const qeal s1 = 1, const qeal s2 = 1, int dim = 3);

__device__ __forceinline__
void getVectorSub(qeal* v1, qeal* v2, qeal* result, const qeal s1 = 1, const qeal s2 = 1, int dim = 3);

__device__ __forceinline__
void getMutilMatrix(qeal* lmat, qeal* rmat, qeal* result, int dim = 3, int rdim = 3, int cdim = 3);

//un-safe & hard code
__device__ __forceinline__
qeal getMatrix3Determinant(qeal* mat3);

//un-safe & hard code
__device__ __forceinline__
void getMatrix3Inverse(qeal* mat3, qeal* inv);

//un-safe & hard code
__device__ __forceinline__
qeal getMatrix4Determinant(qeal* mat4);

//un-safe & hard code
__device__ __forceinline__
void getMatrix4Inverse(qeal* mat4, qeal* inv);

__device__ __forceinline__
void getMutilVVT(qeal* v, qeal* mat, int dim = 3);

__device__ __forceinline__
void getMutilMV(qeal*mat, qeal* v, qeal* result, int rdim = 3, int cdim = 3);

__device__ __forceinline__
void getMatrix3EigenValue(qeal* mat, qeal* result);

__device__ __forceinline__
void addMatrix3(qeal* mat0, qeal* mat1, qeal* result, qeal S0 = 1.0, qeal S1 = 1.0);

__device__ __forceinline__
bool getMatrix3EigenvalueAndEigenvector(const int dim, qeal* matrix, qeal* eigenvalue, qeal* eigenvectors, qeal tol = 1.e-12, int maxIter = 100);


#endif