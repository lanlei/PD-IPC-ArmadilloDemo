#pragma once
#ifndef CUDA_HEADER_H
#define CUDA_HEADER_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "math.h" // CUDA math library
#include "DataCore.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/count.h>

#include "cublas.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#include "cusparse.h"
#include "cusolverSp.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"
#include "cuda_runtime.h"


#define CUDA_CALL(x) { const cudaError_t a = (x); if(a != cudaSuccess){ printf("\nCUDA Error: %s (err_num = %d) \n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}

#define MAX_THREADS_NUM 1024
#define MAX_BLOCKS_NUM 1024
#define THREADS_NUM_768 768
#define THREADS_NUM_576 576
#define THREADS_NUM 512
#define THREADS_NUM_256 256
#define THREADS_NUM_192 192
#define THREADS_NUM_128 128
#define THREADS_NUM_64 64
#define THREADS_NUM_32 32
#define THREADS_NUM_16 16

#define MAX_CELLS_NUM 10000
#define MAX_CELLS_CONTSTRAINT_TRIANGLE_NUM 500

#ifdef USE_DOUBLE_PRECISION
#define CUDA_ZERO 5e-13
#define IS_Numerical_ZERO(d) (CUDA_ABS(d) < CUDA_ZERO)
#define CUDA_SQRT(d) sqrt(d)
#define CUDA_ABS(d) fabs(d)
#define CUDA_LOG(d) log(d)
#define IS_CUDA_ZERO(d) (fabs(d) < CUDA_ZERO)
#define Check_CUDA_ZERO(d) (fabs(d) > CUDA_ZERO ? d : 0)
#define CUDA_QEAL_MAX(a, b) fmax(a, b)
#define CUDA_QEAL_MIN(a, b) fmin(a, b)
#define CUDA_QEAL_MAX3(a, b, c) fmax(fmax(a, b), c)
#define CUDA_QEAL_MIN3(a, b, c) fmin(fmin(a, b), c)
#else
#define CUDA_ZERO 1e-6
#define IS_Numerical_ZERO(d) (CUDA_ABS(d) < CUDA_ZERO)
#define CUDA_SQRT(d) sqrtf(d)
#define CUDA_ABS(d) fabsf(d)
#define CUDA_LOG(d) logf(d)
#define IS_CUDA_ZERO(d) (fabsf(d) < CUDA_ZERO)
#define Check_CUDA_ZERO(d) (fabsf(d) > CUDA_ZERO ? d : 0)
#define CUDA_QEAL_MAX(a, b) fmaxf(a, b)
#define CUDA_QEAL_MIN(a, b) fminf(a, b)
#define CUDA_QEAL_MAX3(a, b, c) fmaxf(fmaxf(a, b), c)
#define CUDA_QEAL_MIN3(a, b, c) fminf(fminf(a, b), c)
#endif


template<class T> inline
void SwapPtr(T &a, T &b)
{
	T c = a; a = b; b = c;
}


#endif