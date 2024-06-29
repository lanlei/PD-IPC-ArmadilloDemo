#include "PdIpcGpuFunc.cuh"
#include "Simulator\Cuda\CudaSVD.cu"
#include "Simulator\Cuda\CudaMatrixOperator.cu"

namespace PD_IPC
{
	__global__ void laplacianDampingKernel
	(
		const int simPointsNum,
		const qeal rate,
		const qeal* V,
		qeal* nextV,
		const int* pointsNeigborsList,
		const int* pointsNeigborsNum,
		const int* pointsNeigborsOffset,
		const int* positionConstraintFlag
	)
	{
		__shared__ qeal sharedNV[1536];
		__shared__ qeal sharedV[1536];
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		int localIndex = 3 * threadIdx.x;

		qeal* v = sharedV + localIndex;
		qeal* nv = sharedNV + localIndex;
		nv[0] = 0;
		nv[1] = 0;
		nv[2] = 0;

		if (tid < simPointsNum)
		{
			v[0] = V[3 * tid];
			v[1] = V[3 * tid + 1];
			v[2] = V[3 * tid + 2];
		}
		__syncthreads();

		if (tid < simPointsNum)
		{
			if (positionConstraintFlag[tid] != -1)	return;
			int num = pointsNeigborsNum[tid];
			int offset = pointsNeigborsOffset[tid];
			int j;
			for (int i = 1; i < num; i++)
			{
				j = pointsNeigborsList[offset + i];
				nv[0] += V[3 * j] - v[0];
				nv[1] += V[3 * j + 1] - v[1];
				nv[2] += V[3 * j + 2] - v[2];
			}
			nextV[3 * tid] = v[0] + nv[0] * rate;
			nextV[3 * tid + 1] = v[1] + nv[1] * rate;
			nextV[3 * tid + 2] = v[2] + nv[2] * rate;
		}

	}


	__host__ void laplacianDampingHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simPointsNum,
		const qeal rate,
		const qeal* V,
		qeal* nextV,
		const int* pointsNeigborsList,
		const int* pointsNeigborsNum,
		const int* pointsNeigborsOffset,
		const int* positionConstraintFlag
	)
	{
		laplacianDampingKernel << <gridSize, blockSize >> >
			(
				simPointsNum,
				rate,
				V, 
				nextV, 
				pointsNeigborsList, 
				pointsNeigborsNum, 
				pointsNeigborsOffset, 
				positionConstraintFlag
				);
		cudaDeviceSynchronize();
	}

	///
	__global__ void predictiveKernel
	(
		const int simPointsNum,
		const qeal dt,
		qeal* X,
		qeal* V,
		qeal* initRhs,
		qeal* Diagonal,
		qeal* gravityForce,
		qeal* mouseForce,
		qeal* penaltyForce,
		qeal* invMass,
		const int* positionConstraintFlag
	)
	{
		__shared__ qeal vx[THREADS_NUM];
		__shared__ qeal vy[THREADS_NUM];
		__shared__ qeal vz[THREADS_NUM];
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		if (tid >= simPointsNum)	return;

		int vidx = 3 * tid;
		int vidy = vidx + 1;
		int vidz = vidx + 2;
		if (positionConstraintFlag[tid] != -1)
		{
			V[vidx] = 0;
			V[vidy] = 0;
			V[vidz] = 0;
			return;
		}
		qeal m = invMass[tid] * dt;

		vx[threadIdx.x] = V[vidx];
		vy[threadIdx.x] = V[vidy];
		vz[threadIdx.x] = V[vidz];

		vx[threadIdx.x] += m * (gravityForce[vidx] + mouseForce[vidx] + penaltyForce[vidx]);
		vy[threadIdx.x] += m * (gravityForce[vidy] + mouseForce[vidy] + penaltyForce[vidy]);
		vz[threadIdx.x] += m * (gravityForce[vidz] + mouseForce[vidz] + penaltyForce[vidz]);

		vx[threadIdx.x] *= dt;
		vy[threadIdx.x] *= dt;
		vz[threadIdx.x] *= dt;

		vx[threadIdx.x] += X[vidx];
		vy[threadIdx.x] += X[vidy];
		vz[threadIdx.x] += X[vidz];

		m = 1.0 / (dt * m);
		Diagonal[tid] = 1.0 / (invMass[tid] * dt * dt);

		initRhs[vidx] = m * vx[threadIdx.x];
		initRhs[vidy] = m * vy[threadIdx.x];
		initRhs[vidz] = m * vz[threadIdx.x];

		X[vidx] = vx[threadIdx.x];
		X[vidy] = vy[threadIdx.x];
		X[vidz] = vz[threadIdx.x];
	}


	__host__ void predictiveHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simPointsNum,
		const qeal dt,
		qeal* X,
		qeal* V,
		qeal* initRhs,
		qeal* Diagonal,
		qeal* gravityForce,
		qeal* mouseForce,
		qeal* penaltyForce,
		qeal* invMass,
		const int* positionConstraintFlag
	)
	{
		predictiveKernel << <gridSize, blockSize >> >
			(
				simPointsNum,
				dt,
				X,
				V,
				initRhs,
				Diagonal,
				gravityForce,
				mouseForce,
				penaltyForce,
				invMass,
				positionConstraintFlag
				);
		cudaDeviceSynchronize();
	}


	__global__ void computeDirectionKernel
	(
		const int simDims,
		qeal* devLastPoints,
		qeal* devCurPoints,
		qeal* devDirection
	)
	{
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		for (; tid < simDims; tid += length)
		{
			devDirection[tid] = devCurPoints[tid] - devLastPoints[tid];
		}
	}

	__host__ void computeDirectionHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simDims,
		qeal* devLastPoints,
		qeal* devCurPoints,
		qeal* devDirection
	)
	{
		computeDirectionKernel << <gridSize, blockSize >> >
			(
				simDims,
				devLastPoints,
				devCurPoints,
				devDirection
				);
		cudaDeviceSynchronize();
	}

	__global__ void clampDirectionKernel
	(
		const int simDims,
		qeal* devLastPoints,
		qeal* devDirection,
		qeal* devCurPoints,
		qeal* devToi
	)
	{
		__shared__ qeal toi;
		if (threadIdx.x == 0)
		{
			toi = *devToi;
		}
		__syncthreads();
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		for (; tid < simDims; tid += length)
		{
			devCurPoints[tid] = devLastPoints[tid] + toi * devDirection[tid];
		}
	}

	__host__ void clampDirectionHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simDims,
		qeal* devLastPoints,
		qeal* devDirection,
		qeal* devCurPoints,
		qeal* devToi
	)
	{
		clampDirectionKernel << <gridSize, blockSize >> >
			(
				simDims,
				devLastPoints,
				devDirection,
				devCurPoints,
				devToi
				);
		cudaDeviceSynchronize();
	}

	__global__ void updateVelocityKernel
	(
		const int simDims,
		const qeal damping,
		const qeal dt,
		qeal* X,
		qeal* X0,
		qeal* V
	)
	{
		__shared__ qeal invt;
		if (threadIdx.x == 0)
			invt = 1.0 / dt * damping;
		__syncthreads();

		int i = blockDim.x * blockIdx.x + threadIdx.x;
		if (i >= simDims)	return;
		V[i] = (X[i] - X0[i]) * invt;
	}


	__host__ void updateVelocityHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simDims,
		const qeal damping,
		const qeal dt,
		qeal* X,
		qeal* X0,
		qeal* V
	)
	{
		updateVelocityKernel << <gridSize, blockSize >> >
			(
				simDims,
				damping,
				dt,
				X,
				X0,
				V
				);
		cudaDeviceSynchronize();
	}

__device__ __forceinline__
		void clamp(qeal* n, qeal* result)
	{
		if ((*n) < 1.001)
			*result = (*n);
		else *result = 1.001;

		if ((*result) < 0.999)
			*result = 0.999;
	}

///////////////////////////////////////////////////////////////////////////////////////////
//  math kernels
///////////////////////////////////////////////////////////////////////////////////////////
__device__ void dev_Matrix_Product_3(const float *A, const float *B, float *R)				//R=A*B
{
	R[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
	R[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
	R[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
	R[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
	R[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
	R[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
	R[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
	R[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
	R[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

__device__ void dev_Matrix_Substract_3(float *A, float *B, float *R)						//R=A-B
{
	for (int i = 0; i < 9; i++)	R[i] = A[i] - B[i];
}

__device__ void dev_Matrix_Product(float *A, float *B, float *R, int nx, int ny, int nz)	//R=A*B
{
	memset(R, 0, sizeof(float)*nx*nz);
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < nz; j++)
			for (int k = 0; k < ny; k++)
				R[i*nz + j] += A[i*ny + k] * B[k*nz + j];
}

__device__ void Get_Rotation(float F[3][3], float R[3][3])
{
	float C[3][3];
	memset(&C[0][0], 0, sizeof(float) * 9);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				C[i][j] += F[k][i] * F[k][j];

	float C2[3][3];
	memset(&C2[0][0], 0, sizeof(float) * 9);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				C2[i][j] += C[i][k] * C[j][k];

	float det = F[0][0] * F[1][1] * F[2][2] +
		F[0][1] * F[1][2] * F[2][0] +
		F[1][0] * F[2][1] * F[0][2] -
		F[0][2] * F[1][1] * F[2][0] -
		F[0][1] * F[1][0] * F[2][2] -
		F[0][0] * F[1][2] * F[2][1];

	float I_c = C[0][0] + C[1][1] + C[2][2];
	float I_c2 = I_c * I_c;
	float II_c = 0.5*(I_c2 - C2[0][0] - C2[1][1] - C2[2][2]);
	float III_c = det * det;
	float k = I_c2 - 3 * II_c;

	float inv_U[3][3];
	if (k < 1e-10f)
	{
		float inv_lambda = 1 / sqrt(I_c / 3);
		memset(inv_U, 0, sizeof(float) * 9);
		inv_U[0][0] = inv_lambda;
		inv_U[1][1] = inv_lambda;
		inv_U[2][2] = inv_lambda;
	}
	else
	{
		float l = I_c * (I_c*I_c - 4.5*II_c) + 13.5*III_c;
		float k_root = sqrt(k);
		float value = l / (k*k_root);
		if (value < -1.0) value = -1.0;
		if (value > 1.0) value = 1.0;
		float phi = acos(value);
		float lambda2 = (I_c + 2 * k_root*cos(phi / 3)) / 3.0;
		float lambda = sqrt(lambda2);

		float III_u = sqrt(III_c);
		if (det < 0)   III_u = -III_u;
		float I_u = lambda + sqrt(-lambda2 + I_c + 2 * III_u / lambda);
		float II_u = (I_u*I_u - I_c)*0.5;

		float U[3][3];
		float inv_rate, factor;

		inv_rate = 1 / (I_u*II_u - III_u);
		factor = I_u * III_u*inv_rate;

		memset(U, 0, sizeof(float) * 9);
		U[0][0] = factor;
		U[1][1] = factor;
		U[2][2] = factor;

		factor = (I_u*I_u - II_u)*inv_rate;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				U[i][j] += factor * C[i][j] - inv_rate * C2[i][j];

		inv_rate = 1 / III_u;
		factor = II_u * inv_rate;
		memset(inv_U, 0, sizeof(float) * 9);
		inv_U[0][0] = factor;
		inv_U[1][1] = factor;
		inv_U[2][2] = factor;

		factor = -I_u * inv_rate;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				inv_U[i][j] += factor * U[i][j] + inv_rate * C[i][j];
	}

	memset(&R[0][0], 0, sizeof(float) * 9);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				R[i][j] += F[i][k] * inv_U[k][j];
}


__global__ void solveVolumetricStrainConstraintKernel
	(
		const int constraintsNum,
		int* devSimElements,
		qeal* devSimElementsRestShape,
		qeal* devVolumetricStrainStiffness,
		qeal* devPoints,
		qeal* devRhs,
		qeal* devVolumetricStrainEnergy
	)
	{
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		for (; tid < constraintsNum; tid += length)
		{
			const int offset9 = tid * 9;
			const int offset4 = tid * 4;

			int idx0 = 3 * devSimElements[offset4];
			int idx1 = 3 * devSimElements[offset4 + 1];
			int idx2 = 3 * devSimElements[offset4 + 2];
			int idx3 = 3 * devSimElements[offset4 + 3];

			const float* idm = &devSimElementsRestShape[offset9];

			qeal x0 = devPoints[idx0]; qeal x1 = devPoints[idx1 + 0]; qeal x2 = devPoints[idx2 + 0]; qeal x3 = devPoints[idx3 + 0];
			qeal y0 = devPoints[idx0 + 1]; qeal y1 = devPoints[idx1 + 1]; qeal y2 = devPoints[idx2 + 1]; qeal y3 = devPoints[idx3 + 1];
			qeal z0 = devPoints[idx0 + 2]; qeal z1 = devPoints[idx1 + 2]; qeal z2 = devPoints[idx2 + 2]; qeal z3 = devPoints[idx3 + 2];

			float Ds[9];
			Ds[0] = x1 - x0;
			Ds[3] = y1 - y0;
			Ds[6] = z1 - z0;
			Ds[1] = x2 - x0;
			Ds[4] = y2 - y0;
			Ds[7] = z2 - z0;
			Ds[2] = x3 - x0;
			Ds[5] = y3 - y0;
			Ds[8] = z3 - z0;

			qeal F[9];
			qeal new_R[9];
			dev_Matrix_Product_3(Ds, idm, F);

			//Get_Rotation((float(*)[3])F, (float(*)[3])new_R);

			qeal U[9];
			qeal S[9];
			qeal V[9];

			CudaSVD::svd(F, U, S, V);
			qeal det_F = S[0] * S[4] * S[8];

			if (det_F < 0)
			{
				S[8] *= -1;
				qeal high_val = S[0];
				qeal mid_val = S[4];
				qeal low_val = S[8];

				if (mid_val < low_val) {
					qeal temp = low_val;
					low_val = mid_val;
					mid_val = temp;
				}

				if (high_val < low_val) {
					qeal temp = low_val;
					low_val = high_val;
					high_val = temp;
				}
				if (high_val < mid_val) {
					qeal temp = mid_val;
					mid_val = high_val;
					high_val = temp;
				}

				S[0] = high_val;
				S[4] = mid_val;
				S[8] = low_val;
			}
			qeal sigma_new[3];

			clamp(S, sigma_new);
			clamp(S + 4, sigma_new + 1);
			clamp(S + 8, sigma_new + 2);

			qeal SVt[9]; // = new_sigma * vt
			SVt[0] = sigma_new[0] * V[0];
			SVt[1] = sigma_new[1] * V[3];
			SVt[2] = sigma_new[2] * V[6];

			SVt[3] = sigma_new[0] * V[1];
			SVt[4] = sigma_new[1] * V[4];
			SVt[5] = sigma_new[2] * V[7];

			SVt[6] = sigma_new[0] * V[2];
			SVt[7] = sigma_new[1] * V[5];
			SVt[8] = sigma_new[2] * V[8];
			//
			// u * s * vt
			new_R[0] = (U[0] * SVt[0] + U[3] * SVt[1] + U[6] * SVt[2]);
			new_R[1] = (U[1] * SVt[0] + U[4] * SVt[1] + U[7] * SVt[2]);
			new_R[2] = (U[2] * SVt[0] + U[5] * SVt[1] + U[8] * SVt[2]);

			new_R[3] = (U[0] * SVt[3] + U[3] * SVt[4] + U[6] * SVt[5]);
			new_R[4] = (U[1] * SVt[3] + U[4] * SVt[4] + U[7] * SVt[5]);
			new_R[5] = (U[2] * SVt[3] + U[5] * SVt[4] + U[8] * SVt[5]);

			new_R[6] = (U[0] * SVt[6] + U[3] * SVt[7] + U[6] * SVt[8]);
			new_R[7] = (U[1] * SVt[6] + U[4] * SVt[7] + U[7] * SVt[8]);
			new_R[8] = (U[2] * SVt[6] + U[5] * SVt[7] + U[8] * SVt[8]);


			float half_matrix[3][4], result_matrix[3][4];
			half_matrix[0][0] = -idm[0] - idm[3] - idm[6];
			half_matrix[0][1] = idm[0];
			half_matrix[0][2] = idm[3];
			half_matrix[0][3] = idm[6];
			half_matrix[1][0] = -idm[1] - idm[4] - idm[7];
			half_matrix[1][1] = idm[1];
			half_matrix[1][2] = idm[4];
			half_matrix[1][3] = idm[7];
			half_matrix[2][0] = -idm[2] - idm[5] - idm[8];
			half_matrix[2][1] = idm[2];
			half_matrix[2][2] = idm[5];
			half_matrix[2][3] = idm[8];

			//dev_Matrix_Substract_3(new_R, F, new_R);
			dev_Matrix_Product(new_R, &half_matrix[0][0], &result_matrix[0][0], 3, 3, 4);

			qeal stiffness = devVolumetricStrainStiffness[tid];

			//qeal e = 0;
			//for (int k = 0; k < 9; k++)
			//	e += (new_R[k] - F[k]) *  (new_R[k] - F[k]);
			//devVolumetricStrainEnergy[tid] = 0.5 * stiffness * e;
	
			for (int r = 0; r < 3; r++)
				for (int c = 0; c < 4; c++)
					result_matrix[r][c] *= stiffness;

			atomicAdd(devRhs + idx0, result_matrix[0][0]);
			atomicAdd(devRhs + idx0 + 1, result_matrix[1][0]);
			atomicAdd(devRhs + idx0 + 2, result_matrix[2][0]);

			atomicAdd(devRhs + idx1, result_matrix[0][1]);
			atomicAdd(devRhs + idx1 + 1, result_matrix[1][1]);
			atomicAdd(devRhs + idx1 + 2, result_matrix[2][1]);

			atomicAdd(devRhs + idx2, result_matrix[0][2]);
			atomicAdd(devRhs + idx2 + 1, result_matrix[1][2]);
			atomicAdd(devRhs + idx2 + 2, result_matrix[2][2]);

			atomicAdd(devRhs + idx3, result_matrix[0][3]);
			atomicAdd(devRhs + idx3 + 1,result_matrix[1][3]);
			atomicAdd(devRhs + idx3 + 2, result_matrix[2][3]);
		}

	}

	__host__ void solveVolumetricStrainConstraintHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int constraintsNum,
		int* devSimElements,
		qeal* devSimElementsRestShape,
		qeal* devVolumetricStrainStiffness,
		qeal* devPoints,
		qeal* devRhs,
		qeal* devVolumetricStrainEnergy
	)
	{
		solveVolumetricStrainConstraintKernel << <gridSize, blockSize >> >
			(
				constraintsNum,
				devSimElements,
				devSimElementsRestShape,
				devVolumetricStrainStiffness,
				devPoints,
				devRhs,
				devVolumetricStrainEnergy
				);
		cudaDeviceSynchronize();
	}

	__global__ void solveClothPointsStrainAndBendingConstraintKernel
	(
		const int simClothPointsNum,
		const qeal* devStrainStiffness,
		const qeal* bendingStiffness,
		const int* devClothPointsNeighborsList,
		const int* devClothPointsNeighborsNum,
		const int* devClothPointsNeighborsOffset,
		const qeal* devClothPointsNeighborsRestLength,
		const qeal* devClothPointsNeighborsLbo,
		const qeal* dveClothPointsMeanCurvatureNorm,
		qeal* devPoints,
		qeal* devRhs
	)
	{
		__shared__ qeal sharedX[1536];
		__shared__ qeal sharedNX[1536];
		__shared__ qeal sharedTX[1536];
		__shared__ qeal sharedAq[1536];

		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		int localIndex = 3 * threadIdx.x;

		qeal* x = sharedX + localIndex;
		qeal* nx = sharedNX + localIndex;
		qeal* tx = sharedTX + localIndex;
		qeal* aq = sharedAq + localIndex;

		if (tid < simClothPointsNum)
		{	
			int num = devClothPointsNeighborsNum[tid];
			int offset = devClothPointsNeighborsOffset[tid];
			
			int pointId = 3 * devClothPointsNeighborsList[offset];
			x[0] = devPoints[pointId];
			x[1] = devPoints[pointId +1];
			x[2] = devPoints[pointId + 2];

			qeal curLen, restLen, rc;
			tx[0] = 0; tx[1] = 0; tx[2] = 0;

			qeal lbo = devClothPointsNeighborsLbo[offset];
			aq[0] = x[0] * lbo; aq[1] = x[1] * lbo; aq[2] = x[2] * lbo;

			qeal dir[3];
			qeal strain = devStrainStiffness[offset];
			
			for (int i = 1; i < num; i++)
			{
				int neighborId = 3 *devClothPointsNeighborsList[offset + i];
				nx[0] = devPoints[neighborId];
				nx[1] = devPoints[neighborId + 1];
				nx[2] = devPoints[neighborId + 2];

				dir[0] = x[0] - nx[0]; dir[1] = x[1] - nx[1]; dir[2] = x[2] - nx[2];

				curLen = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];

				curLen = CUDA_SQRT(curLen);

				restLen = devClothPointsNeighborsRestLength[offset + i];
				rc = strain * (restLen / curLen);

				tx[0] += rc * dir[0];
				tx[1] += rc * dir[1];
				tx[2] += rc * dir[2];

				lbo = devClothPointsNeighborsLbo[offset + i];
				aq[0] += nx[0] * lbo;
				aq[1] += nx[1] * lbo;
				aq[2] += nx[2] * lbo;
			}

			qeal aqnorm = getVectorNorm(aq);
			qeal cur;

			if (aqnorm >= 1e-10)
			{
				cur = dveClothPointsMeanCurvatureNorm[tid] / aqnorm;
				aq[0] *= cur; aq[1] *= cur; aq[2] *= cur;
			}
			qeal bending = bendingStiffness[offset];
			tx[0] += bending * aq[0];
			tx[1] += bending * aq[1];
			tx[2] += bending * aq[2];

			atomicAdd(devRhs + pointId, tx[0]);
			atomicAdd(devRhs + pointId + 1, tx[1]);
			atomicAdd(devRhs + pointId + 2, tx[2]);

			for (int j = 1; j < num; j++)
			{
				bending = bendingStiffness[offset + j];
				int neighborId = 3 * devClothPointsNeighborsList[offset + j];

				atomicAdd(devRhs + neighborId, bending * aq[0]);
				atomicAdd(devRhs + neighborId + 1, bending * aq[1]);
				atomicAdd(devRhs + neighborId + 2, bending * aq[2]);
			}

		}
	}

	__host__ void solveClothPointsStrainAndBendingConstraintHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simClothPointsNum,
		const qeal* devStrainStiffness,
		const qeal* bendingStiffness,
		const int* devClothPointsNeighborsList,
		const int* devClothPointsNeighborsNum,
		const int* devClothPointsNeighborsOffset,
		const qeal* devClothPointsNeighborsRestLength,
		const qeal* devClothPointsNeighborsLbo,
		const qeal* dveClothPointsMeanCurvatureNorm,
		qeal* devPoints,
		qeal* devRhs
	)
	{
		solveClothPointsStrainAndBendingConstraintKernel << <gridSize, blockSize >> >
			(
				simClothPointsNum,
				devStrainStiffness,
				bendingStiffness,
				devClothPointsNeighborsList,
				devClothPointsNeighborsNum,
				devClothPointsNeighborsOffset,
				devClothPointsNeighborsRestLength,
				devClothPointsNeighborsLbo,
				dveClothPointsMeanCurvatureNorm,
				devPoints,
				devRhs
				);
		cudaDeviceSynchronize();
	}

	__global__ void solvePointsConstraintKernel
	(
		const int pointsConstraintNum,
		qeal* devRestPoints,
		const int* pointConstraintIndex,
		const qeal* positionStiffness,
		qeal* devRhs
	)
	{
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		if (tid < pointsConstraintNum)
		{
			int vid = 3 * pointConstraintIndex[tid];
			qeal stiffness = positionStiffness[tid];

			qeal restPos[3];
			restPos[0] = stiffness * devRestPoints[vid];
			restPos[1] = stiffness * devRestPoints[vid + 1];
			restPos[2] = stiffness * devRestPoints[vid + 2];
			atomicAdd(devRhs + vid, restPos[0]);
			atomicAdd(devRhs + vid + 1, restPos[1]);
			atomicAdd(devRhs + vid + 2, restPos[2]);
		}
	
	}

	__host__ void solvePointsConstraintHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int pointsConstraintNum,
		qeal* devRestPoints,
		const int* pointConstraintIndex,
		const qeal* positionStiffness,
		qeal* devRhs
	)
	{
		solvePointsConstraintKernel << <gridSize, blockSize >> >
			(
				pointsConstraintNum,
				devRestPoints,
				pointConstraintIndex,
				positionStiffness,
				devRhs
				);
		cudaDeviceSynchronize();
	}



	__global__ void solvePointsLocalConstraintsKernel
	(
		const int simPointsNum,
		qeal* devPoints,
		qeal* devRestPoints,
		qeal* devSn,
		const qeal* pointsMass,
		const qeal dt,
		const qeal* positionStiffness,
		const int* positionConstraintFlag,
		qeal* devRhs,
		qeal* inertialEnergy,
		qeal* positionEnergy
	)
	{
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		if (tid < simPointsNum)
		{
			qeal x[3];
			x[0] = devPoints[3 * tid];
			x[1] = devPoints[3 * tid + 1];
			x[2] = devPoints[3 * tid + 2];

			qeal q[3];
			q[0] = x[0] - devSn[3 * tid];
			q[1] = x[1] - devSn[3 * tid + 1];
			q[2] = x[2] - devSn[3 * tid + 2];
			//inertialEnergy[tid] = (0.5 * pointsMass[tid] / (dt * dt)) * getVectorSqNorm(q);

			
			qeal restPos[3];
			restPos[0] = devRestPoints[3 * tid];
			restPos[1] = devRestPoints[3 * tid + 1];
			restPos[2] = devRestPoints[3 * tid + 2];

			int fixed = positionConstraintFlag[tid];
			if (fixed != -1)
			{			
				qeal stiffness = positionStiffness[fixed];
				atomicAdd(devRhs + 3 * tid, stiffness * restPos[0]);
				atomicAdd(devRhs + 3 * tid + 1, stiffness * restPos[1]);
				atomicAdd(devRhs + 3 * tid + 2, stiffness * restPos[2]);
				//positionEnergy[fixed] = 0.5 * stiffness * ((x[0] - restPos[0]) * (x[0] - restPos[0]) + (x[1] - restPos[1]) * (x[1] - restPos[1]) + (x[2] - restPos[2]) * (x[2] - restPos[2]));
			}
		}

	}


	__host__ void solvePointsLocalConstraintsHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simPointsNum,
		qeal* devPoints,
		qeal* devRestPoints,
		qeal* devSn,
		const qeal* pointsMass,
		const qeal dt,
		const qeal* positionStiffness,
		const int* positionConstraintFlag,
		qeal* devRhs,
		qeal* inertialEnergy,
		qeal* positionEnergy
	)
	{
		solvePointsLocalConstraintsKernel << <gridSize, blockSize >> >
			(
				simPointsNum,
				devPoints,
				devRestPoints,
				devSn,
				pointsMass,
				dt,
				positionStiffness,
				positionConstraintFlag,
				devRhs,
				inertialEnergy,
				positionEnergy
				);
		cudaDeviceSynchronize();
	}

	__global__ void updateR1OperatorsKernel
	(
		const int hessR1OperatorsNum,
		qeal* devSysDiagonal,
		qeal* devHessR1Value,
		int* devHessR1DiiIndex,
		qeal* devHessR1Operators
	)
	{
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		int index;
		for (; tid < hessR1OperatorsNum; tid += length)
		{
			index = devHessR1DiiIndex[tid];
			devHessR1Operators[tid] = devHessR1Value[tid] / devSysDiagonal[index];
		}
	}

	__global__ void updateR2OperatorsKernel
	(
		const int hessR2OperatorsNum,
		qeal* devSysDiagonal,
		qeal* devHessR2Aijk,
		int* devHessR2DiiIndex,
		int* devHessR2DjjIndex,
		qeal* devHessR2Operators
	)
	{
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		int diiIndex, djjIndex;
		for (; tid < hessR2OperatorsNum; tid += length)
		{
			diiIndex = devHessR2DiiIndex[tid];
			djjIndex = devHessR2DjjIndex[tid];
			devHessR2Operators[tid] = devHessR2Aijk[tid] / (devSysDiagonal[diiIndex] * devSysDiagonal[djjIndex]);
		}
	}

	__global__ void updateR2RowOperatorsKernel
	(
		const int hessR2RowOperatorsTotalNum,
		int* devHessR2RowOperatorsNum,
		int* devHessR2RowOperatorsOffset,
		qeal* devHessR1Operators,
		qeal* devHessR2Operators,
		qeal* devHessR2ConstOrderOperators,
		int* devHessR2FirstOrderOperators,
		qeal* devHessR2RowOperators,
		const qeal relex =2.0 / 3
	)
	{
		__shared__ qeal val[THREADS_NUM];
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		int num;
		for (; tid < hessR2RowOperatorsTotalNum; tid += length)
		{
			num = devHessR2RowOperatorsNum[tid];
			qeal* operators = devHessR2Operators + devHessR2RowOperatorsOffset[tid];
			val[threadIdx.x] = 0;
			for (int i = 0; i < num; i++)
				val[threadIdx.x] += operators[i];

			val[threadIdx.x] *= relex * relex;
			val[threadIdx.x] += (1 - relex) * (1 - relex) * devHessR2ConstOrderOperators[tid];
			int firstOrderIndex = devHessR2FirstOrderOperators[tid];
			if (firstOrderIndex > -1)
			{
				val[threadIdx.x] += 2 * (1 - relex) * relex * devHessR1Operators[firstOrderIndex];
			}

			devHessR2RowOperators[tid] = val[threadIdx.x];
		}
	}

	__host__ void updateSuperJacobiOperatorsHost
	(
		const dim3 R1OperatorsBlockSize,
		const dim3 R1OperatorsGridSize,
		
		const int hessR1OperatorsNum,
		qeal* devSysDiagonal,
		qeal* devHessR1Value,
		int* devHessR1DiiIndex,
		qeal* devHessR1Operators,
		//
		const dim3 R2OperatorsBlockSize,
		const dim3 R2OperatorsGridSize,
		const int hessR2OperatorsNum,
		qeal* devHessR2Aijk,
		int* devHessR2DiiIndex,
		int* devHessR2DjjIndex,
		qeal* devHessR2Operators,
		//
		const dim3 R2RowOperatorsBlockSize,
		const dim3 R2RowOperatorsGridSize,
		const int hessR2RowOperatorsTotalNum,
		int* devHessR2RowOperatorsNum,
		int* devHessR2RowOperatorsOffset,
		qeal* devHessR2ConstOrderOperators,
		int* devHessR2FirstOrderOperators,
		qeal* devHessR2RowOperators
	)
	{
		updateR1OperatorsKernel << <R1OperatorsGridSize, R1OperatorsBlockSize >> >
			(
				hessR1OperatorsNum,
				devSysDiagonal,
				devHessR1Value,
				devHessR1DiiIndex,
				devHessR1Operators
				);
		cudaDeviceSynchronize();

		updateR2OperatorsKernel << <R2OperatorsGridSize, R2OperatorsBlockSize >> >
			(
				hessR2OperatorsNum,
				devSysDiagonal,
				devHessR2Aijk,
				devHessR2DiiIndex,
				devHessR2DjjIndex,
				devHessR2Operators
				);
		cudaDeviceSynchronize();

		updateR2RowOperatorsKernel << <R2RowOperatorsGridSize, R2RowOperatorsBlockSize >> >
			(
				hessR2RowOperatorsTotalNum,
				devHessR2RowOperatorsNum,
				devHessR2RowOperatorsOffset,
				devHessR1Operators,
				devHessR2Operators,
				devHessR2ConstOrderOperators,
				devHessR2FirstOrderOperators,
				devHessR2RowOperators
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}

	////
	__global__ void superJacobiR2FirstTermKernel
	(
		const int simPointsNum,
		qeal* Diagonal,
		qeal* Rhs,
		qeal* hessR1Operators,
		int* hessR1Index,
		int* hessR1IndexNum,
		int* hessR1IndexOffset,
		qeal* R2FirstTerm,
		const qeal relex = 2.0 / 3
	)
	{
		__shared__ qeal val[THREADS_NUM_192];

		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		int pointId = tid / 3;
		int xyz = threadIdx.x % 3;
		if (pointId < simPointsNum)
		{
			int num = hessR1IndexNum[pointId];
			int offset = hessR1IndexOffset[pointId];
			int* indexList = hessR1Index + offset;
			qeal* operators = hessR1Operators + offset;
			int index;
			val[threadIdx.x] = 0;
			for (int i = 1; i < num; i++)
			{
				index = indexList[i];
				val[threadIdx.x] +=  operators[i] * (Rhs[3 * index + xyz] / Diagonal[index]);
			}
			R2FirstTerm[3 * pointId + xyz] = relex * relex * val[threadIdx.x] + (2 - relex) * relex * Rhs[3 * pointId + xyz] / Diagonal[pointId];
		}
	}

	__global__ void superJacobiR2IterationKernel
	(
		const int simPointsNum,
		qeal* curX,
		qeal* nextX,
		qeal* R2FirstTerm,
		qeal* hessR2Operators,
		int* hessR2Index,
		int* hessR2IndexNum,
		int* hessR2IndexOffset
	)
	{
		__shared__ qeal val[THREADS_NUM_192]; //
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		int pointId = tid / R2_DIM_SPAN;
		int localId = threadIdx.x / R2_DIM_SPAN;
		int batchId = threadIdx.x - localId * R2_DIM_SPAN;

		int localBatchId = batchId % R2_SPAN;
		int xyz = batchId / R2_SPAN;

		int num = 0, offset, index;
		int* indexList;
		qeal* operators;

		if (pointId < simPointsNum)
		{
			num = hessR2IndexNum[pointId];
			offset = hessR2IndexOffset[pointId];
			indexList = hessR2Index + offset;
			operators = hessR2Operators + offset;
		}
		val[threadIdx.x] = 0;
		__syncthreads();

#pragma unroll R2_SPAN
		for (; localBatchId < num; localBatchId += R2_SPAN)
		{
			index = 3 * indexList[localBatchId] + xyz;
			val[threadIdx.x] += operators[localBatchId] * curX[index];
		}
		__syncthreads();

		if (threadIdx.x % R2_SPAN == 0 && pointId < simPointsNum)
		{
#pragma unroll R2_SPAN
			for (int i = 1; i < R2_SPAN; i++)
				val[threadIdx.x] += val[threadIdx.x + i];
			nextX[3 * pointId + xyz] = val[threadIdx.x] + R2FirstTerm[3 * pointId + xyz];
		}
		__syncthreads();
	}

	__global__ void superJacobiR2IterationKernel
	(
		const int simPointsNum,
		qeal* curX,
		qeal* nextX,
		qeal* hessR2Operators,
		int* hessR2Index,
		int* hessR2IndexNum,
		int* hessR2IndexOffset
	)
	{
		__shared__ qeal val[THREADS_NUM_192]; //
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		int pointId = tid / R2_DIM_SPAN;
		int localId = threadIdx.x / R2_DIM_SPAN;
		int batchId = threadIdx.x - localId * R2_DIM_SPAN;

		int localBatchId = batchId % R2_SPAN;
		int xyz = batchId / R2_SPAN;

		int num = 0, offset, index;
		int* indexList;
		qeal* operators;

		if (pointId < simPointsNum)
		{
			num = hessR2IndexNum[pointId];
			offset = hessR2IndexOffset[pointId];
			indexList = hessR2Index + offset;
			operators = hessR2Operators + offset;
		}
		val[threadIdx.x] = 0;
		__syncthreads();

#pragma unroll R2_SPAN
		for (; localBatchId < num; localBatchId += R2_SPAN)
		{
			index = 3 * indexList[localBatchId] + xyz;
			val[threadIdx.x] += operators[localBatchId] * curX[index];
		}
		__syncthreads();

		if (threadIdx.x % R2_SPAN == 0 && pointId < simPointsNum)
		{
#pragma unroll R2_SPAN
			for (int i = 1; i < R2_SPAN; i++)
				val[threadIdx.x] += val[threadIdx.x + i];
			nextX[3 * pointId + xyz] = val[threadIdx.x];
		}
		__syncthreads();
	}

	__global__ void chebyshevJacobiR2IterationKernel
	(
		const int simPointsNum,
		const qeal Omega,
		qeal* R2FirstTerm,
		qeal* hessR2Operators,
		int* hessR2Index,
		int* hessR2IndexNum,
		int* hessR2IndexOffset,
		qeal* preX,
		qeal* curX,
		qeal* nextX
	)
	{
		__shared__ qeal val[THREADS_NUM_192];
		int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

		int pointId = tid / R2_DIM_SPAN;
		int localId = threadIdx.x / R2_DIM_SPAN;
		int batchId = threadIdx.x - localId * R2_DIM_SPAN;

		int localBatchId = batchId % R2_SPAN;
		int xyz = batchId / R2_SPAN;

		int num = 0, offset, index;
		int* indexList;
		qeal* operators;
		if (pointId < simPointsNum)
		{
			num = hessR2IndexNum[pointId];
			offset = hessR2IndexOffset[pointId];
			indexList = hessR2Index + offset;
			operators = hessR2Operators + offset;
		}
		val[threadIdx.x] = 0;
		__syncthreads();

#pragma unroll R2_SPAN
		for (; localBatchId < num; localBatchId += R2_SPAN)
		{
			index = 3 * indexList[localBatchId] + xyz;
			val[threadIdx.x] += operators[localBatchId] * curX[index];
		}
		__syncthreads();

		if (threadIdx.x % R2_SPAN == 0 && pointId < simPointsNum)
		{
#pragma unroll R2_SPAN
			for (int i = 1; i < R2_SPAN; i++)
				val[threadIdx.x] += val[threadIdx.x + i];
			index = 3 * pointId + xyz;
			val[threadIdx.x] += R2FirstTerm[index];
			nextX[index] = Omega * (/*0.9 **/ val[threadIdx.x] /*+ 0.1 * curX[index]*/ - preX[index]) + preX[index];
		}
		__syncthreads();

	}


	__host__ void chebyshevJacobiR2SolverHost
	(
		const int maxIter,
		cublasHandle_t& cublasHandle,
		const dim3 R2FirstTermBlockSize,
		const dim3 R2FirstTermGridSize,
		const int simDimsNum,
		const int simPointsNum,
		qeal* Diagonal,
		qeal* Rhs,
		qeal* hessR1Operators,
		int* hessR1Index,
		int* hessR1IndexNum,
		int* hessR1IndexOffset,
		qeal* R2FirstTerm,
		const dim3 R2IteraionBlockSize,
		const dim3 R2IteraionGridSize,
		qeal* preX,
		qeal* curX,
		qeal* nextX,
		qeal* hessR2Operators,
		int* hessR2Index,
		int* hessR2IndexNum,
		int* hessR2IndexOffset
	)
	{
		superJacobiR2FirstTermKernel << <R2FirstTermGridSize, R2FirstTermBlockSize >> >
			(
				simPointsNum,
				Diagonal,
				Rhs,
				hessR1Operators,
				hessR1Index,
				hessR1IndexNum,
				hessR1IndexOffset,
				R2FirstTerm
				);
		cudaDeviceSynchronize();

		//advance a step for omega
		superJacobiR2IterationKernel << <R2IteraionGridSize, R2IteraionBlockSize >> >
			(
				simPointsNum,
				curX,
				nextX,
				R2FirstTerm,
				hessR2Operators,
				hessR2Index,
				hessR2IndexNum,
				hessR2IndexOffset
				);
		cudaDeviceSynchronize();

		SwapPtr(curX, preX);
		SwapPtr(curX, nextX);
		////////////////////////////
		//superJacobiR2IterationKernel << <R2IteraionGridSize, R2IteraionBlockSize >> >
		//	(
		//		simPointsNum,
		//		curX,
		//		nextX,
		//		hessR2Operators,
		//		hessR2Index,
		//		hessR2IndexNum,
		//		hessR2IndexOffset
		//		);
		//cudaDeviceSynchronize();

		//qeal Xsq, RXsq;
		//cublasSdot(cublasHandle, simDimsNum, curX, 1, curX, 1, &Xsq);
		//cublasSdot(cublasHandle, simDimsNum, nextX, 1, nextX, 1, &RXsq);

		//qeal rho = std::sqrtf(RXsq) / std::sqrtf(Xsq);
		qeal rho = 0.9992;
		rho *= rho;
		
		qeal omega = 2.0 / (2 - rho);

		int iter = 1;
		for (; iter < maxIter; iter++)
		{
			chebyshevJacobiR2IterationKernel << <R2IteraionGridSize, R2IteraionBlockSize >> >
				(
					simPointsNum,
					omega,
					R2FirstTerm,
					hessR2Operators,
					hessR2Index,
					hessR2IndexNum,
					hessR2IndexOffset,
					preX,
					curX,
					nextX
					);
			cudaDeviceSynchronize();
			omega = 4.0 / (4 - rho * omega);
			SwapPtr(curX, preX);
			SwapPtr(curX, nextX);
		}

		
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	__device__ __forceinline__
	void getFaceBbox(const int index, qeal* x, qeal* dx, const int* facesIndex, const qeal inflation_radius, qeal* min, qeal* max)
	{
		int face[3];
		face[0] = facesIndex[3 * index]; face[1] = facesIndex[3 * index + 1]; 	face[2] = facesIndex[3 * index + 2];
		qeal p0t0, p1t0, p2t0;
		qeal p0t1, p1t1, p2t1;
		for (int k = 0; k < 3; k++)
		{
			p0t0 = x[3 * face[0] + k];
			p1t0 = x[3 * face[1] + k];
			p2t0 = x[3 * face[2] + k];

			p0t1 = p0t0 + dx[3 * face[0] + k];
			p1t1 = p1t0 + dx[3 * face[1] + k];
			p2t1 = p2t0 + dx[3 * face[2] + k];

			min[k] = min3(p0t0, p1t0, p2t0);
			min[k] = min3(min[k], p0t1, p1t1);
			min[k] = min2(min[k], p2t1);
			min[k] -= inflation_radius;
			max[k] = max3(p0t0, p1t0, p2t0);
			max[k] = max3(max[k], p0t1, p1t1);
			max[k] = max2(max[k], p2t1);
			max[k] += inflation_radius;
		}
	}

	__device__ __forceinline__
		void getFaceBbox(const int index, qeal* x, const int* facesIndex, const qeal inflation_radius, qeal* min, qeal* max)
	{
		int face[3];
		face[0] = facesIndex[3 * index]; face[1] = facesIndex[3 * index + 1]; 	face[2] = facesIndex[3 * index + 2];
		qeal p0t0, p1t0, p2t0;
		for (int k = 0; k < 3; k++)
		{
			p0t0 = x[3 * face[0] + k];
			p1t0 = x[3 * face[1] + k];
			p2t0 = x[3 * face[2] + k];

			min[k] = min3(p0t0, p1t0, p2t0);
			min[k] -= inflation_radius;
			max[k] = max3(p0t0, p1t0, p2t0);
			max[k] += inflation_radius;
		}
	}

	__device__ __forceinline__
		void mergeBBox(qeal* lhs_min, qeal* lhs_max, qeal* rhs_min, qeal* rhs_max)
	{
		for (int k = 0; k < 3; k++)
		{
			lhs_min[k] = min2(lhs_min[k], rhs_min[k]);
			lhs_max[k] = max2(lhs_max[k], rhs_max[k]);
		}
	}

	__global__ void updateCcdPatchBboxesKernel
	(
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		qeal* devX,
		qeal* devDx,
		qeal* devInflationRadius,
		qeal* devFacesBboxes,
		qeal* devPatchsBboxes
	)
	{
		extern __shared__ qeal bbox[];
		__shared__ int facesNum;
		__shared__ int offset;
		__shared__ qeal inflation_radius;
		if (threadIdx.x == 0)
		{
			facesNum = devPatchFacesNum[blockIdx.x];
			offset = devPatchFacesOffset[blockIdx.x];
			inflation_radius = *devInflationRadius;
		}
		__syncthreads();

		if (threadIdx.x < facesNum)
		{
			int address = offset + threadIdx.x;
			int fid = devPatchFacesList[address];
			getFaceBbox(fid, devX, devDx, devFaces, inflation_radius, bbox + 6 * threadIdx.x, bbox + 6 * threadIdx.x + 3);
			devFacesBboxes[6 * address] = bbox[6 * threadIdx.x];
			devFacesBboxes[6 * address + 1] = bbox[6 * threadIdx.x + 1];
			devFacesBboxes[6 * address + 2] = bbox[6 * threadIdx.x + 2];

			devFacesBboxes[6 * address + 3] = bbox[6 * threadIdx.x + 3];
			devFacesBboxes[6 * address + 4] = bbox[6 * threadIdx.x + 4];
			devFacesBboxes[6 * address + 5] = bbox[6 * threadIdx.x + 5];
		}
		else
		{
			bbox[6 * threadIdx.x] = QEAL_MAX;
			bbox[6 * threadIdx.x + 1] = QEAL_MAX;
			bbox[6 * threadIdx.x + 2] = QEAL_MAX;
			bbox[6 * threadIdx.x + 3] = -QEAL_MAX;
			bbox[6 * threadIdx.x + 4] = -QEAL_MAX;
			bbox[6 * threadIdx.x + 5] = -QEAL_MAX;
		}
		__syncthreads();

		for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
		{
			if (threadIdx.x < stride)
			{
				mergeBBox(bbox + 6 * threadIdx.x, bbox + 6 * threadIdx.x + 3, bbox + 6 * (threadIdx.x + stride), bbox + 6 * (threadIdx.x + stride) + 3);
			}
			__syncthreads();
		}

		if (threadIdx.x == 0)
		{
			devPatchsBboxes[6 * blockIdx.x] = bbox[6 * threadIdx.x];
			devPatchsBboxes[6 * blockIdx.x + 1] = bbox[6 * threadIdx.x + 1];
			devPatchsBboxes[6 * blockIdx.x + 2] = bbox[6 * threadIdx.x + 2];
			devPatchsBboxes[6 * blockIdx.x + 3] = bbox[6 * threadIdx.x + 3];
			devPatchsBboxes[6 * blockIdx.x + 4] = bbox[6 * threadIdx.x + 4];
			devPatchsBboxes[6 * blockIdx.x + 5] = bbox[6 * threadIdx.x + 5];

		}
	}

	__host__ void updateCcdPatchBboxesHost
	(
		const dim3 blockDim,
		const dim3 gridDim,
		int* hostPatchMaxThreads, // must be 2^n and <= 1024
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		qeal* devX,
		qeal* devDx,
		qeal* devInflationRadius,
		qeal* devFacesBboxes,
		qeal* devPatchsBboxes
	)
	{
		updateCcdPatchBboxesKernel << <gridDim, blockDim, (*hostPatchMaxThreads) * 6 * sizeof(qeal) >> >
			(
				devPatchFacesList,
				devPatchFacesNum,
				devPatchFacesOffset,
				devFaces,
				devX,
				devDx,
				devInflationRadius,
				devFacesBboxes,
				devPatchsBboxes
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}

	__global__ void updateDcdPatchBboxesKernel
	(
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		qeal* devX,
		qeal* devInflationRadius,
		qeal* devFacesBboxes,
		qeal* devPatchsBboxes
	)
	{
		extern __shared__ qeal bbox[];
		__shared__ int facesNum;
		__shared__ int offset;
		__shared__ qeal inflation_radius;
		if (threadIdx.x == 0)
		{
			facesNum = devPatchFacesNum[blockIdx.x];
			offset = devPatchFacesOffset[blockIdx.x];
			inflation_radius = *devInflationRadius;
		}
		__syncthreads();

		if (threadIdx.x < facesNum)
		{
			int address = offset + threadIdx.x;
			int fid = devPatchFacesList[address];
			getFaceBbox(fid, devX, devFaces, inflation_radius, bbox + 6 * threadIdx.x, bbox + 6 * threadIdx.x + 3);
			devFacesBboxes[6 * address] = bbox[6 * threadIdx.x];
			devFacesBboxes[6 * address + 1] = bbox[6 * threadIdx.x + 1];
			devFacesBboxes[6 * address + 2] = bbox[6 * threadIdx.x + 2];

			devFacesBboxes[6 * address + 3] = bbox[6 * threadIdx.x + 3];
			devFacesBboxes[6 * address + 4] = bbox[6 * threadIdx.x + 4];
			devFacesBboxes[6 * address + 5] = bbox[6 * threadIdx.x + 5];
		}
		else
		{
			bbox[6 * threadIdx.x] = QEAL_MAX;
			bbox[6 * threadIdx.x + 1] = QEAL_MAX;
			bbox[6 * threadIdx.x + 2] = QEAL_MAX;
			bbox[6 * threadIdx.x + 3] = -QEAL_MAX;
			bbox[6 * threadIdx.x + 4] = -QEAL_MAX;
			bbox[6 * threadIdx.x + 5] = -QEAL_MAX;
		}
		__syncthreads();

		for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
		{
			if (threadIdx.x < stride)
			{
				mergeBBox(bbox + 6 * threadIdx.x, bbox + 6 * threadIdx.x + 3, bbox + 6 * (threadIdx.x + stride), bbox + 6 * (threadIdx.x + stride) + 3);
			}
			__syncthreads();
		}

		if (threadIdx.x == 0)
		{
			devPatchsBboxes[6 * blockIdx.x] = bbox[6 * threadIdx.x];
			devPatchsBboxes[6 * blockIdx.x + 1] = bbox[6 * threadIdx.x + 1];
			devPatchsBboxes[6 * blockIdx.x + 2] = bbox[6 * threadIdx.x + 2];
			devPatchsBboxes[6 * blockIdx.x + 3] = bbox[6 * threadIdx.x + 3];
			devPatchsBboxes[6 * blockIdx.x + 4] = bbox[6 * threadIdx.x + 4];
			devPatchsBboxes[6 * blockIdx.x + 5] = bbox[6 * threadIdx.x + 5];

		}
	}

	__host__ void updateDcdPatchBboxesHost
	(
		const dim3 blockDim,
		const dim3 gridDim,
		int* hostPatchMaxThreads, // must be 2^n and <= 1024
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		qeal* devX,
		qeal* devInflationRadius,
		qeal* devFacesBboxes,
		qeal* devPatchsBboxes
	)
	{
		updateDcdPatchBboxesKernel << <gridDim, blockDim, (*hostPatchMaxThreads) * 6 * sizeof(qeal) >> >
			(
				devPatchFacesList,
				devPatchFacesNum,
				devPatchFacesOffset,
				devFaces,
				devX,
				devInflationRadius,
				devFacesBboxes,
				devPatchsBboxes
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}

	__global__ void updateBvhsBboxesKernel
	(
		const int bvhsNodesGroups,
		int* devBvhsNodesGroupList,
		int* devBvhsNodesGroupNum,
		int* devBvhsNodesGroupOffset,
		int* devNodeIdToPatchId,
		int* devPatchIdToNodeId,
		qeal* devPatchsBboxes,
		qeal* devBvhsNodesBbox
	)
	{
		__shared__ int nodeNum;
		__shared__ int offset;

		for (int i = 0; i < bvhsNodesGroups; i++)
		{
			if (threadIdx.x == 0)
			{
				nodeNum = devBvhsNodesGroupNum[i];
				offset = devBvhsNodesGroupOffset[i];
			}
			__syncthreads();

			for (int id = threadIdx.x; id < nodeNum; id += blockDim.x)
			{
				int nodeId = devBvhsNodesGroupList[offset + id];

				int elementId = devNodeIdToPatchId[nodeId];
		
				
				if (elementId >= 0)
				{
					devBvhsNodesBbox[6 * nodeId] = devPatchsBboxes[6 * elementId];
					devBvhsNodesBbox[6 * nodeId + 1] = devPatchsBboxes[6 * elementId + 1];
					devBvhsNodesBbox[6 * nodeId + 2] = devPatchsBboxes[6 * elementId + 2];
					devBvhsNodesBbox[6 * nodeId + 3] = devPatchsBboxes[6 * elementId + 3];
					devBvhsNodesBbox[6 * nodeId + 4] = devPatchsBboxes[6 * elementId + 4];
					devBvhsNodesBbox[6 * nodeId + 5] = devPatchsBboxes[6 * elementId + 5];
				}
				else
				{
					int childL = 2 * nodeId;
					int childR = childL + 1;
					devBvhsNodesBbox[6 * nodeId] = CUDA_QEAL_MIN(devBvhsNodesBbox[6 * childL], devBvhsNodesBbox[6 * childR]);
					devBvhsNodesBbox[6 * nodeId + 1] = CUDA_QEAL_MIN(devBvhsNodesBbox[6 * childL + 1], devBvhsNodesBbox[6 * childR + 1]);
					devBvhsNodesBbox[6 * nodeId + 2] = CUDA_QEAL_MIN(devBvhsNodesBbox[6 * childL + 2], devBvhsNodesBbox[6 * childR + 2]);
					devBvhsNodesBbox[6 * nodeId + 3] = CUDA_QEAL_MAX(devBvhsNodesBbox[6 * childL + 3], devBvhsNodesBbox[6 * childR + 3]);
					devBvhsNodesBbox[6 * nodeId + 4] = CUDA_QEAL_MAX(devBvhsNodesBbox[6 * childL + 4], devBvhsNodesBbox[6 * childR + 4]);
					devBvhsNodesBbox[6 * nodeId + 5] = CUDA_QEAL_MAX(devBvhsNodesBbox[6 * childL + 5], devBvhsNodesBbox[6 * childR + 5]);
				}
			}
			__syncthreads();

		}
	}

	__host__ void updateBvhsBboxesHost
	(		
		const int bvhsNodesGroups,
		int* devBvhsNodesGroupList,
		int* devBvhsNodesGroupNum,
		int* devBvhsNodesGroupOffset,
		int* devNodeIdToPatchId,
		int* devPatchIdToNodeId,
		qeal* devPatchsBboxes,
		qeal* devBvhsNodesBbox
	)
	{
		dim3 blockSize(THREADS_NUM);
		dim3 gridSize(1);
		updateBvhsBboxesKernel << <gridSize, blockSize >> >
			(
				bvhsNodesGroups,
				devBvhsNodesGroupList,
				devBvhsNodesGroupNum,
				devBvhsNodesGroupOffset,
				devNodeIdToPatchId,
				devPatchIdToNodeId,
				devPatchsBboxes,
				devBvhsNodesBbox
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}

	__device__ __forceinline__
		bool isLeaf
		(
			int* devNodeIdToPatchId,
			int node
		)
	{
		if (devNodeIdToPatchId[node] == -1)
			return false;
		return true;
	}

	__device__ __forceinline__
		bool isNull
		(
			int* devNodeIdToPatchId,
			int node
		)
	{
		if (devNodeIdToPatchId[node] == -2)
			return true;
		return false;
	}

	__device__ __forceinline__
		bool checkOverlap
		(
			qeal* lhs,
			qeal* rhs
		)
	{
		if (lhs[3] < rhs[0] || rhs[3] < lhs[0]) { return false; }
		if (lhs[4] < rhs[1] || rhs[4] < lhs[1]) { return false; }
		if (lhs[5] < rhs[2] || rhs[5] < lhs[2]) { return false; }
		return true;
	}

	__device__ __forceinline__
		void traverseIterative
		(
			const int bvhsNodesNum,
			qeal* bvhsNodesBbox,
			int* nodeIdToPatchId,
			int* patchIdToNodeId,
			const int queryNodeId,
			const int queryElementId,
			int* devPatchICPairsList,
			int& eachPatchICPairsNum
		)
	{
		int stack[64];
		int* stackPtr = stack;
		*stackPtr++ = 0;

		int node = 1; // root node
		eachPatchICPairsNum = 0;

		do
		{
			int childL = 2 * node; // 2 * n
			int childR = childL + 1; // 2 * n + 1
			int childLElementId = nodeIdToPatchId[childL];
			int childRElementId = nodeIdToPatchId[childR];

			bool overlapL = false;
			if (childL < bvhsNodesNum && childLElementId != -2)
			{
				overlapL = checkOverlap(bvhsNodesBbox + 6 * queryNodeId, bvhsNodesBbox + 6 * childL);
			}

			bool overlapR = false;
			if (childR < bvhsNodesNum && childRElementId != -2)
			{
				overlapR = checkOverlap(bvhsNodesBbox + 6 * queryNodeId, bvhsNodesBbox + 6 * childR);
			}

			// Query overlaps a leaf node => report collision.

			if (overlapL && childLElementId != -1 && queryElementId < childLElementId)
			{
				devPatchICPairsList[eachPatchICPairsNum] = queryElementId;
				devPatchICPairsList[eachPatchICPairsNum + 1] = childLElementId;
				eachPatchICPairsNum += 2;

			}

			if (overlapR && childRElementId != -1 && queryElementId < childRElementId)
			{
				devPatchICPairsList[eachPatchICPairsNum] = queryElementId;
				devPatchICPairsList[eachPatchICPairsNum + 1] = childRElementId;
				eachPatchICPairsNum += 2;
			}

			bool traverseL = (overlapL && childLElementId == -1);
			bool traverseR = (overlapR && childRElementId == -1);

			if (!traverseL && !traverseR)
				node = *--stackPtr;
			else
			{
				node = (traverseL) ? childL : childR;
				if (traverseL && traverseR)
					*stackPtr++ = childR;
			}

		} while (node != 0);
	}

	__global__ void findPatchsICPairsKernel
	(
		const int leafNodesNum,
		const int bvhsNodesNum,
		qeal* devBvhsNodesBbox,
		int* devNodeIdToPatchId,
		int* devPatchIdToNodeId,
		int* devPatchICPairsList,
		int* devPatchICPairsListOffset,
		int* devEachPatchICPairsNum
	)
	{
		extern __shared__ int eachPatchICPairsNum[];

		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < leafNodesNum)
		{
			int queryElementId = idx;
			int queryNodeId = devPatchIdToNodeId[queryElementId];
			int offset = devPatchICPairsListOffset[queryElementId];
			traverseIterative
			(
				bvhsNodesNum,
				devBvhsNodesBbox,
				devNodeIdToPatchId,
				devPatchIdToNodeId,
				queryNodeId,
				queryElementId,
				devPatchICPairsList + offset,
				eachPatchICPairsNum[threadIdx.x]
			);
			devEachPatchICPairsNum[idx] = eachPatchICPairsNum[threadIdx.x];
		}
	}

	__global__ void fillCompactPotentialPatchPairListKernel
	(
		const int leafNodesNum,
		int* devPatchICPairsList,
		int* devPatchICPairsListOffset,
		int* devEachPatchICPairsNum,
		int* devEachPatchICPairsOffset,
		int* devPatchICPairs
	)
	{
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < leafNodesNum)
		{
			int num = devEachPatchICPairsNum[idx];
			int patchOffset = devPatchICPairsListOffset[idx];
			int listOffset = devEachPatchICPairsOffset[idx];

			for (int i = 0; i < num; i++)
			{
				devPatchICPairs[listOffset + i] = devPatchICPairsList[patchOffset + i];
			}
		}
	}

	__host__ void findPotentialCollisionsHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int leafNodesNum,
		const int bvhsNodesNum,
		qeal* devBvhsNodesBbox,
		int* devNodeIdToPatchId,
		int* devPatchIdToNodeId,
		int* devPatchICPairsList,
		int* devPatchICPairsListOffset,
		thrust::device_vector<int>* devEachPatchICPairsNum,
		thrust::device_vector<int>* devEachPatchICPairsOffset,
		int* hostPatchICPairsNum,
		int* devPatchICPairsNum,
		int* devPatchICPairs
	)
	{
		int* devEachPatchICPairsNumPtr = thrust::raw_pointer_cast(devEachPatchICPairsNum->data());
		int* devEachPatchICPairsOffsetPtr = thrust::raw_pointer_cast(devEachPatchICPairsOffset->data());
		findPatchsICPairsKernel << <gridSize, blockSize, THREADS_NUM * sizeof(int) >> >
			(
				leafNodesNum,
				bvhsNodesNum,
				devBvhsNodesBbox,
				devNodeIdToPatchId,
				devPatchIdToNodeId,
				devPatchICPairsList,
				devPatchICPairsListOffset,
				devEachPatchICPairsNumPtr
				);
		CUDA_CALL(cudaDeviceSynchronize());

		thrust::exclusive_scan(devEachPatchICPairsNum->begin(), devEachPatchICPairsNum->end(), devEachPatchICPairsOffset->begin());
		(*hostPatchICPairsNum) = (*devEachPatchICPairsNum)[leafNodesNum - 1] + (*devEachPatchICPairsOffset)[leafNodesNum - 1];
		if (*hostPatchICPairsNum == 0)
			return;
		(*hostPatchICPairsNum) /= 2;
		cudaMemcpy(devPatchICPairsNum, hostPatchICPairsNum, sizeof(int), cudaMemcpyHostToDevice);

		fillCompactPotentialPatchPairListKernel << <gridSize, blockSize, THREADS_NUM * sizeof(int) >> >
			(
				leafNodesNum,
				devPatchICPairsList,
				devPatchICPairsListOffset,
				devEachPatchICPairsNumPtr,
				devEachPatchICPairsOffsetPtr,
				devPatchICPairs
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}


	__global__ void computeCcdPatchICPairsVfeeNumKernel
	(
		const int patchICPairsNum,
		int* devPatchICPairs,
		int* devPatchPointsNum,
		int* devPatchFacesNum,
		int* devPatchEdgesNum,
		int* devPatchICPairsVfeeNumList
	)
	{
		int tid = threadIdx.x + blockDim.x * blockIdx.x;
		if (tid < patchICPairsNum)
		{
			int patchAid = devPatchICPairs[2 * tid];
			int patchBid = devPatchICPairs[2 * tid + 1];

			int pointsA = devPatchPointsNum[patchAid];
			int facesA = devPatchFacesNum[patchAid];
			int edgesA = devPatchEdgesNum[patchAid];

			int pointsB = devPatchPointsNum[patchBid];
			int facesB = devPatchFacesNum[patchBid];
			int edgesB = devPatchEdgesNum[patchBid];
			devPatchICPairsVfeeNumList[tid] = 5 * (pointsA * facesB + pointsB * facesA + edgesA * edgesB);
		}
	}

	__global__ void getCcdPotentialICPairsListKernel
	(
		int* devPatchICPairs,
		int* devPatchICPairsVfeeNumList,
		int* devPatchICPairsVfeeOffsetList,
		int* devPatchPointsList,
		int* devPatchPointsNum,
		int* devPatchPointsOffset,
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devPatchEdgesList,
		int* devPatchEdgesNum,
		int* devPatchEdgesOffset,
		int* devFaces,
		int* devEdges,
		int* devPotentialICcdCPairsList
	)
	{
		__shared__ int totalVfee;
		__shared__ int* ICPairsList;

		__shared__ int patchid[2];
		__shared__ int pointsNum[2];
		__shared__ int* pointsList[2];
		__shared__ int facesNum[2];
		__shared__ int* facesList[2];
		__shared__ int edgesNum[2];
		__shared__ int* edgesList[2];

		__shared__ int VaFb;
		__shared__ int VbFa;
		__shared__ int EaEb;

#define PATCH_A 0
#define PATCH_B 1

		if (threadIdx.x == 0)
		{
			totalVfee = devPatchICPairsVfeeNumList[blockIdx.x] / 5;
			ICPairsList = devPotentialICcdCPairsList + devPatchICPairsVfeeOffsetList[blockIdx.x];
		}
		__syncthreads();

		if (threadIdx.x < 2)
		{
			patchid[threadIdx.x] = devPatchICPairs[2 * blockIdx.x + threadIdx.x];
			pointsNum[threadIdx.x] = devPatchPointsNum[patchid[threadIdx.x]];
			pointsList[threadIdx.x] = devPatchPointsList + devPatchPointsOffset[patchid[threadIdx.x]];
			facesNum[threadIdx.x] = devPatchFacesNum[patchid[threadIdx.x]];
			facesList[threadIdx.x] = devPatchFacesList + devPatchFacesOffset[patchid[threadIdx.x]];
			edgesNum[threadIdx.x] = devPatchEdgesNum[patchid[threadIdx.x]];
			edgesList[threadIdx.x] = devPatchEdgesList + devPatchEdgesOffset[patchid[threadIdx.x]];
		}
		__syncthreads();
		if (threadIdx.x == 0)
			VaFb = pointsNum[PATCH_A] * facesNum[PATCH_B];
		if (threadIdx.x == 1)
			VbFa = pointsNum[PATCH_B] * facesNum[PATCH_A];
		if (threadIdx.x == 2)
		{
			EaEb = edgesNum[PATCH_A] * edgesNum[PATCH_B];
		}
		__syncthreads();

		int tid = threadIdx.x;
		for (; tid < totalVfee; tid += blockDim.x)
		{
			int id = tid;
			int aId, int bId;
			if (id < VaFb)
			{
				aId = id / facesNum[PATCH_B];
				bId = id % facesNum[PATCH_B];
				int fid = facesList[PATCH_B][bId];
				//

				ICPairsList[5 * tid] = 0;
				ICPairsList[5 * tid + 1] = pointsList[PATCH_A][aId];
				ICPairsList[5 * tid + 2] = devFaces[3 * fid];
				ICPairsList[5 * tid + 3] = devFaces[3 * fid + 1];
				ICPairsList[5 * tid + 4] = devFaces[3 * fid + 2];

			}
			else if ((id = id - VaFb) < VbFa)
			{
				aId = id / facesNum[PATCH_A];
				bId = id % facesNum[PATCH_A];
				int fid = facesList[PATCH_A][bId];
				//
				ICPairsList[5 * tid] = 0;
				ICPairsList[5 * tid + 1] = pointsList[PATCH_B][aId];
				ICPairsList[5 * tid + 2] = devFaces[3 * fid];
				ICPairsList[5 * tid + 3] = devFaces[3 * fid + 1];
				ICPairsList[5 * tid + 4] = devFaces[3 * fid + 2];
			}
			else if ((id = id - VbFa) < EaEb)
			{
				aId = id / edgesNum[PATCH_B];
				bId = id % edgesNum[PATCH_B];

				int ea = edgesList[PATCH_A][aId];
				int eb = edgesList[PATCH_B][bId];
				ICPairsList[5 * tid] = 1;
				ICPairsList[5 * tid + 1] = devEdges[2 * ea];
				ICPairsList[5 * tid + 2] = devEdges[2 * ea + 1];
				ICPairsList[5 * tid + 3] = devEdges[2 * eb];
				ICPairsList[5 * tid + 4] = devEdges[2 * eb + 1];
			}
		}
	}
	
	__host__ void getCcdPotentialICPairsListHost
		(
			const int hostPatchICPairsNum,
			int* devPatchICPairs,
			int* devPatchPointsList,
			int* devPatchPointsNum,
			int* devPatchPointsOffset,
			int* devPatchFacesList,
			int* devPatchFacesNum,
			int* devPatchFacesOffset,
			int* devPatchEdgesList,
			int* devPatchEdgesNum,
			int* devPatchEdgesOffset,
			int* devFaces,
			int* devEdges,
			thrust::device_vector<int>* devPatchICPairsVfeeNumList,
			thrust::device_vector<int>* devPatchICPairsVfeeOffsetList,
			int* hostPotentialCcdICPairsNum,
			int* devPotentialICcdCPairsNum,
			int* devPotentialICcdCPairsList
		)
		{
			if (hostPatchICPairsNum == 0)
				return;
			dim3 blockSize(THREADS_NUM);
			int num_block = (hostPatchICPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
			dim3 gridSize(num_block);
			int* devPatchICPairsVfeeNumListPtr = thrust::raw_pointer_cast(devPatchICPairsVfeeNumList->data());
			int* devPatchICPairsVfeeOffsetListPtr = thrust::raw_pointer_cast(devPatchICPairsVfeeOffsetList->data());
			computeCcdPatchICPairsVfeeNumKernel << <gridSize, blockSize >> >
				(
					hostPatchICPairsNum,
					devPatchICPairs,
					devPatchPointsNum,
					devPatchFacesNum,
					devPatchEdgesNum,
					devPatchICPairsVfeeNumListPtr
					);
			cudaDeviceSynchronize();
			thrust::exclusive_scan(devPatchICPairsVfeeNumList->begin(), devPatchICPairsVfeeNumList->begin() + (hostPatchICPairsNum), devPatchICPairsVfeeOffsetList->begin());
			(*hostPotentialCcdICPairsNum) = (*devPatchICPairsVfeeNumList)[hostPatchICPairsNum - 1] + (*devPatchICPairsVfeeOffsetList)[hostPatchICPairsNum - 1];
			(*hostPotentialCcdICPairsNum) /= 5;
			if (*hostPotentialCcdICPairsNum == 0)
				return;
			cudaMemcpy(devPotentialICcdCPairsNum, hostPotentialCcdICPairsNum, sizeof(int), cudaMemcpyHostToDevice);

			blockSize = dim3(THREADS_NUM_32);
			gridSize = dim3(hostPatchICPairsNum);

			getCcdPotentialICPairsListKernel << <gridSize, blockSize >> >
				(
					devPatchICPairs,
					devPatchICPairsVfeeNumListPtr,
					devPatchICPairsVfeeOffsetListPtr,
					devPatchPointsList,
					devPatchPointsNum,
					devPatchPointsOffset,
					devPatchFacesList,
					devPatchFacesNum,
					devPatchFacesOffset,
					devPatchEdgesList,
					devPatchEdgesNum,
					devPatchEdgesOffset,
					devFaces,
					devEdges,
					devPotentialICcdCPairsList
					);
			cudaDeviceSynchronize();
		}

	__global__ void computeDcdPatchICPairsF2FNumKernel
	(
		const int patchICPairsNum,
		int* devPatchICPairs,
		int* devPatchFacesNum,
		int* devPatchICPairsF2FNumList
	)
	{
		int tid = threadIdx.x + blockDim.x * blockIdx.x;
		if (tid < patchICPairsNum)
		{
			int patchAid = devPatchICPairs[2 * tid];
			int patchBid = devPatchICPairs[2 * tid + 1];

			int facesA = devPatchFacesNum[patchAid];
			int facesB = devPatchFacesNum[patchBid];
			devPatchICPairsF2FNumList[tid] = 7 * (facesA * facesB); // 6 + 1
		}
	}

	__global__ void getDcdPotentialICPairsListKernel
	(
		int* devPatchICPairs,
		int* devPatchICPairsF2FNumList,
		int* devPatchICPairsF2FOffsetList,
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		int* devPotentialDcdICPairsList
	)
	{
		__shared__ int totalF2F;
		__shared__ int* ICPairsList;
		__shared__ int patchid[2];
		__shared__ int facesNum[2];
		__shared__ int* facesList[2];

#define PATCH_A 0
#define PATCH_B 1

		if (threadIdx.x == 0)
		{
			totalF2F = devPatchICPairsF2FNumList[blockIdx.x] / 7;
			ICPairsList = devPotentialDcdICPairsList + devPatchICPairsF2FOffsetList[blockIdx.x];
		}
		__syncthreads();

		if (threadIdx.x < 2)
		{
			patchid[threadIdx.x] = devPatchICPairs[2 * blockIdx.x + threadIdx.x];
			facesNum[threadIdx.x] = devPatchFacesNum[patchid[threadIdx.x]];
			facesList[threadIdx.x] = devPatchFacesList + devPatchFacesOffset[patchid[threadIdx.x]];
		}
		__syncthreads();


		int tid = threadIdx.x;
		for (; tid < totalF2F; tid += blockDim.x)
		{
			int id = tid;
			int aId, int bId;

			aId = id / facesNum[PATCH_B];
			bId = id % facesNum[PATCH_B];
			int fa = 3 * facesList[PATCH_A][aId];
			int fb = 3 * facesList[PATCH_B][bId];

			ICPairsList[7 * tid] = 0;
			ICPairsList[7 * tid + 1] = devFaces[fa];
			ICPairsList[7 * tid + 2] = devFaces[fa + 1];
			ICPairsList[7 * tid + 3] = devFaces[fa + 2];
			ICPairsList[7 * tid + 4] = devFaces[fb];
			ICPairsList[7 * tid + 5] = devFaces[fb + 1];
			ICPairsList[7 * tid + 6] = devFaces[fb + 2];
		}
	}

	__host__ void getDcdPotentialICPairsListHost
	(
		const int hostPatchICPairsNum,
		int* devPatchICPairs,
		int* devPatchFacesList,
		int* devPatchFacesNum,
		int* devPatchFacesOffset,
		int* devFaces,
		thrust::device_vector<int>* devPatchICPairsF2FNumList,
		thrust::device_vector<int>* devPatchICPairsF2FOffsetList,
		int* hostPotentialDcdICPairsNum,
		int* devPotentialDcdICPairsNum,
		int* devPotentialDcdICPairsList
	)
	{
		if (hostPatchICPairsNum == 0)
			return;
		dim3 blockSize(THREADS_NUM);
		int num_block = (hostPatchICPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);
		int* devPatchICPairsF2FNumListPtr = thrust::raw_pointer_cast(devPatchICPairsF2FNumList->data());
		int* devPatchICPairsF2FOffsetListPtr = thrust::raw_pointer_cast(devPatchICPairsF2FOffsetList->data());
		computeDcdPatchICPairsF2FNumKernel << <gridSize, blockSize >> >
			(
				hostPatchICPairsNum,
				devPatchICPairs,
				devPatchFacesNum,
				devPatchICPairsF2FNumListPtr
				);
		cudaDeviceSynchronize();
		thrust::exclusive_scan(devPatchICPairsF2FNumList->begin(), devPatchICPairsF2FNumList->begin() + (hostPatchICPairsNum), devPatchICPairsF2FOffsetList->begin());
		(*hostPotentialDcdICPairsNum) = (*devPatchICPairsF2FNumList)[hostPatchICPairsNum - 1] + (*devPatchICPairsF2FOffsetList)[hostPatchICPairsNum - 1];
		(*hostPotentialDcdICPairsNum) /= 7;
		if (*hostPotentialDcdICPairsNum == 0)
			return;
		cudaMemcpy(devPotentialDcdICPairsNum, hostPotentialDcdICPairsNum, sizeof(int), cudaMemcpyHostToDevice);

		blockSize = dim3(THREADS_NUM_32);
		gridSize = dim3(hostPatchICPairsNum);

		getDcdPotentialICPairsListKernel << <gridSize, blockSize >> >
			(
				devPatchICPairs,
				devPatchICPairsF2FNumListPtr,
				devPatchICPairsF2FOffsetListPtr,
				devPatchFacesList,
				devPatchFacesNum,
				devPatchFacesOffset,
				devFaces,
				devPotentialDcdICPairsList
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}
		//
		namespace GeometryDistance
		{
			__device__ __forceinline__
				qeal VertexEdgeSqDistance(qeal* v, qeal* ea, qeal* eb, qeal* r, qeal* dir)
			{
				qeal va[3], ba[3];
				getVectorSub(v, ea, va);
				getVectorSub(eb, ea, ba);

				qeal va_ba = getVectorDot(va, ba);
				qeal ba_ba = getVectorDot(ba, ba);

				if (va_ba < 0)
					*r = 0;
				else if (va_ba > ba_ba)
					*r = 1;
				else *r = va_ba / ba_ba;

				dir[0] = v[0] - (ea[0] * (1.0 - *r) + eb[0] * *r);
				dir[1] = v[1] - (ea[1] * (1.0 - *r) + eb[1] * *r);
				dir[2] = v[2] - (ea[2] * (1.0 - *r) + eb[2] * *r);
				return getVectorDot(dir, dir);
			}

			__device__ __forceinline__
				qeal EdgeEdgeSqDistance(qeal* xi, qeal* xj, qeal* xa, qeal* xb, qeal* r, qeal* s, qeal* dir)
			{
				qeal xba[3], xji[3], xai[3];
				getVectorSub(xb, xa, xba);
				getVectorSub(xj, xi, xji);
				getVectorSub(xa, xi, xai);

				getVector3Cross(xji, xba, dir);
				qeal nn = getVectorDot(dir, dir);
				qeal temp[3];

				getVector3Cross(xai, xji, temp);
				qeal weight_aiji = getVectorDot(dir, temp);

				getVector3Cross(xai, xba, temp);
				qeal weight_aiba = getVectorDot(dir, temp);

				if (nn > 1e-24f && weight_aiji >= 0 && weight_aiji <= nn && weight_aiba >= 0 && weight_aiba <= nn)
				{
					*r = weight_aiba / nn;
					*s = weight_aiji / nn;
				}
				else
				{
					qeal minDistance = 999999999;
					qeal distance, v;
					qeal Ndir[3];
					if (weight_aiba < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xb, &v, Ndir)) < minDistance))
					{
						minDistance = distance;
						*r = 0;
						*s = v;
					}
					if (weight_aiba > nn && ((distance = VertexEdgeSqDistance(xj, xa, xb, &v, Ndir)) < minDistance))
					{
						minDistance = distance;
						*r = 1;
						*s = v;
					}
					if (weight_aiji < 0 && ((distance = VertexEdgeSqDistance(xa, xi, xj, &v, Ndir)) < minDistance))
					{
						minDistance = distance;
						*r = v;
						*s = 0;
					}
					if (weight_aiji > nn && ((distance = VertexEdgeSqDistance(xb, xi, xj, &v, Ndir)) < minDistance))
					{
						minDistance = distance;
						*r = v;
						*s = 1;
					}


				}

				dir[0] = xi[0] * (1 - *r) + xj[0] * *r - xa[0] * (1.0 - *s) - xb[0] * *s;
				dir[1] = xi[1] * (1 - *r) + xj[1] * *r - xa[1] * (1.0 - *s) - xb[1] * *s;
				dir[2] = xi[2] * (1 - *r) + xj[2] * *r - xa[2] * (1.0 - *s) - xb[2] * *s;
				return getVectorDot(dir, dir);
			}

			__device__ __forceinline__
				qeal VertexTriangleDistance(qeal* xi, qeal* xa, qeal* xb, qeal* xc, qeal* ba, qeal* bb, qeal* bc, qeal* dir)
			{
				qeal xba[3], xca[3], xia[3];
				getVectorSub(xb, xa, xba);
				getVectorSub(xc, xa, xca);
				getVectorSub(xi, xa, xia);

				getVector3Cross(xba, xca, dir);
				qeal nn = getVectorDot(dir, dir);

				qeal temp[3];
				getVector3Cross(xia, xca, temp);
				qeal weight_iaca = getVectorDot(dir, temp);
				getVector3Cross(xba, xia, temp);
				qeal weight_baia = getVectorDot(dir, temp);

				if (nn > 1e-24f && weight_iaca >= 0 && weight_baia >= 0 && nn - weight_iaca - weight_baia >= 0)
				{
					*bb = weight_iaca / nn;
					*bc = weight_baia / nn;
					*ba = 1 - *bb - *bc;
				}
				else
				{
					qeal minDistance = 999999999;
					qeal r, distance, N[3];
					if (nn - weight_iaca - weight_baia < 0 && ((distance = VertexEdgeSqDistance(xi, xb, xc, &r, N)) < minDistance))
					{
						minDistance = distance;
						*bb = 1 - r;
						*bc = r;
						*ba = 0;
					}
					if (weight_iaca < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xc, &r, N)) < minDistance))
					{
						minDistance = distance;
						*bb = 0;
						*bc = r;
						*ba = 1 - *bb - *bc;
					}
					if (weight_baia < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xb, &r, N)) < minDistance))
					{
						minDistance = distance;
						*bb = r;
						*bc = 0;
						*ba = 1 - *bb - *bc;
					}
				}

				dir[0] = xi[0] - xa[0] * *ba - xb[0] * *bb - xc[0] * *bc;
				dir[1] = xi[1] - xa[1] * *ba - xb[1] * *bb - xc[1] * *bc;
				dir[2] = xi[2] - xa[2] * *ba - xb[2] * *bb - xc[2] * *bc;
				return getVectorDot(dir, dir);
			}

			__device__ __forceinline__
				qeal SimpleVertexTriangleDistance(qeal* xi, qeal* xa, qeal* xb, qeal* xc, qeal* ba, qeal* bb, qeal* bc, qeal* dir)
			{
				qeal  xba[3], xca[3], xia[3];
				getVectorSub(xb, xa, xba);
				getVectorSub(xc, xa, xca);
				getVectorSub(xi, xa, xia);

				getVector3Cross(xia, xca, dir);
				qeal nn = getVectorDot(dir, dir);

				qeal temp[3];
				getVector3Cross(xia, xca, temp);
				qeal weight_iaca = getVectorDot(dir, temp);
				getVector3Cross(xba, xia, temp);
				qeal weight_baia = getVectorDot(dir, temp);

				if (nn > 1e-24f && weight_iaca >= 0 && weight_baia >= 0 && nn - weight_iaca - weight_baia >= 0)
				{
					*bb = weight_iaca / nn;
					*bc = weight_baia / nn;
					*ba = 1.0 - *bb - *bc;
				}
				else	return 999999;

				dir[0] = xi[0] - xa[0] * *ba - xb[0] * *bb - xc[0] * *bc;
				dir[1] = xi[1] - xa[1] * *ba - xb[1] * *bb - xc[1] * *bc;
				dir[2] = xi[2] - xa[2] * *ba - xb[2] * *bb - xc[2] * *bc;
				return getVectorDot(dir, dir);
			}
		}


		__device__ __forceinline__
			bool point_triangle_cd_broadphase(
				qeal* p,
				qeal* t0,
				qeal* t1,
				qeal* t2,
				qeal dist)
		{
			for (int k = 0; k < 3; k++)
			{
				qeal tri_max, tri_min;
				tri_min = min3(t0[k], t1[k], t2[k]);
				tri_max = max3(t0[k], t1[k], t2[k]);
				if ((p[k] - tri_max) > dist || (tri_min - p[k]) > dist)
					return false;
			}
			return true;
		}

		__device__ __forceinline__
			bool point_triangle_ccd_broadphase(
				qeal* p,
				qeal* t0,
				qeal* t1,
				qeal* t2,
				qeal* dp,
				qeal* dt0,
				qeal* dt1,
				qeal* dt2,
				qeal dist)
		{
			qeal p_max, p_min, tri_max, tri_min;
			tri_min = min3(t0[0], t1[0], t2[0]);
			tri_min = min3(tri_min, t0[0] + dt0[0], t1[0] + dt1[0]);
			tri_min = min2(tri_min, t2[0] + dt2[0]);

			tri_max = max3(t0[0], t1[0], t2[0]);
			tri_max = max3(tri_max, t0[0] + dt0[0], t1[0] + dt1[0]);
			tri_max = max2(tri_max, t2[0] + dt2[0]);

			p_min = min2(p[0], p[0] + dp[0]);
			p_max = max2(p[0], p[0] + dp[0]);

			if ((p_min - tri_max) > dist || (tri_min - p_max) > dist)
				return false;

			tri_min = min3(t0[1], t1[1], t2[1]);
			tri_min = min3(tri_min, t0[1] + dt0[1], t1[1] + dt1[1]);
			tri_min = min2(tri_min, t2[1] + dt2[1]);

			tri_max = max3(t0[1], t1[1], t2[1]);
			tri_max = max3(tri_max, t0[1] + dt0[1], t1[1] + dt1[1]);
			tri_max = max2(tri_max, t2[1] + dt2[1]);

			p_min = min2(p[1], p[1] + dp[1]);
			p_max = max2(p[1], p[1] + dp[1]);

			if ((p_min - tri_max) > dist || (tri_min - p_max) > dist)
				return false;

			tri_min = min3(t0[2], t1[2], t2[2]);
			tri_min = min3(tri_min, t0[2] + dt0[2], t1[2] + dt1[2]);
			tri_min = min2(tri_min, t2[2] + dt2[2]);

			tri_max = max3(t0[2], t1[2], t2[2]);
			tri_max = max3(tri_max, t0[2] + dt0[2], t1[2] + dt1[2]);
			tri_max = max2(tri_max, t2[2] + dt2[2]);

			p_min = min2(p[2], p[2] + dp[2]);
			p_max = max2(p[2], p[2] + dp[2]);

			if ((p_min - tri_max) > dist || (tri_min - p_max) > dist)
				return false;

			return true;
		}

		__device__ __forceinline__
			bool edge_edge_cd_broadphase(
				qeal* ea0,
				qeal* ea1,
				qeal* eb0,
				qeal* eb1,
				qeal dist)
		{
			qeal max_a, min_a, max_b, min_b;
			for (int k = 0; k < 3; k++)
			{
				min_a = min2(ea0[k], ea1[k]);
				max_a = max2(ea0[k], ea1[k]);
				min_b = min2(eb0[k], eb1[k]);
				max_b = max2(eb0[k], eb1[k]);
				if ((min_a - max_b) > dist || (min_b - max_a) > dist)
					return false;
			}
			return true;
		}

		__device__ __forceinline__
			bool edge_edge_ccd_broadphase(
				qeal* ea0,
				qeal* ea1,
				qeal* eb0,
				qeal* eb1,
				qeal* dea0,
				qeal* dea1,
				qeal* deb0,
				qeal* deb1,
				qeal dist)
		{
			qeal max_a, min_a, max_b, min_b;
			for (int k = 0; k < 3; k++)
			{
				min_a = min2(ea0[k], ea0[k] + dea0[k]);
				min_a = min2(min_a, ea1[k]);
				min_a = min2(min_a, ea1[k] + dea1[k]);
				max_a = max2(ea0[k], ea0[k] + dea0[k]);
				max_a = max2(max_a, ea1[k]);
				max_a = max2(max_a, ea1[k] + dea1[k]);

				min_b = min2(eb0[k], eb0[k] + deb0[k]);
				min_b = min2(min_b, eb1[k]);
				min_b = min2(min_b, eb1[k] + deb1[k]);
				max_b = max2(eb0[k], eb0[k] + deb0[k]);
				max_b = max2(max_b, eb1[k]);
				max_b = max2(max_b, eb1[k] + deb1[k]);

				if ((min_a - max_b) > dist || (min_b - max_a) > dist)
					return false;
			}
			return true;
		}

		__device__ qeal ccdVertexTriangleToi
		(
			int tid,
			qeal* x0,
			qeal* dx0,
			qeal* x1,
			qeal* dx1,
			qeal* x2,
			qeal* dx2,
			qeal* x3,
			qeal* dx3
		)
		{
			qeal mov[3];
			mov[0] = (dx0[0] + dx1[0] + dx2[0] + dx3[0]) / 4.0;
			mov[1] = (dx0[1] + dx1[1] + dx2[1] + dx3[1]) / 4.0;
			mov[2] = (dx0[2] + dx1[2] + dx2[2] + dx3[2]) / 4.0;

			dx0[0] -= mov[0]; dx0[1] -= mov[1]; dx0[2] -= mov[2];
			dx1[0] -= mov[0]; dx1[1] -= mov[1]; dx1[2] -= mov[2];
			dx2[0] -= mov[0]; dx2[1] -= mov[1]; dx2[2] -= mov[2];
			dx3[0] -= mov[0]; dx3[1] -= mov[1]; dx3[2] -= mov[2];

			qeal disp_mag2_vec = CUDA_QEAL_MAX(getVectorSqNorm(dx1), getVectorSqNorm(dx2));
			disp_mag2_vec = CUDA_QEAL_MAX(disp_mag2_vec, getVectorSqNorm(dx3));
			qeal max_disp_mag = getVectorNorm(dx0) + CUDA_SQRT(disp_mag2_vec);
			if (max_disp_mag == 0.0)
				return 1.0;

			qeal a, b, c, dir[3];
			qeal dist2_cur = GeometryDistance::VertexTriangleDistance(x0, x1, x2, x3, &a, &b, &c, dir);
			qeal dist_cur = CUDA_SQRT(dist2_cur);

			qeal gap = 0.1 * dist_cur;
			qeal toc = 0.0;

			int i = 0;
			for (; i < 10; i++)
			{
				qeal toc_lower_bound = 0.9 * (dist2_cur) / (dist_cur * max_disp_mag);
				x0[0] += toc_lower_bound * dx0[0];
				x0[1] += toc_lower_bound * dx0[1];
				x0[2] += toc_lower_bound * dx0[2];
				//
				x1[0] += toc_lower_bound * dx1[0];
				x1[1] += toc_lower_bound * dx1[1];
				x1[2] += toc_lower_bound * dx1[2];
				//
				x2[0] += toc_lower_bound * dx2[0];
				x2[1] += toc_lower_bound * dx2[1];
				x2[2] += toc_lower_bound * dx2[2];
				//
				x3[0] += toc_lower_bound * dx3[0];
				x3[1] += toc_lower_bound * dx3[1];
				x3[2] += toc_lower_bound * dx3[2];
				dist2_cur = GeometryDistance::VertexTriangleDistance(x0, x1, x2, x3, &a, &b, &c, dir);
				dist_cur = CUDA_SQRT(dist2_cur);


				if (toc && dist_cur < gap) {
					break;
				}
				toc += toc_lower_bound;
				if (toc > 1.0) {
					toc = 1.0;
					break;
				}
			}
			return toc;
		}

		__device__ qeal ccdEdgeEdgeToi
		(
			const int tid, qeal* x0, qeal* dx0,
			qeal* x1, qeal* dx1,
			qeal* x2, qeal* dx2,
			qeal* x3, qeal* dx3
		)
		{
			qeal mov[3];
			mov[0] = (dx0[0] + dx1[0] + dx2[0] + dx3[0]) / 4.0;
			mov[1] = (dx0[1] + dx1[1] + dx2[1] + dx3[1]) / 4.0;
			mov[2] = (dx0[2] + dx1[2] + dx2[2] + dx3[2]) / 4.0;

			dx0[0] -= mov[0]; dx0[1] -= mov[1]; dx0[2] -= mov[2];
			dx1[0] -= mov[0]; dx1[1] -= mov[1]; dx1[2] -= mov[2];
			dx2[0] -= mov[0]; dx2[1] -= mov[1]; dx2[2] -= mov[2];
			dx3[0] -= mov[0]; dx3[1] -= mov[1]; dx3[2] -= mov[2];

			qeal max_disp_mag0 = CUDA_QEAL_MAX(getVectorSqNorm(dx0), getVectorSqNorm(dx1));
			qeal max_disp_mag1 = CUDA_QEAL_MAX(getVectorSqNorm(dx2), getVectorSqNorm(dx3));
			qeal max_disp_mag = CUDA_SQRT(max_disp_mag0) + CUDA_SQRT(max_disp_mag1);

			if (max_disp_mag == 0)
				return 1.0;

			qeal r, s, dir[3];
			qeal dist2_cur = GeometryDistance::EdgeEdgeSqDistance(x0, x1, x2, x3, &r, &s, dir);

			qeal dFunc = dist2_cur;
			if (dFunc <= 0)
			{

				dist2_cur = (x0[0] - x2[0]) * (x0[0] - x2[0]) + (x0[1] - x2[1]) * (x0[1] - x2[1]) + (x0[2] - x2[2]) * (x0[2] - x2[2]);
				dist2_cur = CUDA_QEAL_MIN((x0[0] - x3[0]) * (x0[0] - x3[0]) + (x0[1] - x3[1]) * (x0[1] - x3[1]) + (x0[2] - x3[2]) * (x0[2] - x3[2]), dist2_cur);
				dist2_cur = CUDA_QEAL_MIN((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2]), dist2_cur);
				dist2_cur = CUDA_QEAL_MIN((x1[0] - x3[0]) * (x1[0] - x3[0]) + (x1[1] - x3[1]) * (x1[1] - x3[1]) + (x1[2] - x3[2]) * (x1[2] - x3[2]), dist2_cur);
				dFunc = dist2_cur;
			}
			qeal dist_cur = CUDA_SQRT(dist2_cur);
			qeal gap = 0.1 * dFunc / dist_cur;
			qeal toc = 0.0;

			qeal toc_lower_bound;
			int i = 0;
			for (; i < 10; i++)
			{
				toc_lower_bound = 0.9 * dFunc / (dist_cur * max_disp_mag);
				x0[0] += toc_lower_bound * dx0[0]; x0[1] += toc_lower_bound * dx0[1]; x0[2] += toc_lower_bound * dx0[2];
				x1[0] += toc_lower_bound * dx1[0]; x1[1] += toc_lower_bound * dx1[1]; x1[2] += toc_lower_bound * dx1[2];
				x2[0] += toc_lower_bound * dx2[0]; x2[1] += toc_lower_bound * dx2[1]; x2[2] += toc_lower_bound * dx2[2];
				x3[0] += toc_lower_bound * dx3[0]; x3[1] += toc_lower_bound * dx3[1]; x3[2] += toc_lower_bound * dx3[2];

				dist2_cur = GeometryDistance::EdgeEdgeSqDistance(x0, x1, x2, x3, &r, &s, dir);
				dFunc = dist2_cur;
				if (dFunc <= 0)
				{
					dist2_cur = (x0[0] - x2[0]) * (x0[0] - x2[0]) + (x0[1] - x2[1]) * (x0[1] - x2[1]) + (x0[2] - x2[2]) * (x0[2] - x2[2]);
					dist2_cur = CUDA_QEAL_MIN((x0[0] - x3[0]) * (x0[0] - x3[0]) + (x0[1] - x3[1]) * (x0[1] - x3[1]) + (x0[2] - x3[2]) * (x0[2] - x3[2]), dist2_cur);
					dist2_cur = CUDA_QEAL_MIN((x1[0] - x2[0]) * (x1[0] - x2[0]) + (x1[1] - x2[1]) * (x1[1] - x2[1]) + (x1[2] - x2[2]) * (x1[2] - x2[2]), dist2_cur);
					dist2_cur = CUDA_QEAL_MIN((x1[0] - x3[0]) * (x1[0] - x3[0]) + (x1[1] - x3[1]) * (x1[1] - x3[1]) + (x1[2] - x3[2]) * (x1[2] - x3[2]), dist2_cur);
					dFunc = dist2_cur;
				}
				dist_cur = CUDA_SQRT(dist2_cur);

				if (toc && ((dFunc / dist_cur) < gap)) {
					break;
				}
				toc += toc_lower_bound;
				if (toc > 1.0)
				{
					toc = 1.0;
					break;
				}
			}

			return toc;
		}


		__global__ void getImpactOfTimeKernel
		(
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			qeal* devX0,
			qeal* devX1,
			qeal* devToiList
		)
		{
			__shared__ qeal toiList[THREADS_NUM];
			__shared__ int length;
			if (threadIdx.x == 0)
			{
				length = gridDim.x *  blockDim.x;
			}
			toiList[threadIdx.x] = 1.0;
			__syncthreads();
			int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

			for (; tid < potentialCollisionPairsNum; tid += length)
			{
				int type = devPotentialCollisionPairsList[5 * tid];
				int t0Id = 3 * devPotentialCollisionPairsList[5 * tid + 1];
				int t1Id = 3 * devPotentialCollisionPairsList[5 * tid + 2];
				int t2Id = 3 * devPotentialCollisionPairsList[5 * tid + 3];
				int t3Id = 3 * devPotentialCollisionPairsList[5 * tid + 4];

				qeal x0_0[3], x1_0[3], x2_0[3], x3_0[3];
				qeal dx0[3], dx1[3], dx2[3], dx3[3];
				for (int k = 0; k < 3; k++)
				{
					x0_0[k] = devX0[t0Id + k];
					x1_0[k] = devX0[t1Id + k];
					x2_0[k] = devX0[t2Id + k];
					x3_0[k] = devX0[t3Id + k];
					dx0[k] = devX1[t0Id + k] - x0_0[k];
					dx1[k] = devX1[t1Id + k] - x1_0[k];
					dx2[k] = devX1[t2Id + k] - x2_0[k];
					dx3[k] = devX1[t3Id + k] - x3_0[k];
				}

				qeal toi = 1.0;
				if (type == 0)
				{
					if (t0Id == t1Id || t0Id == t2Id || t0Id == t3Id)
						continue;
					if (!point_triangle_ccd_broadphase(x0_0, x1_0, x2_0, x3_0, dx0, dx1, dx2, dx3, 5e-3))
						continue;

					toi = ccdVertexTriangleToi(tid, x0_0, dx0, x1_0, dx1, x2_0, dx2, x3_0, dx3);
			
					//qeal y = x0_0[1] + dx0[1];
					//if (y < 0)
					//{
					//	qeal xtoi = x0_0[1] / abs(dx0[1]);
					//	if (toi > xtoi)
					//		toi = xtoi;
					//}

					if (toi < toiList[threadIdx.x])
						toiList[threadIdx.x] = toi;
				}
				else if (type == 1)
				{
					if (t0Id == t2Id || t0Id == t3Id || t1Id == t2Id || t1Id == t3Id)
						continue;
					if (!edge_edge_ccd_broadphase(x0_0, x1_0, x2_0, x3_0, dx0, dx1, dx2, dx3, 5e-3))
						continue;

					toi = ccdEdgeEdgeToi(tid, x0_0, dx0, x1_0, dx1, x2_0, dx2, x3_0, dx3);

					if (toi < toiList[threadIdx.x])
						toiList[threadIdx.x] = toi;
				}

			}
			__syncthreads();

			for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
			{
				if (threadIdx.x < stride)
					toiList[threadIdx.x] = CUDA_QEAL_MIN(toiList[threadIdx.x], toiList[threadIdx.x + stride]);
				__syncthreads();
			}
			__syncthreads();

			if (threadIdx.x == 0)
			{
				if (toiList[0] < 1.0)
					toiList[0] *= 0.98;
				devToiList[blockIdx.x] = toiList[0];
			}
			__syncthreads();

		}

		__host__ qeal getImpactOfTimeHost
		(
			cublasHandle_t& blasHandle,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			qeal* devX0,
			qeal* devX1,
			qeal* devToiList
		)
		{
			qeal toi = 1.0;
			if (potentialCollisionPairsNum == 0)
			{
				cudaMemcpy(devToiList, &toi, sizeof(qeal), cudaMemcpyHostToDevice);
				return toi;
			}
			dim3 blockSize(THREADS_NUM);
			int num_block = (potentialCollisionPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
			dim3 gridSize(num_block);

			getImpactOfTimeKernel << <gridSize, blockSize >> >
				(
					potentialCollisionPairsNum,
					devPotentialCollisionPairsList,
					devX0,
					devX1,
					devToiList
					);
			cudaDeviceSynchronize();

			int toiIndex;

			cublasIsamin(blasHandle, num_block, devToiList, 1, &toiIndex);
			cudaMemcpy(devToiList, devToiList + (toiIndex - 1), sizeof(qeal), cudaMemcpyDeviceToDevice);
			CUDA_CALL(cudaMemcpy(&toi, devToiList, sizeof(qeal), cudaMemcpyDeviceToHost));
			return toi;
		}


		__device__ __forceinline__
			qeal barrier(const qeal dist2, const qeal dHat2, const qeal kappa)
		{
			qeal e = 0.0;
			if (dist2 < dHat2) {
				qeal t2 = dist2 - dHat2;
				e = -kappa * t2 * t2 * CUDA_LOG(dist2 / dHat2) / (dHat2 * dHat2);
			}
			return e;
		}

		__global__ void setCollisionConstraintsKernel
		(
			const qeal dHat,
			const qeal dHat2,
			const qeal kappa,
			const qeal dt,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			qeal* devX,
			qeal* devPointsMass,
			qeal* devCollisionDiagonal,
			qeal* devCollisionRhs,
			const int totalPointsNum
		)
		{
			__shared__ int length;
			if (threadIdx.x == 0)
			{
				length = gridDim.x *  blockDim.x;
			}
			__syncthreads();
			int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

			for (; tid < potentialCollisionPairsNum; tid += length)
			{
				int type = devPotentialCollisionPairsList[5 * tid];
				int t0Id = devPotentialCollisionPairsList[5 * tid + 1];
				int t1Id = devPotentialCollisionPairsList[5 * tid + 2];
				int t2Id = devPotentialCollisionPairsList[5 * tid + 3];
				int t3Id = devPotentialCollisionPairsList[5 * tid + 4];

				int buffer0Id = 3 * t0Id;
				int buffer1Id = 3 * t1Id;
				int buffer2Id = 3 * t2Id;
				int buffer3Id = 3 * t3Id;

				qeal x0[3], x1[3], x2[3], x3[3];

				x0[0] = devX[buffer0Id]; x1[0] = devX[buffer1Id]; x2[0] = devX[buffer2Id];  x3[0] = devX[buffer3Id];
				x0[1] = devX[1 + buffer0Id]; x1[1] = devX[1 + buffer1Id]; x2[1] = devX[1 + buffer2Id];  x3[1] = devX[1 + buffer3Id];
				x0[2] = devX[2 + buffer0Id]; x1[2] = devX[2 + buffer1Id]; x2[2] = devX[2 + buffer2Id];  x3[2] = devX[2 + buffer3Id];

				qeal dist, dist2, weight; qeal dir[3];
				qeal mp, mq;
				qeal rp, rq;
				qeal dp, dq;

				if (type == 0)
				{
					if (t0Id == t1Id || t0Id == t2Id || t0Id == t3Id)
						continue;
					if (!point_triangle_cd_broadphase(x0, x1, x2, x3, dHat))
						continue;
					qeal a, b, c;
					dist2 = GeometryDistance::VertexTriangleDistance(x0, x1, x2, x3, &a, &b, &c, dir);
				
					if (dist2 >= dHat2)
						continue;

					dist = CUDA_SQRT(dist2);
					mp = devPointsMass[t0Id];
					mq = a * devPointsMass[t1Id] + b * devPointsMass[t2Id] + c * devPointsMass[t3Id];

					weight = barrier(dist2, dHat2, kappa) + CUDA_QEAL_MIN(mp, mq) / (dt * dt) ;

					atomicAdd(devCollisionDiagonal + t0Id, weight);
					atomicAdd(devCollisionDiagonal + t1Id, weight);
					atomicAdd(devCollisionDiagonal + t2Id, weight);
					atomicAdd(devCollisionDiagonal + t3Id, weight);

					rp = mq / (mp + mq);
					rq = mp / (mp + mq);

					dist = dHat - dist;
				
					dp = dist * rp;
					dq = dist * rq;

					getVectorNormalize(dir);
					qeal dir_dist;
					for (int k = 0; k < 3; k++)
					{
						dir_dist = dp * dir[k];
						atomicAdd(devCollisionRhs + buffer0Id + k, weight * (x0[k] + dir_dist));
						dir_dist = dq * dir[k];
						atomicAdd(devCollisionRhs + buffer1Id + k, weight * (x1[k] - dir_dist));
						atomicAdd(devCollisionRhs + buffer2Id + k, weight * (x2[k] - dir_dist));
						atomicAdd(devCollisionRhs + buffer3Id + k, weight * (x3[k] - dir_dist));
					}
				}
				else if (type == 1)
				{
					if (t0Id == t2Id || t0Id == t3Id || t1Id == t2Id || t1Id == t3Id)
						continue;
					if (!edge_edge_cd_broadphase(x0, x1, x2, x3, dHat))
						continue;
					qeal r, s;
					dist2 = GeometryDistance::EdgeEdgeSqDistance(x0, x1, x2, x3, &r, &s, dir);
					if (dist2 >= dHat2)
						continue;
					dist = CUDA_SQRT(dist2);
					mp = (1 - r) * devPointsMass[t0Id] + r * devPointsMass[t1Id];
					mq = (1 - s) * devPointsMass[t2Id] + s * devPointsMass[t3Id];

					weight = barrier(dist2, dHat2, kappa) + CUDA_QEAL_MIN(mp, mq) / (dt * dt);

					atomicAdd(devCollisionDiagonal + t0Id, weight);
					atomicAdd(devCollisionDiagonal + t1Id, weight);
					atomicAdd(devCollisionDiagonal + t2Id, weight);
					atomicAdd(devCollisionDiagonal + t3Id, weight);

					rp = mq / (mp + mq);
					rq = mp / (mp + mq);

					getVectorNormalize(dir);
					dist = dHat - dist;
			
					dp = dist * rp;
					dq = dist * rq;
					qeal dir_dist;
					for (int k = 0; k < 3; k++)
					{
						dir_dist = dp * dir[k];
						atomicAdd(devCollisionRhs + buffer0Id + k, weight * (x0[k] + dir_dist));
						atomicAdd(devCollisionRhs + buffer1Id + k, weight * (x1[k] + dir_dist));
						dir_dist = dq * dir[k];
						atomicAdd(devCollisionRhs + buffer2Id + k, weight * (x2[k] - dir_dist));
						atomicAdd(devCollisionRhs + buffer3Id + k, weight * (x3[k] - dir_dist));
					}
				}
			}
		}

		__host__ void setCollisionConstraintsHost
		(
			const qeal dHat,
			const qeal dHat2,
			const qeal kappa,
			const qeal dt,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			qeal* devX,
			qeal* devPointsMass,
			qeal* devCollisionDiagonal,
			qeal* devCollisionRhs,
			const int totalPointsNum
		)
		{
			if (potentialCollisionPairsNum == 0)
				return;
			dim3 blockSize(THREADS_NUM);
			int num_block = (potentialCollisionPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
			dim3 gridSize(num_block);

			setCollisionConstraintsKernel << <gridSize, blockSize >> >
				(
					dHat,
					dHat2,
					kappa,
					dt,
					potentialCollisionPairsNum,
					devPotentialCollisionPairsList,
					devX,
					devPointsMass,
					devCollisionDiagonal,
					devCollisionRhs,
					totalPointsNum
					);
			cudaDeviceSynchronize();
		}

		///
		__global__ void prediticeDCPositionKernel
		(
			const int simPointsNum,
			const qeal dt,
			qeal* X,
			qeal* V,
			qeal* gravityForce,
			qeal* mouseForce,
			qeal* invMass,
			const int* positionConstraintFlag,
			qeal* devSn
		)
		{
			__shared__ qeal vx[THREADS_NUM];
			__shared__ qeal vy[THREADS_NUM];
			__shared__ qeal vz[THREADS_NUM];
			int tid = blockDim.x * blockIdx.x + threadIdx.x;
			if (tid >= simPointsNum)	return;
			int vidx = 3 * tid;
			int vidy = vidx + 1;
			int vidz = vidx + 2;
			if (positionConstraintFlag[tid] != -1)
			{
				V[vidx] = 0;
				V[vidy] = 0;
				V[vidz] = 0;
				devSn[vidx] = X[vidx];
				devSn[vidy] = X[vidy];
				devSn[vidz] = X[vidz];
				return;
			}
			qeal m = invMass[tid] * dt;

			vx[threadIdx.x] = V[vidx];
			vy[threadIdx.x] = V[vidy];
			vz[threadIdx.x] = V[vidz];

			vx[threadIdx.x] += m * (gravityForce[vidx] + mouseForce[vidx]);
			vy[threadIdx.x] += m * (gravityForce[vidy] + mouseForce[vidy]);
			vz[threadIdx.x] += m * (gravityForce[vidz] + mouseForce[vidz]);

			vx[threadIdx.x] *= dt;
			vy[threadIdx.x] *= dt;
			vz[threadIdx.x] *= dt;

			devSn[vidx] = X[vidx] + vx[threadIdx.x];
			devSn[vidy] = X[vidy] + vy[threadIdx.x];
			devSn[vidz] = X[vidz] + vz[threadIdx.x];
		}

		__host__ void prediticeDCPositionHost
		(
			const dim3 blockSize,
			const dim3 gridSize,
			const int simPointsNum,
			const qeal dt,
			qeal* X,
			qeal* V,
			qeal* gravityForce,
			qeal* mouseForce,
			qeal* invMass,
			const int* positionConstraintFlag,
			qeal* devSn
		)
		{
			prediticeDCPositionKernel << <gridSize, blockSize >> >
				(
					simPointsNum,
					dt,
					X,
					V,
					gravityForce,
					mouseForce,
					invMass,
					positionConstraintFlag,
					devSn
					);
			cudaDeviceSynchronize();
		}

		__global__ void computePenaltyForceKernal
		(
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			const qeal penaltyStiffness,
			qeal* devX0,
			qeal* devDirection,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes
		)
		{
			__shared__ int length;
			if (threadIdx.x == 0)
			{
				length = gridDim.x *  blockDim.x;
			}
			__syncthreads();
			int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

			for (; tid < potentialCollisionPairsNum; tid += length)
			{
				int type = devPotentialCollisionPairsList[5 * tid];
				int t0Id = 3 * devPotentialCollisionPairsList[5 * tid + 1];
				int t1Id = 3 * devPotentialCollisionPairsList[5 * tid + 2];
				int t2Id = 3 * devPotentialCollisionPairsList[5 * tid + 3];
				int t3Id = 3 * devPotentialCollisionPairsList[5 * tid + 4];

				qeal x0_0[3], x1_0[3], x2_0[3], x3_0[3];
				qeal dx0[3], dx1[3], dx2[3], dx3[3];
				for (int k = 0; k < 3; k++)
				{
					x0_0[k] = devX0[t0Id + k];
					x1_0[k] = devX0[t1Id + k];
					x2_0[k] = devX0[t2Id + k];
					x3_0[k] = devX0[t3Id + k];
					dx0[k] = devDirection[t0Id + k];
					dx1[k] = devDirection[t1Id + k];
					dx2[k] = devDirection[t2Id + k];
					dx3[k] = devDirection[t3Id + k];
				}
				qeal toi = 1.0;
				if (type == 0)
				{
					if (t0Id == t1Id || t0Id == t2Id || t0Id == t3Id)
						continue;
					if (!point_triangle_ccd_broadphase(x0_0, x1_0, x2_0, x3_0, dx0, dx1, dx2, dx3, 1e-4))
						continue;
					toi = ccdVertexTriangleToi(tid, x0_0, dx0, x1_0, dx1, x2_0, dx2, x3_0, dx3);
					if (toi < 1.0)
					{
						qeal tx0[3], tx1[3], tx2[3], tx3[3];
						for (int k = 0; k < 3; k++)
						{
							tx0[k] = devX0[t0Id + k] + devDirection[t0Id + k];
							tx1[k] = devX0[t1Id + k] + devDirection[t1Id + k];
							tx2[k] = devX0[t2Id + k] + devDirection[t2Id + k];
							tx3[k] = devX0[t3Id + k] + devDirection[t3Id + k];
						}

						qeal a, b, c, dir[3];
						qeal dist2 = GeometryDistance::VertexTriangleDistance(tx0, tx1, tx2, tx3, &a, &b, &c, dir);
						dist2 = CUDA_SQRT(dist2);

						qeal force[3];
						force[0] = dir[0] * dist2 * penaltyStiffness / 3.0;
						force[1] = dir[1] * dist2 * penaltyStiffness / 3.0;
						force[2] = dir[2] * dist2 * penaltyStiffness / 3.0;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + t0Id + k, -force[k] * 3.0);
							atomicAdd(devPenaltyForce + t1Id + k, force[k]);
							atomicAdd(devPenaltyForce + t2Id + k, force[k]);
							atomicAdd(devPenaltyForce + t3Id + k, force[k]);
						}

						atomicAdd(devPointsPenaltyTimes + (t0Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t1Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t2Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t3Id / 3), 1);
					}
				}
				else if (type == 1)
				{
					if (t0Id == t2Id || t0Id == t3Id || t1Id == t2Id || t1Id == t3Id)
						continue;
					if (!edge_edge_ccd_broadphase(x0_0, x1_0, x2_0, x3_0, dx0, dx1, dx2, dx3, 1e-4))
						continue;

					toi = ccdEdgeEdgeToi(tid, x0_0, dx0, x1_0, dx1, x2_0, dx2, x3_0, dx3);
					if (toi < 1.0)
					{
						qeal tx0[3], tx1[3], tx2[3], tx3[3];
						for (int k = 0; k < 3; k++)
						{
							tx0[k] = devX0[t0Id + k] + devDirection[t0Id + k];
							tx1[k] = devX0[t1Id + k] + devDirection[t1Id + k];
							tx2[k] = devX0[t2Id + k] + devDirection[t2Id + k];
							tx3[k] = devX0[t3Id + k] + devDirection[t3Id + k];
						}

						qeal r, s, dir[3];
						qeal dist2 = GeometryDistance::EdgeEdgeSqDistance(tx0, tx1, tx2, tx3, &r, &s, dir);
						dist2 = CUDA_SQRT(dist2);

						qeal force[3];
						force[0] = dir[0] * dist2 * penaltyStiffness / 2.0;
						force[1] = dir[1] * dist2 * penaltyStiffness / 2.0;
						force[2] = dir[2] * dist2 * penaltyStiffness / 2.0;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + t0Id + k, -force[k]);
							atomicAdd(devPenaltyForce + t1Id + k, -force[k]);
							atomicAdd(devPenaltyForce + t2Id + k, force[k]);
							atomicAdd(devPenaltyForce + t3Id + k, force[k]);
						}
						atomicAdd(devPointsPenaltyTimes + (t0Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t1Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t2Id / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (t3Id / 3), 1);
					}
				}
			}
			__syncthreads();
		}

		__global__ void uniformPenaltyForceKernal
		(
			const int simPointsNum,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes,
			const qeal penaltyDamping,
			qeal* devVelocity
		)
		{
			int tid = blockDim.x * blockIdx.x + threadIdx.x;
			if (tid >= simPointsNum)	return;

			int times = devPointsPenaltyTimes[tid];
			if (times == 0)
				return;
			int vx = 3 * tid;

			qeal force[3];

			force[0] = devPenaltyForce[vx] / times;
			force[1] = devPenaltyForce[vx + 1] / times;
			force[2] = devPenaltyForce[vx + 2] / times;

			qeal velocity[3];
			velocity[0] = devVelocity[vx];
			velocity[1] = devVelocity[vx + 1];
			velocity[2] = devVelocity[vx + 2];

			devPenaltyForce[vx] = force[0] - penaltyDamping * velocity[0];
			devPenaltyForce[vx + 1] = force[1] - penaltyDamping * velocity[1];
			devPenaltyForce[vx + 2] = force[2] - penaltyDamping * velocity[2];

			getVectorNormalize(force);
			qeal cc = getVectorDot(velocity, force);
			qeal vn[3];
			vn[0] = cc * force[0];
			vn[1] = cc * force[1];
			vn[2] = cc * force[2];

			qeal vt[3];
			vt[0] = velocity[0] - vn[0];
			vt[1] = velocity[1] - vn[1];
			vt[2] = velocity[2] - vn[2];

			vn[0] *= -0.4;
			vn[1] *= -0.4;
			vn[2] *= -0.4;

			vt[0] = vt[0] * 0.9;
			vt[1] = vt[1] * 0.9;
			vt[2] = vt[2] * 0.9;

			devVelocity[vx] = vn[0] + vt[0];
			devVelocity[vx + 1] = vn[1] + vt[1];
			devVelocity[vx + 2] = vn[2] + vt[2];
		}


		__host__ void computePenaltyForceHost
		(
			const dim3 simPointsBlockSize,
			const dim3 simPointsGridSize,
			const int simPointsNum,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			const qeal penaltyStiffness,
			const qeal penaltyDamping,
			qeal* devX0,
			qeal* devDirection,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes,
			qeal* devVelocity
		)
		{
			if (potentialCollisionPairsNum == 0)
			{
				return;
			}
			dim3 blockSize(THREADS_NUM);
			int num_block = (potentialCollisionPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
			dim3 gridSize(num_block);

			computePenaltyForceKernal << <gridSize, blockSize >> >
				(
					potentialCollisionPairsNum,
					devPotentialCollisionPairsList,
					penaltyStiffness,
					devX0,
					devDirection,
					devPenaltyForce,
					devPointsPenaltyTimes
					);
			cudaDeviceSynchronize();

			uniformPenaltyForceKernal << <simPointsGridSize, simPointsBlockSize >> >
				(
					simPointsNum,
					devPenaltyForce,
					devPointsPenaltyTimes,
					penaltyDamping,
					devVelocity
					);
			cudaDeviceSynchronize();


		}


		__global__ void applyPenaltyForceForFloorKernal
		(
			const int simPointsNum,
			const qeal penaltyStiffness,
			qeal* devX,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes
		)
		{
			int tid = blockDim.x * blockIdx.x + threadIdx.x;
			if (tid >= simPointsNum)	return;

			qeal y = devX[3 * tid + 1];

			if (y < 0)
			{
				qeal force = -1 * (y - 0) * penaltyStiffness;
				devPenaltyForce[3 * tid] = 0;
				devPenaltyForce[3 * tid + 1] = force;
				devPenaltyForce[3 * tid + 2] = 0;
				devPointsPenaltyTimes[tid] = 1;
			}

		}

		__host__ void applyPenaltyForceForFloorHost
		(
			const dim3 simPointsBlockSize,
			const dim3 simPointsGridSize,
			const int simPointsNum,
			const qeal penaltyStiffness,
			const qeal penaltyDamping,
			qeal* devX,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes,
			qeal* devVelocity
		)
		{
			applyPenaltyForceForFloorKernal << <simPointsGridSize, simPointsBlockSize >> >
				(
					simPointsNum,
					penaltyStiffness,
					devX,
					devPenaltyForce,
					devPointsPenaltyTimes
					);
			cudaDeviceSynchronize();

			uniformPenaltyForceKernal << <simPointsGridSize, simPointsBlockSize >> >
				(
					simPointsNum,
					devPenaltyForce,
					devPointsPenaltyTimes,
					penaltyDamping,
					devVelocity
					);
			cudaDeviceSynchronize();
		}

		__device__ __forceinline__
			void computeTriangleNormal(qeal* x0, qeal* x1, qeal* x2, qeal* n)
		{
			qeal x0x1[3], x0x2[3];
			for (int k = 0; k < 3; k++)
			{
				x0x1[k] = x1[k] - x0[k];
				x0x2[k] = x2[k] - x0[k];
			}
			getVector3Cross(x0x1, x0x2, n);
			getVectorNormalize(n);
		}

		__device__ __forceinline__
			void crossVV(qeal* Vr, qeal* V1, qeal* V2)
		{
			Vr[0] = V1[1] * V2[2] - V1[2] * V2[1];
			Vr[1] = V1[2] * V2[0] - V1[0] * V2[2];
			Vr[2] = V1[0] * V2[1] - V1[1] * V2[0];
		}

		__device__ __forceinline__
			qeal dotVV(qeal* V1, qeal* V2)
		{
			return (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
		}

		__device__ __forceinline__
			qeal max(qeal a, qeal b, qeal c)
		{
			qeal t = a;
			if (b > t) t = b;
			if (c > t) t = c;
			return t;

		}

		__device__ __forceinline__
			qeal min(qeal a, qeal b, qeal c)
		{
			qeal t = a;
			if (b < t) t = b;
			if (c < t) t = c;
			return t;
		}

		__device__ __forceinline__
			qeal abs_max(qeal a, qeal b, qeal c)
		{
			qeal t = a;
			if (b > t) t = b;
			if (c > t) t = c;
			return CUDA_ABS(t);
		}

		__device__ __forceinline__
			qeal abs_min(qeal a, qeal b, qeal c)
		{
			qeal t = a;
			if (b < t) t = b;
			if (c < t) t = c;
			return CUDA_ABS(t);
		}

		__device__ __forceinline__
			bool project6(qeal* ax, qeal* p1, qeal* p2, qeal* p3, qeal* q1, qeal* q2, qeal* q3)
		{
			qeal P1 = dotVV(ax, p1);
			qeal P2 = dotVV(ax, p2);
			qeal P3 = dotVV(ax, p3);
			qeal Q1 = dotVV(ax, q1);
			qeal Q2 = dotVV(ax, q2);
			qeal Q3 = dotVV(ax, q3);

			qeal mx1 = max(P1, P2, P3);
			qeal mn1 = min(P1, P2, P3);
			qeal mx2 = max(Q1, Q2, Q3);
			qeal mn2 = min(Q1, Q2, Q3);

			if (mn1 > mx2) return false;
			if (mn2 > mx1) return false;
			return true;
		}

		__device__ __forceinline__
			bool triContact(qeal * P1, qeal * P2, qeal * P3, qeal * Q1, qeal * Q2, qeal * Q3)
		{
			qeal p1[3], p2[3], p3[3];
			qeal q1[3], q2[3], q3[3];
			qeal e1[3], e2[3], e3[3];
			qeal f1[3], f2[3], f3[3];
			qeal g1[3], g2[3], g3[3];
			qeal h1[3], h2[3], h3[3];
			qeal n1[3], m1[3];
			//qeal z[3];

			qeal ef11[3], ef12[3], ef13[3];
			qeal ef21[3], ef22[3], ef23[3];
			qeal ef31[3], ef32[3], ef33[3];

			//z[0] = 0.0;  z[1] = 0.0;  z[2] = 0.0;

			p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
			p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
			p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];

			q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
			q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
			q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];

			e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
			e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
			e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

			f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
			f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
			f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];

			crossVV(n1, e1, e2);
			crossVV(m1, f1, f2);

			crossVV(g1, e1, n1);
			crossVV(g2, e2, n1);
			crossVV(g3, e3, n1);
			crossVV(h1, f1, m1);
			crossVV(h2, f2, m1);
			crossVV(h3, f3, m1);

			crossVV(ef11, e1, f1);
			crossVV(ef12, e1, f2);
			crossVV(ef13, e1, f3);
			crossVV(ef21, e2, f1);
			crossVV(ef22, e2, f2);
			crossVV(ef23, e2, f3);
			crossVV(ef31, e3, f1);
			crossVV(ef32, e3, f2);
			crossVV(ef33, e3, f3);

			// now begin the series of tests

			if (!project6(n1, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(m1, p1, p2, p3, q1, q2, q3)) return false;

			if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return false;

			if (!project6(g1, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(g2, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(g3, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(h1, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(h2, p1, p2, p3, q1, q2, q3)) return false;
			if (!project6(h3, p1, p2, p3, q1, q2, q3)) return false;

			return true;
		}

		__device__ __forceinline__ qeal projectToTriangle(qeal * p0, qeal * v0, qeal * v1, qeal * v2, qeal * tri_norm, qeal * dir)
		{
			qeal v1p0[3];
			getVectorSub(p0, v1, v1p0);
			qeal c = getVectorDot(v1p0, tri_norm);
			dir[0] = c * tri_norm[0];
			dir[1] = c * tri_norm[1];
			dir[2] = c * tri_norm[2];
			qeal len = getVectorNorm(dir);

			dir[0] /= len;
			dir[1] /= len;
			dir[2] /= len;

			return len;
		}

		__global__ void computeDcdPenaltyForceKernal
		(
			const int simPointsNum,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			const int* devPointsModelsId,
			const qeal penaltyStiffness,
			const qeal penaltySelfStiffness,
			qeal* devX,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes
		)
		{
			__shared__ int length;
			if (threadIdx.x == 0)
			{
				length = gridDim.x *  blockDim.x;
			}
			__syncthreads();
			int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;

			qeal stiffness;

			for (; tid < potentialCollisionPairsNum; tid += length)
			{
				int type = devPotentialCollisionPairsList[7 * tid];

				int fa0 = 3 * devPotentialCollisionPairsList[7 * tid + 1];
				int fa1 = 3 * devPotentialCollisionPairsList[7 * tid + 2];
				int fa2 = 3 * devPotentialCollisionPairsList[7 * tid + 3];

				int fb0 = 3 * devPotentialCollisionPairsList[7 * tid + 4];
				int fb1 = 3 * devPotentialCollisionPairsList[7 * tid + 5];
				int fb2 = 3 * devPotentialCollisionPairsList[7 * tid + 6];

				if (fa0 == fb0 || fa0 == fb1 || fa0 == fb2)
					continue;
				if (fa1 == fb0 || fa1 == fb1 || fa1 == fb2)
					continue;
				if (fa2 == fb0 || fa2 == fb1 || fa2 == fb2)
					continue;

				qeal fax0[3], fax1[3], fax2[3], faN[3];
				qeal fbx0[3], fbx1[3], fbx2[3], fbN[3];
				for (int k = 0; k < 3; k++)
				{
					fax0[k] = devX[fa0 + k]; fax1[k] = devX[fa1 + k]; fax2[k] = devX[fa2 + k];
					fbx0[k] = devX[fb0 + k]; fbx1[k] = devX[fb1 + k]; fbx2[k] = devX[fb2 + k];
				}

				if (!triContact(fax0, fax1, fax2, fbx0, fbx1, fbx2))
					continue;
	
				if (devPointsModelsId[fa0 / 3] == devPointsModelsId[fb0 / 3])
					stiffness = penaltySelfStiffness;
				else stiffness = penaltyStiffness;

				computeTriangleNormal(fax0, fax1, fax2, faN);
				computeTriangleNormal(fbx0, fbx1, fbx2, fbN);
				//
				qeal dir[3], dist;

				if ((fa0 / 3) < simPointsNum)
				{
					dist = projectToTriangle(fax0, fbx0, fbx1, fbx2, fbN, dir);
					if (getVectorDot(dir, fbN) < 0)
					{
						qeal force[3];
						force[0] = fbN[0] * dist * stiffness;
						force[1] = fbN[1] * dist * stiffness;
						force[2] = fbN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fa0 + k, force[k]);
							atomicAdd(devPenaltyForce + fb0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fa0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb2 / 3), 1);
					}

					dist = projectToTriangle(fax1, fbx0, fbx1, fbx2, fbN, dir);
					if (getVectorDot(dir, fbN) < 0)
					{
						qeal force[3];
						force[0] = fbN[0] * dist * stiffness;
						force[1] = fbN[1] * dist * stiffness;
						force[2] = fbN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fa1 + k, force[k]);
							atomicAdd(devPenaltyForce + fb0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fa1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb2 / 3), 1);
					}

					dist = projectToTriangle(fax2, fbx0, fbx1, fbx2, fbN, dir);
					if (getVectorDot(dir, fbN) < 0)
					{
						qeal force[3];
						force[0] = fbN[0] * dist * stiffness;
						force[1] = fbN[1] * dist * stiffness;
						force[2] = fbN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fa2 + k, force[k]);
							atomicAdd(devPenaltyForce + fb0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fb2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fa2 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fb2 / 3), 1);
					}
				}

				if ((fb0 / 3) < simPointsNum)
				{
					dist = projectToTriangle(fbx0, fax0, fax1, fax2, faN, dir);
					if (getVectorDot(dir, faN) < 0)
					{
						qeal force[3];
						force[0] = faN[0] * dist * stiffness;
						force[1] = faN[1] * dist * stiffness;
						force[2] = faN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fb0 + k, force[k]);
							atomicAdd(devPenaltyForce + fa0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fb0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa2 / 3), 1);
					}

					dist = projectToTriangle(fbx1, fax0, fax1, fax2, faN, dir);
					if (getVectorDot(dir, faN) < 0)
					{
						qeal force[3];
						force[0] = faN[0] * dist * stiffness;
						force[1] = faN[1] * dist * stiffness;
						force[2] = faN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fb1 + k, force[k]);
							atomicAdd(devPenaltyForce + fa0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fb1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa2 / 3), 1);
					}

					dist = projectToTriangle(fbx2, fax0, fax1, fax2, faN, dir);
					if (getVectorDot(dir, faN) < 0)
					{
						qeal force[3];
						force[0] = faN[0] * dist * stiffness;
						force[1] = faN[1] * dist * stiffness;
						force[2] = faN[2] * dist * stiffness;

						for (int k = 0; k < 3; k++)
						{
							atomicAdd(devPenaltyForce + fb2 + k, force[k]);
							atomicAdd(devPenaltyForce + fa0 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa1 + k, -force[k] / 3.0);
							atomicAdd(devPenaltyForce + fa2 + k, -force[k] / 3.0);
						}

						atomicAdd(devPointsPenaltyTimes + (fb2 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa0 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa1 / 3), 1);
						atomicAdd(devPointsPenaltyTimes + (fa2 / 3), 1);
					}
				}				
			}
			__syncthreads();
		}

		__host__ void computeDcdPenaltyForceHost
		(
			const dim3 simPointsBlockSize,
			const dim3 simPointsGridSize,
			const int simPointsNum,
			const int potentialCollisionPairsNum,
			int* devPotentialCollisionPairsList,
			const int* devPointsModelsId,
			const qeal penaltyStiffness,
			const qeal penaltySelfStiffness,
			const qeal penaltyDamping,
			qeal* devX,
			qeal* devPenaltyForce,
			int* devPointsPenaltyTimes,
			qeal* devVelocity
		)
		{
			if (potentialCollisionPairsNum == 0)
			{
				return;
			}
			dim3 blockSize(THREADS_NUM);
			int num_block = (potentialCollisionPairsNum + (THREADS_NUM - 1)) / THREADS_NUM;
			dim3 gridSize(num_block);

			computeDcdPenaltyForceKernal << <gridSize, blockSize >> >
				(
					simPointsNum,
					potentialCollisionPairsNum,
					devPotentialCollisionPairsList,
					devPointsModelsId,
					penaltyStiffness,
					penaltySelfStiffness,
					devX,
					devPenaltyForce,
					devPointsPenaltyTimes
					);
			cudaDeviceSynchronize();

			uniformPenaltyForceKernal << <simPointsGridSize, simPointsBlockSize >> >
				(
					simPointsNum,
					devPenaltyForce,
					devPointsPenaltyTimes,
					penaltyDamping,
					devVelocity
					);
			CUDA_CALL(cudaDeviceSynchronize());
		}


		__global__ void updateR3OperatorsKernel
		(
			const int hessR3OperatorsNum,
			qeal* devSysDiagonal,
			qeal* devHessR3Aijk,
			int* devHessR3DiiIndex,
			int* devHessR3DjjIndex,
			qeal* devHessR1Operators,
			qeal* devHessR2RowOperators,
			qeal* devHessR3Operators
		)
		{
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
			int diiIndex, djjIndex;
			for (; tid < hessR3OperatorsNum; tid += length)
			{
				diiIndex = devHessR3DiiIndex[tid];
				djjIndex = devHessR3DjjIndex[tid];
				devHessR3Operators[tid] = devHessR2RowOperators[diiIndex] * (devHessR3Aijk[tid] / devSysDiagonal[djjIndex]);
			}
		}


		__global__ void updateR3RowOperatorsKernel
		(
			const int hessR3RowOperatorsTotalNum,
			int* devHessR3RowOperatorsNum,
			int* devHessR3RowOperatorsOffset,
			qeal* devHessR2RowOperators,
			qeal* devHessR3Operators,
			qeal* devHessR3ConstOrderOperators,
			int* devHessR3FirstOrderOperators,
			qeal* devHessR3RowOperators,
			const qeal relax = 2.0 / 3
		)
		{
			__shared__ qeal val[THREADS_NUM];
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
			int num;
			for (; tid < hessR3RowOperatorsTotalNum; tid += length)
			{
				num = devHessR3RowOperatorsNum[tid];
				qeal* operators = devHessR3Operators + devHessR3RowOperatorsOffset[tid];
				val[threadIdx.x] = 0;
				for (int i = 0; i < num; i++)
					val[threadIdx.x] += operators[i];

				int firstOrderIndex = devHessR3FirstOrderOperators[tid];
				val[threadIdx.x] *= relax;
				if (firstOrderIndex > -1)
				{
					val[threadIdx.x] += (1 - relax) * devHessR2RowOperators[firstOrderIndex];
				}

				devHessR3RowOperators[tid] = val[threadIdx.x];


			}
		}



		__host__ void updateSuperJacobiOperatorsR3Host
		(
			const dim3 R1OperatorsBlockSize,
			const dim3 R1OperatorsGridSize,

			const int hessR1OperatorsNum,
			qeal* devSysDiagonal,
			qeal* devHessR1Value,
			int* devHessR1DiiIndex,
			qeal* devHessR1Operators,
			//
			const dim3 R2OperatorsBlockSize,
			const dim3 R2OperatorsGridSize,
			const int hessR2OperatorsNum,
			qeal* devHessR2Aijk,
			int* devHessR2DiiIndex,
			int* devHessR2DjjIndex,
			qeal* devHessR2Operators,
			//
			const dim3 R2RowOperatorsBlockSize,
			const dim3 R2RowOperatorsGridSize,
			const int hessR2RowOperatorsTotalNum,
			int* devHessR2RowOperatorsNum,
			int* devHessR2RowOperatorsOffset,
			qeal* devHessR2ConstOrderOperators,
			int* devHessR2FirstOrderOperators,
			qeal* devHessR2RowOperators,

			//
			const dim3 R3OperatorsBlockSize,
			const dim3 R3OperatorsGridSize,
			const int hessR3OperatorsNum,
			qeal* devHessR3Aijk,
			int* devHessR3DiiIndex,
			int* devHessR3DjjIndex,
			qeal* devHessR3Operators,
			//
			const dim3 R3RowOperatorsBlockSize,
			const dim3 R3RowOperatorsGridSize,
			const int hessR3RowOperatorsTotalNum,
			int* devHessR3RowOperatorsNum,
			int* devHessR3RowOperatorsOffset,
			qeal* devHessR3ConstOrderOperators,
			int* devHessR3FirstOrderOperators,
			qeal* devHessR3RowOperators,
			const qeal relax
		)
		{
			updateR1OperatorsKernel << <R1OperatorsGridSize, R1OperatorsBlockSize >> >
				(
					hessR1OperatorsNum,
					devSysDiagonal,
					devHessR1Value,
					devHessR1DiiIndex,
					devHessR1Operators
					);
			cudaDeviceSynchronize();

			updateR2OperatorsKernel << <R2OperatorsGridSize, R2OperatorsBlockSize >> >
				(
					hessR2OperatorsNum,
					devSysDiagonal,
					devHessR2Aijk,
					devHessR2DiiIndex,
					devHessR2DjjIndex,
					devHessR2Operators
					);
			cudaDeviceSynchronize();

			updateR2RowOperatorsKernel << <R2RowOperatorsGridSize, R2RowOperatorsBlockSize >> >
				(
					hessR2RowOperatorsTotalNum,
					devHessR2RowOperatorsNum,
					devHessR2RowOperatorsOffset,
					devHessR1Operators,
					devHessR2Operators,
					devHessR2ConstOrderOperators,
					devHessR2FirstOrderOperators,
					devHessR2RowOperators,
					relax
					);
			cudaDeviceSynchronize();
			//////////////////////////////////////////////////
			updateR3OperatorsKernel << <R3OperatorsGridSize, R3OperatorsBlockSize >> >
				(
					hessR3OperatorsNum,
					devSysDiagonal,
					devHessR3Aijk,
					devHessR3DiiIndex,
					devHessR3DjjIndex,
					devHessR1Operators,
					devHessR2RowOperators,
					devHessR3Operators
					);
			cudaDeviceSynchronize();

			updateR3RowOperatorsKernel << <R3RowOperatorsGridSize, R3RowOperatorsBlockSize >> >
				(
					hessR3RowOperatorsTotalNum,
					devHessR3RowOperatorsNum,
					devHessR3RowOperatorsOffset,
					devHessR2RowOperators,
					devHessR3Operators,
					devHessR3ConstOrderOperators,
					devHessR3FirstOrderOperators,
					devHessR3RowOperators,
					relax
					);
			CUDA_CALL(cudaDeviceSynchronize());


		}



		__global__ void superJacobiR3wFirstTermIR1Kernel
		(
			const int simPointsNum,
			qeal* Diagonal,
			qeal* Rhs,
			qeal* hessR1Operators,
			int* hessR1Index,
			int* hessR1IndexNum,
			int* hessR1IndexOffset,
			qeal* devFirstTerm,
			const qeal relax = 2.0 / 3
		)
		{
			__shared__ qeal val[THREADS_NUM_192];

			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
			int pointId = tid / 3;
			int xyz = threadIdx.x % 3;
			if (pointId < simPointsNum)
			{
				int num = hessR1IndexNum[pointId];
				int offset = hessR1IndexOffset[pointId];
				int* indexList = hessR1Index + offset;
				qeal* operators = hessR1Operators + offset;
				int index;
				val[threadIdx.x] = 0;
				for (int i = 1; i < num; i++)
				{
					index = indexList[i];
					val[threadIdx.x] += operators[i] * (Rhs[3 * index + xyz] / Diagonal[index]);
				}
				devFirstTerm[3 * pointId + xyz] = relax * relax * val[threadIdx.x] + (2 - relax) * relax * Rhs[3 * pointId + xyz] / Diagonal[pointId];
			}
		}

		__global__ void superJacobiR3wFirstTermR2Kernel
		(
			const int simPointsNum,
			qeal* Diagonal,
			qeal* Rhs,
			qeal* hessR2Operators,
			int* hessR2Index,
			int* hessR2IndexNum,
			int* hessR2IndexOffset,
			qeal* devFirstTerm,
			const qeal relax = 2.0 / 3
		)
		{
			__shared__ qeal val[THREADS_NUM_192]; //
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

			int pointId = tid / R2_DIM_SPAN;
			int localId = threadIdx.x / R2_DIM_SPAN;
			int batchId = threadIdx.x - localId * R2_DIM_SPAN;

			int localBatchId = batchId % R2_SPAN;
			int xyz = batchId / R2_SPAN;

			int num = 0, offset, index;
			int* indexList;
			qeal* operators;

			if (pointId < simPointsNum)
			{
				num = hessR2IndexNum[pointId];
				offset = hessR2IndexOffset[pointId];
				indexList = hessR2Index + offset;
				operators = hessR2Operators + offset;
			}
			val[threadIdx.x] = 0;
			__syncthreads();

#pragma unroll R2_SPAN
			for (; localBatchId < num; localBatchId += R2_SPAN)
			{
				index = 3 * indexList[localBatchId] + xyz;
				val[threadIdx.x] += operators[localBatchId] * (Rhs[index] / Diagonal[indexList[localBatchId]]);
			}
			__syncthreads();

			if (threadIdx.x % R2_SPAN == 0 && pointId < simPointsNum)
			{
#pragma unroll R2_SPAN
				for (int i = 1; i < R2_SPAN; i++)
					val[threadIdx.x] += val[threadIdx.x + i];
				devFirstTerm[3 * pointId + xyz] += relax * val[threadIdx.x];
			}
			__syncthreads();
		}

		__global__ void superJacobiR3IterationKernel
		(
			const int simPointsNum,
			qeal* curX,
			qeal* nextX,
			qeal* hessR3Operators,
			int* hessR3Index,
			int* hessR3IndexNum,
			int* hessR3IndexOffset
		)
		{
			__shared__ qeal val[THREADS_NUM_192]; //
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

			int pointId = tid / R3_DIM_SPAN;
			int localId = threadIdx.x / R3_DIM_SPAN;
			int batchId = threadIdx.x - localId * R3_DIM_SPAN;

			int localBatchId = batchId % R3_SPAN;
			int xyz = batchId / R3_SPAN;

			int num = 0, offset, index;
			int* indexList;
			qeal* operators;

			if (pointId < simPointsNum)
			{
				num = hessR3IndexNum[pointId];
				offset = hessR3IndexOffset[pointId];
				indexList = hessR3Index + offset;
				operators = hessR3Operators + offset;
			}
			val[threadIdx.x] = 0;
			__syncthreads();

#pragma unroll R3_SPAN
			for (; localBatchId < num; localBatchId += R3_SPAN)
			{
				index = 3 * indexList[localBatchId] + xyz;
				val[threadIdx.x] += operators[localBatchId] * curX[index];
			}
			__syncthreads();

			if (threadIdx.x % R3_SPAN == 0 && pointId < simPointsNum)
			{
#pragma unroll R3_SPAN
				for (int i = 1; i < R2_SPAN; i++)
					val[threadIdx.x] += val[threadIdx.x + i];
				nextX[3 * pointId + xyz] = val[threadIdx.x];
			}
			__syncthreads();
		}

		__global__ void superJacobiR3IterationKernel
		(
			const int simPointsNum,
			qeal* curX,
			qeal* nextX,
			qeal* R3FirstTerm,
			qeal* hessR3Operators,
			int* hessR3Index,
			int* hessR3IndexNum,
			int* hessR3IndexOffset
		)
		{
			__shared__ qeal val[THREADS_NUM_192]; //
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

			int pointId = tid / R3_DIM_SPAN;
			int localId = threadIdx.x / R3_DIM_SPAN;
			int batchId = threadIdx.x - localId * R3_DIM_SPAN;

			int localBatchId = batchId % R3_SPAN;
			int xyz = batchId / R3_SPAN;

			int num = 0, offset, index;
			int* indexList;
			qeal* operators;

			if (pointId < simPointsNum)
			{
				num = hessR3IndexNum[pointId];
				offset = hessR3IndexOffset[pointId];
				indexList = hessR3Index + offset;
				operators = hessR3Operators + offset;
			}
			val[threadIdx.x] = 0;
			__syncthreads();

#pragma unroll R3_SPAN
			for (; localBatchId < num; localBatchId += R3_SPAN)
			{
				index = 3 * indexList[localBatchId] + xyz;
				val[threadIdx.x] += operators[localBatchId] * curX[index];
			}
			__syncthreads();

			if (threadIdx.x % R3_SPAN == 0 && pointId < simPointsNum)
			{
#pragma unroll R3_SPAN
				for (int i = 1; i < R3_SPAN; i++)
					val[threadIdx.x] += val[threadIdx.x + i];
				index = 3 * pointId + xyz;
				nextX[index] = val[threadIdx.x] + R3FirstTerm[index];
			}
			__syncthreads();
		}

		__global__ void chebyshevJacobiR3IterationKernel
		(
			const int simPointsNum,
			const qeal Omega,
			qeal* R3FirstTerm,
			qeal* hessR3Operators,
			int* hessR3Index,
			int* hessR3IndexNum,
			int* hessR3IndexOffset,
			qeal* preX,
			qeal* curX,
			qeal* nextX
		)
		{
			__shared__ qeal val[THREADS_NUM_192];
			int length = gridDim.x * blockDim.x;
			int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

			int pointId = tid / R3_DIM_SPAN;
			int localId = threadIdx.x / R3_DIM_SPAN;
			int batchId = threadIdx.x - localId * R3_DIM_SPAN;

			int localBatchId = batchId % R3_SPAN;
			int xyz = batchId / R3_SPAN;

			int num = 0, offset, index;
			int* indexList;
			qeal* operators;
			if (pointId < simPointsNum)
			{
				num = hessR3IndexNum[pointId];
				offset = hessR3IndexOffset[pointId];
				indexList = hessR3Index + offset;
				operators = hessR3Operators + offset;
			}
			val[threadIdx.x] = 0;
			__syncthreads();

#pragma unroll R3_SPAN
			for (; localBatchId < num; localBatchId += R3_SPAN)
			{
				index = 3 * indexList[localBatchId] + xyz;
				val[threadIdx.x] += operators[localBatchId] * curX[index];
			}
			__syncthreads();

			if (threadIdx.x % R3_SPAN == 0 && pointId < simPointsNum)
			{
#pragma unroll R3_SPAN
				for (int i = 1; i < R3_SPAN; i++)
					val[threadIdx.x] += val[threadIdx.x + i];
				index = 3 * pointId + xyz;
				val[threadIdx.x] += R3FirstTerm[index];
				nextX[index] = Omega * (val[threadIdx.x] - preX[index]) + preX[index];
			}
			__syncthreads();

		}



		__host__ void chebyshevJacobiR3SolverHost
		(
			const int maxIter,
			cublasHandle_t& cublasHandle,
			const dim3 R2FirstTermBlockSize,
			const dim3 R2FirstTermGridSize,
			const int simDimsNum,
			const int simPointsNum,
			qeal* Diagonal,
			qeal* Rhs,
			qeal* hessR1Operators,
			int* hessR1Index,
			int* hessR1IndexNum,
			int* hessR1IndexOffset,
			qeal* devFirstTerm,
			const dim3 R2IteraionBlockSize,
			const dim3 R2IteraionGridSize,
			qeal* preX,
			qeal* curX,
			qeal* nextX,
			qeal* hessR2Operators,
			int* hessR2Index,
			int* hessR2IndexNum,
			int* hessR2IndexOffset,
			const dim3 R3IteraionBlockSize,
			const dim3 R3IteraionGridSize,
			qeal* hessR3Operators,
			int* hessR3Index,
			int* hessR3IndexNum,
			int* hessR3IndexOffset,
			const qeal relax
		)
		{
			superJacobiR3wFirstTermIR1Kernel << <R2FirstTermGridSize, R2FirstTermBlockSize >> >
				(
					simPointsNum,
					Diagonal,
					Rhs,
					hessR1Operators,
					hessR1Index,
					hessR1IndexNum,
					hessR1IndexOffset,
					devFirstTerm,
					relax
					);
			cudaDeviceSynchronize();

			superJacobiR3wFirstTermR2Kernel << <R2IteraionGridSize, R2IteraionBlockSize >> >
				(
					simPointsNum,
					Diagonal,
					Rhs,
					hessR2Operators,
					hessR2Index,
					hessR2IndexNum,
					hessR2IndexOffset,
					devFirstTerm,
					relax
					);
			cudaDeviceSynchronize();

			//advance a step for omega
			superJacobiR3IterationKernel << <R3IteraionGridSize, R3IteraionBlockSize >> >
				(
					simPointsNum,
					curX,
					nextX,
					devFirstTerm,
					hessR3Operators,
					hessR3Index,
					hessR3IndexNum,
					hessR3IndexOffset
					);
			cudaDeviceSynchronize();

			SwapPtr(curX, preX);
			SwapPtr(curX, nextX);

			qeal rho = 0.9992;
			rho *= rho;
			qeal omega = 2.0 / (2 - rho);
			int iter = 1;
			for (; iter < maxIter; iter++)
			{
				chebyshevJacobiR3IterationKernel << <R3IteraionGridSize, R3IteraionBlockSize >> >
					(
						simPointsNum,
						omega,
						devFirstTerm,
						hessR3Operators,
						hessR3Index,
						hessR3IndexNum,
						hessR3IndexOffset,
						preX,
						curX,
						nextX
						);
				cudaDeviceSynchronize();

				omega = 4.0 / (4 - rho * omega);
				SwapPtr(curX, preX);
				SwapPtr(curX, nextX);
			}

		}
}