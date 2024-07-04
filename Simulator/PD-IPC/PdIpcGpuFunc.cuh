#ifndef PD_IPC_GPU_DEVICE_FUNC_CUH
#define PD_IPC_GPU_DEVICE_FUNC_CUH
#include "Simulator\Cuda\CudaHeader.cuh"


namespace PD_IPC
{	
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
	);

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
	);

	__host__ void computeDirectionHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simDims,
		qeal* devLastPoints,
		qeal* devCurPoints,
		qeal* devDirection
	);

	__host__ void clampDirectionHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int simDims,
		qeal* devLastPoints,
		qeal* devDirection,
		qeal* devCurPoints,
		const qeal toi
	);

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
	);

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
	);

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
	);

	__host__ void solvePointsConstraintHost
	(
		const dim3 blockSize,
		const dim3 gridSize,
		const int pointsConstraintNum,
		qeal* devRestPoints,
		const int* pointConstraintIndex,
		const qeal* positionStiffness,
		qeal* devRhs
	);

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
	);


#define R1_SPAN 4
#define R1_DIM_SPAN 12
#define R2_SPAN 4
#define R2_DIM_SPAN 12
#define R3_SPAN 8
#define R3_Dim_SPAN 24

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
	);

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
	);

#define EACH_MAX_IC_PAIRS 5000
#define MAX_CCD_CANDIDATE_NUM 40000000
#define MAX_BLOCK_TOI 1000000

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
		const qeal inflationRadius,
		qeal* devFacesBboxes,
		qeal* devPatchsBboxes
	);

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
	);

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
	);

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
	);

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
	);

	__host__ qeal getImpactOfTimeHost
	(
		cublasHandle_t& blasHandle,
		const int potentialCollisionPairsNum,
		int* devPotentialCollisionPairsList,
		qeal* devX0,
		qeal* devX1,
		qeal* devToiList
	);

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
		const int totalPointsNum = 0
	);


	////////////
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
	);

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
	);

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
	);

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
	);

#define R3_SPAN 8
#define R3_DIM_SPAN 24

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
	);

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
	);
}




#endif