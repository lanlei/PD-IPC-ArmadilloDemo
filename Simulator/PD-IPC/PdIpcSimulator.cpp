#include "PdIpcSimulator.h"

namespace PD_IPC
{
	template <class T> inline
		void Matrix_Self_Product(T* A, T* R, int nx, int ny)			//R=A'*A
	{
		memset(R, 0, sizeof(T) * ny * ny);
		for (int i = 0; i < ny; i++)
			for (int j = i; j < ny; j++)
			{
				for (int k = 0; k < nx; k++)
					R[i * ny + j] += A[k * ny + i] * A[k * ny + j];
				if (i != j)	R[j * ny + i] = R[i * ny + j];
			}
	}
	
	void PdIpcSimulator::run(int frame)
	{
		doTime(frame);
	}

	void PdIpcSimulator::postRun()
	{	
		cudaMemcpy(points.data(), devPoints, points.size() * sizeof(qeal), cudaMemcpyDeviceToHost);
		std::copy(points.data(), points.data() + 3 * simClothPointsNum, pointsBuffer.buffer.data());
		std::copy(points.data() + 3 * simClothPointsNum, points.data() + 3 * simPointsNum, tetPointsBuffer.buffer.data());
		alignAllMesh(tetPointsBuffer.buffer.data());

		//std::copy(points.data(), points.data() + 3 * simPointsNum, tetPointsBuffer.buffer.data());
		//alignAllMesh(tetPointsBuffer.buffer.data());
	}

	void PdIpcSimulator::initialization()
	{
		totalPointsNum = 0;
		totalElementsNum = 0;
		totalFacesNum = 0;
		
		simPointsNum = 0;
		simElementsNum = 0;
		simClothPointsNum = 0;
		simDims = 0;

		totalBoundaryPointsNum = 0;
		totalBoundaryEdgesNum = 0;
		totalBoundaryFacesNum = 0;

		simBoundaryPointsNum = 0;
		simBoundaryEdgesNum = 0;
		simBoundaryFacesNum = 0;

		staticBoundaryPointsNum = 0;
		staticBoundaryEdgesNum = 0;
		staticBoundaryFacesNum = 0;

		std::vector<std::vector<Vector3>> pointsList;

		std::vector<std::vector<Vector4i>> volumetricElementList;


		std::vector<std::vector<qeal>> simPointsMassList;
		std::vector<std::vector<qeal>> simPointsInvMassList;

		std::vector<std::vector<int>> boundaryPointsList;
		std::vector<std::vector<Vector2i>> boundaryEdgesList;
		std::vector<std::vector<Vector3i>> boundaryFacesList;
	
		for (int i = 0; i < models.size(); i++)
		{
			int ds = totalPointsNum;
			if (getModel(i)->getModelType() == Volumetric)
			{
				PdIpcVolumetricModel* m = dynamic_cast<PdIpcVolumetricModel*> (models[i]);
				m->init
				(
					totalPointsNum,
					totalElementsNum,
					simBoundaryPointsNum,
					simBoundaryEdgesNum,
					simBoundaryFacesNum,
					pointsList,
					volumetricElementList,
					simPointsMassList,
					simPointsInvMassList,
					boundaryPointsList,
					boundaryEdgesList,
					boundaryFacesList,
					volumetricStrainConsList,
					volumetricStrainStiffness,
					volumetricVolumeConsList,
					volumetricVolumeStiffness,
					positionConsList,
					positionStiffness
				);
			}
			else if (getModel(i)->getModelType() == Cloth)
			{
				PdIpcClothModel* m = dynamic_cast<PdIpcClothModel*> (models[i]);
				m->init
				(
					totalPointsNum,
					simClothPointsNum,
					simBoundaryPointsNum,
					simBoundaryEdgesNum,
					simBoundaryFacesNum,
					pointsList,
					simPointsMassList,
					simPointsInvMassList,
					boundaryPointsList,
					boundaryEdgesList,
					boundaryFacesList,
					clothPointsBendingArea,
					clothPointsNeighborsList,
					clothPointsStrainStiffness,
					clothPointsBendingStiffness,
					clothPointsNeighborsRestLength,
					clothPointsNeighborsLbo,
					clothPointsMeanCurvatureNorm,
					positionConsList,
					positionStiffness
				);
			}
		}

		simPointsNum = totalPointsNum;
		simElementsNum = totalElementsNum;
		simDims = simPointsNum * 3;

		simPointsMass.resize(simPointsNum);
		simPointsInvMass.resize(simPointsNum);
		int simPointsOffset = 0;
		for (int i = 0; i < simPointsMassList.size(); i++)
		{
			for (int j = 0; j < simPointsMassList[i].size(); j++)
			{
				simPointsMass.data()[j + simPointsOffset] = simPointsMassList[i][j];
				simPointsInvMass.data()[j + simPointsOffset] = simPointsInvMassList[i][j];
			}
			simPointsOffset += simPointsMassList[i].size();
		}

		for (int i = 0; i < staticModels.size(); i++)
		{
			if (getStaticModel(i)->getModelType() == Volumetric)
			{
				PdIpcVolumetricModel* m = dynamic_cast<PdIpcVolumetricModel*> (staticModels[i]);
				m->init
				(
					totalPointsNum,
					staticBoundaryPointsNum,
					staticBoundaryEdgesNum,
					staticBoundaryFacesNum,
					pointsList,
					boundaryPointsList,
					boundaryEdgesList,
					boundaryFacesList
				);
			}
			else if (getStaticModel(i)->getModelType() == Cloth)
			{
				PdIpcClothModel* m = dynamic_cast<PdIpcClothModel*> (staticModels[i]);
				m->init
				(
					totalPointsNum,
					staticBoundaryPointsNum,
					staticBoundaryEdgesNum,
					staticBoundaryFacesNum,
					pointsList,
					boundaryPointsList,
					boundaryEdgesList,
					boundaryFacesList
				);
			}
		}

		totalDims = 3 * totalPointsNum;
		totalBoundaryPointsNum = simBoundaryPointsNum + staticBoundaryPointsNum;
		totalBoundaryEdgesNum = simBoundaryEdgesNum + staticBoundaryEdgesNum;
		totalBoundaryFacesNum = simBoundaryFacesNum + staticBoundaryFacesNum;
		//
		{
			points.resize(3, totalPointsNum);
			int pointsOffset = 0;
			for (int i = 0; i < pointsList.size(); i++)
			{
				for (int j = 0; j < pointsList[i].size(); j++)
				{
					hostPointsModelsId.push_back(i);
					points.col(j + pointsOffset) = pointsList[i][j];
				}					
				pointsOffset += pointsList[i].size();
			}
		}

		{
			simElements.resize(4, totalElementsNum);
			int elementsOffset = 0;
			for (int i = 0; i < volumetricElementList.size(); i++)
			{
				for (int j = 0; j < volumetricElementList[i].size(); j++)
					simElements.col(j + elementsOffset) = volumetricElementList[i][j];
				elementsOffset += volumetricElementList[i].size();
			}
		}

		{
			boundaryPoints.resize(totalBoundaryPointsNum);
			int pointsOffset = 0;
			for (int i = 0; i <models.size(); i++)
			{
				for (int j = 0; j < boundaryPointsList[i].size(); j++)
					boundaryPoints[pointsOffset + j] = boundaryPointsList[i][j];
				pointsOffset += boundaryPointsList[i].size();
			}
		}

		{
			boundaryEdges.resize(2, totalBoundaryEdgesNum);
			int edgesOffset = 0;
			for (int i = 0; i < boundaryEdgesList.size(); i++)
			{
				for (int j = 0; j < boundaryEdgesList[i].size(); j++)
					boundaryEdges.col(j + edgesOffset) = boundaryEdgesList[i][j];
				edgesOffset += boundaryEdgesList[i].size();
			}
		}

		{
			boundaryFaces.resize(3, totalBoundaryFacesNum);
			int facesOffset = 0;
			for (int i = 0; i < boundaryFacesList.size(); i++)
			{
				for (int j = 0; j < boundaryFacesList[i].size(); j++)
					boundaryFaces.col(j + facesOffset) = boundaryFacesList[i][j];
				facesOffset += boundaryFacesList[i].size();
			}
		}

		restPoints = points;
		velocity = points.block(0, 0, 3, simPointsNum);
		velocity.setZero();
		Sn = points.block(0, 0, 3, simPointsNum);
		Sn.setZero();
		externalForce = points.block(0, 0, 3, simPointsNum);
		externalForce.setZero();
		gravityForce = points.block(0, 0, 3, simPointsNum);
		gravityForce.setZero();
		for (int i = 0; i < simPointsNum; i++)
		{
			gravityForce.col(i) = Vector3(0, getGraviry() *simPointsMass[i], 0);
		}
		mouseForce = points.block(0, 0, 3, simPointsNum);
		mouseForce.setZero();

		setupHessMat();
		gpuGeneralInit();
		gpuSolverInit();
		gpuCollisionInit();

		std::cout << "simPointsNum: " << simPointsNum << std::endl;
		std::cout << "simElementNum: " << simElementsNum << std::endl;
		std::cout << "simClothPointsNum: " << simClothPointsNum << std::endl;	
		std::cout << "boundary edges: " << totalBoundaryEdgesNum << std::endl;
		std::cout << "boundary faces: " << totalBoundaryFacesNum << std::endl;
		std::cout << "boundary edges: " << totalBoundaryEdgesNum << std::endl;

	}

	void PdIpcSimulator::setupHessMat()
	{
		std::vector<TripletX> triplets;

		hostSimPointsConstraintFlag = std::vector<int>(simPointsNum, -1);
		positionConstraintNum = 0;
		for (int i = 0; i < positionConsList.size(); i++)
		{
			qeal pointStiffness = positionStiffness[i];
			for (int j = 0; j < positionConsList[i].size(); j++)
			{
				hostSimPointsConstraintFlag[positionConsList[i][j]] = positionConstraintNum;
				hostSimPointsConstraintIndex.push_back(positionConsList[i][j]);
				hostSimPointsConstraintStiffness.push_back(pointStiffness);
				positionConstraintNum++;
			}

		}

		for (int i = 0; i < simPointsNum; i++)
		{
			int fixId = hostSimPointsConstraintFlag[i];
			if (fixId != -1)
			{
				qeal stiffness = hostSimPointsConstraintStiffness[fixId];
				triplets.push_back(TripletX(i, i, stiffness));
			}
		}
		SparseMatrix pointHess(simPointsNum, simPointsNum);
		pointHess.setFromTriplets(triplets.begin(), triplets.end());

		triplets.clear();
		for (int i = 0; i < volumetricStrainConsList.size(); i++)
		{
			for (int j = 0; j < volumetricStrainConsList[i].size(); j++)
			{
				Vector4i ele = volumetricStrainConsList[i][j];

				hostSimElements.push_back(ele[0]);
				hostSimElements.push_back(ele[1]);
				hostSimElements.push_back(ele[2]);
				hostSimElements.push_back(ele[3]);

				Vector3 p0 = points.col(ele[0]);
				Vector3 p1 = points.col(ele[1]);
				Vector3 p2 = points.col(ele[2]);
				Vector3 p3 = points.col(ele[3]);

				Matrix3 edges, rest;
				edges.row(0) = p1 - p0;
				edges.row(1) = p2 - p0;
				edges.row(2) = p3 - p0;
				rest = edges.inverse();
				qeal V = std::abs(edges.determinant() / 6.0);
				qeal stiffness = volumetricStrainStiffness[i] * V;
			
				hostVolumetricStrainStiffness.push_back(stiffness);

				for (int k = 0; k < 9; k++)
					hostSimElementsRestShape.push_back(rest.data()[k]);

				qeal	half_matrix[3][4];
				half_matrix[0][0] = -rest.data()[0] - rest.data()[3] - rest.data()[6];
				half_matrix[0][1] = rest.data()[0];
				half_matrix[0][2] = rest.data()[3];
				half_matrix[0][3] = rest.data()[6];
				half_matrix[1][0] = -rest.data()[1] - rest.data()[4] - rest.data()[7];
				half_matrix[1][1] = rest.data()[1];
				half_matrix[1][2] = rest.data()[4];
				half_matrix[1][3] = rest.data()[7];
				half_matrix[2][0] = -rest.data()[2] - rest.data()[5] - rest.data()[8];
				half_matrix[2][1] = rest.data()[2];
				half_matrix[2][2] = rest.data()[5];
				half_matrix[2][3] = rest.data()[8];

				qeal	full_matrix[4][4];
				Matrix_Self_Product(&half_matrix[0][0], &full_matrix[0][0], 3, 4);

				for (int r = 0; r < 4; r++)
					for (int c = 0; c < 4; c++)
						triplets.push_back(Eigen::Triplet<qeal>(ele[r], ele[c], full_matrix[r][c] * stiffness));
			}
		}
		volumetricStrainHess.resize(simPointsNum, simPointsNum);
		volumetricStrainHess.setFromTriplets(triplets.begin(), triplets.end());

		/////
		triplets.clear();
		for (int i = 0; i < clothPointsNeighborsList.size(); i++)
		{
			int pi = clothPointsNeighborsList[i][0];
			qeal stiffness = clothPointsStrainStiffness[i][0];

			for (int k = 1; k < clothPointsNeighborsList[i].size(); k++)
			{
				int pj = clothPointsNeighborsList[i][k];
				triplets.push_back(TripletX(pi, pi, stiffness));
				triplets.push_back(TripletX(pi, pj, -stiffness));	
			}
		}

		for (int i = 0; i < clothPointsNeighborsList.size(); i++)
		{
			int lbo_length;
			MatrixX ATA;
			std::vector<int> index;
			lbo_length = clothPointsNeighborsLbo[i].size();
			index.resize(lbo_length);
			for (int v = 0; v < clothPointsNeighborsList[i].size(); v++)
			{
				index[v] = clothPointsNeighborsList[i][v];
			}
			qeal lbo2 = 0;
			VectorX lbov(clothPointsNeighborsLbo[i].size());
			//VectorX bendingStiffness(clothPointsNeighborsLbo[i].size());
			for (int k = 0; k < clothPointsNeighborsLbo[i].size(); k++)
			{
				lbov[k] = clothPointsNeighborsLbo[i][k];
				//bendingStiffness[k] = clothPointsBendingStiffness[i][k];
			}

			//ATA = bendingStiffness * lbov.transpose();
			ATA = clothPointsBendingArea[i] * (lbov * lbov.transpose());
			for (int r = 0; r < lbo_length; r++)
				for (int c = 0; c < lbo_length; c++)
					triplets.push_back(TripletX(index[r], index[c], ATA.data()[lbo_length * r + c]));
		}
		


		clothHess.resize(simPointsNum, simPointsNum);
		clothHess.setFromTriplets(triplets.begin(), triplets.end());

		elasticsHess = pointHess + volumetricStrainHess + clothHess;

		Mass.resize(simPointsNum, simPointsNum);
		triplets.clear();
		for (int i = 0; i < simPointsNum; i++)
			triplets.push_back(TripletX(i, i, simPointsMass[i]));
		Mass.setFromTriplets(triplets.begin(), triplets.end());
		
		sysHess = elasticsHess + Mass / (_timeStep * _timeStep);

		std::cout << "SimPointsNum: " << simPointsNum << std::endl;
	}

	void PdIpcSimulator::gpuGeneralInit()
	{
		totalDimsBlockSize = dim3(THREADS_NUM_128);
		totalDimsGridSize = dim3((totalDims + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		totalPointsBlockSize = dim3(THREADS_NUM_128);
		totalPointsGridSize = dim3((simPointsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		simDimsBlockSize = dim3(THREADS_NUM_128);
		simDimsGridSize = dim3((simDims + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		simPointsBlockSize = dim3(THREADS_NUM_128);
		simPointsGridSize = dim3((simPointsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		simElementsBlockSize = dim3(THREADS_NUM_128);
		simElementsGridSize = dim3((simElementsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		simClothPointsBlockSize = dim3(THREADS_NUM_128);
		simClothPointsGridSize = dim3((simClothPointsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);
		simPointsConstraintBlockSize = dim3(THREADS_NUM_128);
		simPointsConstraintGridSize = dim3((positionConstraintNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128);


		//int* devBoundaryFaces, *devBoundaryEdges;
		CUDA_CALL(cudaMalloc((void**)&devBoundaryFaces, sizeof(int) * boundaryFaces.size()));
		CUDA_CALL(cudaMemcpy(devBoundaryFaces, boundaryFaces.data(), sizeof(int) * boundaryFaces.size(), cudaMemcpyHostToDevice));

		cudaMalloc((void**)&devBoundaryEdges, sizeof(int) * boundaryEdges.size());
		cudaMemcpy(devBoundaryEdges, boundaryEdges.data(), sizeof(int) * boundaryEdges.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devPointsModelsId, sizeof(int) * hostPointsModelsId.size());
		cudaMemcpy(devPointsModelsId, hostPointsModelsId.data(), sizeof(int) * hostPointsModelsId.size(), cudaMemcpyHostToDevice);

		//
		VectorX simZeroVec(3 * totalPointsNum);
		simZeroVec.setZero();
		cudaMalloc((void**)&devZeroVector, sizeof(qeal) * 3 * totalPointsNum);
		cudaMemcpy(devZeroVector, simZeroVec.data(), sizeof(qeal) * 3 * totalPointsNum, cudaMemcpyHostToDevice);

		VectorXi simZeroVecInt(3 * totalPointsNum);
		simZeroVecInt.setZero();
		cudaMalloc((void**)&devZeroVectorInt, sizeof(int) * 3 * totalPointsNum);
		cudaMemcpy(devZeroVectorInt, simZeroVecInt.data(), sizeof(int) * 3 * totalPointsNum, cudaMemcpyHostToDevice);


		cudaMalloc((void**)&devPoints, sizeof(qeal) * totalDims);
		cudaMemcpy(devPoints, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devRestPoints, sizeof(qeal) * totalDims);
		cudaMemcpy(devRestPoints, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);
		
		cudaMalloc((void**)&devX0, sizeof(qeal) * totalDims);
		cudaMemcpy(devX0, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devOldPoints, sizeof(qeal) * totalDims);
		cudaMemcpy(devOldPoints, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devDirection, sizeof(qeal) * totalDims);
		
		cudaMalloc((void**)&devPrePoints, sizeof(qeal) * totalDims);
		cudaMemcpy(devPrePoints, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devNextPoints, sizeof(qeal) * totalDims);
		cudaMemcpy(devNextPoints, points.data(), sizeof(qeal) * totalDims, cudaMemcpyHostToDevice);

		//
		cudaMalloc((void**)&devVelocity, sizeof(qeal) * simDims);
		cudaMemcpy(devVelocity, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devNextVelocity, sizeof(qeal) * simDims);
		cudaMemcpy(devNextVelocity, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		std::vector<std::vector<int>> simPointsNeighborsList;
		getSimPointsNeighborsList(simPointsNeighborsList);
		flatten2DArray(simPointsNeighborsList, hostSimPointsNeighborsList, hostSimPointsNeighborsNum, hostSimPointsNeighborsOffset);
		
		cudaMalloc((void**)&devSimPointsNeighborsList, sizeof(int) * hostSimPointsNeighborsList.size());
		cudaMemcpy(devSimPointsNeighborsList, hostSimPointsNeighborsList.data(), sizeof(int) * hostSimPointsNeighborsList.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSimPointsNeighborsNum, sizeof(int) * hostSimPointsNeighborsNum.size());
		cudaMemcpy(devSimPointsNeighborsNum, hostSimPointsNeighborsNum.data(), sizeof(int) * hostSimPointsNeighborsNum.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSimPointsNeighborsOffset, sizeof(int) * hostSimPointsNeighborsOffset.size());
		cudaMemcpy(devSimPointsNeighborsOffset, hostSimPointsNeighborsOffset.data(), sizeof(int) * hostSimPointsNeighborsOffset.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSimPointsConstraintFlag, sizeof(int) * hostSimPointsConstraintFlag.size());
		cudaMemcpy(devSimPointsConstraintFlag, hostSimPointsConstraintFlag.data(), sizeof(int) * hostSimPointsConstraintFlag.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSimPointsConstraintIndex, sizeof(int) * hostSimPointsConstraintIndex.size());
		cudaMemcpy(devSimPointsConstraintIndex, hostSimPointsConstraintIndex.data(), sizeof(int) * hostSimPointsConstraintIndex.size(), cudaMemcpyHostToDevice);

		

		cudaMalloc((void**)&devSimPointsConstraintStiffness, sizeof(qeal) * hostSimPointsConstraintStiffness.size());
		cudaMemcpy(devSimPointsConstraintStiffness, hostSimPointsConstraintStiffness.data(), sizeof(qeal) * hostSimPointsConstraintStiffness.size(), cudaMemcpyHostToDevice);
		//
		cudaMalloc((void**)&devSimElements, sizeof(int) * hostSimElements.size());
		cudaMemcpy(devSimElements, hostSimElements.data(), sizeof(int) * hostSimElements.size(), cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devSimElementsRestShape, sizeof(qeal) * hostSimElementsRestShape.size());
		cudaMemcpy(devSimElementsRestShape, hostSimElementsRestShape.data(), sizeof(qeal) * hostSimElementsRestShape.size(), cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devVolumetricStrainStiffness, sizeof(qeal) * hostVolumetricStrainStiffness.size());
		cudaMemcpy(devVolumetricStrainStiffness, hostVolumetricStrainStiffness.data(), sizeof(qeal) * hostVolumetricStrainStiffness.size(), cudaMemcpyHostToDevice);
		//
		cudaMalloc((void**)&devInertialEnergy, sizeof(qeal) * simPointsNum);
		cudaMemcpy(devInertialEnergy, simZeroVec.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);
		cudaMalloc((void**)&devPositionEnergy, sizeof(qeal) * positionConstraintNum);
		cudaMemcpy(devPositionEnergy, simZeroVec.data(), sizeof(qeal) * positionConstraintNum, cudaMemcpyHostToDevice);
		
		cudaMalloc((void**)&devVolumetricStrainEnergy, sizeof(qeal) * simElementsNum);
		VectorX simElementZeroVec(simElementsNum);
		simElementZeroVec.setZero();
		cudaMemcpy(devVolumetricStrainEnergy, simElementZeroVec.data(), sizeof(qeal) * simElementsNum, cudaMemcpyHostToDevice);

		std::vector<int> hostClothPointsNeighborsList;
		std::vector<int> hostClothPointsNeighborsNum;
		std::vector<int> hostClothPointsNeighborsOffset;
		flatten2DArray(clothPointsNeighborsList, hostClothPointsNeighborsList, hostClothPointsNeighborsNum, hostClothPointsNeighborsOffset);

		cudaMalloc((void**)&devClothPointsNeighborsList, sizeof(int)* hostClothPointsNeighborsList.size());
		cudaMemcpy(devClothPointsNeighborsList, hostClothPointsNeighborsList.data(), sizeof(int)* hostClothPointsNeighborsList.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devClothPointsNeighborsNum, sizeof(int)* hostClothPointsNeighborsNum.size());
		cudaMemcpy(devClothPointsNeighborsNum, hostClothPointsNeighborsNum.data(), sizeof(int)* hostClothPointsNeighborsNum.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devClothPointsNeighborsOffset, sizeof(int)* hostClothPointsNeighborsOffset.size());
		cudaMemcpy(devClothPointsNeighborsOffset, hostClothPointsNeighborsOffset.data(), sizeof(int)* hostClothPointsNeighborsOffset.size(), cudaMemcpyHostToDevice);

		std::vector<qeal> hostClothPointStrainStiffness;
		std::vector<qeal> hostClothPointBendingStiffness;

		flatten2DArray(clothPointsStrainStiffness, hostClothPointStrainStiffness, hostClothPointsNeighborsNum, hostClothPointsNeighborsOffset);
		cudaMalloc((void**)&devClothPointsStrainStiffness, sizeof(qeal)* hostClothPointStrainStiffness.size());
		cudaMemcpy(devClothPointsStrainStiffness, hostClothPointStrainStiffness.data(), sizeof(qeal)* hostClothPointStrainStiffness.size(), cudaMemcpyHostToDevice);

		flatten2DArray(clothPointsBendingStiffness, hostClothPointBendingStiffness, hostClothPointsNeighborsNum, hostClothPointsNeighborsOffset);
		cudaMalloc((void**)&devClothPointsBendingStiffness, sizeof(qeal)* hostClothPointBendingStiffness.size());
		cudaMemcpy(devClothPointsBendingStiffness, hostClothPointBendingStiffness.data(), sizeof(qeal)* hostClothPointBendingStiffness.size(), cudaMemcpyHostToDevice);

		std::vector<qeal> hostClothPointNeighborsLbo;
		std::vector<qeal> hostClothPointNeighborsRestLength;

		flatten2DArray(clothPointsNeighborsLbo, hostClothPointNeighborsLbo, hostClothPointsNeighborsNum, hostClothPointsNeighborsOffset);
		cudaMalloc((void**)&devClothPointsNeighborsLbo, sizeof(qeal)* hostClothPointNeighborsLbo.size());
		cudaMemcpy(devClothPointsNeighborsLbo, hostClothPointNeighborsLbo.data(), sizeof(qeal)* hostClothPointNeighborsLbo.size(), cudaMemcpyHostToDevice);

		flatten2DArray(clothPointsNeighborsRestLength, hostClothPointNeighborsRestLength, hostClothPointsNeighborsNum, hostClothPointsNeighborsOffset);
		cudaMalloc((void**)&devClothPointsNeighborsRestLength, sizeof(qeal)* hostClothPointNeighborsRestLength.size());
		cudaMemcpy(devClothPointsNeighborsRestLength, hostClothPointNeighborsRestLength.data(), sizeof(qeal)* hostClothPointNeighborsRestLength.size(), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devClothPointsMeanCurvatureNorm, sizeof(qeal)* clothPointsMeanCurvatureNorm.size());
		cudaMemcpy(devClothPointsMeanCurvatureNorm, clothPointsMeanCurvatureNorm.data(), sizeof(qeal)* clothPointsMeanCurvatureNorm.size(), cudaMemcpyHostToDevice);



	}

	void PdIpcSimulator::gpuSolverInit()
	{
		VectorX simZeroVec(3 * totalPointsNum);
		simZeroVec.setZero();

		VectorXi simZeroIntVec(3 * totalPointsNum);
		simZeroIntVec.setZero();

		cudaMalloc((void**)&devInitRhs, sizeof(qeal) * simDims);
		cudaMemcpy(devInitRhs, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSysRhs, sizeof(qeal) * simDims);
		cudaMemcpy(devSysRhs, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devCollisionRhs, sizeof(qeal) * totalPointsNum * 3);
		cudaMemcpy(devCollisionRhs, simZeroVec.data(), sizeof(qeal) * totalPointsNum * 3, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSn, sizeof(qeal) * totalPointsNum * 3);
		cudaMemcpy(devSn, points.data(), sizeof(qeal) * totalPointsNum * 3, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devInertialDiagonal, sizeof(qeal) * simPointsNum);
		cudaMemcpy(devInertialDiagonal, simZeroVec.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devMassDiag, sizeof(qeal) * simPointsNum);
		cudaMemcpy(devMassDiag, simPointsMass.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devInvMassDiag, sizeof(qeal) * simPointsNum);
		cudaMemcpy(devInvMassDiag, simPointsInvMass.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devGravityForce, sizeof(qeal) * simDims);
		cudaMemcpy(devGravityForce, gravityForce.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);


		cudaMalloc((void**)&devMouseForce, sizeof(qeal) * simDims);
		cudaMemcpy(devMouseForce, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devPenaltyForce, sizeof(qeal) * 3 * totalPointsNum);
		cudaMemcpy(devPenaltyForce, simZeroVec.data(), sizeof(qeal) * 3 * totalPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devPointsPenaltyTimes, sizeof(int) * totalPointsNum);
		cudaMemcpy(devPointsPenaltyTimes, simZeroIntVec.data(), sizeof(int) * totalPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devElasticsDiagonal, sizeof(qeal) * simPointsNum);
		VectorX elasticsDiagonal = elasticsHess.diagonal();
		cudaMemcpy(devElasticsDiagonal, elasticsDiagonal.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devCollisionDiagonal, sizeof(qeal) * totalPointsNum);
		cudaMemcpy(devCollisionDiagonal, simZeroVec.data(), sizeof(qeal) * totalPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devSysDiagonal, sizeof(qeal) * simPointsNum);
		cudaMemcpy(devSysDiagonal, simZeroVec.data(), sizeof(qeal) * simPointsNum, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devFirstTerm, sizeof(qeal) * simDims);
		cudaMemcpy(devFirstTerm, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devTempTerm, sizeof(qeal) * simDims);
		cudaMemcpy(devTempTerm, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devDinvRhs, sizeof(qeal) * simDims);
		cudaMemcpy(devDinvRhs, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devXk, sizeof(qeal) * simDims);
		cudaMemcpy(devXk, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devPreXk, sizeof(qeal) * simDims);
		cudaMemcpy(devPreXk, simZeroVec.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devOmega, sizeof(qeal));
		cudaMemcpy(devOmega, simZeroVec.data(), sizeof(qeal), cudaMemcpyHostToDevice);
		///////////////////////////////////////////////////////////////////////////////////////
		SparseMatrix R1, R2, R3, R4, R5;
		std::vector<std::vector<int>> hessR1Index(simPointsNum);
		std::vector<std::vector<qeal>> hessR1Value(simPointsNum);
		std::vector<std::vector<int>> hessR1DiiIndex(simPointsNum);
		std::vector<TripletX> tripletR1;


		for (int i = 0; i < simPointsNum; i++)
		{
			hessR1Index[i].push_back(i);
			hessR1Value[i].push_back(0);
			hessR1DiiIndex[i].push_back(i);
		}
#define IG 1e-12
		for (int k = 0; k < sysHess.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(sysHess, k); it; ++it)
			{
				qeal value = it.value();
				int row_id = it.row();
				int col_id = it.col();
				if (std::abs(value) < IG)
					continue;
				if (row_id != col_id)
				{
					hessR1Index[row_id].push_back(col_id);
					hessR1Value[row_id].push_back(-value);
					hessR1DiiIndex[row_id].push_back(row_id);
					tripletR1.push_back(TripletX(row_id, col_id, -value));
				}
				else
				{
					tripletR1.push_back(TripletX(row_id, col_id, 0.0));
				}

			}
		}

		R1.resize(simPointsNum, simPointsNum);
		R1.setFromTriplets(tripletR1.begin(), tripletR1.end());

		std::vector<int> hostHessR1Index;
		std::vector<int> hostHessR1DiiIndex;
		std::vector<qeal> hostHessR1Value;
		std::vector<int> hostHessR1IndexNum;
		std::vector<int> hostHessR1IndexOffset;

		flatten2DArray(hessR1Index, hostHessR1Index, hostHessR1IndexNum, hostHessR1IndexOffset);
		flatten2DArray(hessR1Value, hostHessR1Value, hostHessR1IndexNum, hostHessR1IndexOffset);
		flatten2DArray(hessR1DiiIndex, hostHessR1DiiIndex, hostHessR1IndexNum, hostHessR1IndexOffset);

		CUDA_CALL(cudaMalloc((void**)&devHessR1Index, hostHessR1Index.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR1Index, hostHessR1Index.data(), sizeof(int) * hostHessR1Index.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR1Value, hostHessR1Value.size() * sizeof(qeal)));
		CUDA_CALL(cudaMemcpy(devHessR1Value, hostHessR1Value.data(), sizeof(qeal) * hostHessR1Value.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR1DiiIndex, hostHessR1DiiIndex.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR1DiiIndex, hostHessR1DiiIndex.data(), sizeof(int) * hostHessR1DiiIndex.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR1IndexNum, hostHessR1IndexNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR1IndexNum, hostHessR1IndexNum.data(), sizeof(int) * hostHessR1IndexNum.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR1IndexOffset, hostHessR1IndexOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR1IndexOffset, hostHessR1IndexOffset.data(), sizeof(int) * hostHessR1IndexOffset.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR1Operators, hostHessR1Value.size() * sizeof(qeal)));

		hostHessR1OperatorsNum = hostHessR1Value.size();
		CUDA_CALL(cudaMalloc((void**)&devHessR1OperatorsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR1OperatorsNum, &hostHessR1OperatorsNum, sizeof(int), cudaMemcpyHostToDevice));


		//////////////////////////////////////////////////////////

		R2 = R1 * R1;
		std::vector<std::vector<int>> hessR2Index(simPointsNum);
		std::vector<std::vector<qeal>> hessR2Value(simPointsNum);
		std::vector<std::vector<bool>> hessR2IndexValid(simPointsNum);

		for (int k = 0; k < R2.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(R2, k); it; ++it)
			{
				qeal value = it.value();
				int row_id = it.row();
				int col_id = it.col();
				if (std::abs(value) < IG)
					continue;
				hessR2Index[row_id].push_back(col_id);
				hessR2IndexValid[row_id].push_back(true);
				hessR2Value[row_id].push_back(value);
			}
		}


		std::vector<std::vector<int>> R2DiiList(simPointsNum);
		std::vector<std::vector<int>> R2DjjList(simPointsNum);
		std::vector<std::vector<qeal>> R2AijkList(simPointsNum);
		std::vector<std::vector<int>> hessR2OperatorsNumList(simPointsNum);
		std::vector<std::vector<qeal>> hessR2ConstOrderOperators(simPointsNum);
		std::vector<std::vector<int>> hessR2FirstOrderOperators(simPointsNum);

		for (int i = 0; i < simPointsNum; i++)
		{
			std::vector<int>& refIndex = hessR2Index[i];
			std::vector<int>& rowIndex = hessR1Index[i];
			for (int refId = 0; refId < refIndex.size(); refId++)
			{
				std::vector<int>& colIndex = hessR1Index[refIndex[refId]];

				int r2num = 0;
				for (int r = 0; r < rowIndex.size(); r++)
				{
					int rid = rowIndex[r];
					for (int c = 0; c < colIndex.size(); c++)
					{
						int cid = colIndex[c];
						if (rid != cid) // The diag value is zero
							continue;
						//if (R1blendingFlag[i][r].second && R1blendingFlag[refIndex[refId]][c].second)
						//	continue;
						qeal ajk = hessR1Value[i][r] * hessR1Value[refIndex[refId]][c];
						if (std::abs(ajk) < IG)
							continue;
						R2DiiList[i].push_back(i);
						R2DjjList[i].push_back(cid);
						R2AijkList[i].push_back(ajk);

						r2num++;
					}
				}
				if (r2num == 0)
					hessR2IndexValid[i][refId] = false;
				else hessR2OperatorsNumList[i].push_back(r2num);
			}
		}

		for (int i = 0; i < simPointsNum; i++)
		{
			std::vector<int> oldIndex = hessR2Index[i];
			std::vector<int> newIndex;
			for (int j = 0; j < oldIndex.size(); j++)
			{
				if (hessR2IndexValid[i][j])
					newIndex.push_back(oldIndex[j]);
			}
			hessR2Index[i] = newIndex;
			hessR2ConstOrderOperators[i].resize(newIndex.size(), 0.0);
			hessR2FirstOrderOperators[i].resize(newIndex.size(), -1);

			int* R1OperatorsIndex = hostHessR1Index.data() + hostHessR1IndexOffset[i];
			int R1Num = hostHessR1IndexNum[i];

			for (int j = 0; j < newIndex.size(); j++)
			{
				if (i == newIndex[j])
					hessR2ConstOrderOperators[i][j] = 1.0;
				for (int k = 0; k < R1Num; k++)
				{
					if (newIndex[j] == R1OperatorsIndex[k])
					{
						hessR2FirstOrderOperators[i][j] = hostHessR1IndexOffset[i] + k;
						break;
					}
				}
			}
		}



		std::vector<int> hostHessR2Index;
		std::vector<int> hostHessR2IndexNum;
		std::vector<int> hostHessR2IndexOffset;
		flatten2DArray(hessR2Index, hostHessR2Index, hostHessR2IndexNum, hostHessR2IndexOffset);

		CUDA_CALL(cudaMalloc((void**)&devHessR2Index, hostHessR2Index.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2Index, hostHessR2Index.data(), sizeof(int) * hostHessR2Index.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR2IndexNum, hostHessR2IndexNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2IndexNum, hostHessR2IndexNum.data(), sizeof(int) * hostHessR2IndexNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR2IndexOffset, hostHessR2IndexOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2IndexOffset, hostHessR2IndexOffset.data(), sizeof(int) * hostHessR2IndexOffset.size(), cudaMemcpyHostToDevice));

		std::vector<qeal> hostHessR2ConstOrderOperators;
		flatten2DArray(hessR2ConstOrderOperators, hostHessR2ConstOrderOperators, hostHessR2IndexNum, hostHessR2IndexOffset);

		CUDA_CALL(cudaMalloc((void**)&devHessR2ConstOrderOperators, hostHessR2ConstOrderOperators.size() * sizeof(qeal)));
		CUDA_CALL(cudaMemcpy(devHessR2ConstOrderOperators, hostHessR2ConstOrderOperators.data(), sizeof(qeal) * hostHessR2ConstOrderOperators.size(), cudaMemcpyHostToDevice));


		std::vector<int> hostHessR2FirstOrderOperators;
		flatten2DArray(hessR2FirstOrderOperators, hostHessR2FirstOrderOperators, hostHessR2IndexNum, hostHessR2IndexOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR2FirstOrderOperators, hostHessR2FirstOrderOperators.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2FirstOrderOperators, hostHessR2FirstOrderOperators.data(), sizeof(int) * hostHessR2FirstOrderOperators.size(), cudaMemcpyHostToDevice));

		//

		std::vector<int> hostHessR2DiiIndex;
		std::vector<int> hostHessR2DjjIndex;
		std::vector<qeal> hostHessR2Aijk;
		std::vector<int> hostHessR2AijkNum;
		std::vector<int> hostHessR2AijkOffset;

		flatten2DArray(R2DiiList, hostHessR2DiiIndex, hostHessR2AijkNum, hostHessR2AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR2DiiIndex, hostHessR2DiiIndex.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2DiiIndex, hostHessR2DiiIndex.data(), sizeof(int) * hostHessR2DiiIndex.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR2AijkNum, hostHessR2AijkNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2AijkNum, hostHessR2AijkNum.data(), sizeof(int) * hostHessR2AijkNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR2AijkOffset, hostHessR2AijkOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2AijkOffset, hostHessR2AijkOffset.data(), sizeof(int) * hostHessR2AijkOffset.size(), cudaMemcpyHostToDevice));

		flatten2DArray(R2DjjList, hostHessR2DjjIndex, hostHessR2AijkNum, hostHessR2AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR2DjjIndex, hostHessR2DjjIndex.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2DjjIndex, hostHessR2DjjIndex.data(), sizeof(int) * hostHessR2DjjIndex.size(), cudaMemcpyHostToDevice));

		flatten2DArray(R2AijkList, hostHessR2Aijk, hostHessR2AijkNum, hostHessR2AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR2Aijk, hostHessR2Aijk.size() * sizeof(qeal)));
		CUDA_CALL(cudaMemcpy(devHessR2Aijk, hostHessR2Aijk.data(), sizeof(qeal) * hostHessR2Aijk.size(), cudaMemcpyHostToDevice));

		hostHessR2OperatorsNum = hostHessR2Aijk.size();
		CUDA_CALL(cudaMalloc((void**)&devHessR2OperatorsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2OperatorsNum, &hostHessR2OperatorsNum, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR2Operators, sizeof(qeal) * hostHessR2OperatorsNum));

		std::vector<int> hostHessR2RowOperatorsNum;
		std::vector<int> hostHessR2RowOperatorsOffset;
		int r2num = 0;
		int r2Offset = 0;
		for (int i = 0; i < hessR2OperatorsNumList.size(); i++)
		{
			for (int j = 0; j < hessR2OperatorsNumList[i].size(); j++)
			{
				hostHessR2RowOperatorsOffset.push_back(r2Offset);
				int num = hessR2OperatorsNumList[i][j];
				hostHessR2RowOperatorsNum.push_back(num);
				r2Offset += num;
			}
		}


		CUDA_CALL(cudaMalloc((void**)&devHessR2RowOperatorsNum, hostHessR2RowOperatorsNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2RowOperatorsNum, hostHessR2RowOperatorsNum.data(), sizeof(int) * hostHessR2RowOperatorsNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR2RowOperatorsOffset, hostHessR2RowOperatorsOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2RowOperatorsOffset, hostHessR2RowOperatorsOffset.data(), sizeof(int) * hostHessR2RowOperatorsOffset.size(), cudaMemcpyHostToDevice));

		hostHessR2RowOperatorsTotal = hostHessR2RowOperatorsOffset.size();
		CUDA_CALL(cudaMalloc((void**)&devHessR2RowOperatorsTotal, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR2RowOperatorsTotal, &hostHessR2RowOperatorsTotal, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR2RowOperators, hostHessR2RowOperatorsOffset.size() * sizeof(qeal)));

		///////////////////////////////////////////////////////


		hostR1MaxOperatorIndexNum = 0;
		for (int i = 0; i < hostHessR1IndexNum.size(); i++)
			if (hostR1MaxOperatorIndexNum < hostHessR1IndexNum[i])
				hostR1MaxOperatorIndexNum = hostHessR1IndexNum[i];
		CUDA_CALL(cudaMalloc((void**)&devR1MaxOperatorIndexNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR1MaxOperatorIndexNum, &hostR1MaxOperatorIndexNum, sizeof(int), cudaMemcpyHostToDevice));

		hostR2MaxOperatorIndexNum = 0;
		for (int i = 0; i < hostHessR2IndexNum.size(); i++)
			if (hostR2MaxOperatorIndexNum < hostHessR2IndexNum[i])
				hostR2MaxOperatorIndexNum = hostHessR2IndexNum[i];
		CUDA_CALL(cudaMalloc((void**)&devR2MaxOperatorIndexNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR2MaxOperatorIndexNum, &hostR2MaxOperatorIndexNum, sizeof(int), cudaMemcpyHostToDevice));

		////////////////////////////////////////////
		hostR2AlignedPointsSharedMemSize = (THREADS_NUM * hostR2MaxOperatorIndexNum) / R2_SPAN;
		CUDA_CALL(cudaMalloc((void**)&devR2AlignedPointsSharedMemSize, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR2AlignedPointsSharedMemSize, &hostR2AlignedPointsSharedMemSize, sizeof(int), cudaMemcpyHostToDevice));

		hostR1AlignedDimsSharedMemSize = (THREADS_NUM_576 * hostR1MaxOperatorIndexNum) / (R1_DIM_SPAN);
		CUDA_CALL(cudaMalloc((void**)&devR1AlignedDimsSharedMemSize, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR1AlignedDimsSharedMemSize, &hostR1AlignedDimsSharedMemSize, sizeof(int), cudaMemcpyHostToDevice));

		hostR2AlignedDimsSharedMemSize = (THREADS_NUM_576 * hostR2MaxOperatorIndexNum) / (R2_DIM_SPAN);
		CUDA_CALL(cudaMalloc((void**)&devR2AlignedDimsSharedMemSize, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR2AlignedDimsSharedMemSize, &hostR2AlignedDimsSharedMemSize, sizeof(int), cudaMemcpyHostToDevice));
		////////////////////////////////////////////////////////////////////

		//R3
		R3 = R2 * R1;
		std::vector<std::vector<int>> hessR3Index(simPointsNum);
		std::vector<std::vector<bool>> hessR3IndexValid(simPointsNum);

		for (int k = 0; k < R3.outerSize(); ++k)
		{
			for (SparseMatrix::InnerIterator it(R3, k); it; ++it)
			{
				qeal value = it.value();
				int row_id = it.row();
				int col_id = it.col();
				if (std::abs(value) < IG)
					continue;
				hessR3Index[row_id].push_back(col_id);
				hessR3IndexValid[row_id].push_back(true);
			}
		}

		std::vector<std::vector<int>> R3DiiList(simPointsNum);
		std::vector<std::vector<int>> R3DjjList(simPointsNum);
		std::vector<std::vector<qeal>> R3AijkList(simPointsNum);

		std::vector<std::vector<int>> R3RowR2OperatorIndexList(simPointsNum);
		std::vector<std::vector<int>> R3ColR1OperatorIndexList(simPointsNum);

		std::vector<std::vector<int>> hessR3OperatorsNumList(simPointsNum);
		std::vector<std::vector<qeal>> hessR3ConstOrderOperators(simPointsNum);
		std::vector<std::vector<int>> hessR3FirstOrderOperators(simPointsNum);

		for (int i = 0; i < simPointsNum; i++)
		{
			std::vector<int>& refIndex = hessR3Index[i];
			std::vector<int>& rowIndex = hessR2Index[i];
			for (int refId = 0; refId < refIndex.size(); refId++)
			{
				int refColId = refIndex[refId];
				std::vector<int>& colIndex = hessR1Index[refColId];
				int r3num = 0;
				for (int r = 0; r < rowIndex.size(); r++)
				{
					int rid = rowIndex[r];
					for (int c = 0; c < colIndex.size(); c++)
					{
						int cid = colIndex[c];
						if (rid != cid) // The diag value is zero
							continue;
						qeal ajk = hessR2Value[i][r] * hessR1Value[refColId][c];
						if (std::abs(ajk) < IG)
							continue;
						R3DiiList[i].push_back(hostHessR2IndexOffset[i] + r);
						R3DjjList[i].push_back(cid);
						R3AijkList[i].push_back(hessR1Value[refColId][c]);

						R3RowR2OperatorIndexList[i].push_back(hostHessR2IndexOffset[i] + r);
						R3ColR1OperatorIndexList[i].push_back(hostHessR1IndexOffset[refColId] + c);

						r3num++;
					}
				}
				if (r3num == 0)
					hessR3IndexValid[i][refId] = false;
				else hessR3OperatorsNumList[i].push_back(r3num);
			}
		}



		for (int i = 0; i < simPointsNum; i++)
		{
			std::vector<int> oldIndex = hessR3Index[i];
			std::vector<int> newIndex;
			for (int j = 0; j < oldIndex.size(); j++)
			{
				if (hessR3IndexValid[i][j])
					newIndex.push_back(oldIndex[j]);
			}
			hessR3Index[i] = newIndex;
			hessR3ConstOrderOperators[i].resize(newIndex.size(), 0.0);
			hessR3FirstOrderOperators[i].resize(newIndex.size(), -1);
			//
			int* R2OperatorsIndex = hostHessR2Index.data() + hostHessR2IndexOffset[i];
			int R2Num = hostHessR2IndexNum[i];

			for (int j = 0; j < newIndex.size(); j++)
			{
				if (i == newIndex[j])
					hessR3ConstOrderOperators[i][j] = 1.0;
				for (int k = 0; k < R2Num; k++)
				{
					if (newIndex[j] == R2OperatorsIndex[k])
					{
						hessR3FirstOrderOperators[i][j] = hostHessR2IndexOffset[i] + k;
						break;
					}
				}
			}
		}

		std::vector<int> hostHessR3Index;
		std::vector<int> hostHessR3IndexNum;
		std::vector<int> hostHessR3IndexOffset;
		flatten2DArray(hessR3Index, hostHessR3Index, hostHessR3IndexNum, hostHessR3IndexOffset);

		CUDA_CALL(cudaMalloc((void**)&devHessR3Index, hostHessR3Index.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3Index, hostHessR3Index.data(), sizeof(int) * hostHessR3Index.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR3IndexNum, hostHessR3IndexNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3IndexNum, hostHessR3IndexNum.data(), sizeof(int) * hostHessR3IndexNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR3IndexOffset, hostHessR3IndexOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3IndexOffset, hostHessR3IndexOffset.data(), sizeof(int) * hostHessR3IndexOffset.size(), cudaMemcpyHostToDevice));

		std::vector<qeal> hostHessR3ConstOrderOperators;
		flatten2DArray(hessR3ConstOrderOperators, hostHessR3ConstOrderOperators, hostHessR3IndexNum, hostHessR3IndexOffset);

		CUDA_CALL(cudaMalloc((void**)&devHessR3ConstOrderOperators, hostHessR3ConstOrderOperators.size() * sizeof(qeal)));
		CUDA_CALL(cudaMemcpy(devHessR3ConstOrderOperators, hostHessR3ConstOrderOperators.data(), sizeof(qeal) * hostHessR3ConstOrderOperators.size(), cudaMemcpyHostToDevice));


		std::vector<int> hostHessR3FirstOrderOperators;
		flatten2DArray(hessR3FirstOrderOperators, hostHessR3FirstOrderOperators, hostHessR3IndexNum, hostHessR3IndexOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR3FirstOrderOperators, hostHessR3FirstOrderOperators.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3FirstOrderOperators, hostHessR3FirstOrderOperators.data(), sizeof(int) * hostHessR3FirstOrderOperators.size(), cudaMemcpyHostToDevice));

		//

		std::vector<int> hostHessR3DiiIndex;
		std::vector<int> hostHessR3DjjIndex;
		std::vector<qeal> hostHessR3Aijk;
		std::vector<int> hostHessR3AijkNum;
		std::vector<int> hostHessR3AijkOffset;

		flatten2DArray(R3DiiList, hostHessR3DiiIndex, hostHessR3AijkNum, hostHessR3AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR3DiiIndex, hostHessR3DiiIndex.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3DiiIndex, hostHessR3DiiIndex.data(), sizeof(int) * hostHessR3DiiIndex.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR3AijkNum, hostHessR3AijkNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3AijkNum, hostHessR3AijkNum.data(), sizeof(int) * hostHessR3AijkNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR3AijkOffset, hostHessR3AijkOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3AijkOffset, hostHessR3AijkOffset.data(), sizeof(int) * hostHessR3AijkOffset.size(), cudaMemcpyHostToDevice));

		flatten2DArray(R3DjjList, hostHessR3DjjIndex, hostHessR3AijkNum, hostHessR3AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR3DjjIndex, hostHessR3DjjIndex.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3DjjIndex, hostHessR3DjjIndex.data(), sizeof(int) * hostHessR3DjjIndex.size(), cudaMemcpyHostToDevice));

		flatten2DArray(R3AijkList, hostHessR3Aijk, hostHessR3AijkNum, hostHessR3AijkOffset);
		CUDA_CALL(cudaMalloc((void**)&devHessR3Aijk, hostHessR3Aijk.size() * sizeof(qeal)));
		CUDA_CALL(cudaMemcpy(devHessR3Aijk, hostHessR3Aijk.data(), sizeof(qeal) * hostHessR3Aijk.size(), cudaMemcpyHostToDevice));

		hostHessR3OperatorsNum = hostHessR3Aijk.size();
		CUDA_CALL(cudaMalloc((void**)&devHessR3OperatorsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3OperatorsNum, &hostHessR3OperatorsNum, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR3Operators, sizeof(qeal) * hostHessR3OperatorsNum));

		std::vector<int> hostHessR3RowOperatorsNum;
		std::vector<int> hostHessR3RowOperatorsOffset;
		int r3num = 0;
		int r3Offset = 0;
		for (int i = 0; i < hessR3OperatorsNumList.size(); i++)
		{
			for (int j = 0; j < hessR3OperatorsNumList[i].size(); j++)
			{
				hostHessR3RowOperatorsOffset.push_back(r3Offset);
				int num = hessR3OperatorsNumList[i][j];
				hostHessR3RowOperatorsNum.push_back(num);
				r3Offset += num;
			}
		}


		CUDA_CALL(cudaMalloc((void**)&devHessR3RowOperatorsNum, hostHessR3RowOperatorsNum.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3RowOperatorsNum, hostHessR3RowOperatorsNum.data(), sizeof(int) * hostHessR3RowOperatorsNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devHessR3RowOperatorsOffset, hostHessR3RowOperatorsOffset.size() * sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3RowOperatorsOffset, hostHessR3RowOperatorsOffset.data(), sizeof(int) * hostHessR3RowOperatorsOffset.size(), cudaMemcpyHostToDevice));

		hostHessR3RowOperatorsTotal = hostHessR3RowOperatorsOffset.size();
		CUDA_CALL(cudaMalloc((void**)&devHessR3RowOperatorsTotal, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devHessR3RowOperatorsTotal, &hostHessR3RowOperatorsTotal, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devHessR3RowOperators, hostHessR3RowOperatorsOffset.size() * sizeof(qeal)));

		CUDA_CALL(cudaMalloc((void**)&devHessR21RowOperators, hostHessR2RowOperatorsOffset.size() * sizeof(qeal)));

		///////////////////////////////////////////////////////



		hostR3MaxOperatorIndexNum = 0;
		for (int i = 0; i < hostHessR3IndexNum.size(); i++)
			if (hostR3MaxOperatorIndexNum < hostHessR3IndexNum[i])
				hostR3MaxOperatorIndexNum = hostHessR3IndexNum[i];
		CUDA_CALL(cudaMalloc((void**)&devR3MaxOperatorIndexNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR3MaxOperatorIndexNum, &hostR3MaxOperatorIndexNum, sizeof(int), cudaMemcpyHostToDevice));

		////////////////////////////////////////////
		hostR3AlignedPointsSharedMemSize = (THREADS_NUM * hostR3MaxOperatorIndexNum) / R3_SPAN;
		CUDA_CALL(cudaMalloc((void**)&devR3AlignedPointsSharedMemSize, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR3AlignedPointsSharedMemSize, &hostR3AlignedPointsSharedMemSize, sizeof(int), cudaMemcpyHostToDevice));

		hostR3AlignedDimsSharedMemSize = (THREADS_NUM_576 * hostR3MaxOperatorIndexNum) / (R3_DIM_SPAN);
		CUDA_CALL(cudaMalloc((void**)&devR3AlignedDimsSharedMemSize, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devR3AlignedDimsSharedMemSize, &hostR3AlignedDimsSharedMemSize, sizeof(int), cudaMemcpyHostToDevice));

		//////////////////////////////////////////
		R1OperatorsBlockSize = THREADS_NUM_128;
		R1OperatorsGridSize = (hostHessR1OperatorsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128;

		R1IteraionBlockSize = THREADS_NUM_192;
		R1IteraionGridSize = (simDims + (THREADS_NUM_192 - 1)) / THREADS_NUM_192;


		R2OperatorsBlockSize = THREADS_NUM_128;
		R2OperatorsGridSize = (hostHessR2OperatorsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128;
		R2RowOperatorsBlockSize = THREADS_NUM_128;
		R2RowOperatorsGridSize = (hostHessR2RowOperatorsTotal + (THREADS_NUM_128 - 1)) / THREADS_NUM_128;

		R2FirstTermBlockSize = THREADS_NUM_192;
		R2FirstTermGridSize = (simDims * R2_SPAN + (THREADS_NUM_192 - 1)) / THREADS_NUM_192;

		R2IteraionBlockSize = THREADS_NUM_192;
		R2IteraionGridSize = (simDims * R2_SPAN + (THREADS_NUM_192 - 1)) / THREADS_NUM_192;

		R3OperatorsBlockSize = THREADS_NUM_128;
		R3OperatorsGridSize = (hostHessR3OperatorsNum + (THREADS_NUM_128 - 1)) / THREADS_NUM_128;
		R3RowOperatorsBlockSize = THREADS_NUM_128;
		R3RowOperatorsGridSize = (hostHessR3RowOperatorsTotal + (THREADS_NUM_128 - 1)) / THREADS_NUM_128;

		R3FirstTermBlockSize = THREADS_NUM_192;
		R3FirstTermGridSize = (simDims * R3_SPAN + (THREADS_NUM_192 - 1)) / THREADS_NUM_192;

		R3IteraionBlockSize = THREADS_NUM_192;
		R3IteraionGridSize = (simDims * R3_SPAN + (THREADS_NUM_192 - 1)) / THREADS_NUM_192;
	}

	void PdIpcSimulator::gpuCollisionInit()
	{
		enableAllEC = true;
		qeal inflation_radius = 0.5 * dHat;
		//m_patches
		std::vector<std::vector<int>> facesSharedEdges(totalBoundaryFacesNum);
		std::vector<std::vector<int>> facesSharedPoints(totalBoundaryFacesNum);
		std::vector<bool> boundaryEdgesFlag(totalBoundaryEdgesNum, false);
		int pointOffset = 0;
		int faceOffset = 0;
		int edgeOffset = 0;
		for (int i = 0; i < models.size() + staticModels.size(); i++)
		{
			PdIpcModel* m;
			if (i < models.size())
				m = getModel(i);
			else m = getStaticModel(i - models.size());
			//m_patches
			for (int j = 0; j < m->collisionPatches.size(); j++)
			{
				std::vector<int>& patch = m->collisionPatches[j];
				std::vector<int> temp;
				for (int k = 0; k < patch.size(); k++)
					temp.push_back(int(patch[k] + faceOffset));
				if (temp.size() > 0)
				{
					m_patchesFaces.push_back(temp);
					QColor color;
					color = QColor::colorNames()[j % 34];
					m_patchesFacesColor.push_back(color);
					m->collisionPatches[j].clear();
				}
			}
			
			for (int j = 0; j < m->boundaryFacesNum; j++)
			{
				std::vector<int>& list = m->boundaryFacesSharedEdges[j];
				std::vector<int> temp;
				for (int k = 0; k < list.size(); k++)
					temp.push_back(int(list[k] + edgeOffset));
				facesSharedEdges[faceOffset + j] = temp;

				list = m->boundaryFacesSharedPoints[j];
				temp.clear();
				for (int k = 0; k < list.size(); k++)
					temp.push_back(int(list[k] + pointOffset));
				facesSharedPoints[faceOffset + j] = temp;
			}

			for (int j = 0; j < m->boundaryEdgesNum; j++)
				boundaryEdgesFlag[j + edgeOffset] = m->boundaryEdgesFlag[j];

			pointOffset += m->usePointsNum;
			faceOffset += m->boundaryFacesNum;
			edgeOffset += m->boundaryEdgesNum;
		}

		m_patchesEdges.resize(m_patchesFaces.size());
		m_patchesPoints.resize(m_patchesFaces.size());
		m_patchesExPoints.resize(m_patchesFaces.size());

		for (int i = 0; i < m_patchesFaces.size(); i++)
		{
			std::set<int> edgeSet, pointSet, pointExSet;
			for (int j = 0; j < m_patchesFaces[i].size(); j++)
			{
				int fid = m_patchesFaces[i][j];

				Vector3i face = boundaryFaces.col(fid);
				for (int k = 0; k < facesSharedEdges[fid].size(); k++)
				{
					int eid = facesSharedEdges[fid][k];
					if (enableAllEC)
						edgeSet.insert(eid);
					else
					{
						if (boundaryEdgesFlag[eid])
							edgeSet.insert(eid);
					}
				}

				for (int k = 0; k < facesSharedPoints[fid].size(); k++)
					pointSet.insert(facesSharedPoints[fid][k]);
				pointExSet.insert(face[0]);
				pointExSet.insert(face[1]);
				pointExSet.insert(face[2]);
			}

			m_patchesEdges[i].resize(edgeSet.size());
			m_patchesPoints[i].resize(pointSet.size());
			m_patchesExPoints[i].resize(pointExSet.size());


			std::copy(edgeSet.begin(), edgeSet.end(), m_patchesEdges[i].begin());
			std::copy(pointSet.begin(), pointSet.end(), m_patchesPoints[i].begin());
			std::copy(pointExSet.begin(), pointExSet.end(), m_patchesExPoints[i].begin());
		}


		hostPatchNum = m_patchesFaces.size();
		cudaMalloc((void**)&devPatchNum, sizeof(int));
		cudaMemcpy(devPatchNum, &hostPatchNum, sizeof(int), cudaMemcpyHostToDevice);

		flatten2DArray(m_patchesFaces, hostPatchFacesList, hostPatchFacesNum, hostPatchFacesOffset);
		int d = totalFacesNum;
		hostFacesToPatchIndex.resize(totalBoundaryFacesNum);
		for (int i = 0; i < hostPatchFacesList.size(); i++)
		{
			hostFacesToPatchIndex[hostPatchFacesList[i]] = i;
		}

		CUDA_CALL(cudaMalloc((void**)&devFacesToPatchIndex, sizeof(int) * hostFacesToPatchIndex.size()));
		CUDA_CALL(cudaMemcpy(devFacesToPatchIndex, hostFacesToPatchIndex.data(), sizeof(int) * hostFacesToPatchIndex.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchFacesList, sizeof(int) * hostPatchFacesList.size()));
		CUDA_CALL(cudaMemcpy(devPatchFacesList, hostPatchFacesList.data(), sizeof(int) * hostPatchFacesList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchFacesNum, sizeof(int) * hostPatchFacesNum.size()));
		CUDA_CALL(cudaMemcpy(devPatchFacesNum, hostPatchFacesNum.data(), sizeof(int) * hostPatchFacesNum.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchFacesOffset, sizeof(int) * hostPatchFacesOffset.size()));
		CUDA_CALL(cudaMemcpy(devPatchFacesOffset, hostPatchFacesOffset.data(), sizeof(int) * hostPatchFacesOffset.size(), cudaMemcpyHostToDevice));

		flatten2DArray(m_patchesEdges, hostPatchEdgesList, hostPatchEdgesNum, hostPatchEdgesOffset);
		CUDA_CALL(cudaMalloc((void**)&devPatchEdgesList, sizeof(int) * hostPatchEdgesList.size()));
		CUDA_CALL(cudaMemcpy(devPatchEdgesList, hostPatchEdgesList.data(), sizeof(int) * hostPatchEdgesList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchEdgesNum, sizeof(int) * hostPatchEdgesNum.size()));
		CUDA_CALL(cudaMemcpy(devPatchEdgesNum, hostPatchEdgesNum.data(), sizeof(int) * hostPatchEdgesNum.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchEdgesOffset, sizeof(int) * hostPatchEdgesOffset.size()));
		CUDA_CALL(cudaMemcpy(devPatchEdgesOffset, hostPatchEdgesOffset.data(), sizeof(int) * hostPatchEdgesOffset.size(), cudaMemcpyHostToDevice));

		flatten2DArray(m_patchesPoints, hostPatchPointsList, hostPatchPointsNum, hostPatchPointsOffset);
		CUDA_CALL(cudaMalloc((void**)&devPatchPointsList, sizeof(int) * hostPatchPointsList.size()));
		CUDA_CALL(cudaMemcpy(devPatchPointsList, hostPatchPointsList.data(), sizeof(int) * hostPatchPointsList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchPointsNum, sizeof(int) * hostPatchPointsNum.size()));
		CUDA_CALL(cudaMemcpy(devPatchPointsNum, hostPatchPointsNum.data(), sizeof(int) * hostPatchPointsNum.size(), cudaMemcpyHostToDevice));

		CUDA_CALL(cudaMalloc((void**)&devPatchPointsOffset, sizeof(int) * hostPatchPointsOffset.size()));
		CUDA_CALL(cudaMemcpy(devPatchPointsOffset, hostPatchPointsOffset.data(), sizeof(int) * hostPatchPointsOffset.size(), cudaMemcpyHostToDevice));

		flatten2DArray(m_patchesExPoints, hostPatchExPointsList, hostPatchExPointsNum, hostPatchExPointsOffset);
		CUDA_CALL(cudaMalloc((void**)&devPatchExPointsList, sizeof(int) * hostPatchExPointsList.size()));
		CUDA_CALL(cudaMemcpy(devPatchExPointsList, hostPatchExPointsList.data(), sizeof(int) * hostPatchExPointsList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchExPointsNum, sizeof(int) * hostPatchExPointsNum.size()));
		CUDA_CALL(cudaMemcpy(devPatchExPointsNum, hostPatchExPointsNum.data(), sizeof(int) * hostPatchExPointsNum.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchExPointsOffset, sizeof(int) * hostPatchExPointsOffset.size()));
		CUDA_CALL(cudaMemcpy(devPatchExPointsOffset, hostPatchExPointsOffset.data(), sizeof(int) * hostPatchExPointsOffset.size(), cudaMemcpyHostToDevice));

		//
		hostFacesBboxes.resize(6 * totalBoundaryFacesNum, 0);
		CUDA_CALL(cudaMalloc((void**)&devFacesBboxes, sizeof(qeal) * 6 * totalBoundaryFacesNum));
		CUDA_CALL(cudaMemcpy(devFacesBboxes, hostFacesBboxes.data(), sizeof(qeal) * 6 * totalBoundaryFacesNum, cudaMemcpyHostToDevice));

		hostPatchMaxThreads = 0;
		for (int i = 0; i < m_patchesFaces.size(); i++)
			if (m_patchesFaces[i].size() > hostPatchMaxThreads)
				hostPatchMaxThreads = m_patchesFaces[i].size();
		for (int i = 0; 1; i++)
		{
			int base = 2;
			for (int j = 0; j < i; j++)
				base *= 2;
			if (hostPatchMaxThreads <= base)
			{
				hostPatchMaxThreads = base;
				break;
			}
		}

		hostPatchBboxes.resize(6 * hostPatchNum, 0);
		cudaMalloc((void**)&devPatchBboxes, sizeof(qeal) * 6 * hostPatchNum);
		cudaMemcpy(devPatchBboxes, hostPatchBboxes.data(), sizeof(qeal) * 6 * hostPatchNum, cudaMemcpyHostToDevice);

		m_patchBbox.resize(hostPatchNum); 
		tbb::parallel_for(size_t(0), size_t(hostPatchNum), [&](size_t i)
		{
			int id = i;
			m_patchBbox[id][0] = Vector3(DBL_MAX, DBL_MAX, DBL_MAX);
			m_patchBbox[id][1] = -Vector3(DBL_MAX, DBL_MAX, DBL_MAX);
			for (int j = 0; j < m_patchesExPoints[id].size(); j++)
			{
				Vector3 p = points.col(m_patchesExPoints[id][j]);
				for (int k = 0; k < 3; k++)
				{
					m_patchBbox[id][0][k] = min2(m_patchBbox[id][0][k], p[k] - inflation_radius);
					m_patchBbox[id][1][k] = max2(m_patchBbox[id][1][k], p[k] + inflation_radius);
				}
			}
		});

		AipcBvhs patchBvhs;
		patchBvhs.init(m_patchBbox);

		hostBvhsNodesNum = patchBvhs.boxlist.size();
		hostBvhsNodesBbox.resize(patchBvhs.boxlist.size() * 6);
		for (int i = 0; i < patchBvhs.boxlist.size(); i++)
		{
			hostBvhsNodesBbox[6 * i] = patchBvhs.boxlist[i][0][0];
			hostBvhsNodesBbox[6 * i + 1] = patchBvhs.boxlist[i][0][1];
			hostBvhsNodesBbox[6 * i + 2] = patchBvhs.boxlist[i][0][2];
			hostBvhsNodesBbox[6 * i + 3] = patchBvhs.boxlist[i][1][0];
			hostBvhsNodesBbox[6 * i + 4] = patchBvhs.boxlist[i][1][1];
			hostBvhsNodesBbox[6 * i + 5] = patchBvhs.boxlist[i][1][2];
		}

		hostPatchIdToNodeId.resize(patchBvhs.new2old.size());
		for (int i = 0; i < patchBvhs.new2old.size(); i++)
			hostPatchIdToNodeId[i] = patchBvhs.oldboxId[i];

		hostNodeIdToPatchId = patchBvhs.nodeIdToElementId;
		cudaMalloc((void**)&devBvhsNodesNum, sizeof(int));
		cudaMemcpy(devBvhsNodesNum, &hostBvhsNodesNum, sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devBvhsNodesBbox, hostBvhsNodesBbox.size() * sizeof(qeal));
		cudaMemcpy(devBvhsNodesBbox, hostBvhsNodesBbox.data(), hostBvhsNodesBbox.size() * sizeof(qeal), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devNodeIdToPatchId, hostNodeIdToPatchId.size() * sizeof(int));
		cudaMemcpy(devNodeIdToPatchId, hostNodeIdToPatchId.data(), hostNodeIdToPatchId.size() * sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc((void**)&devPatchIdToNodeId, hostPatchIdToNodeId.size() * sizeof(int));
		cudaMemcpy(devPatchIdToNodeId, hostPatchIdToNodeId.data(), hostPatchIdToNodeId.size() * sizeof(int), cudaMemcpyHostToDevice);
		//------------------------------------------------------------------------------------------

		std::vector<std::vector<int>> patchICPairsList;
		patchICPairsList.resize(hostPatchNum);
		for (int i = 0; i < hostPatchNum; i++)
			patchICPairsList[i].resize(EACH_MAX_IC_PAIRS, 0);

		std::vector<int> hostPatchICPairsListNum;
		flatten2DArray(patchICPairsList, hostPatchICPairsList, hostPatchICPairsListNum, hostPatchICPairsListOffset);
		CUDA_CALL(cudaMalloc((void**)&devPatchICPairsList, sizeof(int) * hostPatchICPairsList.size()));
		CUDA_CALL(cudaMemcpy(devPatchICPairsList, hostPatchICPairsList.data(), sizeof(int) * hostPatchICPairsList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchICPairsListOffset, sizeof(int) * hostPatchICPairsListOffset.size()));
		CUDA_CALL(cudaMemcpy(devPatchICPairsListOffset, hostPatchICPairsListOffset.data(), sizeof(int) * hostPatchICPairsListOffset.size(), cudaMemcpyHostToDevice));

		hostPatchICPairsNum = 0;
		CUDA_CALL(cudaMalloc((void**)&devPatchICPairsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devPatchICPairsNum, &hostPatchICPairsNum, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPatchICPairs, sizeof(int) * hostPatchICPairsList.size()));

		hostEachPatchICPairsNum.resize(hostPatchNum, 0);
		devEachPatchICPairsNum = hostEachPatchICPairsNum;
		hostEachPatchICPairsOffset.resize(hostPatchNum, 0);
		devEachPatchICPairsOffset = hostEachPatchICPairsOffset;

		//------------------------------------------------------------------------------------
		//////////////////////////////////////////////////////////////////
		std::vector<std::vector<int>> bvhsNodeLevelGroup(patchBvhs.nodeLevelGroup.size());
		for (int i = 0; i < bvhsNodeLevelGroup.size(); i++)
		{
			bvhsNodeLevelGroup[i].resize(patchBvhs.nodeLevelGroup[i].size());
			for (int j = 0; j < patchBvhs.nodeLevelGroup[i].size(); j++)
				bvhsNodeLevelGroup[i][j] = patchBvhs.nodeLevelGroup[i][j];
		}

		hostBvhsNodesGroups = bvhsNodeLevelGroup.size();
		CUDA_CALL(cudaMalloc((void**)&devBvhsNodesGroups, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devBvhsNodesGroups, &hostBvhsNodesGroups, sizeof(int), cudaMemcpyHostToDevice));

		flatten2DArray(bvhsNodeLevelGroup, hostBvhsNodesGroupList, hostBvhsNodesGroupNum, hostBvhsNodesGroupOffset);
		CUDA_CALL(cudaMalloc((void**)&devBvhsNodesGroupList, sizeof(int) * hostBvhsNodesGroupList.size()));
		CUDA_CALL(cudaMemcpy(devBvhsNodesGroupList, hostBvhsNodesGroupList.data(), sizeof(int) * hostBvhsNodesGroupList.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devBvhsNodesGroupNum, sizeof(int) * hostBvhsNodesGroupNum.size()));
		CUDA_CALL(cudaMemcpy(devBvhsNodesGroupNum, hostBvhsNodesGroupNum.data(), sizeof(int) * hostBvhsNodesGroupNum.size(), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devBvhsNodesGroupOffset, sizeof(int) * hostBvhsNodesGroupOffset.size()));
		CUDA_CALL(cudaMemcpy(devBvhsNodesGroupOffset, hostBvhsNodesGroupOffset.data(), sizeof(int) * hostBvhsNodesGroupOffset.size(), cudaMemcpyHostToDevice));

		int tolListNum = 0;
		tolListNum = hostPatchNum;
		CUDA_CALL(cudaMalloc((void**)&devToiList, MAX_BLOCK_TOI * sizeof(qeal)));

		//Ccd
		hostPotentialCcdICPairsNum = 0;
		CUDA_CALL(cudaMalloc((void**)&devPotentialCcdICPairsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devPotentialCcdICPairsNum, &hostPotentialCcdICPairsNum, sizeof(int), cudaMemcpyHostToDevice));

		std::vector<int> hostPotentialCcdSCPairsList;
		for (int i = 0; i < hostPatchNum; i++)
		{
			//std::vector<int>& faces = m_patchesFaces[i];
			//std::vector<int>& points = m_patchesPoints[i];
			//std::vector<int>& edges = m_patchesEdges[i];
			//continue; // remove
			//for (int j = 0; j < points.size(); j++)
			//{
			//	int vid = points[j];
			//	if (vid >= simPointsNum)
			//		continue;
			//	for (int k = 0; k < faces.size(); k++)
			//	{
			//		int fid = faces[k];
			//		int f0 = boundaryFaces.col(fid)[0];
			//		int f1 = boundaryFaces.col(fid)[1];
			//		int f2 = boundaryFaces.col(fid)[2];

			//		if (vid == f0 || vid == f1 || vid == f2)
			//			continue;
			//		hostPotentialCcdSCPairsList.push_back(0);
			//		hostPotentialCcdSCPairsList.push_back(vid);
			//		hostPotentialCcdSCPairsList.push_back(f0);
			//		hostPotentialCcdSCPairsList.push_back(f1);
			//		hostPotentialCcdSCPairsList.push_back(f2);
			//	}
			//}

			//for (int j = 0; j < edges.size(); j++)
			//{
			//	int eaId = edges[j];
			//	int ea0 = boundaryEdges.col(eaId)[0];
			//	int ea1 = boundaryEdges.col(eaId)[1];
			//	if (ea0 >= simPointsNum)
			//		continue;
			//	for (int k = j + 1; k < edges.size(); k++)
			//	{
			//		int ebId = edges[k];
			//		int eb0 = boundaryEdges.col(ebId)[0];
			//		int eb1 = boundaryEdges.col(ebId)[1];
			//		if (ea0 == eb0 || ea0 == eb1)
			//			continue;
			//		if (ea1 == eb0 || ea1 == eb1)
			//			continue;
			//		hostPotentialCcdSCPairsList.push_back(1);
			//		hostPotentialCcdSCPairsList.push_back(ea0);
			//		hostPotentialCcdSCPairsList.push_back(ea1);
			//		hostPotentialCcdSCPairsList.push_back(eb0);
			//		hostPotentialCcdSCPairsList.push_back(eb1);
			//	}
			//}
		}
		hostPotentialCcdSCPairsNum = hostPotentialCcdSCPairsList.size() / 5;
		CUDA_CALL(cudaMalloc((void**)&devPotentialCcdSCPairsNum, sizeof(int)));
		CUDA_CALL(cudaMemcpy(devPotentialCcdSCPairsNum, &hostPotentialCcdSCPairsNum, sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMalloc((void**)&devPotentialCcdPairsList, (5 * MAX_CCD_CANDIDATE_NUM + hostPotentialCcdSCPairsList.size()) * sizeof(int)));

		CUDA_CALL(cudaMemcpy(devPotentialCcdPairsList, hostPotentialCcdSCPairsList.data(), hostPotentialCcdSCPairsList.size() * sizeof(int), cudaMemcpyHostToDevice));

		devPotentialCcdICPairsList = devPotentialCcdPairsList + hostPotentialCcdSCPairsList.size();

		hostPatchICPairsVfeeNumList.resize(hostPatchICPairsList.size());
		devPatchICPairsVfeeNumList = hostPatchICPairsVfeeNumList;
		hostPatchICPairsVfeeOffsetList.resize(hostPatchICPairsList.size());
		devPatchICPairsVfeeOffsetList = hostPatchICPairsVfeeOffsetList;
		hostPatchICPairsVfeeNumList.clear();
		hostPatchICPairsVfeeOffsetList.clear();


		hostCollisionPointsMass.resize(totalPointsNum, 10000.0);
		std::copy(simPointsMass.data(), simPointsMass.data() + simPointsNum, hostCollisionPointsMass.data());
		cudaMalloc((void**)&devCollisionPointsMass, totalPointsNum * sizeof(qeal));
		cudaMemcpy(devCollisionPointsMass, hostCollisionPointsMass.data(), totalPointsNum * sizeof(qeal), cudaMemcpyHostToDevice);
		///
		updatePatchBboxBlockSize = hostPatchMaxThreads;
		updatePatchBboxGridSize = hostPatchNum;

		findPotentialCollisionsBlockSize = THREADS_NUM_128;
		findPotentialCollisionsGridSize = ((hostPatchNum)+(THREADS_NUM_128 - 1)) / THREADS_NUM_128;
	}

	void PdIpcSimulator::getSimPointsNeighborsList(std::vector<std::vector<int>>& simPointsNeighborsList)
	{
		std::vector<std::set<int>> simPointsNeighborsSet(simPointsNum);
		for (int i = 0; i < simElements.cols(); i++)
		{
			Vector4i ele = simElements.col(i);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
					if (j != k)
						simPointsNeighborsSet[ele[j]].insert(ele[k]);
			}
		}

		//clothPointsNeighborsList
		for (int i = 0; i < simClothPointsNum; i++)
		{
			for (int j = 1; j < clothPointsNeighborsList[i].size(); j++)
				simPointsNeighborsSet[clothPointsNeighborsList[i][0]].insert(clothPointsNeighborsList[i][j]);
		}

		simPointsNeighborsList.resize(simPointsNum);
		for (int i = 0; i < simPointsNum; i++)
		{
			std::vector<int> list; 
			list.push_back(i);
			for (auto it : simPointsNeighborsSet[i])
				list.push_back(it);
			simPointsNeighborsList[i] = list;
		}
	}

	void PdIpcSimulator::doTime(int frame)
	{
		//if (frame == 280)
		//{		
		//	cudaMemcpy(devGravityForce, devTempGravityForce, sizeof(qeal) * simDims, cudaMemcpyDeviceToDevice);
		//}
		
		subStep = 1;
		qeal dt = _timeStep / subStep;
		for (int i = 0; i < subStep; i++)
		{
			advance(frame, dt);
		}
	
	}

	int PdIpcSimulator::LocalGlobalIteration(const int frame, const qeal dt)
	{
		// local-global
		cudaMemcpy(devSysDiagonal, devElasticsDiagonal, sizeof(qeal) * simPointsNum, cudaMemcpyDeviceToDevice);
		cublasSaxpy(blasHandle, simPointsNum, &posOne, devInertialDiagonal, 1, devSysDiagonal, 1);
		if (activeCollisionList)
			cublasSaxpy(blasHandle, simPointsNum, &posOne, devCollisionDiagonal, 1, devSysDiagonal, 1);
			
		updateSuperJacobiOperatorsHost
		(
			R1OperatorsBlockSize,
			R1OperatorsGridSize,
			hostHessR1OperatorsNum,
			devSysDiagonal,
			devHessR1Value,
			devHessR1DiiIndex,
			devHessR1Operators,
			//
			R2OperatorsBlockSize,
			R2OperatorsGridSize,
			hostHessR2OperatorsNum,
			devHessR2Aijk,
			devHessR2DiiIndex,
			devHessR2DjjIndex,
			devHessR2Operators,
			//
			R2RowOperatorsBlockSize,
			R2RowOperatorsGridSize,
			hostHessR2RowOperatorsTotal,
			devHessR2RowOperatorsNum,
			devHessR2RowOperatorsOffset,
			devHessR2ConstOrderOperators,
			devHessR2FirstOrderOperators,
			devHessR2RowOperators
		);

		ConstraintEnergyInfo e0, e1;	
		int lgIter = 0;
		for (; lgIter < 10; lgIter++)
		{
			doLocalProjection(frame, dt, e0);
			chebyshevJacobiR2SolverHost
			(
				15,
				blasHandle,
				R2FirstTermBlockSize,
				R2FirstTermGridSize,
				simDims,
				simPointsNum,
				devSysDiagonal,
				devSysRhs,
				devHessR1Operators,
				devHessR1Index,
				devHessR1IndexNum,
				devHessR1IndexOffset,
				devFirstTerm,
				R2IteraionBlockSize,
				R2IteraionGridSize,
				devPrePoints,
				devPoints,
				devNextPoints,
				devHessR2RowOperators,
				devHessR2Index,
				devHessR2IndexNum,
				devHessR2IndexOffset
			);
		}
		return lgIter;
	}

	void PdIpcSimulator::doLocalProjection(int frame, const qeal dt, ConstraintEnergyInfo& energy)
	{	
		cudaMemcpy(devSysRhs, devInitRhs, sizeof(qeal) * simDims, cudaMemcpyDeviceToDevice); // done
		if (simElementsNum > 0)
		{
			solveVolumetricStrainConstraintHost
			(
				simElementsBlockSize, // done
				simElementsGridSize, // done
				simElementsNum, // done
				devSimElements, // done
				devSimElementsRestShape, //done
				devVolumetricStrainStiffness, //done
				devPoints, // done
				devSysRhs, //done
				devVolumetricStrainEnergy//done
			);
		}

		if (simClothPointsNum > 0)
		{
			solveClothPointsStrainAndBendingConstraintHost
			(
				simClothPointsBlockSize,
				simClothPointsGridSize,
				simClothPointsNum,
				devClothPointsStrainStiffness,
				devClothPointsBendingStiffness,
				devClothPointsNeighborsList,
				devClothPointsNeighborsNum,
				devClothPointsNeighborsOffset,
				devClothPointsNeighborsRestLength,
				devClothPointsNeighborsLbo,
				devClothPointsMeanCurvatureNorm,
				devPoints,
				devSysRhs
			);
		}
		

		if (activeCollisionList)
			cublasSaxpy(blasHandle, simDims, &posOne, devCollisionRhs, 1, devSysRhs, 1);

		if (positionConstraintNum > 0)
		{
			solvePointsConstraintHost
			(
				simPointsConstraintBlockSize,
				simPointsConstraintGridSize,
				positionConstraintNum,
				devRestPoints,
				devSimPointsConstraintIndex,
				devSimPointsConstraintStiffness,
				devSysRhs
			);
		}

	}

	void PdIpcSimulator::advance(const int frame, const qeal dt)
	{
		const int maxOutterIters = 30;
		int outterIter = 0;
		qeal toi = 1.0;
		activeCollisionList = false;

		cudaMemcpy(devX0, devPoints, sizeof(qeal) * simDims, cudaMemcpyDeviceToDevice);
		cudaMemcpy(devOldPoints, devPoints, sizeof(qeal) * simDims, cudaMemcpyDeviceToDevice);

		predictiveHost
		(
			simPointsBlockSize,
			simPointsGridSize,
			simPointsNum,
			dt,
			devPoints,
			devVelocity,
			devInitRhs,
			devInertialDiagonal,
			devGravityForce,
			devMouseForce,
			devPenaltyForce,
			devInvMassDiag,
			devSimPointsConstraintFlag
		);


		for (; outterIter < maxOutterIters; outterIter++)
		{
			LocalGlobalIteration(frame, dt);

			computeDirectionHost //
			(
				simDimsBlockSize,
				simDimsGridSize,
				simDims,
				devOldPoints,
				devPoints,
				devDirection
			);

			if (outterIter == 0)
			{
				gpuCcdCulling(devOldPoints, devDirection);
			}

			toi = getImpactOfTimeHost
			(
				blasHandle,
				hostPotentialCcdICPairsNum,
				devPotentialCcdPairsList,
				devOldPoints,
				devPoints,
				devToiList
			);

			clampDirectionHost
			(
				simDimsBlockSize,
				simDimsGridSize,
				simDims,
				devOldPoints,
				devDirection,
				devPoints,
				toi
			);

			cudaMemcpy(devCollisionDiagonal, devZeroVector, simPointsNum * sizeof(qeal), cudaMemcpyDeviceToDevice);
			CUDA_CALL(cudaMemcpy(devCollisionRhs, devZeroVector, simDims * sizeof(qeal), cudaMemcpyDeviceToDevice));
				
			// test projective at the state of free collision
			// remove tangent
			setCollisionConstraintsHost
			(
				dHat,
				dHat2,
				kappa,
				dt,
				hostPotentialCcdICPairsNum,
				devPotentialCcdPairsList,
				devPoints,
				devCollisionPointsMass,
				devCollisionDiagonal,
				devCollisionRhs,
				totalPointsNum
			);

			activeCollisionList = true;

			qeal residual;
			cublasSdot(blasHandle, simDims, devDirection, 1, devDirection, 1, &residual);

			if ((residual < 1e-5) && outterIter > 0)
			{
				break;
			}
			cudaMemcpy(devOldPoints, devPoints, sizeof(qeal) * simDims, cudaMemcpyDeviceToDevice);
		}


		updateVelocityHost
		(
			simDimsBlockSize,
			simDimsGridSize,
			simDims,
			1. - dt,
			dt,
			devPoints,
			devX0,
			devVelocity
		);

	}


	void PdIpcSimulator::gpuCcdCulling(qeal* devX, qeal* devDir)
	{
		hostPotentialCcdICPairsNum = 0;
		updateCcdPatchBboxesHost
		(
			updatePatchBboxBlockSize,
			updatePatchBboxGridSize,
			&hostPatchMaxThreads, // must be 2^n and <= 1024
			devPatchFacesList,
			devPatchFacesNum,
			devPatchFacesOffset,
			devBoundaryFaces,
			devX,
			devDir,
			inflation_radius,
			devFacesBboxes,
			devPatchBboxes
		);

		updateBvhsBboxesHost
		(
			hostBvhsNodesGroups,
			devBvhsNodesGroupList,
			devBvhsNodesGroupNum,
			devBvhsNodesGroupOffset,
			devNodeIdToPatchId,
			devPatchIdToNodeId,
			devPatchBboxes,
			devBvhsNodesBbox
		);

		findPotentialCollisionsHost
		(
			findPotentialCollisionsBlockSize,
			findPotentialCollisionsGridSize,
			hostPatchNum,
			hostBvhsNodesNum,
			devBvhsNodesBbox,
			devNodeIdToPatchId,
			devPatchIdToNodeId,
			devPatchICPairsList,
			devPatchICPairsListOffset,
			&devEachPatchICPairsNum,
			&devEachPatchICPairsOffset,
			&hostPatchICPairsNum,
			devPatchICPairsNum,
			devPatchICPairs
		);

		getCcdPotentialICPairsListHost
		(
			hostPatchICPairsNum,
			devPatchICPairs,
			devPatchPointsList,
			devPatchPointsNum,
			devPatchPointsOffset,
			devPatchFacesList,
			devPatchFacesNum,
			devPatchFacesOffset,
			devPatchEdgesList,
			devPatchEdgesNum,
			devPatchEdgesOffset,
			devBoundaryFaces,
			devBoundaryEdges,
			&devPatchICPairsVfeeNumList,
			&devPatchICPairsVfeeOffsetList,
			&hostPotentialCcdICPairsNum,
			devPotentialCcdICPairsNum,
			devPotentialCcdICPairsList
		);
	}



}