#pragma once
#ifndef PD_IPC_SIMULATOR
#define PD_IPC_SIMULATOR
#include "Simulator\BaseSimulator.h"
#include <direct.h>
#include <oneapi/tbb.h>
#include "PdIpcVolumetricModel.h"
#include "PdIpcClothModel.h"
#include "Commom\FileIO.h"
#include "Commom\GeometryComputation.h"
#include "PdIpcGpuFunc.cuh"
#include "TetConstraint.h"
#include "Bow/Utils/Timer.h"
#include "Simulator/Bvhs/AipcBvhs.h"


namespace PD_IPC
{
	using namespace Bow;
	using namespace AIPC;
	using T = qeal;
	using _StorageIndex = int;
	using StorageIndex = int;

	class PdIpcSimulator : public BaseSimulator
	{
	public:
		std::function<void(void)> timestep_callback = []() {};

		struct ConstraintEnergyInfo
		{
			qeal inertialEnergy = 0;
			qeal volumetricStrainEnergy = 0;
			qeal clothStrainEnergy = 0;
			qeal clothBendingEnergy = 0;
			qeal positionEnergy = 0;
			qeal collisionEnergy = 0;

			qeal energy()
			{
				return inertialEnergy + volumetricStrainEnergy + clothStrainEnergy + clothBendingEnergy + positionEnergy + collisionEnergy;
			}
			void reset()
			{
				inertialEnergy = 0, volumetricStrainEnergy = 0, clothStrainEnergy = 0, clothBendingEnergy = 0, positionEnergy = 0, collisionEnergy = 0;
			}
		};

		qeal damping, dHat, dHat2, kappa, tol;
		int subStep = 1;
		int lapDamping = 4;
		qeal penaltyStiffness = 500, penaltySelfStiffness = 500, penaltyDamping = 0.1;

		qeal fps = 0;
		int totalJacobiNum = 0;
		int totalLGNum = 0;
		int totalOutterNum = 0;
		int maxCollisionNum = 0;

		qeal dcdTiming = 0;
		qeal ccdTiming = 0;
		qeal collisionTiming = 0;
		qeal pdSolverTiming = 0;
		qeal constraintsTiming = 0;

		int totalPointsNum;
		int totalElementsNum;
		int totalFacesNum;
		int totalDims;

		int simPointsNum;
		int simElementsNum;
		int simClothPointsNum;
		int simDims;
		

		int totalBoundaryPointsNum;
		int totalBoundaryEdgesNum;
		int totalBoundaryFacesNum;

		int simBoundaryPointsNum;
		int simBoundaryEdgesNum;
		int simBoundaryFacesNum;

		int staticBoundaryPointsNum;
		int staticBoundaryEdgesNum;
		int staticBoundaryFacesNum;

		MatrixX points;
		MatrixX restPoints;
		
		MatrixXi simElements;

		std::vector<qeal> clothPointsBendingArea;
		std::vector<std::vector<qeal>> clothPointsStrainStiffness;
		std::vector<std::vector<qeal>> clothPointsBendingStiffness;
		std::vector<std::vector<int>> clothPointsNeighborsList;
		std::vector<std::vector<qeal>> clothPointsNeighborsLbo;
		std::vector<std::vector<qeal>> clothPointsNeighborsRestLength;
		std::vector<qeal> clothPointsMeanCurvatureNorm;



		MatrixXi boundaryEdges;
		MatrixXi boundaryFaces;
		VectorXi boundaryPoints;


		VectorX simPointsMass;
		VectorX simPointsInvMass;

		MatrixX velocity;
		MatrixX Sn;
		MatrixX externalForce;
		MatrixX gravityForce;
		MatrixX mouseForce;

		SparseMatrix Mass;
		SparseMatrix volumetricStrainHess;
		SparseMatrix clothHess;
		SparseMatrix elasticsHess;
		SparseMatrix collisionHess;
		SparseMatrix sysHess;

		qeal posOne = 1.0, negOne = -1.0;
		std::vector<std::vector<Vector4i>> volumetricStrainConsList;
		std::vector<qeal> volumetricStrainStiffness;
		std::vector<std::vector<Vector4i>> volumetricVolumeConsList;
		std::vector<qeal> volumetricVolumeStiffness;
		std::vector<std::vector<int>> positionConsList;
		std::vector<qeal> positionStiffness;

		dim3 totalDimsBlockSize;
		dim3 totalDimsGridSize;

		dim3 totalPointsBlockSize;
		dim3 totalPointsGridSize;

		dim3 simDimsBlockSize;
		dim3 simDimsGridSize;

		dim3 simPointsBlockSize;
		dim3 simPointsGridSize;

		dim3 simElementsBlockSize;
		dim3 simElementsGridSize;

		dim3 simClothPointsBlockSize;
		dim3 simClothPointsGridSize;

		dim3 simPointsConstraintBlockSize;
		dim3 simPointsConstraintGridSize;

		qeal* devZeroVector;
		qeal* devZeroVectorInt;
		int* devBoundaryFaces, *devBoundaryEdges;

		qeal* devPoints, *devRestPoints, *devX0, *devOldPoints, *devDirection;
		qeal* devPrePoints, *devNextPoints;
		qeal* devVelocity;
		qeal* devNextVelocity;

		std::vector<int> hostPointsModelsId;
		std::vector<int> hostSimPointsNeighborsList;
		std::vector<int> hostSimPointsNeighborsNum;
		std::vector<int> hostSimPointsNeighborsOffset;
		int* devPointsModelsId;
		int* devSimPointsNeighborsList;
		int* devSimPointsNeighborsNum;
		int* devSimPointsNeighborsOffset;

		// constraints
		qeal* devInertialEnergy, *devPositionEnergy, *devVolumetricStrainEnergy;

		int positionConstraintNum;
		std::vector<int> hostSimPointsConstraintFlag;
		std::vector<int> hostSimPointsConstraintIndex;
		std::vector<qeal> hostSimPointsConstraintStiffness;
		int* devSimPointsConstraintFlag;
		int* devSimPointsConstraintIndex;
		qeal* devSimPointsConstraintStiffness;

		std::vector<int> hostSimElements;
		std::vector<qeal> hostVolumetricStrainStiffness;
		std::vector<qeal> hostSimElementsRestShape;
		int* devSimElements;
		qeal* devSimElementsRestShape, *devVolumetricStrainStiffness;
		//
		qeal* devClothPointsStrainStiffness;
		qeal* devClothPointsBendingStiffness;
		int* devClothPointsNeighborsList;
		int* devClothPointsNeighborsNum;
		int* devClothPointsNeighborsOffset;
		qeal* devClothPointsNeighborsLbo;
		qeal* devClothPointsNeighborsRestLength;
		qeal* devClothPointsMeanCurvatureNorm;

		//
		qeal* devInitRhs, *devSn, *devInertialDiagonal;
		qeal* devGravityForce, *devMouseForce, *devPenaltyForce, *devMassDiag, *devInvMassDiag;
		qeal* devSysRhs, *devCollisionRhs;
		int* devPointsPenaltyTimes;

		qeal* devTempGravityForce;

		qeal *devSysDiagonal, *devElasticsDiagonal, *devCollisionDiagonal;
		qeal* devFirstTerm, *devTempTerm, *devDinvRhs, *devXk, *devPreXk, *devOmega;

		//solver
		dim3 R1OperatorsBlockSize;
		dim3 R1OperatorsGridSize;
		dim3 R1IteraionBlockSize;
		dim3 R1IteraionGridSize;

		dim3 R2OperatorsBlockSize;
		dim3 R2OperatorsGridSize;
		dim3 R2RowOperatorsBlockSize;
		dim3 R2RowOperatorsGridSize;
		dim3 R2FirstTermBlockSize;
		dim3 R2FirstTermGridSize;
		dim3 R2IteraionBlockSize;
		dim3 R2IteraionGridSize;

		dim3 R3OperatorsBlockSize;
		dim3 R3OperatorsGridSize;
		dim3 R3RowOperatorsBlockSize;
		dim3 R3RowOperatorsGridSize;
		dim3 R3FirstTermBlockSize;
		dim3 R3FirstTermGridSize;
		dim3 R3IteraionBlockSize;
		dim3 R3IteraionGridSize;



		//R3
		int hostHessR3OperatorsNum;
		int* devHessR3Index, * devHessR3IndexNum, * devHessR3IndexOffset;
		int* devHessR3DiiIndex, * devHessR3DjjIndex, * devHessR3AijkNum, * devHessR3AijkOffset, * devHessR3OperatorsNum;
		qeal* devHessR3Aijk, * devHessR3Operators;




		qeal* devHessR3ConstOrderOperators;
		int* devHessR3FirstOrderOperators;

		// R3 row
		int hostHessR3RowOperatorsTotal;
		int* devHessR3RowOperatorsNum, * devHessR3RowOperatorsOffset, * devHessR3RowOperatorsTotal;
		qeal* devHessR3RowOperators;

		//R21
		qeal* devHessR21RowOperators;


		int hostR3MaxOperatorIndexNum;
		int* devR3MaxOperatorIndexNum;

		int hostR3AlignedPointsSharedMemSize;
		int* devR3AlignedPointsSharedMemSize;

		int hostR3AlignedDimsSharedMemSize;
		int* devR3AlignedDimsSharedMemSize;
		//R1
		int hostHessR1OperatorsNum;
		int* devHessR1Index, *devHessR1DiiIndex;
		qeal* devHessR1Value, *devHessR1Operators;
		int* devHessR1IndexNum, *devHessR1IndexOffset, *devHessR1OperatorsNum;
		
		int hostR1MaxOperatorIndexNum;
		int* devR1MaxOperatorIndexNum;

		//R2
		int hostHessR2OperatorsNum;
		int* devHessR2Index, *devHessR2IndexNum, *devHessR2IndexOffset;
		int* devHessR2DiiIndex, *devHessR2DjjIndex, *devHessR2AijkNum, *devHessR2AijkOffset, *devHessR2OperatorsNum;
		qeal *devHessR2Aijk, *devHessR2Operators;
	
		qeal* devHessR2ConstOrderOperators;
		int* devHessR2FirstOrderOperators;

		// R2 row
		int hostHessR2RowOperatorsTotal;
		int* devHessR2RowOperatorsNum, *devHessR2RowOperatorsOffset, *devHessR2RowOperatorsTotal;
		qeal* devHessR2RowOperators;

		int hostR2MaxOperatorIndexNum;
		int* devR2MaxOperatorIndexNum;

		int hostR1AlignedDimsSharedMemSize;
		int* devR1AlignedDimsSharedMemSize;

		int hostR2AlignedPointsSharedMemSize;
		int* devR2AlignedPointsSharedMemSize;

		int hostR2AlignedDimsSharedMemSize;
		int* devR2AlignedDimsSharedMemSize;

		
		// collision detection
		bool enableAllEC = true; // enable all edge edge constraints
		qeal* devInflationRadius;
		std::vector<std::vector<int>> m_patchesFaces;
		std::vector<std::vector<int>> m_patchesEdges;
		std::vector<std::vector<int>> m_patchesPoints;
		std::vector<std::vector<int>> m_patchesExPoints;
		std::vector<QColor> m_patchesFacesColor;

		int hostPatchNum, *devPatchNum;
		std::vector<int> hostFacesToPatchIndex;
		std::vector<int> hostPatchFacesList;
		std::vector<int> hostPatchFacesNum;
		std::vector<int> hostPatchFacesOffset;
		int* devFacesToPatchIndex;
		int* devPatchFacesList;
		int* devPatchFacesNum;
		int* devPatchFacesOffset;

		std::vector<int> hostPatchEdgesList;
		std::vector<int> hostPatchEdgesNum;
		std::vector<int> hostPatchEdgesOffset;
		int* devPatchEdgesList;
		int* devPatchEdgesNum;
		int* devPatchEdgesOffset;

		std::vector<int> hostPatchPointsList;
		std::vector<int> hostPatchPointsNum;
		std::vector<int> hostPatchPointsOffset;
		int* devPatchPointsList;
		int* devPatchPointsNum;
		int* devPatchPointsOffset;

		std::vector<int> hostPatchExPointsList;
		std::vector<int> hostPatchExPointsNum;
		std::vector<int> hostPatchExPointsOffset;
		int* devPatchExPointsList;
		int* devPatchExPointsNum;
		int* devPatchExPointsOffset;

		std::vector<std::array<Vector3, 2>> m_patchBbox;
		int hostPatchMaxThreads;
		std::vector<qeal> hostPatchBboxes;
		qeal* devPatchBboxes;
		std::vector<qeal> hostFacesBboxes;
		qeal* devFacesBboxes;

		int hostBvhsNodesNum;
		int* devBvhsNodesNum;
		std::vector<qeal> hostBvhsNodesBbox;
		qeal* devBvhsNodesBbox;
		std::vector<int> hostNodeIdToPatchId;
		int* devNodeIdToPatchId;
		std::vector<int> hostPatchIdToNodeId;
		int* devPatchIdToNodeId;
		//------------------------------------------
		//---------------------------------------------------------
		std::vector<int> hostPatchICPairsList;
		std::vector<int> hostPatchICPairsListOffset;
		int* devPatchICPairsList;
		int* devPatchICPairsListOffset;
		//
		int hostPatchICPairsNum;
		int* devPatchICPairsNum;
		int* devPatchICPairs;
		thrust::host_vector<int> hostEachPatchICPairsNum;
		thrust::device_vector<int> devEachPatchICPairsNum;
		thrust::host_vector<int> hostEachPatchICPairsOffset;
		thrust::device_vector<int> devEachPatchICPairsOffset;

		thrust::host_vector<int> hostPatchICPairsVfeeNumList;
		thrust::device_vector<int> devPatchICPairsVfeeNumList;
		thrust::host_vector<int> hostPatchICPairsVfeeOffsetList;
		thrust::device_vector<int> devPatchICPairsVfeeOffsetList;

		int hostPotentialCcdICPairsNum;
		int* devPotentialCcdICPairsNum;
		int* devPotentialCcdICPairsList;

		int hostPotentialCcdSCPairsNum;
		int* devPotentialCcdSCPairsNum;
		int* devPotentialCcdPairsList;

		thrust::host_vector<int> hostPatchICPairsF2FNumList;
		thrust::device_vector<int> devPatchICPairsF2FNumList;
		thrust::host_vector<int> hostPatchICPairsF2FOffsetList;
		thrust::device_vector<int> devPatchICPairsF2FOffsetList;

		int hostPotentialDcdICPairsNum;
		int* devPotentialDcdICPairsNum;
		int* devPotentialDcdICPairsList;

		int hostPotentialDcdSCPairsNum;
		int* devPotentialDcdSCPairsNum;
		int* devPotentialDcdPairsList;


		int* devCollisionConstraintsFlag;
		qeal* devCollisionConstraintsEnergy;
		qeal* devCollisionConstraintsTargetPos;
		bool activeCollisionList;

		qeal* devBarrier;

		///////////////////////////////////////
		int hostBvhsNodesGroups;
		int* devBvhsNodesGroups;
		std::vector<int> hostBvhsNodesGroupList;
		std::vector<int> hostBvhsNodesGroupNum;
		std::vector<int> hostBvhsNodesGroupOffset;
		int* devBvhsNodesGroupList;
		int* devBvhsNodesGroupNum;
		int* devBvhsNodesGroupOffset;
		//

		int hostPatchPairsMaxElementsPairNum;
		int* devPatchPairsMaxElementsPairNum;
		int hostPatchPairsMaxElementsNum;
		int* devPatchPairsMaxElementsNum;

		qeal* devToiList;

		int hostCandidateNum;
		int *devCandidateNum;
		int* devCandidateList;

		std::vector<qeal> hostCollisionPointsMass;
		qeal* devCollisionPointsMass;


		dim3 updatePatchBboxBlockSize;
		dim3 updatePatchBboxGridSize;

		dim3 findPotentialCollisionsBlockSize;
		dim3 findPotentialCollisionsGridSize;

		std::vector<Vector3> mousePointsList;
		std::vector<int> mousePointsIds;
		Vector3 mouseForceLine;

		//
		PdIpcSimulator(std::string simName = "pd_ipc_simulator", RunPlatform runPlatform = RunPlatform::CUDA) :
			BaseSimulator(simName, runPlatform),
			damping(0.999),
			dHat(1e-3),
			kappa(1e4),
			tol(1e-2)
		{
			dHat2 = dHat * dHat;
			_simulatorName = std::string("pd_ipc_simulator");
		}

		PdIpcModel* getModel(const int id)
		{
			return dynamic_cast<PdIpcModel*> (models[id]);
		}

		PdIpcModel* getStaticModel(const int id)
		{
			return dynamic_cast<PdIpcModel*> (staticModels[id]);
		}

		virtual bool addModelFromConfigFile(const std::string filename, TiXmlElement* item)
		{
			TiXmlElement* childItem = item->FirstChildElement();
			qeal sx = 1.0, sy = 1.0, sz = 1.0;
			qeal tx = 0.0, ty = 0.0, tz = 0.0;
			qeal rx = 0.0, ry = 0.0, rz = 0.0, rsita = 0.0;
			qeal density = 10.0;
			qeal stiffness0, stiffness1, stiffness2;
			int patchDim = 10; 

			int modelType = 1;
			std::strstream ss;
			while (childItem)
			{
				ss.clear();
				std::string itemName = childItem->Value();
				if (itemName == std::string("type"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> modelType >> stiffness0 >> stiffness1 >> stiffness2;
				}
				else if (itemName == std::string("patchDim"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> patchDim;
				}
				else if (itemName == std::string("scale"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> sx >> sy >> sz;
				}
				else if (itemName == std::string("translation"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> tx >> ty >> tz;
				}
				else if (itemName == std::string("rotation"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> rx >> ry >> rz >> rsita;
				}
				else if (itemName == std::string("density"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> density;
					//m->density = density;
				}
				childItem = childItem->NextSiblingElement();
			}

			PdIpcModel *m;
			if (modelType == 1)
				m = new PdIpcVolumetricModel(stiffness0, stiffness1, stiffness2);
			else m = new PdIpcClothModel(stiffness0, stiffness1, stiffness2);
			bool isReadMesh = BaseSimulator::addModelFromConfigFile(filename, item, m);
			if (!isReadMesh)
				return false;

			m->scaleModel(sx, sy, sz);
			m->translateModel(tx, ty, tz);
			m->rotateModel(rx, ry, rz, rsita);
			m->density = density;
			m->setCollisionPatchDim(patchDim);
			m->initMeshesHandel();

			return true;
		}

		virtual bool addStaticModelFromConfigFile(const std::string filename, TiXmlElement* item)
		{
			TiXmlElement* childItem = item->FirstChildElement();
			qeal sx = 1.0, sy = 1.0, sz = 1.0;
			qeal tx = 0.0, ty = 0.0, tz = 0.0;
			qeal rx = 0.0, ry = 0.0, rz = 0.0, rsita = 0.0;
			int modelType = 1;
			qeal stiffness0, stiffness1, stiffness2;
			int patchDim = 10;
			std::strstream ss;
			while (childItem)
			{
				ss.clear();
				std::string itemName = childItem->Value();
				if (itemName == std::string("type"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> modelType >> stiffness0 >> stiffness1 >> stiffness2;
				}else if (itemName == std::string("patchDim"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> patchDim;
				}
				else if (itemName == std::string("scale"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> sx >> sy >> sz;
				}
				else if (itemName == std::string("translation"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> tx >> ty >> tz;
				}
				else if (itemName == std::string("rotation"))
				{
					std::string str = childItem->GetText();
					ss << str;
					ss >> rx >> ry >> rz >> rsita;
				}
				childItem = childItem->NextSiblingElement();
			}

			PdIpcModel *m;
			if (modelType == 1)
				m = new PdIpcVolumetricModel();
			else m = new PdIpcClothModel();
			bool isReadMesh = BaseSimulator::addStaticModelFromConfigFile(filename, item, m);
			if (!isReadMesh)
				return false;

			m->scaleModel(sx, sy, sz);
			m->translateModel(tx, ty, tz);
			m->rotateModel(rx, ry, rz, rsita);
			m->setCollisionPatchDim(patchDim);
			m->initMeshesHandel();
			m->density = 100000;

			QColor color("#00a2e9");
			m->setSurfaceNormalColor(color);
			return true;
		}

		virtual void readExtraAttributeFromConfigFile(TiXmlElement* item)
		{
			if (!item) return;
			BaseSimulator::readExtraAttributeFromConfigFile(item);
			std::string subItemName = item->Value();

			if (subItemName == std::string("tol"))
			{
				QString strVal(item->GetText());
				tol = strVal.toDouble();
			}
			if (subItemName == std::string("damping"))
			{
				QString strVal(item->GetText());
				damping = strVal.toDouble();
			}
			else if (subItemName == std::string("dhat"))
			{
				QString strVal(item->GetText());
				dHat = strVal.toDouble();
				dHat2 = dHat * dHat;
			}
			else if (subItemName == std::string("kappa"))
			{
				QString strVal(item->GetText());
				kappa = strVal.toDouble();
			}
			else if (subItemName == std::string("substep"))
			{
				QString strVal(item->GetText());
				subStep = strVal.toInt();
			}
			else if (subItemName == std::string("penalty"))
			{
				std::string str = item->GetText();
				std::strstream ss;
				ss << str;
				ss >> penaltyStiffness >> penaltySelfStiffness >> penaltyDamping;
			}

		}

		virtual void handleMouseForce(int nid, qeal& x, qeal& y, qeal& z)
		{
			qeal scale = 40000;
			mouseForce.setZero();
			if (nid > 0 && nid < simPointsNum)
			{
				mouseForce.col(nid) = Vector3(x, y, z) * scale;
				//std::cout << "mouseForce: " << nid <<" " << Vector3(x, y, z).transpose() << std::endl;
			}
			cudaMemcpy(devMouseForce, mouseForce.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);
		}

		virtual void handleMouseForce(std::vector<int>&list, qeal& x, qeal& y, qeal& z)
		{
			qeal scale = 1000;
			mouseForce.setZero();
			mousePointsList.clear();
			mousePointsList.resize(list.size());
			mousePointsIds = list;
			for (int i = 0; i < list.size(); i++)
			{
				mouseForce.col(list[i]) = Vector3(x, y, z) * scale;
				mousePointsList[i] = Vector3(x, y, z);
			}
			mouseForceLine = Vector3(x, y, z) * scale;

			//std::cout << list.size() << std::endl;

			//if (nid > 0 && nid < simPointsNum)
			//{
			//	mouseForce.col(nid) = Vector3(x, y, z) * scale;
			//	//std::cout << "mouseForce: " << nid <<" " << Vector3(x, y, z).transpose() << std::endl;
			//}

			std::cout << "mouseForce: " << list.size() << " " << Vector3(x, y, z).transpose() << std::endl;
			cudaMemcpy(devMouseForce, mouseForce.data(), sizeof(qeal) * simDims, cudaMemcpyHostToDevice);
		}
		
		/*
		virtual void saveSimulator(int frame)
		{
			std::string filename = _sceneDir + "saveDat/";
			int code = mkdir(filename.c_str());
			if (frame < 10)
				filename += "frame_000" + QString().setNum(frame).toStdString() + ".dat";
			else if (frame < 100)
				filename += "frame_00" + QString().setNum(frame).toStdString() + ".dat";
			else if (frame < 1000)
				filename += "frame_0" + QString().setNum(frame).toStdString() + ".dat";
			else
				filename += "frame_" + QString().setNum(frame).toStdString() + ".dat";

			std::ofstream fout;
			fout.open(filename, std::ios::binary);
			EigenMatrixIO::write_binary(fout, points);
			fout.close();


		//	std::vector<Vector3> mousePointsList;
		//	std::vector<int> mousePointsIds;
		//	Vector3 mouseForceLine;
			filename = _sceneDir + "saveMouseForce/";
			code = mkdir(filename.c_str());
			if (frame < 10)
				filename += "mouse_000" + QString().setNum(frame).toStdString() + ".txt";
			else if (frame < 100)
				filename += "mouse_00" + QString().setNum(frame).toStdString() + ".txt";
			else if (frame < 1000)
				filename += "mouse_0" + QString().setNum(frame).toStdString() + ".txt";
			else
				filename += "mouse_" + QString().setNum(frame).toStdString() + ".txt";
			fout.open(filename);
			fout << mousePointsIds.size() << std::endl;
			fout << mouseForceLine[0] << " " << mouseForceLine[1] <<" " << mouseForceLine[2] << std::endl;
			for (int i = 0; i < mousePointsIds.size(); i++)
				fout << mousePointsIds[i] << " " << mousePointsList[i][0] << " " << mousePointsList[i][1] << " " << mousePointsList[i][2] << std::endl;
			fout.close();

			return;
			if (frame % 3 != 0)
				return;
			frame /= 3;


			filename = _sceneDir + "saveObj/";
			code = mkdir(filename.c_str());
			if (frame < 10)
				filename += "frame_000" + QString().setNum(frame).toStdString() + ".obj";
			else if (frame < 100)
				filename += "frame_00" + QString().setNum(frame).toStdString() + ".obj";
			else if (frame < 1000)
				filename += "frame_0" + QString().setNum(frame).toStdString() + ".obj";
			else
				filename += "frame_" + QString().setNum(frame).toStdString() + ".obj";

	
			fout.open(filename);
			for (int i = 0; i < models[0]->pointsNum; i++)
			{
				Vector3 p = models[0]->getSurfacePoint(i);
				fout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			}
			for (int i = 0; i < models[0]->facesNum; i++)
			{
				Vector3i p = models[0]->getSurfaceFace(i);
				fout << "f " << p[0] + 1<< " " << p[1] + 1 << " " << p[2] + 1 << std::endl;
			}

			fout.close();




		}
		*/

		void getSphereObj(Vector3 cent, qeal r, std::vector<Vector3>& spherePoints, std::vector<Vector3i>& sphereFaces)
		{
			std::string sphereFile = "C:\\Users\\lan6\\Desktop\\video\\sphere.obj";
			std::ifstream fin;
			fin.open(sphereFile.c_str());
			while (!fin.eof())
			{
				std::string line;
				getline(fin, line);
				if (line.c_str()[0] == 'v' && line.c_str()[1] == ' ')
				{
					std::stringstream ss;
					ss << line;
					char ch;
					qeal x, y, z;
					ss >> ch >> x >> y >> z;
					spherePoints.push_back(Vector3(x, y, z));
				}
				else if (line.c_str()[0] == 'f')
				{
					std::stringstream ss;
					ss << line;
					char ch;
					int x, y, z;
					ss >> ch >> x >> y >> z;
					sphereFaces.push_back(Vector3i(x - 1, y - 1, z - 1));
				}
			}
			fin.close();

			Vector3 oriCent = Vector3::Zero();
			qeal maxLen = -1;
			for (int i = 0; i < spherePoints.size(); i++)
				oriCent += spherePoints[i];
			oriCent /= spherePoints.size();

			for (int i = 0; i < spherePoints.size(); i++)
			{
				qeal len = (spherePoints[i] - oriCent).norm();
				if (len > maxLen)
					maxLen = len;
			}

			for (int i = 0; i < spherePoints.size(); i++)
			{
				spherePoints[i] -= oriCent;
				spherePoints[i] /= maxLen;
				spherePoints[i] *= r;
				spherePoints[i] += cent;
			}
		}


		void getLineObj(Vector3 spos, Vector3 epos, qeal r, std::vector<Vector3>& pList, std::vector<Vector3i>& fList)
		{
			const int cone_stacks = 30;
			const int cone_slices = 30;

			std::vector<Vector3> cone_pos;
			std::vector<int> indices;

			qeal height = (spos - epos).norm();

			for (int i = 0; i <= cone_stacks; ++i)
			{
				double t = (double)i / (double)cone_stacks;
				double phi = t * M_PI;
				double dr = r;
				double dh = t * height;
				for (int j = 0; j <= cone_slices; ++j)
				{
					double U = (double)j / (double)cone_slices;
					double theta = 2.0 * U * M_PI;

					// use spherical coordinates to calculate the positions.
					double x = dr * cos(theta);
					double y = dh;
					double z = dr * sin(theta);
					Vector3 p(x, y, z);
					cone_pos.push_back(p);
				}

			}

			int ver_start = 0;
			// Calc The Index Positions
			for (int i = 0; i < cone_slices * cone_stacks + 1; ++i) {
				indices.push_back(int(i) + ver_start);
				indices.push_back(int(i + cone_slices + 1) + ver_start);
				indices.push_back(int(i + cone_slices) + ver_start);

				indices.push_back(int(i + cone_slices + 1) + ver_start);
				indices.push_back(int(i) + ver_start);
				indices.push_back(int(i + 1) + ver_start);
			}

			Vector3 tl = (epos);
			Vector3 dir = (epos - spos);
			dir.normalize();
			qglviewer::Vec axis = qglviewer::Vec(dir[0], dir[1], dir[2]);
			qglviewer::Vec y_axis(0, 1, 0);
			qglviewer::Quaternion quat(y_axis, -axis);

			BaseFrame frame;
			frame.setPosition(tl[0], tl[1], tl[2]);
			frame.rotate(quat);

			for (int i = 0; i < cone_pos.size(); i++)
			{
				qglviewer::Vec p = qglviewer::Vec(cone_pos[i][0], cone_pos[i][1], cone_pos[i][2]);
				qglviewer::Vec wp = frame.inverseCoordinatesOf(p);
				pList.push_back(Vector3(wp.x, wp.y, wp.z));
			}
			for (int i = 0; i < indices.size() / 3; i++)
			{
				fList.push_back(Vector3i(indices[3 * i], indices[3 * i + 2], indices[3 * i + 1]));
			}
		}

		virtual void saveSimulator(int frame)
		{
			std::ofstream fout;
			std::string filename = _sceneDir + "saveDat/";
			int code = mkdir(filename.c_str());
			if (frame < 10)
				filename += "frame_000" + QString().setNum(frame).toStdString() + ".dat";
			else if (frame < 100)
				filename += "frame_00" + QString().setNum(frame).toStdString() + ".dat";
			else if (frame < 1000)
				filename += "frame_0" + QString().setNum(frame).toStdString() + ".dat";
			else
				filename += "frame_" + QString().setNum(frame).toStdString() + ".dat";

			fout.open(filename, std::ios::binary);
			EigenMatrixIO::write_binary(fout, points);
			fout.close();


			filename = _sceneDir + "saveMouseForce/";
			code = mkdir(filename.c_str());
			if (frame < 10)
				filename += "mouse_000" + QString().setNum(frame).toStdString() + ".txt";
			else if (frame < 100)
				filename += "mouse_00" + QString().setNum(frame).toStdString() + ".txt";
			else if (frame < 1000)
				filename += "mouse_0" + QString().setNum(frame).toStdString() + ".txt";
			else
				filename += "mouse_" + QString().setNum(frame).toStdString() + ".txt";
			fout.open(filename);
			fout << mousePointsIds.size() << std::endl;
			fout << mouseForceLine[0] << " " << mouseForceLine[1] <<" " << mouseForceLine[2] << std::endl;
			for (int i = 0; i < mousePointsIds.size(); i++)
				fout << mousePointsIds[i] << " " << mousePointsList[i][0] << " " << mousePointsList[i][1] << " " << mousePointsList[i][2] << std::endl;
			fout.close();


		}

		virtual void run(int frame = 0);
		virtual void postRun();
		virtual void initialization();
		virtual void setupHessMat();
		virtual void gpuGeneralInit();
		virtual void gpuSolverInit();
		virtual void gpuCollisionInit();
		void getSimPointsNeighborsList(std::vector<std::vector<int>>& simPointsNeighborsList);
		void doTime(int frame = 0);
		virtual void advance(const int frame, const qeal dt);
		virtual int LocalGlobalIteration(const int frame, const qeal dt);
		virtual void doLocalProjection(int frame, const qeal dt, ConstraintEnergyInfo& energy);
		virtual void gpuCcdCulling(qeal* devX, qeal* devDir);
	};



	template <class T1, class T2>
	void flatten2DArray(std::vector<std::vector<T1>>& array2D, std::vector<T2>& outArray, std::vector<int>& num, std::vector<int>& offset)
	{
		num.resize(array2D.size());
		offset.resize(array2D.size());

		int bufferOffset = 0;
		for (int i = 0; i < array2D.size(); i++)
		{
			num[i] = array2D[i].size();
			offset[i] = bufferOffset;
			for (int j = 0; j < array2D[i].size(); j++)
				outArray.push_back(T2(array2D[i][j]));
			bufferOffset += array2D[i].size();
		}
	}
}


#endif