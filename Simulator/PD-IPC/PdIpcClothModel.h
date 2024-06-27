#pragma once
#ifndef PD_IPC_CLOTH_MODEL
#define PD_IPC_CLOTH_MODEL

#include "PdIpcModel.h"

namespace PD_IPC
{
	class PdIpcClothModel : public PdIpcModel
	{
	public:
		qeal strainStiffness = 1e3;
		qeal bendingStiffness = 1e-3;
		qeal positionStiffness = 1e6;

		PdIpcClothModel(){ modelType = Cloth; }

		PdIpcClothModel(const qeal s0, const qeal s1, const qeal s2) :
			strainStiffness(s0), bendingStiffness(s1), positionStiffness(s2)
		{
			modelType = Cloth;
		}

		virtual void registerBoundaryInfo()
		{
			if (loadBoundaryInfo())
				return;
			boundaryPoints.resize(pointsNum);
			boundaryEdges.resize(edgesNum);
			boundaryFaces.resize(facesNum);

			for (int i = 0; i < pointsNum; i++)
				boundaryPoints[i] = i;
			for (int i = 0; i < edgesNum; i++)
			{
				boundaryEdges[i] = getSurfaceEdge(i);
			}
			for (int i = 0; i < facesNum; i++)
			{
				boundaryFaces[i] = getSurfaceFace(i);
			}
			saveBoundaryInfo();
		}

		virtual void setCollisionPatches()
		{
			if (loadCollisionPatchInfo())
				return;

			Vector3 min = Vector3(DBL_MAX, DBL_MAX, DBL_MAX);
			Vector3 max = -Vector3(DBL_MAX, DBL_MAX, DBL_MAX);
			for (int i = 0; i < boundaryPoints.size(); i++)
			{
				Vector3 p = getSurfacePoint(boundaryPoints[i]);
				for (int k = 0; k < 3; k++)
				{
					min[k] = min2(min[k], p[k]);
					max[k] = max2(max[k], p[k]);
				}
			}

			int len = std::sqrt((qeal)pointsNum) < collisionPatchDim ? std::sqrt((qeal)pointsNum) : collisionPatchDim;
			qeal patchSize[3];
			for (int k = 0; k < 3; k++)
				patchSize[k] = max[k] - min[k];

			int u, v;
			if (patchSize[0] < patchSize[1] && patchSize[0] < patchSize[2])
			{
				u = 1;
				v = 2;
			}
			if (patchSize[1] < patchSize[0] && patchSize[1] < patchSize[2])
			{
				u = 0;
				v = 2;
			}
			if (patchSize[2] < patchSize[0] && patchSize[2] < patchSize[1])
			{
				u = 0;
				v = 1;
			}

			//std::ofstream fout("clotex.txt");
			//for (int i = 0; i < pointsNum; i++)
			//{
			//	Vector3 p = getSurfacePoint(i);
			//	qeal tu = (p[u] - min[u]) / (max[u] - min[u]);
			//	qeal tv = (p[v] - min[v]) / (max[v] - min[v]);
			//	fout << "vt " << tu << " " << tv << std::endl;
			//}
			//fout.close();
			//system("pause");
			int gridSize[3];
			for (int k = 0; k < 3; k++)
			{
				if (k == u || k == v)
					gridSize[k] = len;
				else gridSize[k] = 1;
			}

			std::vector<std::array<Vector2, 2>> patch_bbox;
			for (int r = 0; r < gridSize[u]; r++)
				for (int c = 0; c < gridSize[v]; c++)
				{
					Vector2 bbox_min = Vector2(min[u], min[v]) + Vector2(r * (patchSize[u] / len), c *  (patchSize[v] / len));
					Vector2 bbox_max = bbox_min + Vector2(patchSize[u] / len, patchSize[v] / len);
					std::array<Vector2, 2> pbbox;
					pbbox[0] = bbox_min;
					pbbox[1] = bbox_max;
					patch_bbox.push_back(pbbox);
				}

			std::vector<std::vector<int>> facePocket(patch_bbox.size());
			for (int i = 0; i < boundaryFaces.size(); i++)
			{
				Vector3i face = boundaryFaces[i];
				Vector3 cent = Vector3(0, 0, 0);
				for (int k = 0; k < 3; k++)
					cent += getSurfacePoint(face[k]);
				cent /= 3;
				for (int j = 0; j < patch_bbox.size(); j++)
					if (cent[u] >= patch_bbox[j][0][0] && cent[u] <= patch_bbox[j][1][0] && cent[v] >= patch_bbox[j][0][1] && cent[v] <= patch_bbox[j][1][1])
					{
						facePocket[j].push_back(i);
						break;
					}
			}

			for (int i = 0; i < facePocket.size(); i++)
			{
				if (facePocket[i].size() > 0)
				{
					collisionPatches.push_back(facePocket[i]);
				}
					
			}

			saveCollisionPatchInfo();
		}

		virtual void init
		(
			int& simPointsNum,
			int& simClothPointsNum,
			int& boundaryPointsNum,
			int& boundaryEdgesNum,
			int& boundaryFacesNum,
			std::vector<std::vector<Vector3>>& simPoints,
			std::vector<std::vector<qeal>>& simPointsMass,
			std::vector<std::vector<qeal>>& simPointsInvMass,
			std::vector<std::vector<int>>& boundaryPointsList,
			std::vector<std::vector<Vector2i>>& boundaryEdgesList,
			std::vector<std::vector<Vector3i>>& boundaryFacesList,
			std::vector<qeal>& clothPointsBendingArea,
			std::vector<std::vector<int>>& clothPointsNeighborsList,
			std::vector<std::vector<qeal>>& clothPointsStrainStiffness,
			std::vector<std::vector<qeal>>& clothPointsBendingStiffness,
			std::vector<std::vector<qeal>>& clothPointsNeighborsRestLength,
			std::vector<std::vector<qeal>>& clothPointsNeighborsLbo,
			std::vector<qeal>& clothPointsMeanCurvatureNorm,
			std::vector<std::vector<int>>& positionConsList,
			std::vector<qeal>& positionStiffnessList
		)
		{
			PdIpcModel::init();
			usePointsNum = pointsNum;

			std::vector<Vector3> tempPoints(pointsNum);
			for (int i = 0; i < pointsNum; i++)
			{
				Vector3 p = getSurfacePoint(i);
				tempPoints[i] = p;
			}
			simPoints.push_back(tempPoints);

			std::vector<qeal> faceArea(facesNum, 0), pointsArea(pointsNum, 0);
			for (int f = 0; f < facesNum; f++)
			{
				Vector3i face = getSurfaceFace(f);
				qeal area = TriangleArea(getSurfacePoint(face[0]), getSurfacePoint(face[1]), getSurfacePoint(face[2]));
				faceArea[f] = area;
				for (int k = 0; k < 3; k++)
					pointsArea[face[k]] += area / 3;
			}
			std::vector<qeal> pointsMass(pointsNum);
			std::vector<qeal> pointsInvMass(pointsNum);
			for (int v = 0; v < pointsNum; v++)
			{
				pointsMass[v] = pointsArea[v] * density;
				pointsInvMass[v] = 1.0 / pointsMass[v];
				clothPointsBendingArea.push_back(pointsArea[v] * bendingStiffness);
			}
			simPointsMass.push_back(pointsMass);
			simPointsInvMass.push_back(pointsInvMass);

			//
			std::vector<int> tempBoundaryPoints(boundaryPoints.size());
			for (int i = 0; i < boundaryPoints.size(); i++)
			{
				tempBoundaryPoints[i] = boundaryPoints[i] + simPointsNum;
			}
			boundaryPointsList.push_back(tempBoundaryPoints);

			std::vector<Vector2i> tempBoundaryEdges(boundaryEdges.size());
			for (int i = 0; i < boundaryEdges.size(); i++)
			{
				tempBoundaryEdges[i] = boundaryEdges[i] + Vector2i(simPointsNum, simPointsNum);
			}
			boundaryEdgesList.push_back(tempBoundaryEdges);

			std::vector<Vector3i> tempBoundaryFaces(boundaryFaces.size());
			for (int i = 0; i < tempBoundaryFaces.size(); i++)
			{
				tempBoundaryFaces[i] = boundaryFaces[i] + Vector3i(simPointsNum, simPointsNum, simPointsNum);
			}
			boundaryFacesList.push_back(tempBoundaryFaces);
			//
			///////////////////////////////////////////////////////////////////////////////
			std::vector<qeal> restEdgeLenth(edgesNum);
			std::vector<qeal> restEdgeCotWeight(edgesNum);
			std::vector<std::vector<int>> edgeOppositePoints(edgesNum);

			for (int i = 0; i < edgesNum; i++)
			{
				Vector2i edge = getSurfaceEdge(i);
				for (int j = 0; j < edgeFaceIndices[i].span; j++)
				{
					int fid = edgeFaceIndices[i].buffer[j];
					Vector3i face = getSurfaceFace(fid);
					int v = -1;
					for (int k = 0; k < 3; k++)
					{
						if (face[k] != edge[0] && face[k] != edge[1])
							v = face[k];
						if (v != -1) break;
					}
					edgeOppositePoints[i].push_back(v);
				}
			}

			for (int e = 0; e < edgesNum; e++)
			{
				Vector2i edge = getSurfaceEdge(e);
				Vector3 p1 = getSurfacePoint(edge[0]);
				Vector3 p2 = getSurfacePoint(edge[1]);
				restEdgeLenth[e] = (p1 - p2).norm();

				qeal cotan0 = 0, cotan1 = 0;
				for (int k = 0; k < edgeOppositePoints[e].size(); k++)
				{
					int opp = edgeOppositePoints[e][k];
					Vector3 p = getSurfacePoint(opp);
					Vector3 x10 = p1 - p;
					Vector3 x20 = p2 - p;
					qeal len10 = x10.norm();
					qeal len20 = x20.norm();
					qeal theta = acos(x10.dot(x20) / (len10 * len20));
					qeal cotan = 1.0 / tan(theta);
					restEdgeCotWeight[e] += -0.5 * cotan;
				}
			}

			std::vector<std::vector<int>> pointsNeighbors(pointsNum);
			std::vector<std::vector<qeal>> pointsLBO(pointsNum);
			std::vector < std::vector <qeal>> pointsNeghborsRestLength(pointsNum);
			std::vector<qeal> restPointMeanCurvatureNorm(pointsNum);
			std::vector<std::vector<qeal>> clothStrainStiffness(pointsNum);
			std::vector<std::vector<qeal>> clothBendingStiffness(pointsNum);

			for (int v = 0; v < pointsNum; v++)
			{
				Vector3 pv = getSurfacePoint(v);
				int neighbor = pointEdgeIndices[v].span;
				pointsLBO[v].resize(neighbor + 1);
				pointsNeighbors[v].resize(neighbor + 1);
				pointsNeighbors[v][0] = v + simPointsNum;
				pointsNeghborsRestLength[v].resize(neighbor + 1, 0);
				pointsLBO[v].resize(neighbor + 1, 0);

				clothStrainStiffness[v].resize(neighbor + 1, 0);
				clothStrainStiffness[v][0] = strainStiffness;
				clothBendingStiffness[v].resize(neighbor + 1, 0);

				qeal lbo = 0;
				VectorX q[3];
				q[0].resize(neighbor + 1); q[0][0] = pv[0];
				q[1].resize(neighbor + 1); q[1][0] = pv[1];
				q[2].resize(neighbor + 1); q[2][0] = pv[2];

				for (int k = 0; k < neighbor; k++)
				{
					int e = pointEdgeIndices[v].buffer[k];
					Vector2i edge = getSurfaceEdge(e);
					int ve = (v == edge[0] ? edge[1] : edge[0]);
					pointsNeighbors[v][k + 1] = ve + simPointsNum;
					Vector3 pe =getSurfacePoint(ve);
					qeal restLength = (pv - pe).norm();
					pointsNeghborsRestLength[v][k + 1] = restLength;
					clothStrainStiffness[v][k + 1] = strainStiffness;
					qeal ew = restEdgeCotWeight[e] / pointsArea[v];
					lbo += ew;
					pointsLBO[v][1 + k] = ew;
					q[0][k + 1] = pe[0]; q[1][k + 1] = pe[1]; q[2][k + 1] = pe[2];
				}
				pointsLBO[v][0] = -lbo;
				Vector3 pointCurvature;
				pointCurvature.setZero();

				for (int k = 0; k < neighbor + 1; k++)
				{
					pointCurvature[0] += pointsLBO[v][k] * q[0][k];
					pointCurvature[1] += pointsLBO[v][k] * q[1][k];
					pointCurvature[2] += pointsLBO[v][k] * q[2][k];
					clothBendingStiffness[v][k] = pointsLBO[v][k] * bendingStiffness * pointsArea[v];
				}
				restPointMeanCurvatureNorm[v] = pointCurvature.norm();
			}

			for (int v = 0; v < pointsNum; v++)
			{
				clothPointsNeighborsList.push_back(pointsNeighbors[v]);
				clothPointsStrainStiffness.push_back(clothStrainStiffness[v]);
				clothPointsBendingStiffness.push_back(clothBendingStiffness[v]);
				clothPointsNeighborsRestLength.push_back(pointsNeghborsRestLength[v]);
				clothPointsNeighborsLbo.push_back(pointsLBO[v]);
				clothPointsMeanCurvatureNorm.push_back(restPointMeanCurvatureNorm[v]);
			}
			
			std::vector<int> tempPositionCons(fixedIds.size());
			for (int i = 0; i < fixedIds.size(); i++)
			{
				tempPositionCons[i] = fixedIds[i] + simPointsNum;
			}
			positionConsList.push_back(tempPositionCons);
			positionStiffnessList.push_back(positionStiffness);

			simPointsNum += usePointsNum;
			simClothPointsNum += usePointsNum;
			boundaryPointsNum += boundaryPoints.size();
			boundaryEdgesNum += boundaryEdges.size();
			boundaryFacesNum += boundaryFaces.size();
		}

		virtual void init
		(
			int& simPointsNum,
			int& boundaryPointsNum,
			int& boundaryEdgesNum,
			int& boundaryFacesNum,
			std::vector<std::vector<Vector3>>& simPoints,
			std::vector<std::vector<int>>& boundaryPointsList,
			std::vector<std::vector<Vector2i>>& boundaryEdgesList,
			std::vector<std::vector<Vector3i>>& boundaryFacesList
		)
		{
			PdIpcModel::init();
			std::vector<Vector3> tempPoints(pointsNum);
			for (int i = 0; i < pointsNum; i++)
			{
				Vector3 p = getSurfacePoint(i);
				tempPoints[i] = p;
			}
			simPoints.push_back(tempPoints);

			std::vector<int> tempBoundaryPoints(boundaryPoints.size());
			for (int i = 0; i < boundaryPoints.size(); i++)
			{
				tempBoundaryPoints[i] = boundaryPoints[i] + simPointsNum;
			}
			boundaryPointsList.push_back(tempBoundaryPoints);

			std::vector<Vector2i> tempBoundaryEdges(boundaryEdges.size());
			for (int i = 0; i < boundaryEdges.size(); i++)
			{
				tempBoundaryEdges[i] = boundaryEdges[i] + Vector2i(simPointsNum, simPointsNum);
			}
			boundaryEdgesList.push_back(tempBoundaryEdges);

			std::vector<Vector3i> tempBoundaryFaces(boundaryFaces.size());
			for (int i = 0; i < tempBoundaryFaces.size(); i++)
			{
				tempBoundaryFaces[i] = boundaryFaces[i] + Vector3i(simPointsNum, simPointsNum, simPointsNum);
			}
			boundaryFacesList.push_back(tempBoundaryFaces);
			simPointsNum += pointsNum;
			boundaryPointsNum += boundaryPoints.size();
			boundaryEdgesNum += boundaryEdges.size();
			boundaryFacesNum += boundaryFaces.size();
		}


	};


}

#endif