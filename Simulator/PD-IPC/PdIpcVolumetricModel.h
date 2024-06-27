#pragma once
#ifndef PD_IPC_VOLUMETRIC_MODEL
#define PD_IPC_VOLUMETRIC_MODEL

#include "PdIpcModel.h"
#include "Commom\Utils.h"

namespace PD_IPC
{
	class PdIpcVolumetricModel : public PdIpcModel
	{
	public:
		qeal strainStiffness = 1e3;
		qeal volumeStiffness = 1e3;
		qeal positionStiffness = 1e6;

		PdIpcVolumetricModel(){
			modelType = Volumetric;
		}

		PdIpcVolumetricModel(const qeal s0, const qeal s1, const qeal s2):
			strainStiffness(s0), volumeStiffness(s1), positionStiffness(s2)
		{
			modelType = Volumetric;
		}

		virtual void registerBoundaryInfo()
		{
			if (loadBoundaryInfo())
				return;

			std::vector<Vector3i> elementBoundary(4);
			elementBoundary[0] = Vector3i(3, 1, 0);
			elementBoundary[1] = Vector3i(3, 2, 1);
			elementBoundary[2] = Vector3i(3, 0, 2);
			elementBoundary[3] = Vector3i(2, 1, 0);

			std::unordered_map<Vector3i, std::vector<Vector4i>, VectorHash3i> boundaryFaceMap;
			for (int ele = 0; ele < tetElementNum; ele++)
			{
				Vector4i indices = getTetElement(ele);
				for (int k = 0; k < 4; k++)
				{
					Vector3i face;
					for (int j = 0; j < 3; j++)
						face[j] = indices[elementBoundary[k][j]];
					Vector3i sort_face = face;
					std::sort(sort_face.data(), sort_face.data() + 3);
					boundaryFaceMap[sort_face].push_back(Vector4i(face[0], face[1], face[2], ele));
				}
			}

			for (auto it : boundaryFaceMap)
				if (it.second.size() == 1)
				{
					Vector4i data = it.second[0];
					Vector3i face = Vector3i(data[0], data[1], data[2]);
					int eleId = data[3];

					Vector4i indices = getTetElement(eleId);
					Vector3 p0 = getTetPoint(indices[0]);
					Vector3 p1 = getTetPoint(indices[1]);
					Vector3 p2 = getTetPoint(indices[2]);
					Vector3 p3 = getTetPoint(indices[3]);
					Vector3 cent = (p0 + p1 + p2 + p3) / 4;

					Vector3 fp0 = getTetPoint(face[0]);
					Vector3 fp1 = getTetPoint(face[1]);
					Vector3 fp2 = getTetPoint(face[2]);
					Vector3 fcent = (fp0 + fp1 + fp2) / 3;

					Vector3 normal = (fp1 - fp0).cross(fp2 - fp0);
					normal.normalize();

					Vector3 centDir = cent - fcent;

					centDir.normalize();
					if (centDir.dot(normal) > 0)
					{
						face = Vector3i(data[0], data[2], data[1]);
					}
					boundaryFaces.push_back(face);
				}
					
			///
			std::unordered_set<int> boundaryPointsSet;
			std::unordered_set<Vector2i, VectorHash2i> boundaryEdgeSet;

			for (int f = 0; f < boundaryFaces.size(); f++)
			{
				Vector3i face = boundaryFaces[f];
				boundaryPointsSet.insert(face[0]);
				boundaryPointsSet.insert(face[1]);
				boundaryPointsSet.insert(face[2]);

				Vector2i edge01(face[0], face[1]);
				std::sort(edge01.data(), edge01.data() + 2);
				Vector2i edge12(face[1], face[2]);
				std::sort(edge12.data(), edge12.data() + 2);
				Vector2i edge20(face[2], face[0]);
				std::sort(edge20.data(), edge20.data() + 2);

				boundaryEdgeSet.insert(edge01);
				boundaryEdgeSet.insert(edge12);
				boundaryEdgeSet.insert(edge20);
			}

			boundaryPoints = std::vector<int>(boundaryPointsSet.begin(), boundaryPointsSet.end());
			for (auto it : boundaryEdgeSet)
				boundaryEdges.push_back(it);

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
				Vector3 p = getTetPoint(boundaryPoints[i]);
				for (int k = 0; k < 3; k++)
				{
					min[k] = min2(min[k], p[k]);
					max[k] = max2(max[k], p[k]);
				}
			}
			qeal patchSize[3];
			int len = collisionPatchDim;
			for (int k = 0; k < 3; k++)
				patchSize[k] = max[k] - min[k];
			int gridSize[3];
			for (int k = 0; k < 3; k++)
				gridSize[k] = len;
			std::vector<std::array<Vector3, 2>> patch_bbox;
			for (int s = 0; s < gridSize[0]; s++)
				for (int r = 0; r < gridSize[1]; r++)
					for (int c = 0; c < gridSize[2]; c++)
					{
						Vector3 bbox_min = min + Vector3(s * (patchSize[0] / len), r * (patchSize[1] / len), c *  (patchSize[2] / len));
						Vector3 bbox_max = bbox_min + Vector3((patchSize[0] / len), (patchSize[1] / len), (patchSize[2] / len));
						std::array<Vector3, 2> pbbox;
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
					cent += getTetPoint(face[k]);
				cent /= 3;
				for (int j = 0; j < patch_bbox.size(); j++)
					if (cent[0] >= patch_bbox[j][0][0] &&
						cent[0] <= patch_bbox[j][1][0] &&
						cent[1] >= patch_bbox[j][0][1] &&
						cent[1] <= patch_bbox[j][1][1] &&
						cent[2] >= patch_bbox[j][0][2] &&
						cent[2] <= patch_bbox[j][1][2])
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
			int& simElementsNum,
			int& boundaryPointsNum,
			int& boundaryEdgesNum,
			int& boundaryFacesNum,
			std::vector<std::vector<Vector3>>& simPoints,
			std::vector<std::vector<Vector4i>>& volumetricElements,
			std::vector<std::vector<qeal>>& simPointsMass,
			std::vector<std::vector<qeal>>& simPointsInvMass,
			std::vector<std::vector<int>>& boundaryPointsList,
			std::vector<std::vector<Vector2i>>& boundaryEdgesList,
			std::vector<std::vector<Vector3i>>& boundaryFacesList,
			std::vector<std::vector<Vector4i>>& volumetricStrainConsList,
			std::vector<qeal>& volumetricStrainStiffness,
			std::vector<std::vector<Vector4i>>& volumetricVolumeConsList,
			std::vector<qeal>& volumetricVolumeStiffness,
			std::vector<std::vector<int>>& positionConsList,
			std::vector<qeal>& positionStiffnessList
		)
		{
			PdIpcModel::init();

			usePointsNum = tetPointsNum;

			std::vector<Vector3> tempPoints(tetPointsNum);
			for (int i = 0; i < tetPointsNum; i++)
			{
				Vector3 p = getTetPoint(i);
				tempPoints[i] = p;
			}
			simPoints.push_back(tempPoints);

			std::vector<Vector4i> tempElements(tetElementNum);
			for (int i = 0; i < tetElementNum; i++)
			{
				Vector4i ele = getTetElement(i);
				tempElements[i] = ele + Vector4i(simPointsNum, simPointsNum, simPointsNum, simPointsNum);
			}
			volumetricElements.push_back(tempElements);
			//
			std::vector<qeal> tempPointsMass(tetPointsNum);
			std::vector<qeal> tempPointsInvMass(tetPointsNum);
			for (int i = 0; i < tetPointsNum; i++)
			{
				qeal vol = getTetMeshHandle()->getTetNodeParam(i)->volume;
				tempPointsMass[i] = vol * density;
				tempPointsInvMass[i] = 1.0 / (vol * density);	
			}
			simPointsMass.push_back(tempPointsMass);
			simPointsInvMass.push_back(tempPointsInvMass);
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

			std::vector<int> tempPositionCons(fixedIds.size());
			for (int i = 0; i < fixedIds.size(); i++)
			{
				tempPositionCons[i] = fixedIds[i] + simPointsNum;
			}
			positionConsList.push_back(tempPositionCons);
			positionStiffnessList.push_back(positionStiffness);

			volumetricStrainConsList.push_back(tempElements);
			volumetricStrainStiffness.push_back(strainStiffness);

			simPointsNum += usePointsNum;
			simElementsNum += tetElementNum;
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
			usePointsNum = tetPointsNum;
			std::vector<Vector3> tempPoints(tetPointsNum);
			for (int i = 0; i < tetPointsNum; i++)
			{
				Vector3 p = getTetPoint(i);
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

			simPointsNum += tetPointsNum;
			boundaryPointsNum += boundaryPoints.size();
			boundaryEdgesNum += boundaryEdges.size();
			boundaryFacesNum += boundaryFaces.size();

		}


	};



}


#endif