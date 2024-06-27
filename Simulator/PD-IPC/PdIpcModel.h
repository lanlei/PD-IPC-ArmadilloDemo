#pragma once
#ifndef PD_IPC_MODEL
#define PD_IPC_MODEL
#include "Model\BaseModel.h"
#include <queue>

namespace PD_IPC
{
	enum ModelType
	{
		Volumetric = 0,
		Cloth = 1
	};
	
	class PdIpcModel : public BaseModel
	{
	protected:
		int collisionPatchDim = 10;
		ModelType modelType;
	public:
		qeal density;
		std::vector<int> fixedIds;
		int usePointsNum;
		int boundaryPointsNum;
		int boundaryEdgesNum;
		int boundaryFacesNum;
		std::vector<int> boundaryPoints;
		std::vector<Vector2i> boundaryEdges;
		std::vector<Vector3i> boundaryFaces;	

		std::vector<std::vector<int>> boundaryPointsSharedFaces;
		std::vector<std::vector<int>> boundaryEdgesSharedFaces;
		std::vector<std::vector<int>> boundaryFacesSharedEdges;
		std::vector<std::vector<int>> boundaryFacesSharedPoints;
		std::vector<bool> boundaryEdgesFlag;

		std::vector<std::vector<int>> collisionPatches;


		ModelType getModelType()
		{
			return modelType;
		}


		void setCollisionPatchDim(int dim)
		{
			collisionPatchDim = dim;
		}

		PdIpcModel() : density(1)
		{

		}

		virtual void readFixedInfo()
		{
			std::string filename = dir + "fixed_points.txt";
			std::ifstream fin(filename.c_str());
			if (fin.is_open())
			{
				int num = 0;
				char ch;
				fin >> ch >> num;
				if (ch == 'c')
				{
					for (int j = 0; j < num; j++)
					{
						int vid;
						fin >> vid;
						fixedIds.push_back(vid);
					}
				}
			}
			fin.close();
			std::unique(fixedIds.begin(), fixedIds.end());
		}

		virtual void registerBoundaryInfo()
		{

		}

		virtual void postProcessBoundaryInfo()
		{
			boundaryPointsNum = boundaryPoints.size();
			boundaryEdgesNum = boundaryEdges.size();
			boundaryFacesNum = boundaryFaces.size();
			boundaryPointsSharedFaces.resize(boundaryPointsNum);

			std::map<int, int> pointsIds;
			for (int i = 0; i < boundaryPointsNum; i++)
			{
				pointsIds[boundaryPoints[i]] = i;
			}

			for (int i = 0; i < boundaryFacesNum; i++)
			{
				Vector3i face = boundaryFaces[i];
				for (int k = 0; k < 3; k++)
				{
					int p = pointsIds[face[k]];					
					boundaryPointsSharedFaces[p].push_back(i);
				}					
			}

			std::map<Vector2i, int, Vec2IndComp> edgesIds;
			for (int i = 0; i < boundaryEdgesNum; i++)
			{
				edgesIds[boundaryEdges[i]] = i;
			}

			boundaryEdgesSharedFaces.resize(boundaryEdgesNum);
			for (int i = 0; i < boundaryFacesNum; i++)
			{
				Vector3i face = boundaryFaces[i];
				Vector2i edge01(face[0], face[1]);
				std::sort(edge01.data(), edge01.data() + 2);
				Vector2i edge12(face[1], face[2]);
				std::sort(edge12.data(), edge12.data() + 2);
				Vector2i edge20(face[2], face[0]);
				std::sort(edge20.data(), edge20.data() + 2);

				int e0 = edgesIds[edge01];
				int e1 = edgesIds[edge12];
				int e2 = edgesIds[edge20];
				boundaryEdgesSharedFaces[e0].push_back(i);
				boundaryEdgesSharedFaces[e1].push_back(i);
				boundaryEdgesSharedFaces[e2].push_back(i);
			}
			boundaryEdgesFlag.resize(boundaryEdgesNum, false);
			for (int i = 0; i < boundaryEdgesNum; i++)
				if (boundaryEdgesSharedFaces[i].size() == 1)
					boundaryEdgesFlag[i] = true;

			boundaryFacesSharedEdges.resize(boundaryFacesNum);
			for (int i = 0; i < boundaryEdgesNum; i++)
			{
				int min_fid = boundaryEdgesSharedFaces[i][0];
				int min_size = boundaryFacesSharedEdges[min_fid].size();

				for (int k = 1; k < boundaryEdgesSharedFaces[i].size(); k++)
				{
					int fid = boundaryEdgesSharedFaces[i][k];
					int size = boundaryFacesSharedEdges[fid].size();
					if (size < min_size)
					{
						min_fid = fid;
						min_size = size;
					}
				}
				boundaryFacesSharedEdges[min_fid].push_back(i);
			}

			boundaryFacesSharedPoints.resize(boundaryFacesNum);
			for (int i = 0; i < boundaryPointsNum; i++) 
			{				
				int min_fid = boundaryPointsSharedFaces[i][0];
				int min_size = boundaryFacesSharedPoints[min_fid].size();
				for (int k = 1; k < boundaryPointsSharedFaces[i].size(); k++)
				{
					int fid = boundaryPointsSharedFaces[i][k];
					int size = boundaryFacesSharedPoints[fid].size();
					if (size < min_size)
					{
						min_fid = fid;
						min_size = size;
					}
				}
				boundaryFacesSharedPoints[min_fid].push_back(boundaryPoints[i]);
			}
		}

		virtual void setCollisionPatches()
		{

		}

		bool loadBoundaryInfo()
		{
			std::string boundaryFile = dir + "boundary.dat";
			std::ifstream fin(boundaryFile.c_str(), std::ios::binary);
			if (!fin.is_open())
				return false;
			int bpNum;
			fin.read((char*)(&bpNum), sizeof(int));
			boundaryPoints.resize(bpNum);
			for (int j = 0; j < boundaryPoints.size(); j++)
			{
				int vid;
				fin.read((char*)(&vid), sizeof(int));
				boundaryPoints[j] = vid;
			}
			int beNum;
			fin.read((char*)(&beNum), sizeof(int));
			boundaryEdges.resize(beNum);
			for (int j = 0; j < boundaryEdges.size(); j++)
			{
				int vi, vj;
				fin.read((char*)(&vi), sizeof(int));
				fin.read((char*)(&vj), sizeof(int));
				boundaryEdges[j] = Vector2i(vi, vj);
			}
			int bfNum;
			fin.read((char*)(&bfNum), sizeof(int));
			boundaryFaces.resize(bfNum);
			for (int j = 0; j < boundaryFaces.size(); j++)
			{
				int vi, vj, vk;
				fin.read((char*)(&vi), sizeof(int));
				fin.read((char*)(&vj), sizeof(int));
				fin.read((char*)(&vk), sizeof(int));
				boundaryFaces[j] = Vector3i(vi, vj, vk);
			}
			fin.close();
			return true;
		}

		void saveBoundaryInfo()
		{
			std::string boundaryFile = dir + "boundary.dat";
			std::ofstream fout(boundaryFile.c_str(), std::ios::binary);
			int bpNum = boundaryPoints.size();
			int beNum = boundaryEdges.size();
			int bfNum = boundaryFaces.size();
			fout.write((char*)&bpNum, sizeof(int));
			for (int j = 0; j < boundaryPoints.size(); j++)
			{
				int vid = boundaryPoints[j];
				fout.write((char*)&vid, sizeof(int));
			}
			fout.write((char*)&beNum, sizeof(int));
			for (int j = 0; j < boundaryEdges.size(); j++)
			{
				int vi = boundaryEdges[j][0];
				int vj = boundaryEdges[j][1];
				fout.write((char*)&vi, sizeof(int));
				fout.write((char*)&vj, sizeof(int));
			}
			fout.write((char*)&bfNum, sizeof(int));
			for (int j = 0; j < boundaryFaces.size(); j++)
			{
				int vi = boundaryFaces[j][0];
				int vj = boundaryFaces[j][1];
				int vk = boundaryFaces[j][2];
				fout.write((char*)&vi, sizeof(int));
				fout.write((char*)&vj, sizeof(int));
				fout.write((char*)&vk, sizeof(int));
			}
			fout.close();
		}

		bool loadCollisionPatchInfo()
		{
			std::string patchFile = dir + "patchInfo.dat";
			std::ifstream fin;
			fin.open(patchFile.c_str(), std::ios::binary);
			if (!fin.is_open())
				return false;

			int bpNum, beNum, bfNum, patchNum;
			fin.read((char*)(&bpNum), sizeof(int));
			fin.read((char*)(&beNum), sizeof(int));
			fin.read((char*)(&bfNum), sizeof(int));

			if (!(bpNum == boundaryPoints.size() && beNum == boundaryEdges.size() && bfNum == boundaryFaces.size()))
			{
				fin.close();
				return false;
			}

			fin.read((char*)(&patchNum), sizeof(int));
			collisionPatches.resize(patchNum);
			for (int j = 0; j < patchNum; j++)
			{
				int fn;
				fin.read((char*)(&fn), sizeof(int));
				collisionPatches[j].resize(fn);
				for (int k = 0; k < fn; k++)
				{
					int fid;
					fin.read((char*)(&fid), sizeof(int));
					collisionPatches[j][k] = fid;
				}
			}
			fin.close();

			return true;
		}

		void saveCollisionPatchInfo()
		{
			std::string patchFile = dir + "patchInfo.dat";
			int n = 0;
			std::ofstream fout(patchFile.c_str(), std::ios::binary);

			int bpNum = boundaryPoints.size();
			int beNum = boundaryEdges.size();
			int bfNum = boundaryFaces.size();
			int patchNum = collisionPatches.size();

			fout.write((char*)&bpNum, sizeof(int));
			fout.write((char*)&beNum, sizeof(int));
			fout.write((char*)&bfNum, sizeof(int));
			fout.write((char*)&patchNum, sizeof(int));
			for (int j = 0; j < collisionPatches.size(); j++)
			{
				int fn = collisionPatches[j].size();
				n += fn;
				fout.write((char*)&fn, sizeof(int));
				for (int k = 0; k < fn; k++)
				{
					int fid = collisionPatches[j][k];
					fout.write((char*)&fid, sizeof(int));
				}
			}
			fout.close();
		}

		virtual void init
		()
		{
			readFixedInfo();
			registerBoundaryInfo();
			setCollisionPatches();
			postProcessBoundaryInfo();

			//std::vector<Vector3> ps;
			//for (int i = 0; i < boundaryFaces.size(); i++)
			//{
			//	Vector3i face = boundaryFaces[i];
			//	for (int k = 0; k < 3; k++)
			//	{
			//		Vector3 p = getTetPoint(face[k]);
			//		ps.push_back(p);
			//	}
			//}

			//std::ofstream fout("test.off");
			//fout << "OFF" << std::endl;
			//fout << ps.size() << " " << boundaryFaces.size() << " 0" << std::endl;
			//for (int i = 0; i < ps.size(); i++)
			//{
			//	fout << ps[i][0] << " " << ps[i][1] << " " << ps[i][2] << std::endl;
			//}
			//for (int i = 0; i < boundaryFaces.size(); i++)
			//{
			//	fout << "3 " << 3 * i << " " << 3 * i + 1 << " " << 3 * i + 2 << std::endl;
			//}
			//fout.close();
			//system("pause");

				//
			//std::vector<std::set<int>> vIds;
			//std::vector<std::set<int>> fIds;
			//std::vector<bool> vFlag(pointsNum, false);



			//while (true)
			//{
			//	int start = -1;
			//	for (int i = 0; i < pointsNum; i++)
			//	{
			//		if (!vFlag[i])
			//		{
			//			start = i;
			//			break;
			//		}
			//	}
			//	if (start == -1)
			//		break;

			//	std::queue<int> vidQueue;
			//	vidQueue.push(start);
			//	std::set<int> vidList;
			//	std::set<int> faceList;
			//	while (!vidQueue.empty())
			//	{
			//		int vid = vidQueue.front();
			//		vidQueue.pop();

			//		if (vFlag[vid])
			//			continue;
			//		vFlag[vid] = true;
			//		vidList.insert(vid);
			//		for (int k = 0; k < pointFaceIndices[vid].span; k++)
			//		{
			//			int fid = pointFaceIndices[vid].buffer[k];
			//			faceList.insert(fid);
			//			Vector3i face = getSurfaceFace(fid);
			//			for (int f = 0; f < 3; f++)
			//			{
			//				if (!vFlag[face[f]])
			//					vidQueue.push(face[f]);
			//			}

			//		}

			//	}

			//	vIds.push_back(vidList);
			//	fIds.push_back(faceList);
			//}


			//int patch = vIds.size();
			//for (int i = 0; i < patch; i++)
			//{
			//	int v = 0;
			//	std::vector<int> newVid(pointsNum, -1);
			//	std::vector<std::pair<int, int>> newOld;
			//	for (auto it = vIds[i].begin(); it != vIds[i].end(); it++, v++)
			//	{
			//		newOld.push_back(std::pair<int, int>(*it, v));
			//		newVid[*it] = v;
			//	}
			//	std::vector<Vector3i> newFacess;
			//	for (auto it = fIds[i].begin(); it != fIds[i].end(); it++)
			//	{
			//		int fid = *it;
			//		Vector3i face = getSurfaceFace(fid);
			//		int newV0 = newVid[face[0]];
			//		int newV1 = newVid[face[1]];
			//		int newV2 = newVid[face[2]];
			//		face = Vector3i(newV0, newV1, newV2);
			//		newFacess.push_back(face);
			//	}
			//				
			//	std::ofstream fout;
			//	std::string filename = "pathc_" + QString().setNum(i).toStdString() + ".obj";
			//	fout.open(filename);
			////	fout << "OFF" << std::endl;
			////	fout << vIds[i].size() << " " << newFacess.size() << " 0" << std::endl;
			//	for (auto it = vIds[i].begin(); it != vIds[i].end(); it++)
			//	{
			//		Vector3 points = getSurfacePoint(*it);
			//		fout << "v " << points[0] << " " << points[1] << " " << points[2] << std::endl;
			//	}

			//	for (int k = 0; k < newFacess.size(); k++)
			//		fout << "f " << newFacess[k][0] +1 << " " << newFacess[k][1] + 1<< " " << newFacess[k][2] + 1<< std::endl;
			//	fout.close();

			//	/*filename = "pathc_" + QString().setNum(i).toStdString() + "_info.txt";
			//	fout.open(filename);
			//	fout << newOld.size() << std::endl;
			//	for (int k = 0; k < newOld.size(); k++)
			//		fout << newOld[i].first << " " << newOld[i].second << std::endl;*/

			//	fout.close();

			//}
			//system("pause");
			

		

		}
	};

	



}

#endif // !PD_IPC_MODELs
