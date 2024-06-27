#include "BaseModel.h"
#include <oneapi/tbb.h>
BaseModel::BaseModel():_valid(false),
	_tx(0),_ty(0),_tz(0),_sx(1.0), _sy(1.0), _sz(1.0), _rx(0),_ry(0), _rz(0), _rsita(0)
{

}

BaseModel::~BaseModel()
{

}

void BaseModel::refreshBuffer(BaseSurfaceMeshBufferPool* sPool, BaseTetMeshBufferPool* tPool, BaseMedialMeshBufferPool* mPool)
{
	((BaseSurfaceMesh*)this)->refresh(sPool);
	if(isTetMeshValid()) ((BaseTetMesh*)this)->refresh(tPool);
	if(isMedialMeshValid()) ((BaseMedialMesh*)this)->refresh(mPool);
}

void BaseModel::init(bool setUniform, bool setOrigin)
{	
	if (setOrigin)
		placeOrigin();
	if(setUniform)
		uniform();
	initVBO();
}

void BaseModel::initMeshesHandel()
{
	_smHandle = new BaseSurfaceHandle(this);
	_smHandle->init();
	if (isTetMeshValid())
	{
		_tmHandle = new BaseTetMeshHandle(this);
		_tmHandle->init();
		computeTetMeshInterpolation();
	}
	if (isMedialMeshValid())
	{
		_mmHandle = new BaseMedialMeshHandle(this);
		_mmHandle->init();
	}
}

qeal BaseModel::uniform()
{
	qeal div = BaseSurfaceMesh::uniform();
	if(isTetMeshValid()) BaseTetMesh::uniform(div);
	if (isMedialMeshValid()) BaseMedialMesh::uniform(div);
	return div;
}

void BaseModel::placeOrigin()
{
	computeBBox();
	qeal cx = (bbox[0] + bbox[3]) / 2.0;
	qeal cy = (bbox[1] + bbox[4]) / 2.0;
	qeal cz = (bbox[2] + bbox[5]) / 2.0;

	BaseSurfaceMesh::placeOrigin(cx, cy, cz);
	if (isTetMeshValid()) BaseTetMesh::placeOrigin(cx, cy, cz);
	if (isMedialMeshValid()) BaseMedialMesh::placeOrigin(cx, cy, cz);
	computeBBox();
}

void BaseModel::conbineTetAndMat(BaseTetMeshBufferPool * tetPool, BaseMedialMeshBufferPool* matPool, const std::string newTetNodeFilename, const std::string newTetElementFilename, const std::string newMedialMeshFilename)
{
	if (!(isTetMeshValid() && isMedialMeshValid()))
		return;
	if (isBindedTetMesh())
	{
		assert(bindedTM.size() == medialPointsNum);
		for (size_t i = 0; i < medialPointsNum; i++)
		{
			if (bindedTM[i] < tetPointsNum)
			{
				Vector3 c = getMedialPoint(i);
				Vector3 p = getTetPoint(bindedTM[i]);
				if ((c - p).norm() > 1e-12)
				{
					setBindedTetMesh(false);
					break;
				}
			}
			else
			{
				setBindedTetMesh(false);
				break;
			}
		}
	}

	if (isBindedTetMesh())
	{
		bindedInverseTM.resize(tetPointsNum, -1);
		for (size_t i = 0; i < bindedTM.size(); i++)
			bindedInverseTM[bindedTM[i]] = i;
		return;
	}

	for (size_t i = 0; i < medialPointsNum; i++)
	{
		Vector3 c = getMedialPoint(i);
		qeal dist = QEAL_MAX;
		Vector3 np;
		int index = -1;
		for (int j = 0; j < tetPointsNum; j++)
		{
			Vector3 v = getTetPoint(j);
			qeal temp = (c - v).norm();
			if (temp < dist)
			{
				dist = temp;
				np = v;
				index = j;
				
			}
		}

		setMedialPoint(i, np.data());
		bindedTM[i] = index;
	}

	//std::vector<int> newTetIndices(4 * tetElementNum);
	//std::vector<qeal> newPointPos(3 * (tetPointsNum + medialPointsNum));

	//for (int i = 0; i < tetPointsNum; i++)
	//{
	//	Vector3 p = getTetPoint(i);
	//	newPointPos[3 * i] = p.data()[0];
	//	newPointPos[3 * i + 1] = p.data()[1];
	//	newPointPos[3 * i + 2] = p.data()[2];
	//}

	//for (int i = 0; i < medialPointsNum; i++)
	//{
	//	bindedTM[i] = tetPointsNum + i;
	//	Vector3 p = getCenterofMedialPoint(i);
	//	newPointPos[3 * (tetPointsNum + i)] = p.data()[0];
	//	newPointPos[3 * (tetPointsNum + i) + 1] = p.data()[1];
	//	newPointPos[3 * (tetPointsNum + i) + 2] = p.data()[2];
	//}

	//for (int i = 0; i < tetElementNum; i++)
	//{
	//	Vector4i indices = getTetElement(i);
	//	for (int j = 0; j < 4; j++)
	//		newTetIndices[4 * i + j] = indices.data()[j];
	//}

	//for (int i = 0; i < medialPointsNum; i++)
	//{
	//	int eleId = insertIntoElement[i];
	//	Vector4i indices = getTetElement(eleId);

	//	int p0 = indices[0];
	//	int p1 = indices[1];
	//	int p2 = indices[2];
	//	int p3 = indices[3];
	//	int p4 = tetPointsNum + i;

	//	// 0-1-2-4
	//	Vector4i indices0 = Vector4i(p0, p1, p2, p4);
	//	// 0-1-3-4
	//	Vector4i indices1 = Vector4i(p0, p1, p3, p4);
	//	// 1-2-3-4
	//	Vector4i indices2 = Vector4i(p1, p2, p3, p4);
	//	// 0-2-3-4
	//	Vector4i indices3 = Vector4i(p0, p2, p3, p4);

	//	for (int j = 0; j < 4; j++)
	//		newTetIndices[4 * eleId + j] = indices0.data()[j];
	//	for (int j = 0; j < 4; j++)
	//		newTetIndices.push_back(indices1.data()[j]);
	//	for (int j = 0; j < 4; j++)
	//		newTetIndices.push_back(indices2.data()[j]);
	//	for (int j = 0; j < 4; j++)
	//		newTetIndices.push_back(indices3.data()[j]);
	//}

	//std::vector<std::vector<int>> newPointTetList(tetPointsNum + medialPointsNum);
	//std::vector<std::vector<int>> newPointNeighborList(tetPointsNum + medialPointsNum);
	//for (int i = 0; i < tetElementNum; i++)
	//{
	//	int n0 = newTetIndices.data()[4 * i];
	//	int n1 = newTetIndices.data()[4 * i + 1];
	//	int n2 = newTetIndices.data()[4 * i + 2];
	//	int n3 = newTetIndices.data()[4 * i + 3];

	//	newPointTetList[n0].push_back(i);
	//	newPointTetList[n1].push_back(i);
	//	newPointTetList[n2].push_back(i);
	//	newPointTetList[n3].push_back(i);

	//	newPointNeighborList[n0].push_back(n1);
	//	newPointNeighborList[n0].push_back(n2);
	//	newPointNeighborList[n0].push_back(n3);

	//	newPointNeighborList[n1].push_back(n0);
	//	newPointNeighborList[n1].push_back(n2);
	//	newPointNeighborList[n1].push_back(n3);

	//	newPointNeighborList[n2].push_back(n0);
	//	newPointNeighborList[n2].push_back(n1);
	//	newPointNeighborList[n2].push_back(n3);

	//	newPointNeighborList[n3].push_back(n0);
	//	newPointNeighborList[n3].push_back(n1);
	//	newPointNeighborList[n3].push_back(n2);
	//}

	//for (size_t i = 0; i < newPointNeighborList.size(); i++)
	//{
	//	std::sort(newPointNeighborList[i].begin(), newPointNeighborList[i].end());
	//	newPointNeighborList[i].erase(std::unique(newPointNeighborList[i].begin(), newPointNeighborList[i].end()), newPointNeighborList[i].end());
	//}

	//tetPool->totalTetPointsNum += medialPointsNum;
	//tetPointsNum += medialPointsNum;

	//updateOverallBufferLastFrag(newPointPos, tetPool->tetPointsBuffer, this->tetPoints);

	//tetPool->totalTetElementNum += newTetIndices.size() / 4 - tetElementNum;
	//tetElementNum += newTetIndices.size() / 4 - tetElementNum;

	//updateOverallBufferLastFrag(newTetIndices, tetPool->tetElementIndicesBuffer, this->tetElementIndices);

	//updateOverallBufferLastFrag(newPointTetList, tetPool->tetPointElementListBuffer, this->tetPointElementList);

	//updateOverallBufferLastFrag(newPointNeighborList, tetPool->tetPointNeighborListBuffer, this->tetPointNeighborList);

	//writeMeshToTetFormat(newTetNodeFilename);
	//writeMeshToTetFormat(newTetElementFilename);
	writeMeshToMatFormat(newMedialMeshFilename);

	bindedInverseTM.resize(tetPointsNum, -1);
	for (size_t i = 0; i < bindedTM.size(); i++)
		bindedInverseTM[bindedTM[i]] = i;
}

void BaseModel::alignAllMesh(qeal* tetPointsPtr)
{
	updateSurfaceFromTetMesh(tetPointsPtr + tetPoints.offset);
	updateMedialMeshFromTetMesh(tetPointsPtr + tetPoints.offset);
}

void BaseModel::computeTetMeshInterpolation()
{
	if (!isTetMeshValid())
		return;
	std::string filename = dir + "tetMeshInterpolation.txt";
	tetMeshInterpolationIndices.resize(pointsNum);
	tetMeshInterpolationWeight.resize(4 * pointsNum);
	// read from file
	std::ifstream fin(filename.c_str());
	if (fin.is_open())
	{
		int num;
		fin >> num;
		assert(num == pointsNum);
		for (int i = 0; i < pointsNum; i++)
		{
			int index;
			int eid;
			qeal w0, w1, w2, w3;
			int temp;
			fin >> index >> eid >> w0 >> w1 >> w2 >> w3 >> temp;
			assert(index == i);
			tetMeshInterpolationIndices[i] = eid;
			tetMeshInterpolationWeight[4 * i] = w0;
			tetMeshInterpolationWeight[4 * i + 1] = w1;
			tetMeshInterpolationWeight[4 * i + 2] = w2;
			tetMeshInterpolationWeight[4 * i + 3] = w3;
		}
		fin.close();
		return;
	}

	for (int i = 0; i < pointsNum; i++)
	{
		Vector3 p = Vector3(points.buffer[3 * i], points.buffer[3 * i + 1], points.buffer[3 * i + 2]);
		int index = -1;
		Vector4 w;
		for (int j = 0; j < tetElementNum; j++)
		{
			w.setZero();
			if (!isInsideTetElement(j, p.data(), w.data()))
				continue;
			index = j;

			qeal sum = w.sum();
			w.data()[0] = w.data()[0] / sum;
			w.data()[1] = w.data()[1] / sum;
			w.data()[2] = w.data()[2] / sum;
			w.data()[3] = w.data()[3] / sum;

			break;
		}

		if (index == -1)
		{
			int tetPointId = searchCloseTetNode(p.data());
			points.buffer[3 * i] = tetPoints.buffer[3 * tetPointId];
			points.buffer[3 * i + 1] = tetPoints.buffer[3 * tetPointId + 1];
			points.buffer[3 * i + 2] = tetPoints.buffer[3 * tetPointId + 2];
			p = Vector3(points.buffer[3 * i], points.buffer[3 * i + 1], points.buffer[3 * i + 2]);
			int eid = tetPointElementList[tetPointId].buffer[0];
			Vector4i ele = getTetElement(eid);
			w.setZero();
			for (int k = 0; k < 4; k++)
				if (tetPointId == ele[k]) w[k] = 1.0;
			index = eid;
		}

		/*
		if (index == -1)
		{
			Vector3 fp = searchCloseTetElement(p.data());
			for (int j = 0; j < tetElementNum; j++)
			{
				w.setZero();
				if (!isInsideTetElement(j, fp.data(), w.data()))
					continue;
				index = j;

				qeal sum = w.sum();
				w.data()[0] = w.data()[0] / sum;
				w.data()[1] = w.data()[1] / sum;
				w.data()[2] = w.data()[2] / sum;
				w.data()[3] = w.data()[3] / sum;

				break;
			}
		}
		*/
		Vector4i ele = getTetElement(index);
		Vector3 v0 = getTetPoint(ele[0]);
		Vector3 v1 = getTetPoint(ele[1]);
		Vector3 v2 = getTetPoint(ele[2]);
		Vector3 v3 = getTetPoint(ele[3]);
		p = w.data()[0] * v0 +w.data()[1] * v1 +w.data()[2] * v2 +w.data()[3] * v3;
		setSurfacePoint(i, p.data());
		tetMeshInterpolationIndices[i] = index;
		tetMeshInterpolationWeight[4 * i] = w.data()[0];
		tetMeshInterpolationWeight[4 * i + 1] = w.data()[1];
		tetMeshInterpolationWeight[4 * i + 2] = w.data()[2];
		tetMeshInterpolationWeight[4 * i + 3] = w.data()[3];
	}
	// write file
	std::ofstream fout(filename.c_str());
	fout << pointsNum << std::endl;
	for (int i = 0; i < pointsNum; i++)
	{
		fout << i << " " << tetMeshInterpolationIndices[i] << " " << tetMeshInterpolationWeight[4 * i] << " " << tetMeshInterpolationWeight[4 * i + 1] << " " << tetMeshInterpolationWeight[4 * i + 2] << " " << tetMeshInterpolationWeight[4 * i + 3] << " " << -1 << std::endl;
	}
	fout.close();
}

void BaseModel::updateSurfaceFromTetMesh()
{
	if (!isTetMeshValid())
		return;
	for (int i = 0; i < pointsNum; i++)
	{
		int eid = tetMeshInterpolationIndices[i];
		for (int j = 0; j < 3; j++)
		{
			points.buffer[3 * i + j] = tetMeshInterpolationWeight[4 * i] * tetPoints.buffer[3 * tetElementIndices.buffer[4 * eid] + j];
		}
		for (int j = 0; j < 3; j++)
		{
			for(int k = 1; k < 4; k++)
				points.buffer[3 * i + j] += tetMeshInterpolationWeight[4 * i + k] * tetPoints.buffer[3 * tetElementIndices.buffer[4 * eid + k] + j];
		}
	}
}

void BaseModel::updateSurfaceFromTetMesh(qeal * currentTetPoints)
{
	if (!isTetMeshValid())
		return;

	tbb::parallel_for(0, pointsNum, [&](int i) {
		int eid = tetMeshInterpolationIndices[i];
		for (int j = 0; j < 3; j++)
		{
			points.buffer[3 * i + j] = tetMeshInterpolationWeight[4 * i] * currentTetPoints[3 * tetElementIndices.buffer[4 * eid] + j];
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 1; k < 4; k++)
				points.buffer[3 * i + j] += tetMeshInterpolationWeight[4 * i + k] * currentTetPoints[3 * tetElementIndices.buffer[4 * eid + k] + j];
		}
	});
}

void BaseModel::updateMedialMeshFromTetMesh(qeal * currentTetPoints)
{
	if (!isTetMeshValid() || !isMedialMeshValid() || !isBindedTetMesh())
		return;
	for (int i = 0; i < bindedTM.size(); i++)
	{
		int nid = bindedTM[i];
		setMedialPoint(i, currentTetPoints + 3 * nid);
	}

}

void BaseModel::translateModel(const qeal x, const qeal y, const qeal z)
{	
	BaseSurfaceMesh::translateMesh(x - _tx, y - _ty, z - _tz);
	if (isTetMeshValid()) BaseTetMesh::translateMesh(x - _tx, y - _ty, z - _tz);
	if (isMedialMeshValid()) BaseMedialMesh::translateMesh(x - _tx, y - _ty, z - _tz);

	_tx = x;
	_ty = y;
	_tz = z;
}

void BaseModel::scaleModel(const qeal sx, const qeal sy, const qeal sz)
{
	qeal cx, cy, cz;
	getCenter(cx, cy, cz);
	BaseSurfaceMesh::scaleMesh(sx /_sx, sy/_sy, sz /_sz, cx, cy, cz);
	if (isTetMeshValid()) BaseTetMesh::scaleMesh(sx / _sx, sy / _sy, sz / _sz, cx, cy, cz);
	if (isMedialMeshValid()) BaseMedialMesh::scaleMesh(sx / _sx, sy / _sy, sz / _sz, cx, cy, cz);

	_sx = sx;
	_sy = sy;
	_sz = sz;
}

static Eigen::Matrix3d rotationMatrix(const qeal rx, const qeal ry, const qeal rz)
{
	if (Check_QEAL_ZERO(rx) && Check_QEAL_ZERO(ry) && Check_QEAL_ZERO(rz))
		return Eigen::Matrix3d::Identity();
	Eigen::Vector3d rotation(rx, ry, rz);
	rotation *= M_PI / 180.0;
	return (Eigen::AngleAxisd(rotation.z(), Eigen::Vector3d::UnitZ())
		* Eigen::AngleAxisd(rotation.y(), Eigen::Vector3d::UnitY())
		* Eigen::AngleAxisd(rotation.x(), Eigen::Vector3d::UnitX()))
		.toRotationMatrix();
}

void BaseModel::rotateModel(const qeal rx, const qeal ry, const qeal rz, const qeal sita, bool oldClear)
{
	oldClear = false;
	BaseFrame oldFrame; BaseFrame newFrame;
	Eigen::Matrix3d r_mat;
	Eigen::Matrix4d t_mat;

	t_mat.setIdentity();
	t_mat.block(0, 0, 3, 3) = rotationMatrix(rx, ry, rz);
	newFrame.setFromMatrix(t_mat.transpose().data());
	qglviewer::Quaternion rm = newFrame.rotation();
	t_mat.setIdentity();
	t_mat.block(0, 0, 3, 3) = rotationMatrix(_rx, _ry, _rz);
	oldFrame.setFromMatrix(t_mat.data());
	qglviewer::Quaternion old_rm = oldFrame.rotation();
	qeal cx, cy, cz;
	getCenter(cx, cy, cz);
	BaseSurfaceMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (isTetMeshValid()) BaseTetMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (isMedialMeshValid()) BaseMedialMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (!oldClear)
	{
		_rx = rx; _ry = ry; _rz = rz; _rsita = sita;
	}
}

