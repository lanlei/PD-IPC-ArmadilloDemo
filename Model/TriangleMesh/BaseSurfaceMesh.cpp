#include "BaseSurfaceMesh.h"
#include "Commom\FileIO.h"
#include "Commom\DataConversion.h"
#include "Model\tiny_obj_loader.h"
#include <stdlib.h>


void BaseSurfaceMesh::setNameAndDir(const std::string filename)
{
	getFilenameInfo(filename, dir, name, format);
	if (nickName == std::string(""))
		nickName = name;
}

bool BaseSurfaceMesh::readMeshFromObjFormat(const std::string filename, BaseSurfaceMeshBufferPool* pool)
{
	std::ifstream fin(filename.c_str());
	if (!fin.is_open())
		return false;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> groupMaterials;

	std::string err;
	std::string base_dir = getPathDir(filename);
	if (base_dir.empty())
	{
		base_dir = ".";
	}

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &groupMaterials, &err, filename.c_str(), base_dir.c_str(), true);
	if (!err.empty())
	{
		std::cerr << err << std::endl;
	}
	if (!ret) {
		std::cerr << "Error: Failed to load " << filename << " !" << std::endl;
		fin.close();
		return false;
	}

	hasNormals = false;
	hasTextureCoords = false;
	bool hasMaterials = false;

	if ((int)(attrib.normals.size()) / 3 > 0)
		hasNormals = true;
	if ((int)(attrib.texcoords.size()) / 3 > 0)
		hasTextureCoords = true;
	if ((int)groupMaterials.size() > 0)
		hasMaterials = true;

	int pointsNum = (int)(attrib.vertices.size()) / 3;
	int textureCoordNum = (int)(attrib.texcoords.size()) / 2;

	std::vector<qeal> points(attrib.vertices.size());
	for (int i = 0; i < attrib.vertices.size(); i++)
		points[i] = attrib.vertices[i];

	Vector3 cent = Vector3::Zero();
	for (int i = 0; i < points.size() / 3; i++)
	{
		cent[0] += points[3 * i];
		cent[1] += points[3 * i + 1];
		cent[2] += points[3 * i + 2];
	}
	cent /= pointsNum;

	std::vector<qeal> pointNormals(3 * pointsNum);

	std::vector<qeal> texCoords(attrib.texcoords.size());
	for (int i = 0; i < attrib.texcoords.size(); i++)
		texCoords[i] = attrib.texcoords[i];

	std::vector<std::vector<int>> pointFaceList(pointsNum);

	int renderGroupNum = shapes.size();
	renderMaterials.resize(renderGroupNum);
	for (int i = 0; i < renderGroupNum; i++)
	{
		renderMaterials[i] = new BaseRenderMaterial();
	}

	std::vector<int> groupFacesNum(renderGroupNum);

	std::vector<int> grounpMaterialIndices(renderGroupNum);

	std::vector<int> faceTexCoord;
	std::vector<int> faceIndices;

	for (int i = 0; i < renderGroupNum; i++)
	{
		groupFacesNum[i] = shapes[i].mesh.indices.size() / 3;

		for (size_t f = 0; f < shapes[i].mesh.indices.size() / 3; f++)
		{
			tinyobj::index_t idx0 = shapes[i].mesh.indices[3 * f + 0];
			tinyobj::index_t idx1 = shapes[i].mesh.indices[3 * f + 1];
			tinyobj::index_t idx2 = shapes[i].mesh.indices[3 * f + 2];

			int vid0 = idx0.vertex_index;
			int vid1 = idx1.vertex_index;
			int vid2 = idx2.vertex_index;

			faceIndices.push_back(vid0);
			faceIndices.push_back(vid1);
			faceIndices.push_back(vid2);

			int face_id = faceIndices.size() / 3 - 1;
			pointFaceList[vid0].push_back(face_id);
			pointFaceList[vid1].push_back(face_id);
			pointFaceList[vid2].push_back(face_id);

			int tid0 = idx0.texcoord_index;
			int tid1 = idx1.texcoord_index;
			int tid2 = idx2.texcoord_index;

			faceTexCoord.push_back(tid0);
			faceTexCoord.push_back(tid1);
			faceTexCoord.push_back(tid2);
		}
		if (hasMaterials)
		{
			int meterial_id = shapes[i].mesh.material_ids[0];
			tinyobj::material_t m = groupMaterials[meterial_id];
			renderMaterials[i]->ambient = getQColorRGB(QVector3D(m.ambient[0], m.ambient[1], m.ambient[2]));
			renderMaterials[i]->diffuse = getQColorRGB(QVector3D(m.diffuse[0], m.diffuse[1], m.diffuse[2]));
			renderMaterials[i]->specular = getQColorRGB(QVector3D(m.specular[0], m.specular[1], m.specular[2]));
			renderMaterials[i]->shinness = m.shininess;

			if (m.ambient_texname.length() > 0)
			{
				QString imgFile = QString(base_dir.c_str()) + QString(m.ambient_texname.c_str());
				renderMaterials[i]->readTextureMap(imgFile, AmbientMapIndex);
			}
			if (m.diffuse_texname.length() > 0)
			{
				QString imgFile = QString(base_dir.c_str()) + QString(m.diffuse_texname.c_str());
				renderMaterials[i]->readTextureMap(imgFile, DiffuseMapIndex);
			}
			if (m.specular_texname.length() > 0)
			{
				QString imgFile = QString(base_dir.c_str()) + QString(m.specular_texname.c_str());
				renderMaterials[i]->readTextureMap(imgFile, SpecularMapIndex);
			}
			if (m.bump_texname.length() > 0)
			{
				QString imgFile = QString(base_dir.c_str()) + QString(m.bump_texname.c_str());
				renderMaterials[i]->readTextureMap(imgFile, BumpMapIndex);
			}
		}
	}

	int facesNum = (int)faceIndices.size() / 3;

	std::vector<qeal> faceNormals(3 * facesNum);
	//alignBuffer
	std::vector<std::vector<int>> pointTexCoords(pointsNum);

	for (int i = 0; i < facesNum; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int vid = faceIndices[3 * i + j];
			int tid = faceTexCoord[3 * i + j];
			pointTexCoords[vid].push_back(tid);
		}
	}

	if (hasTextureCoords)
	{
		std::vector<qeal> tempCoords;
		for (int i = 0; i < pointsNum; i++)
		{
			int texId = pointTexCoords[i][0];
			tempCoords.push_back(texCoords[2 * texId]);
			tempCoords.push_back(texCoords[2 * texId + 1]);
			/*if (pointTexCoords[i].size() == 1)
			{
				int texId = pointTexCoords[i][0];
				tempCoords.push_back(texCoords[2 * i]);
				tempCoords.push_back(texCoords[2 * i + 1]);
			}else if (pointTexCoords[i].size() > 1)
			{
				int texId = pointTexCoords[i][0];
				tempCoords.push_back(texCoords[2 * i]);
				tempCoords.push_back(texCoords[2 * i + 1]);
			}*/
		}
		texCoords = tempCoords;
	}

	////////////////////////////////////

	this->pointsNum = pointsNum;
	this->facesNum = facesNum;
	this->renderGroupNum = renderGroupNum;

	pool->totalPointsNum += this->pointsNum;
	pool->totalFacesNum += this->facesNum;
	pool->totalRenderGroupNum += this->renderGroupNum;

	registerToOverallBuffer<qeal>(points, pool->pointsBuffer, this->points);
	registerToOverallBuffer<qeal>(texCoords, pool->texCoordsBuffer, this->texCoords);

	registerToOverallBuffer<qeal>(pointNormals, pool->pointsNormalBuffer, this->pointNormals);
	registerToOverallBuffer<qeal>(faceNormals, pool->facesNormalBuffer, this->faceNormals);
	registerToOverallBuffer<int>(faceIndices, pool->faceIndicesBuffer, this->faceIndices);
	registerToOverallBuffer<int>(pointFaceList, pool->pointFaceListBuffer, this->pointFaceIndices);
	registerToOverallBuffer<int>(groupFacesNum, pool->groupFacesNumBuffer, this->renderGroupFaceNum);
	//////////////////////////////////////
	setNameAndDir(filename);
	pool->meshNum++;

	pointsOffset = this->points.offset / 3;
	faceOffset = this->faceIndices.offset / 3;

	// edge Info
	edgesNum = 0;
	std::vector<Vector2i> edgeList;
	for (int j = 0; j < facesNum; j++)
	{
		Vector3i face = Vector3i(faceIndices[3 * j], faceIndices[3 * j + 1], faceIndices[3 * j + 2]);
		for (int k = 0; k < 3; k++)
		{
			Vector2i e;
			if (face[k % 3] < face[(k + 1) % 3]) e = Vector2i(face[k % 3], face[(k + 1) % 3]);
			else e = Vector2i(face[(k + 1) % 3], face[k % 3]);
			edgeList.push_back(e);
		}
	}
	std::sort(edgeList.begin(), edgeList.end(), v2IndComp);
	edgeList.erase(std::unique(edgeList.begin(), edgeList.end()), edgeList.end());
	std::vector<int> edgeIndices(2 * edgeList.size());
	
	for (int i = 0; i < edgeList.size(); i++)
	{
		edgeIndices[2 * i] = edgeList[i][0];
		edgeIndices[2 * i + 1] = edgeList[i][1];
	}
	
	std::vector<std::vector<int>> pointSharedEdges;
	std::vector<std::vector<int>> faceSharedEdges;
	std::vector<std::vector<int>> edgeSharedFaces;

	pointSharedEdges.resize(pointsNum);
	edgeSharedFaces.resize(edgeList.size());
	faceSharedEdges.resize(facesNum);
	std::map<Vector2i, int, Vec2IndComp> edgesIds;
	for (int i = 0; i < edgeList.size(); i++)
	{
		pointSharedEdges[edgeList[i][0]].push_back(i);
		pointSharedEdges[edgeList[i][1]].push_back(i);
		edgesIds[edgeList[i]] = i;
	}
	for (int i = 0; i < facesNum; i++)
	{
		Vector3i face = Vector3i(faceIndices[3 * i], faceIndices[3 * i + 1], faceIndices[3 * i + 2]);
		for (int j = 0; j < 3; j++)
		{
			Vector2i e;
			if (face[j % 3] < face[(j + 1) % 3]) e = Vector2i(face[j % 3], face[(j + 1) % 3]);
			else e = Vector2i(face[(j + 1) % 3], face[j % 3]);
			int eid = edgesIds[e];
			edgeSharedFaces[eid].push_back(i);
			faceSharedEdges[i].push_back(eid);
		}
	}

	this->edgesNum = edgeList.size();
	registerToOverallBuffer<int>(edgeIndices, pool->edgeIndicesBuffer, this->edgeIndices);
	registerToOverallBuffer<int>(pointSharedEdges, pool->pointEdgeListBuffer, this->pointEdgeIndices);
	registerToOverallBuffer<int>(edgeSharedFaces, pool->edgeFaceListBuffer, this->edgeFaceIndices);
	registerToOverallBuffer<int>(faceSharedEdges, pool->faceEdgeListBuffer, this->faceEdgeIndices);
	edgeOffset = this->edgeIndices.offset / 2;



	return true;
}

bool BaseSurfaceMesh::writeMeshToObjFormat(const std::string filename)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (format != std::string("obj"))
		return false;

	std::ofstream fout(filename.c_str());
	if (!fout.is_open())
		return false;
	
	for (size_t i = 0; i < pointsNum; i++)
		fout << "v " << points.buffer[3 * i] << " " << points.buffer[3 * i + 1] << " " << points.buffer[3 * i + 2] << std::endl;
	if (hasTextureCoords)
	{
		for (size_t i = 0; i < (texCoords.span) / 2; i++)
		{
			fout << "vt " << texCoords.buffer[2 * i] << " " << texCoords.buffer[2 * i + 1] << std::endl;
		}
	}

	int offset = 0;
	for (size_t i = 0; i < renderGroupNum; i++)
	{
		fout << "g group" << i << std::endl;
		int* ptr = faceIndices.buffer + offset;

		for (size_t f = 0; f < renderGroupFaceNum.buffer[i]; f++)
		{			
			fout << "f";
			int v0 = ptr[3 * f];
			fout << " " << v0 + 1;
		//	if(hasTextureCoords)
		//		fout <<"/"<< 
			int v1 = ptr[3 * f + 1];
			fout << " " << v1 + 1;
		//	if(hasTextureCoords)
		//		fout <<"/"<< 
			int v2 = ptr[3 * f + 2];
			fout << " " << v2 + 1;
		//	if(hasTextureCoords)
		//		fout <<"/"<< 
			fout << std::endl;
		}
		offset += renderGroupFaceNum.buffer[i] * 3;
	}

	return true;
}

void BaseSurfaceMesh::render(QOpenGLShaderProgram * program, QOpenGLFunctions* f, bool drawEdge)
{
	if (_hide)
		return;
	updateVBO();
	_vertexArrayBuf->bind();

	int vertexLocation;
	int textureLocation;
	vertexLocation = program->attributeLocation("a_position");
	textureLocation = program->attributeLocation("a_texcoord");
	quintptr offset = 0;
	program->enableAttributeArray(vertexLocation);
	program->setAttributeBuffer(vertexLocation, GL_QEAL, offset, 3, 3 * sizeof(qeal));

	if (hasTextureCoords)
	{
		offset += points.span * sizeof(qeal);
		program->enableAttributeArray(textureLocation);
		program->setAttributeBuffer(textureLocation, GL_QEAL, offset, 2, 2 * sizeof(qeal));
	}

	_vertexArrayBuf->release();

	if(drawEdge)
		program->setUniformValue("enableLineMode", 1);
	else program->setUniformValue("enableLineMode", 0);

	program->setUniformValue("useShadowMap", shadow);

	offset = 0;
	for (int i = 0; i < renderGroupNum; i++)
	{
		renderMaterials[i]->transferToShader(program);
		f->glDrawElements(GL_TRIANGLES, renderGroupFaceNum.buffer[i] * 3, GL_UNSIGNED_INT, faceIndices.buffer + offset);
		offset += renderGroupFaceNum.buffer[i] * 3;
	}
	program->disableAttributeArray(vertexLocation);
	program->disableAttributeArray(textureLocation);
}

void BaseSurfaceMesh::updateVBO()
{
	if (_type == STATICS_MESH)
		return;
	std::copy(points.buffer, points.buffer + points.span, _renderVertexBuffer.begin());
	_vertexArrayBuf->bind();
	_vertexArrayBuf->allocate(_renderVertexBuffer.data(), _renderVertexBuffer.size() * sizeof(qeal));
	_vertexArrayBuf->release();
}

qeal BaseSurfaceMesh::uniform()
{
	qeal div = bbox[3] - bbox[0];
	if (div < (bbox[4] - bbox[1]))
		div = (bbox[4] - bbox[1]);
	if (div < (bbox[5] - bbox[2]))
		div = (bbox[5] - bbox[2]);
	if (fabs(div) < 1e-6) div = 1.0;

	for (int i = 0; i < pointsNum; i++)
	{
		points.buffer[3 * i] /= div;
		points.buffer[3 * i + 1] /= div;
		points.buffer[3 * i + 2] /= div;
	}


	return div;
}

void BaseSurfaceMesh::computeBBox()
{
	bbox[0] = QEAL_MAX;
	bbox[1] = QEAL_MAX;
	bbox[2] = QEAL_MAX;
	bbox[3] = -QEAL_MAX;
	bbox[4] = -QEAL_MAX;
	bbox[5] = -QEAL_MAX;
	for (int i = 0; i < pointsNum; i++)
	{
		qeal x = points.buffer[3 * i];
		qeal y = points.buffer[3 * i + 1];
		qeal z = points.buffer[3 * i + 2];
		if (x < bbox[0])
			bbox[0] = x;
		else if (x > bbox[3])
			bbox[3] = x;

		if (y < bbox[1])
			bbox[1] = y;
		else if (y > bbox[4])
			bbox[4] = y;

		if (z < bbox[2])
			bbox[2] = z;
		else if (z > bbox[5])
			bbox[5] = z;
	}
}

void BaseSurfaceMesh::getCenter(qeal & cx, qeal & cy, qeal & cz)
{
	computeBBox();
	cx = (bbox[0] + bbox[3]) / 2.0;
	cy = (bbox[1] + bbox[4]) / 2.0;
	cz = (bbox[2] + bbox[5]) / 2.0;
}

void BaseSurfaceMesh::placeOrigin(const qeal lx, const qeal ly, const qeal lz)
{
	for (int i = 0; i < pointsNum; i++)
	{
		points.buffer[3 * i] -= lx;
		points.buffer[3 * i + 1] -= ly;
		points.buffer[3 * i + 2] -= lz;
	}
}

void BaseSurfaceMesh::translateMesh(const qeal x, const qeal y, const qeal z)
{
	for (int i = 0; i < pointsNum; i++)
	{
		points.buffer[3 * i] += x;
		points.buffer[3 * i + 1] += y;
		points.buffer[3 * i + 2] += z;
	}
	updateVBO();
}

void BaseSurfaceMesh::scaleMesh(const qeal sx, const qeal sy, const qeal sz, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < pointsNum; i++)
	{
		points.buffer[3 * i] = (points.buffer[3 * i] - cx) * sx + cx;
		points.buffer[3 * i + 1] = (points.buffer[3 * i + 1] - cy) * sy + cy;
		points.buffer[3 * i + 2] = (points.buffer[3 * i + 2] - cz) * sz + cz;
	}
	updateVBO();
}

void BaseSurfaceMesh::rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < pointsNum; i++)
	{
		qglviewer::Vec pos(points.buffer[3 * i] - cx, points.buffer[3 * i + 1] - cy, points.buffer[3 * i + 2] - cz);
		pos = oldR.inverseRotate(pos);
		pos = R.rotate(pos);

		points.buffer[3 * i] = pos.x + cx;
		points.buffer[3 * i + 1] = pos.y + cy;
		points.buffer[3 * i + 2] = pos.z + cz;
	}
	updateVBO();
}

void BaseSurfaceMesh::initVBO()
{
	_renderVertexBuffer.resize(points.span + texCoords.span);

	std::copy(points.buffer, points.buffer + points.span, _renderVertexBuffer.begin());
	if(hasTextureCoords)
		std::copy(texCoords.buffer, texCoords.buffer + texCoords.span, _renderVertexBuffer.begin() + points.span);
	_vertexArrayBuf->bind();
	_vertexArrayBuf->allocate(_renderVertexBuffer.data(), _renderVertexBuffer.size() * sizeof(qeal));
	_vertexArrayBuf->release();
}

void BaseSurfaceMeshBuffer::copySurfacePointsToBuffer(qeal* buffer, int size)
{
	if(size == 0)
		std::copy(points.buffer, points.buffer + points.span, buffer);	
	else std::copy(points.buffer, points.buffer + size, buffer);
}

void BaseSurfaceMeshBuffer::copyBufferToSurfacePoints(qeal* buffer, int size)
{
	if (size == 0)
		std::copy(buffer, buffer + points.span, points.buffer);
	else std::copy(buffer, buffer + size, points.buffer);
}

Vector3 BaseSurfaceMeshBuffer::getSurfacePoint(int nid)
{
	return Vector3(points.buffer[3 * nid], points.buffer[3 * nid + 1], points.buffer[3 * nid + 2]);
}

void BaseSurfaceMeshBuffer::setSurfacePoint(int nid, qeal* p)
{
	points.buffer[3 * nid] = p[0];
	points.buffer[3 * nid + 1] = p[1];
	points.buffer[3 * nid + 2] = p[2];
}

Vector3 BaseSurfaceMeshBuffer::getSurfacePointNormal(int nid)
{
	return Vector3(pointNormals.buffer[3 * nid], pointNormals.buffer[3 * nid + 1], pointNormals.buffer[3 * nid + 2]);
}

void BaseSurfaceMeshBuffer::setSurfacePointNormal(int nid, qeal* n)
{
	pointNormals.buffer[3 * nid] = n[0];
	pointNormals.buffer[3 * nid + 1] = n[1];
	pointNormals.buffer[3 * nid + 2] = n[2];
}

Vector3 BaseSurfaceMeshBuffer::getSurfaceFaceNormal(int nid)
{
	return Vector3(faceNormals.buffer[3 * nid], faceNormals.buffer[3 * nid + 1], faceNormals.buffer[3 * nid + 2]);
}

void BaseSurfaceMeshBuffer::setSurfaceFaceNormal(int nid, qeal* n)
{
	faceNormals.buffer[3 * nid] = n[0];
	faceNormals.buffer[3 * nid + 1] = n[1];
	faceNormals.buffer[3 * nid + 2] = n[2];
}

Vector3i BaseSurfaceMeshBuffer::getSurfaceFace(int fid)
{
	return Vector3i(faceIndices.buffer[3 * fid], faceIndices.buffer[3 * fid + 1], faceIndices.buffer[3 * fid + 2]);
}

Vector2i BaseSurfaceMeshBuffer::getSurfaceEdge(int eid)
{
	return Vector2i(edgeIndices.buffer[2 * eid], edgeIndices.buffer[2 * eid + 1]);
}

int BaseSurfaceMeshBuffer::getSurfacePointOverallId(int nid)
{
	return nid + pointsOffset;
}

int BaseSurfaceMeshBuffer::getSurfaceFaceOverallId(int fid)
{
	return fid + faceOffset;
}

int BaseSurfaceMeshBuffer::getSurfaceEdgeOverallId(int eid)
{
	return eid + edgeOffset;
}