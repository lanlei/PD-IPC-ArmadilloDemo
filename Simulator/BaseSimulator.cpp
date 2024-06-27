#include "BaseSimulator.h"
#include <direct.h>

bool BaseSimulator::readSimulatorFromConfigFile(const std::string filename, TiXmlElement * item)
{
	std::string itemName = item->Value();
	if (itemName != std::string("simulator"))
		return false;

	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (_sceneDir == std::string())
		_sceneDir = dir;

	TiXmlElement* subItem = item->FirstChildElement();
	while (subItem)
	{
		std::string subItemName = subItem->Value();
		if (subItemName == std::string("object"))
		{
			bool flag = addModelFromConfigFile(filename, subItem);
			if (!flag)
			{
				std::cout << "Error: can't read a object !!!" << std::endl;
			}
		}
		else if (subItemName == std::string("static_object"))
		{
			bool flag = addStaticModelFromConfigFile(filename, subItem);
			if (!flag)
			{
				std::cout << "Error: can't read a static object !!!" << std::endl;
			}
		}
		else if (subItemName == std::string("decorative_object"))
		{
			bool flag = addDecorativeModelFromConfigFile(filename, subItem);
			if (!flag)
			{
				std::cout << "Error: can't read a static object !!!" << std::endl;
			}
		}
		else if (subItemName == std::string("timeStep"))
		{
			std::string text = subItem->GetText();
			std::strstream ss;
			ss << text;
			qeal ts;
			ss >> ts;
			setTimeStep(ts);
		}
		else if (subItemName == std::string("gravity"))
		{
			std::string text = subItem->GetText();
			std::strstream ss;
			ss << text;
			qeal g;
			ss >> g;
			setGraviry(g);
		}
		else readExtraAttributeFromConfigFile(subItem);
		subItem = subItem->NextSiblingElement();
	}
	return true;
}

void BaseSimulator::readExtraAttributeFromConfigFile(TiXmlElement * item)
{
	if (!item)
		return;
}

bool BaseSimulator::addModelFromFile(std::vector<std::string> files, BaseModel* model)
{
	if (model == NULL)
		return false;

	std::string tetNodeFilename;
	std::string tetElementFilename;
	std::string medialMeshFilename;
	bool setUniform = true;
	bool setOrigin = true;
	for (int i = 0; i < files.size(); i++)
	{
		std::string dir, name, format;
		getFilenameInfo(files[i], dir, name, format);
		if (format == std::string("obj"))
		{
			model->setValid(model->readMeshFromObjFormat(files[i], this));
		}
		else if (format == std::string("node") || format == std::string("ele"))
		{

			if (!model->readMeshFromTetFormat(files[i], this))
			{
				std::cout << "Warning: can't read tet file ! " << std::endl;
			}
			if(format == std::string("node")) tetNodeFilename = files[i];
			if (format == std::string("ele")) tetElementFilename = files[i];
		}
		else if (format == std::string("mat") || format == std::string("ma"))
		{
			medialMeshFilename = files[i];
			if (!model->readMeshFromMatFormat(files[i], this))
			{
				std::cout << "Warning: can't read mat file ! " << std::endl;
			}
		}
		else if (format == std::string("nick_name"))
		{
			model->nickName = name;
		}
		else if (format == std::string("uniform"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setUniform = false;
		}
		else if (format == std::string("origin"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setOrigin = false;
		}
	}

	for (int i = 0; i < models.size(); i++)
	{
		models[i]->refreshBuffer(this, this, this);
	}

	if (model->isValid())
	{
		model->init(setUniform, setOrigin);
		model->initMeshesHandel();
		models.push_back(model);
	}

	return model->isValid();
}

bool BaseSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item)
{
	BaseModel* m = new BaseModel();
	return BaseSimulator::addModelFromConfigFile(filename, item, m);
}

bool BaseSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item, BaseModel* model)
{
	std::vector<std::string> files;
	getModelFilename(filename, item, files);
	return addModelFromFile(files, model);
}

bool BaseSimulator::getModelFilename(const std::string filename, TiXmlElement * item, std::vector<std::string>& files)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	TiXmlAttribute* attri = item->FirstAttribute();
	files.clear();
	std::string objName;
	if (attri->Name() == std::string("name"))
	{
		objName = attri->Value();
		std::string objFilename = dir + objName + "/" + objName + ".obj";
		files.push_back(objFilename);
	}
	else return false;

	attri = attri->Next();
	while (attri)
	{
		if (attri->Name() == std::string("volumetric_mesh") && std::string(attri->Value()) == std::string("tet"))
		{
			std::string tetNodeFilename = dir + objName + "/" + objName + ".node";
			std::string tetElementFilename = dir + objName + "/" + objName + ".ele";
			files.push_back(tetNodeFilename);
			files.push_back(tetElementFilename);
		}
		if (attri->Name() == std::string("medial_mesh") && (std::string(attri->Value()) == std::string("mat") || std::string(attri->Value()) == std::string("ma")))
		{
			std::string medialMeshFilename = dir + objName + "/" + objName + "." + std::string(attri->Value());
			files.push_back(medialMeshFilename);
		}
		if (attri->Name() == std::string("nick_name"))
		{
			std::string nickName = std::string(attri->Value()) + ".nick_name";
			files.push_back(nickName);
		}
		if (attri->Name() == std::string("uniform"))
		{
			std::string setUniform = std::string(attri->Value()) + ".uniform";
			files.push_back(setUniform);
		}
		if (attri->Name() == std::string("origin"))
		{
			std::string setOrigin = std::string(attri->Value()) + ".origin";
			files.push_back(setOrigin);
		}
		attri = attri->Next();
	}
	return true;
}

bool BaseSimulator::addStaticModelFromFile(std::vector<std::string> files, BaseModel * model)
{
	if (model == NULL)
		return false;

	std::string tetNodeFilename;
	std::string tetElementFilename;
	std::string medialMeshFilename;
	bool setUniform = true;
	bool setOrigin = true;
	for (int i = 0; i < files.size(); i++)
	{
		std::string dir, name, format;
		getFilenameInfo(files[i], dir, name, format);
		if (format == std::string("obj"))
		{
			model->setValid(model->readMeshFromObjFormat(files[i], &staticModelPool));
		}
		else if (format == std::string("node") || format == std::string("ele"))
		{

			if (!model->readMeshFromTetFormat(files[i], &staticModelPool))
			{
				std::cout << "Warning: can't read tet file ! " << std::endl;
			}
			if (format == std::string("node")) tetNodeFilename = files[i];
			if (format == std::string("ele")) tetElementFilename = files[i];
		}
		else if (format == std::string("mat") || format == std::string("ma"))
		{
			medialMeshFilename = files[i];
			if (!model->readMeshFromMatFormat(files[i], &staticModelPool))
			{
				std::cout << "Warning: can't read mat file ! " << std::endl;
			}
		}
		else if (format == std::string("nick_name"))
		{
			model->nickName = name;
		}
		else if (format == std::string("uniform"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setUniform = false;
		}
		else if (format == std::string("origin"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setOrigin = false;
		}
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		staticModels[i]->refreshBuffer(&staticModelPool, &staticModelPool, &staticModelPool);
	}
	if (model->isValid())
	{
		model->init(setUniform, setOrigin);
		staticModels.push_back(model);
	}

	return model->isValid();
}

bool BaseSimulator::addStaticModelFromConfigFile(const std::string filename, TiXmlElement * item)
{
	BaseModel* m = new BaseModel();
	return BaseSimulator::addStaticModelFromConfigFile(filename, item, m);
}

bool BaseSimulator::addStaticModelFromConfigFile(const std::string filename, TiXmlElement * item, BaseModel * model)
{
	std::vector<std::string> files;
	getModelFilename(filename, item, files);
	return addStaticModelFromFile(files, model);
}

bool BaseSimulator::addDecorativeModelFromFile(std::vector<std::string> files, BaseModel* model)
{
	if (model == NULL)
		return false;

	bool setUniform = true;
	bool setOrigin = true;
	for (int i = 0; i < files.size(); i++)
	{
		std::string dir, name, format;
		getFilenameInfo(files[i], dir, name, format);
		if (format == std::string("obj"))
		{
			model->setValid(model->readMeshFromObjFormat(files[i], &decorativeModelPool));
		}
		else if (format == std::string("nick_name"))
		{
			model->nickName = name;
		}
		else if (format == std::string("uniform"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setUniform = false;
		}
		else if (format == std::string("origin"))
		{
			if (name == std::string("off") || name == std::string("false") || name == std::string("0"))
				setOrigin = false;
		}
	}

	for (int i = 0; i < decorativeModels.size(); i++)
	{
		decorativeModels[i]->refreshBuffer(&decorativeModelPool, &decorativeModelPool, &decorativeModelPool);
	}
	if (model->isValid())
	{
		model->init(setUniform, setOrigin);
		decorativeModels.push_back(model);
	}

	return model->isValid();
}

bool BaseSimulator::addDecorativeModelFromConfigFile(const std::string filename, TiXmlElement* item)
{
	BaseModel* m = new BaseModel();
	bool isReadMesh = BaseSimulator::addDecorativeModelFromConfigFile(filename, item, m);

	if (!isReadMesh)
		return false;
	m = decorativeModels[decorativeModels.size() - 1];

	int hide = 0;
	qeal sx = 1.0, sy = 1.0, sz = 1.0;
	qeal tx = 0.0, ty = 0.0, tz = 0.0;
	qeal rx = 0.0, ry = 0.0, rz = 0.0, rsita = 0.0;
	TiXmlElement* childItem = item->FirstChildElement();
	std::strstream ss;

	while (childItem)
	{
		ss.clear();
		std::string itemName = childItem->Value();
		if (itemName == std::string("scale"))
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
		else if (itemName == std::string("hide"))
		{
			std::string str = childItem->GetText();
			ss << str;
			ss >> hide;
		}
		childItem = childItem->NextSiblingElement();
	}

	m->scaleModel(sx, sy, sz);
	m->translateModel(tx, ty, tz);
	m->rotateModel(rx, ry, rz, rsita);
	m->initMeshesHandel();
	m->getSurfaceHandle()->getMesh()->enableHide(hide);
	//	m->initialize();
		//#00a2e9
	QColor color("#00a2e9");
	m->setSurfaceNormalColor(color);

	return isReadMesh;
}

bool BaseSimulator::addDecorativeModelFromConfigFile(const std::string filename, TiXmlElement* item, BaseModel* model)
{
	std::vector<std::string> files;
	getModelFilename(filename, item, files);
	return addDecorativeModelFromFile(files, model);
}


void BaseSimulator::saveFile()
{
	//to do 
}

void BaseSimulator::saveSimulator(int frame)
{	
	QString snapshotFileName = QString(_sceneDir.c_str()) + QString("output/");
	int code = mkdir(snapshotFileName.toStdString().c_str());

	QString filename = snapshotFileName;
	if (frame < 10)
		filename = snapshotFileName + QString("//frame_000") + QString().setNum(frame) + QString(".obj");
	else if(frame < 100)
		filename = snapshotFileName + QString("//frame_00") + QString().setNum(frame) + QString(".obj");
	else if (frame < 1000)
		filename = snapshotFileName + QString("//frame_0") + QString().setNum(frame) + QString(".obj");
	else 
		filename = snapshotFileName + QString("//frame_") + QString().setNum(frame) + QString(".obj");

	std::ofstream fout(filename.toStdString().c_str());
	int offset = 1;
	for (int i = 0; i < models.size(); i++)
	{
		fout << "g " << models[i]->nickName << std::endl;
		for (int j = 0; j < models[i]->pointsNum; j++)
		{
			Vector3 p = models[i]->getSurfacePoint(j);
			fout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		}
		for (int j = 0; j < models[i]->facesNum; j++)
		{
			Vector3i f = models[i]->getSurfaceFace(j);
			fout << "f " << f[0] + offset << " " << f[1] + offset << " " << f[2] + offset << std::endl;
		}
		offset += models[i]->pointsNum;
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		for (int j = 0; j < staticModels[i]->pointsNum; j++)
		{
			Vector3 p = staticModels[i]->getSurfacePoint(j);
			fout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		}
		for (int j = 0; j < staticModels[i]->facesNum; j++)
		{
			Vector3i f = staticModels[i]->getSurfaceFace(j);
			fout << "f " << f[0] + offset << " " << f[1] + offset << " " << f[2] + offset << std::endl;
		}
		offset += staticModels[i]->pointsNum;
	}
	fout.close();
}

void BaseSimulator::render(QOpenGLShaderProgram * program, QOpenGLFunctions * f, bool drawEdge)
{
	for (int i = 0; i < models.size(); i++)
	{
		models[i]->render(program, f, drawEdge);
	}

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int i = 0; i < staticModels.size(); i++)
	{
		staticModels[i]->render(program, f, drawEdge);
	}
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < decorativeModels.size(); i++)
	{
		decorativeModels[i]->render(program, f, drawEdge);
	}
}

void BaseSimulator::renderExtraElementOnCPU(QOpenGLFunctions * f)
{	
	for (int i = 0; i < models.size(); i++)
	{
		if (!((BaseTetMesh*)models[i])->isHide() && ((BaseTetMesh*)models[i])->isTetMeshValid())
			models[i]->getTetMeshHandle()->renderTetMesh(f);

		if (!((BaseMedialMesh*)models[i])->isHide() && ((BaseMedialMesh*)models[i])->isMedialMeshValid())
			models[i]->getMedialMeshHandle()->renderMedialMesh(f);
	}

	for (int i = 0; i < staticModels.size(); i++)
	{
		//if (!((BaseTetMesh*)staticModels[i])->isHide() && ((BaseTetMesh*)staticModels[i])->isTetMeshValid())
		//	staticModels[i]->getTetMeshHandle()->renderTetMesh(f);	
		if (!((BaseMedialMesh*)staticModels[i])->isHide() && ((BaseMedialMesh*)staticModels[i])->isMedialMeshValid())
			staticModels[i]->getMedialMeshHandle()->renderMedialMesh(f);
	}
}

void BaseSimulator::handleTetPointsConstraint()
{
}

SparseMatrix BaseSimulator::computeMassMatrix()
{
	SparseMatrix M = SparseMatrix(3 * totalTetPointsNum, 3 * totalTetPointsNum);
	std::vector<TripletX> triplet;
	int index = 0;
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetPointsNum; j++)
		{
			qeal density = models[i]->getTetMeshHandle()->getNodeMaterial(j)->getDensity();
			qeal volume = models[i]->getTetMeshHandle()->getTetNodeParam(j)->volume;

			qeal v = density * volume;

			triplet.push_back(TripletX(3 * index, 3 * index, v));
			triplet.push_back(TripletX(3 * index + 1, 3 * index + 1, v));
			triplet.push_back(TripletX(3 * index + 2, 3 * index + 2, v));
			index++;
		}
	}
	M.setFromTriplets(triplet.begin(), triplet.end());
	return M;
}

SparseMatrix BaseSimulator::computeInverseMassMatrix()
{
	SparseMatrix M = SparseMatrix(3 * totalTetPointsNum, 3 * totalTetPointsNum);
	std::vector<TripletX> triplet;
	int index = 0;
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetPointsNum; j++)
		{
			qeal density = models[i]->getTetMeshHandle()->getNodeMaterial(j)->getDensity();
			qeal volume = models[i]->getTetMeshHandle()->getTetNodeParam(j)->volume;

			qeal v = 1.0 / (density * volume);

			triplet.push_back(TripletX(3 * index, 3 * index, v));
			triplet.push_back(TripletX(3 * index + 1, 3 * index + 1, v));
			triplet.push_back(TripletX(3 * index + 2, 3 * index + 2, v));
			index++;
		}
	}
	M.setFromTriplets(triplet.begin(), triplet.end());
	return M;
}

void BaseSimulator::alignAllMesh(qeal* tetPointsPtr)
{
	for (size_t mid = 0; mid < models.size(); mid++)
		models[mid]->alignAllMesh(tetPointsPtr);
}
