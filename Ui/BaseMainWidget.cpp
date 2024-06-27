 #include "BaseMainWidget.h"
#include <process.h>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QFile>
#include <QJsonObject>
#include <QDebug>
#include <QJsonArray>
#include "Simulator\SimulatorFactor.h"

void BaseMainWidget::bindWidget()
{
	bindTopToolBar();
	bindRightAndMainWidget();
	bindBottomAndScene();
	bindRightAndBottomWidget();
}

void BaseMainWidget::bindTopToolBar()
{
	connect(_topToolBar->openFileAction, SIGNAL(triggered()), this, SLOT(handleOpenFileAction()));
	connect(_topToolBar->saveFileAction, SIGNAL(triggered()), this, SLOT(handleSaveFileAction()));
	connect(_topToolBar->tetGenAction, SIGNAL(triggered()), this, SLOT(handleTetgenAction()));
	connect(_topToolBar->renderLineModeAction, SIGNAL(triggered()), this, SLOT(handleRenderLineModeAction()));
	connect(_topToolBar->selectSurfacePointsIAction, SIGNAL(triggered()), this, SLOT(handleSelectSurfacePointsAction()));
	connect(_topToolBar->selectTetPointsAction, SIGNAL(triggered()), this, SLOT(handleSelectTetPointsAction()));
}

void BaseMainWidget::bindRightAndMainWidget()
{
	connect(this, SIGNAL(addFileSuccess(bool)), _right, SLOT(refreshModelList()));
}

void BaseMainWidget::bindBottomAndScene()
{
	connect(_bottom, SIGNAL(recordAnimation(bool)), this, SLOT(handleRecorAnimationBtn(bool)));
	connect(_bottom, SIGNAL(playAnimation(bool)), this, SLOT(handlePlayAnimationBtn(bool)));
	connect(_bottom, SIGNAL(resetAnimation()), this, SLOT(handleResetAnimationBtn()));
	connect(_bottom, SIGNAL(changeViewerNum(int)), this, SLOT(handleMultiViewerButtonGroup(int)));
	connect(_scene, SIGNAL(needInitSim()), _bottom, SLOT(resetInitSimulatorStatus()));
}

void BaseMainWidget::bindRightAndBottomWidget()
{
	connect(_right, SIGNAL(needInitSimulator()), _bottom, SLOT(resetInitSimulatorStatus()));
}

void BaseMainWidget::handleOpenFileAction()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a Surface to import"), "./example/", tr(""));
	if (filename.isEmpty())
		return;
	std::string dir, name, format;
	getFilenameInfo(filename.toStdString(), dir, name, format);
	if (format == std::string("xml"))
	{
		bool flag = importFromConfigFile(filename.toStdString());
		Q_EMIT addFileSuccess(flag);
	}
	else
	{
		// read binary file

		/*if (_scene == nullptr)
			return;
		std::vector<std::string> list;
		list.push_back(filename.toStdString());
		bool flag = _scene->getCurrentSimulator()->addModelFromFile(list);
		Q_EMIT addFileSuccess(flag);*/
	}
}

void BaseMainWidget::handleSaveFileAction()
{
	// debug
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a JSON to import"), "./example/", tr(""));

	std::vector<Vector3> nodes;
	std::vector<Vector4i> elements;

	int type = 0;

	std::ifstream fin(filename.toStdString().c_str());
	while (!fin.eof())
	{
		std::string line;
		getline(fin, line);
		if (line == "$Nodes")
		{
			std::string node_num_str;
			int nodesNum;
			getline(fin, node_num_str);
			stringstream ss;
			ss << node_num_str;
			ss >> nodesNum;
			nodes.resize(nodesNum);
			for (int i = 0; i < nodesNum; i++)
			{
				std::string node_pos_str;
				getline(fin, node_pos_str);
				int index;
				qeal x, y, z;
				stringstream sss;
				sss << node_pos_str;
				sss >> index >> x >> y >> z;
				nodes[index - 1] = Vector3(x, y, z);
			}
		}
		if (line == "$Elements")
		{
			std::string ele_num_str;
			int eleNum;
			getline(fin, ele_num_str);
			stringstream ss;
			ss << ele_num_str;
			ss >> eleNum;
			elements.resize(eleNum);
			for (int i = 0; i < eleNum; i++)
			{
				std::string node_pos_str;
				getline(fin, node_pos_str);
				int index, index0, index1;
				int x, y, z, w;
				stringstream sss;
				sss << node_pos_str;
				sss >> index >> index0 >> index1 >> x >> y >> z >> w;
				elements[index - 1] = Vector4i(x - 1, y - 1, z - 1, w - 1);
			}
		}
	}
	fin.close();

	std::ofstream fout("result.ele");
	fout << elements.size() << " 4 0" << std::endl;
	for (int i = 0; i < elements.size(); i++)
	{
		fout << i << " " << elements[i][0] << " " << elements[i][1] << " " << elements[i][2] << " " << elements[i][3]  << std::endl;
	}
	fout.close();

	Vector3 cent = Vector3::Zero();
	for (int i = 0; i < nodes.size(); i++)
		cent += nodes[i];
	cent /= nodes.size();
	//for (int i = 0; i < nodes.size(); i++)
	//{
	//	Vector3 p = nodes[i] - cent;
	//	p *= 1.01;
	//	p += cent;
	//	nodes[i] = p;
	//}

	fout.open("result.node");
	fout << nodes.size() << " 3 0 0" << std::endl;
	for (int i = 0; i < nodes.size(); i++)
	{
		fout << i << " " << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2]  << std::endl;
	}
	fout.close();
}

void BaseMainWidget::handleTetgenAction()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a obj, off, node or batch txt"), "./example/", tr(""));
	std::string xdir, xname, xformat;
	getFilenameInfo(filename.toStdString(), xdir, xname, xformat);

	std::ifstream fin(filename.toStdString().c_str());
	for (int i = 0; i <= 500 && !fin.eof(); i++)
	{
		int num;
		fin >> num;
		if (num != 77940)
		{
			std::cout << "num " << num << std::endl;
			break;
		}
		int n0 = 75554;
		int n1 = 2386;
		std::string filename0 = xdir + xname + +"_0" + ".txt";
		std::string filename1 = xdir + xname + +"_1" + ".txt";

		std::ofstream fout;
		fout.open(filename0.c_str(), std::ios::app);
		fout << n0 << std::endl;
		for (int k = 0; k < n0; k++)
		{
			qeal x, y, z;
			fin >> x >> y >> z;
			fout << x << " " << y<< " " << z << std::endl;
		}

		fout.close();
		fout.open(filename1.c_str(), std::ios::app);
		fout << n1 << std::endl;
		for (int k = 0; k < n1; k++)
		{
			qeal x, y, z;
			fin >> x >> y >> z;
			fout << x << " " << y << " " << z << std::endl;
		}

		fout.close();
	}

	fin.close();

	//std::ifstream fin(filename.toStdString().c_str());
	//std::vector<std::string> strs;
	//int id = 0;
	//while (!fin.eof())
	//{
	//	std::string line;
	//	getline(fin, line);
	//	if (!(line.c_str()[0] == 'v' || line.c_str()[0] == 'f'))
	//	{
	//		if (line.c_str()[0] == 'g')
	//			id = 0;
	//		if (line.c_str()[0] == 'u')
	//		{
	//			QString n;
	//			n.setNum(id);
	//			line += " " + n.toStdString();
	//		}
	//		strs.push_back(line);
	//	}
	//	id++;
	//}
	//fin.close();
	//for (int i = 0; i < strs.size(); i++)
	//	std::cout << strs[i] << std::endl;

	//std::ofstream fout(xdir + "sd.txt");
	//for (int i = 0; i < strs.size(); i++)
	//{
	//	fout << strs[i] << std::endl;
	//}
	//fout.close();
	//system("pause");

	return;
}

void BaseMainWidget::handleRenderLineModeAction()
{
	if (_scene == nullptr)
		return;
	_scene->enableRenderLineMode(_topToolBar->renderLineModeAction->isChecked());
}

void BaseMainWidget::handleSelectSurfacePointsAction()
{
	if (_scene == nullptr)
		return;
	_scene->enableSelectSurfacePointsMode(_topToolBar->selectSurfacePointsIAction->isChecked());
	if (_topToolBar->selectSurfacePointsIAction->isChecked())
		_topToolBar->selectTetPointsAction->setChecked(false);
}

void BaseMainWidget::handleSelectTetPointsAction()
{
	if (_scene == nullptr)
		return;
	if (_scene == nullptr)
		return;
	_scene->enableSelectTetNodesMode(_topToolBar->selectTetPointsAction->isChecked());
	if (_topToolBar->selectTetPointsAction->isChecked())
		_topToolBar->selectSurfacePointsIAction->setChecked(false);
}

void BaseMainWidget::handlePlayAnimationBtn(bool flag)
{
	if (flag)
	{
		_scene->startAnimation();
	}		
	else
	{
		_scene->stopAnimation();
	}
}

void BaseMainWidget::handleRecorAnimationBtn(bool flag)
{
	_scene->getSceneStaus()->recordAnimation = flag;
}

void BaseMainWidget::handleResetAnimationBtn()
{
	_scene->getCurrentSimulator()->reset();
}

void BaseMainWidget::handleMultiViewerButtonGroup(int id)
{
	// to do
}

bool BaseMainWidget::importFromConfigFile(const std::string filename)
{
	TiXmlDocument doc(filename.c_str());
	doc.LoadFile();
	if (doc.Error() && doc.ErrorId() == TiXmlBase::TIXML_ERROR_OPENING_FILE) {
		std::cout << "Error: can't read config file !" << std::endl;
		return false;
	}
	TiXmlElement* headerItem = doc.FirstChildElement();
	if (!headerItem)
		return false;
	std::string checkItemName = headerItem->Value();
	if (checkItemName != std::string("SimFramework"))
		return false;

	// readSimulator
	if (!readSimulatorConfigFile(headerItem)) return false;

	// read OtherInfo
	TiXmlElement* subItem = headerItem->FirstChildElement();

	while (subItem)
	{
		std::string subItemName = subItem->Value();
		if (subItemName == std::string("simulator"))
		{
			if (!readSimulatorFromConfigFile(getCenteralScene()->getCurrentSimulator(), filename, subItem))
			{
				std::cout << "Error: can't read config file !!!" << std::endl;
				return false;
			}
		}
		else if (subItemName == std::string("scene"))
		{		
			getCenteralScene()->setFromConfigFile(filename, subItem);
		}
		subItem = subItem->NextSiblingElement();
	}

	return true;
}

bool BaseMainWidget::readSimulatorConfigFile(TiXmlElement* subItem)
{
	TiXmlAttribute * attri = subItem->FirstAttribute();
	std::string simName;
	std::string runPlaformName;
	while (attri)
	{
		std::string attriName = attri->Name();
		if (attriName == std::string("name"))
		{
			simName = attri->Value();
		}
		else if (attriName == std::string("run"))
		{
			runPlaformName = attri->Value();
		}
		attri = attri->Next();
	}

	RunPlatform run = RunPlatform::CPU;
	if (runPlaformName == "OPEMMP")
		run = RunPlatform::OPEMMP;
	else if (runPlaformName == "CUDA")
		run = RunPlatform::CUDA;

	getCenteralScene()->bindSimulator(createSimulator(simName, run));
	return true;
}
