#include "BaseScene.h"
#include<direct.h>
//#include "Simulator\ProjectiveDynamics\ProjectiveDynamicsSimulator.h"
//#include "Simulator\ReducedIpc\ReducedIPCSimulator.h"
#include "Simulator\PD-IPC\PdIpcSimulator.h"

BaseScene::BaseScene() :_sim(nullptr)
{
	QGLFormat glFormat;
	glFormat.setVersion(3, 2);
	glFormat.setProfile(QGLFormat::CoreProfile);

	QSurfaceFormat glSurfaceFormat;
	glSurfaceFormat.setSamples(32);
	setFormat(glSurfaceFormat);

	_selectedTetPointId = -1;
	_mouseForce.setZero();

}

BaseScene::BaseScene(QWidget* parent) :
	QGLViewer(parent)
{

}

BaseScene::~BaseScene()
{
	if (_camera == nullptr)
		delete _camera;
}

bool BaseScene::setFromConfigFile(const std::string filename, TiXmlElement* item)
{
	std::string itemName = item->Value();
	if (itemName != std::string("scene"))
		return false;

	TiXmlElement* subItem = item->FirstChildElement();
	std::strstream ss;
	while (subItem)
	{
		ss.clear();
		std::string subItemName = subItem->Value();
		if (subItemName == std::string("cameraPos"))
		{
			std::string str = subItem->GetText();
			ss << str;
			qeal x, y, z;
			ss >> x >> y >> z;
			camera()->setPosition(qglviewer::Vec(x, y, z));
		}
		else if (subItemName == std::string("cameraDir"))
		{
			std::string str = subItem->GetText();
			ss << str;
			qeal x, y, z;
			ss >> x >> y >> z;
			camera()->setViewDirection(qglviewer::Vec(x, y, z));
		}
		else if (subItemName == std::string("startFrame"))
		{
			std::string str = subItem->GetText();
			ss << str;
			ss >> _status.animationLifeTime;
		}
		else if (subItemName == std::string("maxFrame"))
		{
			std::string str = subItem->GetText();
			ss << str;
			ss >> _status.animationMaxLifeTime;
		}
		else if (subItemName == std::string("scale_floor"))
		{
			std::string str = subItem->GetText();
			ss << str;
			int scale = 1.0;
			ss >> scale;
			_floor->scaleXZ(scale);
		}
		else if (subItemName == std::string("translate_floor"))
		{
			std::string str = subItem->GetText();
			ss << str;
			qeal x = 0.0;
			qeal y = 0.0;
			qeal z = 0.0;
			ss >> x >> y >> z;
			_floor->setTranslation(x, y, z);
		}
		subItem = subItem->NextSiblingElement();
	}

	return true;
}

void BaseScene::bindSimulator(BaseSimulator* sim)
{
	if (_sim)
		free(_sim);
	_sim = sim;
	Q_EMIT changeSimulator();
}

void BaseScene::draw()
{
	renderDepthToBuffer();

	glViewport(0, 0, width(), height());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	setBackgroundColor(QColor(255, 255, 255));

	_phongProgram.bind();
	_phongProgram.setUniformValue("projectMatrix", ((StandardCamera*)camera())->getProjectionQMatrix());
	_phongProgram.setUniformValue("viewMatrix", ((StandardCamera*)camera())->getViewQMatrix());
	QVector3D lightPos = _lights.getMainLightPosition();

	_lightCamera.setPosition(qglviewer::Vec(lightPos[0], lightPos[1], lightPos[2]));
	_lightCamera.lookAt(qglviewer::Vec(0.0, 0.0, 0.0));

	_phongProgram.setUniformValue("lightSpaceMatrix", _lightCamera.getProjectionViewQMatrix());
	QVector3D viewPos = QVector3D(camera()->position().x, camera()->position().y, camera()->position().z);
	_phongProgram.setUniformValue("viewPos", viewPos);

	_lights.transferToShader(&_phongProgram);

	_phongProgram.setUniformValue("shadowMap", 0);
	_phongProgram.setUniformValue("ambientMap", 1);
	_phongProgram.setUniformValue("diffuseMap", 2);
	_phongProgram.setUniformValue("specularMap", 3);
	_phongProgram.setUniformValue("bumpMap", 4);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, _depthMapId);

	renderScene(&_phongProgram);

	glBindTexture(GL_TEXTURE_2D, 0);
	_phongProgram.release();

	renderScene();
	update();
}

void BaseScene::animate()
{
	if (!_sim)
		return;

	//	std::ofstream fout("timing.txt ", std::ios::app);
	_status.animationTimer.start();
	computeMouseForce();
	if (_status.saveAnimation)
		_sim->saveSimulator(_status.animationLifeTime);

	_sim->animate(_status.animationLifeTime);

	qeal time = _status.animationTimer.nsecsElapsed() / 1e6;

	int costPtr = _status.animationLifeTime % 8;
	_status.animationCost[costPtr] = time;


	if (_sim->getTimeStep() == 0.01 && _status.animationLifeTime % 3 == 0)
	{
		if (_status.recordAnimation)
			Q_EMIT recordAnimation(true);
	}
	else if (_sim->getTimeStep() == 0.005 && _status.animationLifeTime % 6 == 0)
	{
		if (_status.recordAnimation)
			Q_EMIT recordAnimation(true);
	}
	else if (_sim->getTimeStep() == 0.00125 && _status.animationLifeTime % 24 == 0)
	{
		if (_status.recordAnimation)
			Q_EMIT recordAnimation(true);
	}
	else if (_sim->getTimeStep() == 0.004 && _status.animationLifeTime % 10 == 0)
	{
		if (_status.recordAnimation)
			Q_EMIT recordAnimation(true);
	}
	else if (_status.animationLifeTime % 3 == 0)
	{
		if (_status.recordAnimation)
			Q_EMIT recordAnimation(true);
	}


	_status.animationLifeTime++;

	/*if (_status.animationLifeTime > _status.animationMaxLifeTime && _status.animationMaxLifeTime > 0)
	{
		int goOn = 0;
		printf("continue or not? ");
		scanf("%d", &goOn);
		printf("\n");
		if (goOn <= 0)
		{
			stopAnimation();
		}
		else
		{
			int maxFrame = goOn;
			if (maxFrame > _status.animationMaxLifeTime)_status.animationMaxLifeTime = maxFrame;
			else _status.animationMaxLifeTime += maxFrame;
		}
	}*/
}

void BaseScene::postDraw()
{
	renderSelectionRect();
	renderSelectionPoints();
	renderForceLine();
	renderCornerAxis();
	renderInfoOnScene();
	update();
}

void BaseScene::init()
{
	initializeOpenGLFunctions();
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	initShaderProgram();
	initCamera();
	initShadowFramebuffer();
	initAnimation();

	_floor = new Floor();
	//getFloor()->setYPlane(0.0);
	Q_EMIT initSceneSignal();

	select_tet_group = 1;

	//
}

void BaseScene::mousePressEvent(QMouseEvent* e)
{
	_rectangle = QRect(e->pos(), e->pos());
	if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier))
		_selectionMode = SHIFT_SELECT;
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::ShiftModifier))
		_selectionMode = SHIFT_REMOVE;
	else if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
		_selectionMode = ALT_SELECT;
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::AltModifier))
		_selectionMode = ALT_REMOVE;
	else if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ControlModifier))
	{
		_selectionMode = CTRL_SELECT;
		setSelectRegionWidth(80);
		setSelectRegionHeight(80);
		select(e->pos());
		_forceLine.setP1(e->pos());
	}
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::ControlModifier))
		_selectionMode = CTRL_REMOVE;
	else
		QGLViewer::mousePressEvent(e);
	update();
}

void BaseScene::mouseMoveEvent(QMouseEvent* e)
{
	if (_selectionMode != NONE_SELECT && _status.enableSelectSurfacePoints)
	{
		_rectangle.setBottomRight(e->pos());
		update();
	}
	else if (_selectionMode != NONE_SELECT && _status.enableSelectTetNodes)
	{
		_rectangle.setBottomRight(e->pos());
		update();
	}
	else if (_selectionMode == CTRL_SELECT)
	{
		_forceLine.setP2(e->pos());
	}
	else QGLViewer::mouseMoveEvent(e);
}

void BaseScene::mouseReleaseEvent(QMouseEvent* e)
{
	if (_selectionMode != NONE_SELECT && _status.enableSelectSurfacePoints)
	{
		// to do
		_rectangle = _rectangle.normalized();
		setSelectRegionWidth(_rectangle.width());
		setSelectRegionHeight(_rectangle.height());
		select(_rectangle.center());
		_selectionMode = NONE_SELECT;
		_rectangle = QRect();
		update();
	}
	else if (_selectionMode != NONE_SELECT && _status.enableSelectTetNodes)
	{
		_rectangle = _rectangle.normalized();
		setSelectRegionWidth(_rectangle.width());
		setSelectRegionHeight(_rectangle.height());
		select(_rectangle.center());
		_selectionMode = NONE_SELECT;
		_rectangle = QRect();
		update();
	}
	else if (_selectionMode == CTRL_SELECT)
	{
		_forceLine = QLine();
		_selectionMode = NONE_SELECT;
		_selectedTetPointId = -1;
		_mouseForce.setZero();
		selectedTetPointList.clear();
	}
	else QGLViewer::mouseReleaseEvent(e);
}

void BaseScene::keyPressEvent(QKeyEvent* e)
{
	if (e->key() == Qt::Key_D)
		startAnimation();
	else if (e->key() == Qt::Key_C)
	{
		std::cout << "-----Scen Info-----" << std::endl;
		std::cout << "animationFilePath: " << _status.animationFilePath << std::endl;
		std::cout << "animationLifeTime: " << _status.animationLifeTime << std::endl;
		std::cout << "animationMaxLifeTime: " << _status.animationMaxLifeTime << std::endl;
		std::cout << "cameraPos: " << camera()->position().x << " " << camera()->position().y << " " << camera()->position().z << std::endl;
		std::cout << "cameraDir: " << camera()->viewDirection().x << " " << camera()->viewDirection().y << " " << camera()->viewDirection().z << std::endl;
		std::cout << "-------------------" << std::endl;
	}
	else if (e->key() == Qt::Key_Z)
	{
		debugCount++;
	}
	else if (e->key() == Qt::Key_X)
	{
		debugCount--;
	}
	else if (e->key() == Qt::Key_H)
	{
		int num = _sim->selectSurfacePoints.size();
		std::ofstream fout("surface_points.txt");
		fout << "c " << num << std::endl;
		cout << "------------------------------" << std::endl;
		cout << "c " << num << std::endl;
		auto it = _sim->selectSurfacePoints.begin();
		for (; it != _sim->selectSurfacePoints.end(); ++it)
		{
			fout << *it << std::endl;
			cout << *it << std::endl;
		}

		fout.close();
	}
	else if (e->key() == Qt::Key_F)
	{
		int num = _sim->selectTetPoints.size();
		std::ofstream fout("cc.txt");
		fout << "c " << num << std::endl;
		cout << "------------------------------" << std::endl;
		cout << "c " << num << std::endl;
		auto it = _sim->selectTetPoints.begin();
		for (; it != _sim->selectTetPoints.end(); ++it)
		{
			fout << *it << std::endl;
			cout << *it << std::endl;
		}

		fout.close();
	}
}

void BaseScene::wheelEvent(QWheelEvent* e)
{
	QGLViewer::wheelEvent(e);
}

void BaseScene::drawWithNames()
{
	if (!_sim)
		return;
	if (_sim->getModelsNum() == 0)return;

	PD_IPC::PdIpcSimulator* sim = dynamic_cast<PD_IPC::PdIpcSimulator*>(_sim);

	for(int i = 0; i < sim->simPointsNum; i++)
	{
		glPushName(i);
		glBegin(GL_POINTS);
		qeal x, y, z;
		Vector3 p = sim->points.col(i);
		glVertex3f(p[0], p[1], p[2]);
		glEnd();
		glPopName();
	}
}

void BaseScene::endSelection(const QPoint& point)
{
	selectedTetPointList.clear();
	selectedTetPointListCent = Vector3(-100, -100, -100);
	update();
	if (_sim == nullptr)
		return;
	GLint nbHits = glRenderMode(GL_RENDER);
	if (nbHits <= 0)
		return;
	std::cout << "Select : "<< nbHits << std::endl;
	for (int i = 0; i < nbHits; ++i)
	{
		int id = (selectBuffer())[4 * i + 3];
		//std::cout << "tet " << id << std::endl;
		if (_status.enableSelectSurfacePoints)
		{
			if (id >= _sim->totalPointsNum) continue;
			if (_selectionMode == ALT_SELECT)
			{
				_sim->selectSurfacePoints.insert(id);
				//std::cout << id << std::endl;
			}
			else if (_selectionMode == ALT_REMOVE)
			{
				_sim->selectSurfacePoints.erase(id);
			}
		}
		else if (_status.enableSelectTetNodes)
		{
			if (id >= _sim->totalTetPointsNum) continue;
			if (_selectionMode == ALT_SELECT)
			{
				_sim->selectTetPoints.insert(id);
				//std::cout << id << std::endl;
			}
			else if (_selectionMode == ALT_REMOVE)
			{
				_sim->selectTetPoints.erase(id);
			}
		}
	}

	PD_IPC::PdIpcSimulator* sim = dynamic_cast<PD_IPC::PdIpcSimulator*>(_sim);

	if (_selectionMode == CTRL_SELECT)
	{
		_selectedTetPointId = -1;
		_mouseForce.setZero();

		std::vector<Vector3> selectList;
		std::vector<int> select_index_list;
		Vector3 cent = Vector3(0, 0, 0);
		for (int i = 0; i < nbHits; i++)
		{
			int buf_id = selectBuffer()[4 * i + 3];
			if (buf_id < sim->simPointsNum)
			{	
				qeal x, y, z;
				int boundaryId = buf_id;
				Vector3 pos = sim->points.col(boundaryId);
				selectList.push_back(pos);
				select_index_list.push_back(buf_id);
				cent += pos;
			}
		}

		cent /= selectList.size();
		selectedTetPointListCent = cent;
		qeal min_dist = QEAL_MAX;
		for (int i = 0; i < selectList.size(); i++)
		{
			qeal len = (selectList[i] - cent).norm();
			if (len < min_dist)
			{
				min_dist = len;
				_selectedTetPointId = select_index_list[i];
			}
		}

		//hostPointsModelsId
		if (_selectedTetPointId > 0)
		{
			int selectModel = sim->hostPointsModelsId[_selectedTetPointId];

			//selectedTetPointList.push_back(buf_id);
			for (int i = 0; i < select_index_list.size(); i++)
			{
				if(sim->hostPointsModelsId[select_index_list[i]] == selectModel)
					selectedTetPointList.push_back(select_index_list[i]);
			}
		}
	
	}

}

void BaseScene::initCamera()
{
	_camera = new StandardCamera();
	qglviewer::Camera* c = camera();
	setCamera(_camera);
	delete c;

	qglviewer::Vec world_origin = qglviewer::Vec(0.f, 0.f, 0.f);
	setSceneCenter(world_origin);
	camera()->setType(qglviewer::Camera::PERSPECTIVE);
	camera()->setPosition(qglviewer::Vec(0.0, 3.0, 10.0));
	camera()->setViewDirection(-qglviewer::Vec(0.0, 3.0, 10.0));
	_lightCamera.setType(qglviewer::Camera::ORTHOGRAPHIC);
}

void BaseScene::initShaderProgram()
{
	std::string vertexShader = "Shader/Program/phong_shader_program.vs";
	std::string fragShader = "Shader/Program/phong_shader_program.fs";
	_phongProgram.initProgram(vertexShader.data(), fragShader.data());
	vertexShader = "Shader/Program/phong_inv_shader_program.vs";
	fragShader = "Shader/Program/phong_inv_shader_program.fs";
	_phongInvProgram.initProgram(vertexShader.data(), fragShader.data());
	vertexShader = "Shader/Program/depth_map.vs";
	fragShader = "Shader/Program/depth_map.fs";
	_shadowProgram.initProgram(vertexShader.data(), fragShader.data());
	vertexShader = "Shader/Program/text_shader_program.vs";
	fragShader = "Shader/Program/text_shader_program.fs";
	_textProgram.initProgram(vertexShader.data(), fragShader.data());
	_textPainter.generateFont(QOpenGLContext::currentContext()->functions());
}

void BaseScene::initShadowFramebuffer()
{
	glGenFramebuffers(1, &_depthFramebuffer);
	glGenTextures(1, &_depthMapId);
	glBindTexture(GL_TEXTURE_2D, _depthMapId);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 1024, 1024, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	GLfloat borderColor[] = { 1.0, 1.0, 1.0, 1.0 };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

	glBindFramebuffer(GL_FRAMEBUFFER, _depthFramebuffer);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, _depthMapId, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void BaseScene::initAnimation()
{
	stopAnimation();
	setAnimationPeriod(0);
	QDateTime current_date_time = QDateTime::currentDateTime();
	QString snapshotFileName = QString(_status.animationFilePath.c_str()) + QString("//") + current_date_time.toString("yyyy-MM-dd-hh-mm");
	int code = mkdir(snapshotFileName.toStdString().c_str());
	snapshotFileName += QString("/frame");
	setSnapshotFileName(snapshotFileName);
	setSnapshotFormat(QString("PNG"));
	setSnapshotQuality(100);
	connect(this, SIGNAL(recordAnimation(bool)), SLOT(saveSnapshot(bool)));
}

void BaseScene::renderDepthToBuffer()
{
	if (!_lights.useShadow())
		return;
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	_shadowProgram.bind();
	QVector3D lightPos = _lights.getMainLightPosition();
	_lightCamera.setPosition(qglviewer::Vec(lightPos[0], lightPos[1], lightPos[2]));
	_lightCamera.lookAt(qglviewer::Vec(0.0, 0.0, 0.0));

	_shadowProgram.setUniformValue("lightSpaceMatrix", _lightCamera.getProjectionViewQMatrix());
	glBindFramebuffer(GL_FRAMEBUFFER, _depthFramebuffer);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, 1024, 1024);
	glCullFace(GL_FRONT);
	renderScene(&_shadowProgram);
	glCullFace(GL_BACK);
	_shadowProgram.release();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	//glDisable (GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

void BaseScene::renderScene(QOpenGLShaderProgram* program)
{
	_floor->render(program, QOpenGLContext::currentContext()->functions());
	if (_sim != nullptr)
	{
		//
		//_sim->render(&_phongInvProgram, QOpenGLContext::currentContext()->functions());
		_sim->render(program, QOpenGLContext::currentContext()->functions());
		if (_status.renderLineMode)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			_sim->render(program, QOpenGLContext::currentContext()->functions(), true);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
	}
}

void BaseScene::renderScene()
{
	if (_sim != nullptr)
	{
		_sim->renderExtraElementOnCPU(QOpenGLContext::currentContext()->functions());
		drawDebug();
	}
}

void BaseScene::renderCornerAxis()
{

	int viewport[4];
	int scissor[4];

	// The viewport and the scissor are changed to fit the lower left
	// corner. Original values are saved.
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetIntegerv(GL_SCISSOR_BOX, scissor);

	// Axis viewport size, in pixels
	const int size = 150;
	glViewport(0, 0, size, size);
	glScissor(0, 0, size, size);

	// The Z-buffer is cleared to make the axis appear over the
	// original image.
	glClear(GL_DEPTH_BUFFER_BIT);

	// Tune for best line rendering
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd(camera()->orientation().inverse().matrix());

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(1.0, 0.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 1.0, 0.0);

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 1.0);
	glEnd();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glEnable(GL_LIGHTING);

	// The viewport and the scissor are restored.
	glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
	glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
}

void BaseScene::renderInfoOnScene()
{
	double sw = camera()->screenWidth();
	double sh = camera()->screenHeight();
	QMatrix4x4 ortho_matrix = QMatrix4x4(2.0 / (sw), 0, 0, -1,
		0, 2.0 / (sh), 0, -1,
		0, 0, 2, -1,
		0, 0, 0, 1);

	_textProgram.bind();
	_textProgram.setUniformValue("mvp_matrix", ortho_matrix);

	const int w_step = 150;
	const int h_step = 20;

	int fps_w = sw - w_step;
	int fps_h = sh - h_step;

	QVector3D fontColor = QVector3D(0.3, 0.8, 0.3);
	qeal avg_cost = _status.animationCost[0] + _status.animationCost[1] + _status.animationCost[2] + _status.animationCost[3] + _status.animationCost[4] + _status.animationCost[5] + _status.animationCost[6] + _status.animationCost[7];
	avg_cost /= 8.0;
	if (avg_cost != 0.0 && _status.animationLifeTime > 0)
		_status.avg_fps = 1000.0 / (avg_cost);
	else _status.avg_fps = 1000.0;
	if (_status.avg_fps > 1000.0) _status.avg_fps = 1000.0;
	QString fps_str = QString::number(_status.avg_fps, 'f', 2);
	_textPainter.renderText(&_textProgram, "FPS: " + fps_str.toStdString() + " Hz", fps_w, fps_h, 0.5, fontColor);

	QString frame_str = QString::number(_status.animationLifeTime);
	_textPainter.renderText(&_textProgram, "Frames:   " + frame_str.toStdString(), fps_w, fps_h - h_step, 0.5, fontColor);

	_textProgram.release();
}

void BaseScene::renderSelectionRect()
{
	glBlendFunc(GL_ONE, GL_ONE);
	startScreenCoordinatesSystem();
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

	glColor4f(0.0, 0.0, 0.3f, 0.3f);
	glBegin(GL_QUADS);
	glVertex2i(_rectangle.left(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.bottom());
	glVertex2i(_rectangle.left(), _rectangle.bottom());
	glEnd();

	glLineWidth(2.0);
	glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
	glBegin(GL_LINE_LOOP);
	glVertex2i(_rectangle.left(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.bottom());
	glVertex2i(_rectangle.left(), _rectangle.bottom());
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	stopScreenCoordinatesSystem();
	update();
}

void BaseScene::renderSelectionPoints()
{
	if (_sim == nullptr)
		return;
	glPointSize(10.0);
	if (_status.enableSelectSurfacePoints)
	{
		for (auto it = _sim->selectSurfacePoints.begin(); it != _sim->selectSurfacePoints.end(); it++)
		{
			int i = *it;
			glColor3f(0.8, 0.2, 0.2);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getSurfacePoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
	else if (_status.enableSelectTetNodes)
	{
		for (auto it = _sim->selectTetPoints.begin(); it != _sim->selectTetPoints.end(); it++)
		{
			int i = *it;
			glColor3f(0.2, 0.8, 0.2);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getTetPoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
}

void BaseScene::renderForceLine()
{
	//if (_selectedTetPointId == -1 || _selectionMode != CTRL_SELECT || _sim == nullptr)
		//return;
	if (selectedTetPointList.size() == 0 || _selectionMode != CTRL_SELECT || _sim == nullptr)
		return;
	if (_mouseForce.norm() < 1e-8)
		return;

	PD_IPC::PdIpcSimulator* sim = dynamic_cast<PD_IPC::PdIpcSimulator*>(_sim);

	qeal x, y, z;
	//int boundaryId = _selectedTetPointId;
	//Vector3 selectedPos = sim->points.col(boundaryId);
	selectedTetPointListCent = Vector3::Zero();
	for (int i = 0; i < selectedTetPointList.size(); i++)
	{
		selectedTetPointListCent += sim->points.col(selectedTetPointList[i]);
	}
	selectedTetPointListCent /= selectedTetPointList.size();

	Vector3 selectedPos = selectedTetPointListCent;
	x = selectedPos[0]; y = selectedPos[1]; z = selectedPos[2];
	qglviewer::Vec pro_pos = camera()->projectedCoordinatesOf(qglviewer::Vec(x, y, z));

	qglviewer::Vec norm = qglviewer::Vec(_mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
	norm.normalize();
	qeal scale = qglviewer::Vec(_mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]).norm();
	if (scale > 2.0)
		scale = 2.0;

	qglviewer::Vec epos = qglviewer::Vec(selectedPos.data()[0], selectedPos.data()[1], selectedPos.data()[2]) + norm * scale;// qglviewer::Vec(_mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
	if (epos.y < 0.001) epos.y = 0.001;
	qglviewer::Vec spos = qglviewer::Vec(selectedPos.data()[0], selectedPos.data()[1], selectedPos.data()[2]);

	glBegin(GL_LINES);
	glColor3f(240.0 / 255.0, 240.0 / 255.0, 240.0 / 255.0);
	glVertex3d(spos.x, spos.y, spos.z);
	glColor3f(240.0 / 255.0, 240.0 / 255.0, 240.0 / 255.0);
	glVertex3d(epos.x, epos.y, epos.z);
	glEnd();

	double r = 0.1;

	glPushMatrix();
	glTranslated(spos.x, spos.y, spos.z);
	glColor3f(0.8, 0.2, 0.2);
	gluSphere(gluNewQuadric(), r, 50, 50);
	glPopMatrix();
	update();

	glPushMatrix();
	glTranslated(epos.x, epos.y, epos.z);
	glColor3f(0.2, 0.8, 0.2);
	gluSphere(gluNewQuadric(), r, 50, 50);
	glPopMatrix();
}

void BaseScene::computeMouseForce()
{
	//if (_selectedTetPointId == -1 || _sim == nullptr)
	//{
	//	_mouseForce = Vector3(0, 0, 0);
	//	getCurrentSimulator()->handleMouseForce(-1, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
	//	return;
	//}

	if (selectedTetPointList.size() == 0 || _sim == nullptr)
	{
		_mouseForce = Vector3(0, 0, 0);
		lastMouseForce = _mouseForce;
		selectedTetPointList.clear();
		getCurrentSimulator()->handleMouseForce(selectedTetPointList, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
		return;
	}

	if (abs(_forceLine.x1() - _forceLine.x2()) < 0.001&& abs(_forceLine.y1() == _forceLine.y2()) < 0.001)
	{
		_mouseForce = Vector3(0, 0, 0);
		lastMouseForce = _mouseForce;
		selectedTetPointList.clear();
		getCurrentSimulator()->handleMouseForce(selectedTetPointList, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
		return;
	}

	PD_IPC::PdIpcSimulator* sim = dynamic_cast<PD_IPC::PdIpcSimulator*>(_sim);

	qeal x, y, z;

	//int boundaryId = _selectedTetPointId;
	//Vector3 selectedPos = sim->points.col(boundaryId);
	selectedTetPointListCent = Vector3::Zero();
	for (int i = 0; i < selectedTetPointList.size(); i++)
	{
		selectedTetPointListCent += sim->points.col(selectedTetPointList[i]);
	}
	selectedTetPointListCent /= selectedTetPointList.size();

	Vector3 selectedPos = selectedTetPointListCent;
	x = selectedPos[0]; y = selectedPos[1]; z = selectedPos[2];
	qglviewer::Vec pro_pos = camera()->projectedCoordinatesOf(qglviewer::Vec(x, y, z));

	qglviewer::Vec endPos = qglviewer::Vec(_forceLine.x2(), _forceLine.y2(), pro_pos.z);
	qglviewer::Vec spos = qglviewer::Vec(x, y, z);
	qglviewer::Vec epos = camera()->unprojectedCoordinatesOf(endPos);

	//
	//std::cout << "mouse " << std::endl;
	//std::cout << " --- epos" << epos.x << " " << epos.y << " " << epos.z << std::endl;
	//std::cout << " --- spos" << spos.x << " " << spos.y << " " << spos.z << std::endl;
	//
	//
	qglviewer::Vec dir = epos - spos;
	_mouseForce = Vector3(dir.x, dir.y, dir.z);

	if ((lastMouseForce - _mouseForce).norm() > 8)
	{
		_mouseForce.setZero();
		selectedTetPointList.clear();
	}

	//qeal scale = 5000;

	//_mouseForce *= scale;

//	getCurrentSimulator()->handleMouseForce(boundaryId, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);

	getCurrentSimulator()->handleMouseForce(selectedTetPointList, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
	
	lastMouseForce = _mouseForce;
}

void BaseScene::drawDebug()
{
	QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
	if (!_sim)
		return;
	return;

	if (_sim->getModelsNum() == 0) return;
	PD_IPC::PdIpcSimulator* sim = dynamic_cast<PD_IPC::PdIpcSimulator*>(_sim);
	if (sim->simPointsNum <= 0) return;

	for (int i = 0; i < sim->boundaryFaces.cols(); i++)
	{
		Vector3i face = sim->boundaryFaces.col(i);
		Vector3 p[3];
		Vector3 normal;
		for (int k = 0; k < 3; k++)
		{
			p[k] = sim->points.col(face[k]);
		}
		normal = (p[1] - p[0]).cross(p[2] - p[0]);
		normal.normalize();

		Vector3 cent = (p[0] + p[1] + p[2]) / 3;

		Vector3 cent1 = cent + normal * 0.03;

		glBegin(GL_LINES);
		glColor3d(1, 0, 0);
		glVertex3d(cent[0], cent[1], cent[2]);
		glColor3d(0, 0, 1);
		glVertex3d(cent1[0], cent1[1], cent1[2]);
		glEnd();
	}


}

void BaseScene::removeTet()
{
	BaseTetElementSet ele_set = _sim->models[0]->getTetMeshHandle()->_elementSet[select_tet_group];
	std::vector<int> elelist;
	std::set<int> list;
	ele_set.getElements(list);
	std::set<int>::iterator cit = list.begin();
	for (; cit != list.end(); ++cit)
	{
		elelist.push_back(*cit);
	}

	std::set<int>::iterator it = select_tet_element.begin();
	for (; it != select_tet_element.end(); ++it)
	{
		elelist[*it] = -1;
	}

	ele_set.clear();
	for (int i = 0; i < elelist.size(); i++)
	{
		if (elelist[i] != -1)
		{
			ele_set.insert(elelist[i]);
		}
	}

	_sim->models[0]->getTetMeshHandle()->_elementSet[select_tet_group] = ele_set;
	removeSelect();
}

void BaseScene::removeSelect()
{
	select_tet_element.clear();
}
