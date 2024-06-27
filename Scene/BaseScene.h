#pragma once
#ifndef BASE_Scene_H
#define BASE_Scene_H
#include <windows.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions>
#include <QOpenGLFramebufferObject>
#include <QtWidgets>

#include <string>
#include <iostream>
#include <strstream>
#include <vector>
#include <set>

#include <QMatrix4x4>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QElapsedTimer>

#include "StandardCamera.h"
#include "BaseLighting.h"
#include "Model\GeometryElement.h"
#include "Shader\ShaderProgram.h"
#include "Shader\ScreenFontType.h"
#include "Commom\tinyxml\tinyxml.h"
#include "Simulator\BaseSimulator.h"

//select
enum SelectionMode {
	NONE_SELECT,
	ALT_SELECT,
	ALT_REMOVE,
	CTRL_SELECT,
	CTRL_REMOVE,
	SHIFT_SELECT,
	SHIFT_REMOVE
};

struct SceneStaus
{
	std::string animationFilePath = "./animation/";
	bool recordAnimation = false;
	bool saveAnimation = false;
	qeal animationCost[8] = {0,0,0,0,0,0,0,0};
	unsigned int animationLifeTime = 0;
	unsigned int animationMaxLifeTime = 2500;
	QElapsedTimer animationTimer;
	qeal avg_fps = 0;
	bool renderLineMode = false;
	bool enableSelectSurfacePoints = false;
	bool enableSelectTetNodes = false;

};

class BaseScene : public QGLViewer, protected QOpenGLFunctions
{
	Q_OBJECT
public:
	BaseScene();
	BaseScene(QWidget* parent);
	~BaseScene();

	virtual bool setFromConfigFile(const std::string filename, TiXmlElement* item);
	void bindSimulator(BaseSimulator* sim);
	BaseSimulator* getCurrentSimulator() { return _sim; }

	void enableRenderLineMode(bool flag) { _status.renderLineMode = flag; }
	void enableSelectSurfacePointsMode(bool flag) { _status.enableSelectSurfacePoints = flag; if (flag)  _status.enableSelectTetNodes = false; }
	void enableSelectTetNodesMode(bool flag) { _status.enableSelectTetNodes = flag;  if (flag)  _status.enableSelectSurfacePoints = false; }
Q_SIGNALS:
	void initSceneSignal();
	void changeSimulator();
	void recordAnimation(bool);
	void needInitSim();

protected:
	virtual void draw() override;
	virtual void animate() override;
	virtual void postDraw() override;
	virtual void init() override;
	virtual void mousePressEvent(QMouseEvent* e) override;
	virtual void mouseMoveEvent(QMouseEvent* e) override;
	virtual void mouseReleaseEvent(QMouseEvent *e) override;
	virtual void keyPressEvent(QKeyEvent *e) override;
	virtual void wheelEvent(QWheelEvent *e) override;
	virtual void drawWithNames() override;
	virtual void endSelection(const QPoint& point) override;

	virtual void initCamera();
	virtual void initShaderProgram();
	virtual void initShadowFramebuffer();
	virtual void initAnimation();
	//
	virtual void renderDepthToBuffer();
	virtual void renderScene(QOpenGLShaderProgram* program);
	virtual void renderScene();
	// post draw
	void renderCornerAxis();
	void renderInfoOnScene();
	void renderSelectionRect();
	void renderSelectionPoints();
	virtual void renderForceLine();

	virtual void computeMouseForce();


	void drawDebug();
	void removeTet();
	void removeSelect();
public:
	SceneLighting* getLightSetting() { return &_lights; }
	Floor* getFloor() { return _floor; }

	SceneStaus* getSceneStaus() { return &_status; }
protected:
	StandardCamera* _camera;
	StandardCamera _lightCamera;	

	SceneLighting _lights;
	Floor* _floor;

	// Shaders
	ShaderProgram _phongInvProgram;
	ShaderProgram _phongProgram;
	ShaderProgram _shadowProgram;
	ShaderProgram _textProgram;
	ScreenTextPainter::TextPainter _textPainter;

	GLuint _depthFramebuffer;
	GLuint _depthMapId;

	BaseSimulator* _sim;

	// status
	SceneStaus _status;
	SelectionMode _selectionMode;
	QRect _rectangle;

	//mouse force
	QLine _forceLine;
	int _selectedTetPointId;
	std::vector<int> selectedTetPointList;
	Vector3 selectedTetPointListCent;
	Vector3 _mouseForce;
	Vector3 lastMouseForce;

	//
	int select_tet_group;
	std::set<int> select_tet_element;

	int debugCount = 0;
};



#endif