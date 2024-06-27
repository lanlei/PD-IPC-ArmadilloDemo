#pragma on
#ifndef BASE_BOTTOM_WIDGET_H
#define BASE_BOTTOM_WIDGET_H

#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include "ui_BottomFormWidget.h"
#include "Scene\BaseScene.h"

class BaseBottomWidget : public QWidget, public Ui::BottomWidget
{
	Q_OBJECT
public:
	BaseBottomWidget();
	~BaseBottomWidget();

	void resizeEvent(QResizeEvent* size);

	void setBaseScene(BaseScene* scene);
Q_SIGNALS:
	void changeViewerNum(int);
	void recordAnimation(bool);
	void saveAnimation(bool);
	void playAnimation(bool);
	void resetAnimation();
public Q_SLOTS:
	void bindMultiViewerGroup();
	void bindAnimationGroup();
	void handleMultiViewerButtonGroup(QAbstractButton *);
	void handleInitSimulator();
	void handleRecordAnimation();
	void handleSaveAnimation();
	void handlePlayAnimation();
	void handleResetAnimation();
	void resetInitSimulatorStatus();
protected:
	BaseScene* _scene;
	QButtonGroup* multiViewerButtonGroup;
	QIcon* playAnimationIcon;
	QIcon* pauseAnimationIcon;
	QIcon* recordAnimationIcon;
	QIcon* saveAnimationIcon;
	QIcon* resetAnimationIcon;
	QIcon* initSimulatorIcon;

};


#endif
