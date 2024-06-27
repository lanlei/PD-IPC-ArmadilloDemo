#include "BaseBottomWidget.h"

BaseBottomWidget::BaseBottomWidget()
{
	setupUi(this);
}

BaseBottomWidget::~BaseBottomWidget()
{
}

void BaseBottomWidget::resizeEvent(QResizeEvent* size)
{
	int h = height();
	int w = width();
	resize(QSize(w, h));
}

void BaseBottomWidget::setBaseScene(BaseScene * scene)
{
	_scene = scene; 
	if (_scene == nullptr)
		return;
	connect(_scene, SIGNAL(initSceneSignal()), this, SLOT(bindMultiViewerGroup()));
	connect(_scene, SIGNAL(initSceneSignal()), this, SLOT(bindAnimationGroup()));
}

void BaseBottomWidget::bindMultiViewerGroup()
{
	multiViewerButtonGroup = new QButtonGroup(this);
	multiViewerButtonGroup->setExclusive(true);
	multiViewerButtonGroup->addButton(single_viewer_radioButton);
	multiViewerButtonGroup->addButton(two_viewer_radioButton);
	multiViewerButtonGroup->addButton(three_viewer_radioButton);
	multiViewerButtonGroup->addButton(four_viewer_radioButton);

	multiViewerButtonGroup->setId(single_viewer_radioButton, 0);
	multiViewerButtonGroup->setId(two_viewer_radioButton, 1);
	multiViewerButtonGroup->setId(three_viewer_radioButton, 2);
	multiViewerButtonGroup->setId(four_viewer_radioButton, 3);
	single_viewer_radioButton->setChecked(true);

	connect(multiViewerButtonGroup, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(handleMultiViewerButtonGroup(QAbstractButton*)));
}

void BaseBottomWidget::bindAnimationGroup()
{
	playAnimationIcon = new QIcon("./icon/play_animation.png");
	pauseAnimationIcon = new QIcon("./icon/pause_animation.png");
	recordAnimationIcon = new QIcon("./icon/record_animation.png");
	resetAnimationIcon = new QIcon("./icon/reset_animation.png");
	initSimulatorIcon = new QIcon("./icon/init_simulator.png");
	saveAnimationIcon = new QIcon("./icon/save_obj.png");
	
	play_animation_Button->setIcon(*playAnimationIcon);
	play_animation_Button->setIconSize(play_animation_Button->size());

	init_simulator_Button->setIcon(*initSimulatorIcon);
	init_simulator_Button->setIconSize(init_simulator_Button->size());

	reset_Button->setIcon(*resetAnimationIcon);
	reset_Button->setIconSize(reset_Button->size());

	record_animation_Button->setIcon(*recordAnimationIcon);
	record_animation_Button->setIconSize(record_animation_Button->size());

	save_animation_Button->setIcon(*saveAnimationIcon);
	save_animation_Button->setIconSize(save_animation_Button->size());

	connect(play_animation_Button, SIGNAL(clicked(bool)), this, SLOT(handlePlayAnimation()));
	connect(reset_Button, SIGNAL(clicked(bool)), this, SLOT(handleResetAnimation()));
	connect(record_animation_Button, SIGNAL(clicked(bool)), this, SLOT(handleRecordAnimation()));
	connect(save_animation_Button, SIGNAL(clicked(bool)), this, SLOT(handleSaveAnimation()));
	connect(init_simulator_Button, SIGNAL(clicked(bool)), this, SLOT(handleInitSimulator()));
}

void BaseBottomWidget::handleRecordAnimation()
{
	Q_EMIT recordAnimation(record_animation_Button->isChecked());
}

void BaseBottomWidget::handleSaveAnimation()
{
	_scene->getSceneStaus()->saveAnimation = save_animation_Button->isChecked();
}

void BaseBottomWidget::handlePlayAnimation()
{
	Q_EMIT playAnimation(play_animation_Button->isChecked());

	if (play_animation_Button->isChecked())
	{
		play_animation_Button->setIcon(*pauseAnimationIcon);
		reset_Button->setEnabled(false);
	}
	else
	{
		play_animation_Button->setIcon(*playAnimationIcon);
		reset_Button->setEnabled(true);
	}
}

void BaseBottomWidget::handleResetAnimation()
{
	Q_EMIT resetAnimation();
	resetInitSimulatorStatus();
}

void BaseBottomWidget::resetInitSimulatorStatus()
{
	play_animation_Button->setChecked(false);
	play_animation_Button->setEnabled(false);
	init_simulator_Button->setChecked(false);
}

void BaseBottomWidget::handleMultiViewerButtonGroup(QAbstractButton * btn)
{
	quint16 id = multiViewerButtonGroup->checkedId();
	switch (id)
	{
	case 0:
		Q_EMIT changeViewerNum(0);
		break;
	case 1:
		Q_EMIT changeViewerNum(1);
		break;
	case 2:
		Q_EMIT changeViewerNum(2);
		break;
	case 3:
		Q_EMIT changeViewerNum(3);
		break;
	default:
		Q_EMIT changeViewerNum(0);
		break;
	}
}

void BaseBottomWidget::handleInitSimulator()
{
	_scene->getCurrentSimulator()->initialization();
	play_animation_Button->setEnabled(true);
	std::cout << "Simulator have been initialized !" << std::endl;
}
