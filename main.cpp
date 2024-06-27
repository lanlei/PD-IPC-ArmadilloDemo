#include "SimFramework.h"
#include <QtWidgets/QApplication>
#include <igl/Timer.h>

int main(int argc, char *argv[])
{	
	
	QApplication a(argc, argv);
	SimFramework w;
	BaseMainWidget* simWidget = new BaseMainWidget();
	simWidget->connectUI(&w);

	BaseSimulator* sim = new BaseSimulator();
	simWidget->getCenteralScene()->bindSimulator(sim);

	w.setWindowTitle("Sim Framework");
	w.showMaximized();

	return a.exec();
}
