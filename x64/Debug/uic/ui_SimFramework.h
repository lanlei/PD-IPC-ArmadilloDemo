/********************************************************************************
** Form generated from reading UI file 'SimFramework.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SIMFRAMEWORK_H
#define UI_SIMFRAMEWORK_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SimFrameworkClass
{
public:
    QWidget *centralWidget;

    void setupUi(QMainWindow *SimFrameworkClass)
    {
        if (SimFrameworkClass->objectName().isEmpty())
            SimFrameworkClass->setObjectName(QString::fromUtf8("SimFrameworkClass"));
        SimFrameworkClass->resize(1030, 855);
        centralWidget = new QWidget(SimFrameworkClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        SimFrameworkClass->setCentralWidget(centralWidget);

        retranslateUi(SimFrameworkClass);

        QMetaObject::connectSlotsByName(SimFrameworkClass);
    } // setupUi

    void retranslateUi(QMainWindow *SimFrameworkClass)
    {
        SimFrameworkClass->setWindowTitle(QCoreApplication::translate("SimFrameworkClass", "SimFramework", nullptr));
    } // retranslateUi

};

namespace Ui {
    class SimFrameworkClass: public Ui_SimFrameworkClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SIMFRAMEWORK_H
