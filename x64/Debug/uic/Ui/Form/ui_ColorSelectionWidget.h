/********************************************************************************
** Form generated from reading UI file 'ColorSelectionWidget.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_COLORSELECTIONWIDGET_H
#define UI_COLORSELECTIONWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ColorSelectionWidget
{
public:
    QLabel *r_color_label;
    QLabel *g_color_label;
    QLabel *b_color_label;
    QLineEdit *r_color_lineEdit;
    QLineEdit *b_color_lineEdit;
    QLineEdit *g_color_lineEdit;
    QLineEdit *alpha_color_lineEdit;
    QLabel *alpha_color_label;

    void setupUi(QWidget *ColorSelectionWidget)
    {
        if (ColorSelectionWidget->objectName().isEmpty())
            ColorSelectionWidget->setObjectName(QString::fromUtf8("ColorSelectionWidget"));
        ColorSelectionWidget->resize(380, 200);
        QFont font;
        font.setFamily(QString::fromUtf8("Times New Roman"));
        font.setPointSize(10);
        ColorSelectionWidget->setFont(font);
        r_color_label = new QLabel(ColorSelectionWidget);
        r_color_label->setObjectName(QString::fromUtf8("r_color_label"));
        r_color_label->setGeometry(QRect(220, 30, 21, 41));
        g_color_label = new QLabel(ColorSelectionWidget);
        g_color_label->setObjectName(QString::fromUtf8("g_color_label"));
        g_color_label->setGeometry(QRect(250, 100, 21, 41));
        b_color_label = new QLabel(ColorSelectionWidget);
        b_color_label->setObjectName(QString::fromUtf8("b_color_label"));
        b_color_label->setGeometry(QRect(180, 110, 21, 41));
        r_color_lineEdit = new QLineEdit(ColorSelectionWidget);
        r_color_lineEdit->setObjectName(QString::fromUtf8("r_color_lineEdit"));
        r_color_lineEdit->setGeometry(QRect(30, 80, 21, 20));
        b_color_lineEdit = new QLineEdit(ColorSelectionWidget);
        b_color_lineEdit->setObjectName(QString::fromUtf8("b_color_lineEdit"));
        b_color_lineEdit->setGeometry(QRect(70, 130, 21, 20));
        g_color_lineEdit = new QLineEdit(ColorSelectionWidget);
        g_color_lineEdit->setObjectName(QString::fromUtf8("g_color_lineEdit"));
        g_color_lineEdit->setGeometry(QRect(100, 110, 21, 20));
        alpha_color_lineEdit = new QLineEdit(ColorSelectionWidget);
        alpha_color_lineEdit->setObjectName(QString::fromUtf8("alpha_color_lineEdit"));
        alpha_color_lineEdit->setGeometry(QRect(140, 150, 21, 20));
        alpha_color_label = new QLabel(ColorSelectionWidget);
        alpha_color_label->setObjectName(QString::fromUtf8("alpha_color_label"));
        alpha_color_label->setGeometry(QRect(250, 170, 21, 41));

        retranslateUi(ColorSelectionWidget);

        QMetaObject::connectSlotsByName(ColorSelectionWidget);
    } // setupUi

    void retranslateUi(QWidget *ColorSelectionWidget)
    {
        ColorSelectionWidget->setWindowTitle(QCoreApplication::translate("ColorSelectionWidget", "Form", nullptr));
        r_color_label->setText(QCoreApplication::translate("ColorSelectionWidget", "R", nullptr));
        g_color_label->setText(QCoreApplication::translate("ColorSelectionWidget", "G", nullptr));
        b_color_label->setText(QCoreApplication::translate("ColorSelectionWidget", "B", nullptr));
        alpha_color_label->setText(QCoreApplication::translate("ColorSelectionWidget", "A", nullptr));
    } // retranslateUi

};

namespace Ui {
    class ColorSelectionWidget: public Ui_ColorSelectionWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_COLORSELECTIONWIDGET_H
