/********************************************************************************
** Form generated from reading UI file 'BottomFormWidget.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BOTTOMFORMWIDGET_H
#define UI_BOTTOMFORMWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_BottomWidget
{
public:
    QHBoxLayout *horizontalLayout_2;
    QGroupBox *multi_viewer_groupBox;
    QGridLayout *gridLayout;
    QRadioButton *single_viewer_radioButton;
    QRadioButton *two_viewer_radioButton;
    QRadioButton *three_viewer_radioButton;
    QRadioButton *four_viewer_radioButton;
    QGroupBox *animation_groupBox;
    QGridLayout *gridLayout_2;
    QPushButton *reset_Button;
    QPushButton *record_animation_Button;
    QPushButton *init_simulator_Button;
    QPushButton *play_animation_Button;
    QPushButton *save_animation_Button;

    void setupUi(QWidget *BottomWidget)
    {
        if (BottomWidget->objectName().isEmpty())
            BottomWidget->setObjectName(QString::fromUtf8("BottomWidget"));
        BottomWidget->resize(1025, 91);
        horizontalLayout_2 = new QHBoxLayout(BottomWidget);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        multi_viewer_groupBox = new QGroupBox(BottomWidget);
        multi_viewer_groupBox->setObjectName(QString::fromUtf8("multi_viewer_groupBox"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(multi_viewer_groupBox->sizePolicy().hasHeightForWidth());
        multi_viewer_groupBox->setSizePolicy(sizePolicy);
        multi_viewer_groupBox->setMinimumSize(QSize(200, 70));
        multi_viewer_groupBox->setMaximumSize(QSize(450, 70));
        QFont font;
        font.setFamily(QString::fromUtf8("Times New Roman"));
        font.setPointSize(12);
        multi_viewer_groupBox->setFont(font);
        multi_viewer_groupBox->setLayoutDirection(Qt::LeftToRight);
        multi_viewer_groupBox->setAlignment(Qt::AlignCenter);
        gridLayout = new QGridLayout(multi_viewer_groupBox);
        gridLayout->setSpacing(0);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, -1, 0, 0);
        single_viewer_radioButton = new QRadioButton(multi_viewer_groupBox);
        single_viewer_radioButton->setObjectName(QString::fromUtf8("single_viewer_radioButton"));
        QFont font1;
        font1.setPointSize(10);
        single_viewer_radioButton->setFont(font1);

        gridLayout->addWidget(single_viewer_radioButton, 0, 0, 1, 1);

        two_viewer_radioButton = new QRadioButton(multi_viewer_groupBox);
        two_viewer_radioButton->setObjectName(QString::fromUtf8("two_viewer_radioButton"));
        two_viewer_radioButton->setFont(font1);

        gridLayout->addWidget(two_viewer_radioButton, 0, 1, 1, 1);

        three_viewer_radioButton = new QRadioButton(multi_viewer_groupBox);
        three_viewer_radioButton->setObjectName(QString::fromUtf8("three_viewer_radioButton"));
        three_viewer_radioButton->setFont(font1);

        gridLayout->addWidget(three_viewer_radioButton, 1, 0, 1, 1);

        four_viewer_radioButton = new QRadioButton(multi_viewer_groupBox);
        four_viewer_radioButton->setObjectName(QString::fromUtf8("four_viewer_radioButton"));
        four_viewer_radioButton->setFont(font1);

        gridLayout->addWidget(four_viewer_radioButton, 1, 1, 1, 1);


        horizontalLayout_2->addWidget(multi_viewer_groupBox);

        animation_groupBox = new QGroupBox(BottomWidget);
        animation_groupBox->setObjectName(QString::fromUtf8("animation_groupBox"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(animation_groupBox->sizePolicy().hasHeightForWidth());
        animation_groupBox->setSizePolicy(sizePolicy1);
        animation_groupBox->setMinimumSize(QSize(0, 70));
        animation_groupBox->setMaximumSize(QSize(16777215, 70));
        animation_groupBox->setFont(font);
        animation_groupBox->setLayoutDirection(Qt::LeftToRight);
        animation_groupBox->setAlignment(Qt::AlignCenter);
        gridLayout_2 = new QGridLayout(animation_groupBox);
        gridLayout_2->setSpacing(3);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(-1, 0, -1, 3);
        reset_Button = new QPushButton(animation_groupBox);
        reset_Button->setObjectName(QString::fromUtf8("reset_Button"));
        reset_Button->setMinimumSize(QSize(45, 45));
        reset_Button->setMaximumSize(QSize(45, 45));
        reset_Button->setFont(font1);
        reset_Button->setIconSize(QSize(96, 96));
        reset_Button->setFlat(true);

        gridLayout_2->addWidget(reset_Button, 0, 7, 1, 1);

        record_animation_Button = new QPushButton(animation_groupBox);
        record_animation_Button->setObjectName(QString::fromUtf8("record_animation_Button"));
        record_animation_Button->setMinimumSize(QSize(45, 45));
        record_animation_Button->setMaximumSize(QSize(45, 45));
        record_animation_Button->setFont(font1);
        record_animation_Button->setIconSize(QSize(96, 96));
        record_animation_Button->setCheckable(true);
        record_animation_Button->setAutoExclusive(false);
        record_animation_Button->setFlat(true);

        gridLayout_2->addWidget(record_animation_Button, 0, 0, 1, 1);

        init_simulator_Button = new QPushButton(animation_groupBox);
        init_simulator_Button->setObjectName(QString::fromUtf8("init_simulator_Button"));
        init_simulator_Button->setMinimumSize(QSize(45, 45));
        init_simulator_Button->setMaximumSize(QSize(45, 45));
        init_simulator_Button->setFont(font1);
        init_simulator_Button->setIconSize(QSize(96, 96));
        init_simulator_Button->setCheckable(true);
        init_simulator_Button->setAutoExclusive(false);
        init_simulator_Button->setFlat(true);

        gridLayout_2->addWidget(init_simulator_Button, 0, 5, 1, 1);

        play_animation_Button = new QPushButton(animation_groupBox);
        play_animation_Button->setObjectName(QString::fromUtf8("play_animation_Button"));
        play_animation_Button->setEnabled(false);
        play_animation_Button->setMinimumSize(QSize(45, 45));
        play_animation_Button->setMaximumSize(QSize(45, 45));
        play_animation_Button->setFont(font1);
        play_animation_Button->setIconSize(QSize(96, 96));
        play_animation_Button->setCheckable(true);
        play_animation_Button->setAutoExclusive(false);
        play_animation_Button->setFlat(true);

        gridLayout_2->addWidget(play_animation_Button, 0, 6, 1, 1);

        save_animation_Button = new QPushButton(animation_groupBox);
        save_animation_Button->setObjectName(QString::fromUtf8("save_animation_Button"));
        save_animation_Button->setMinimumSize(QSize(45, 45));
        save_animation_Button->setMaximumSize(QSize(45, 45));
        save_animation_Button->setIconSize(QSize(96, 96));
        save_animation_Button->setCheckable(true);
        save_animation_Button->setFlat(true);

        gridLayout_2->addWidget(save_animation_Button, 0, 1, 1, 1);


        horizontalLayout_2->addWidget(animation_groupBox);


        retranslateUi(BottomWidget);

        QMetaObject::connectSlotsByName(BottomWidget);
    } // setupUi

    void retranslateUi(QWidget *BottomWidget)
    {
        BottomWidget->setWindowTitle(QCoreApplication::translate("BottomWidget", "Form", nullptr));
        multi_viewer_groupBox->setTitle(QCoreApplication::translate("BottomWidget", "MultiViewer", nullptr));
        single_viewer_radioButton->setText(QCoreApplication::translate("BottomWidget", "Single Viewer", nullptr));
        two_viewer_radioButton->setText(QCoreApplication::translate("BottomWidget", "Two Viewer", nullptr));
        three_viewer_radioButton->setText(QCoreApplication::translate("BottomWidget", "Three Viewer", nullptr));
        four_viewer_radioButton->setText(QCoreApplication::translate("BottomWidget", "Four Viewer", nullptr));
        animation_groupBox->setTitle(QCoreApplication::translate("BottomWidget", "Animation", nullptr));
        reset_Button->setText(QString());
        record_animation_Button->setText(QString());
        init_simulator_Button->setText(QString());
        play_animation_Button->setText(QString());
        save_animation_Button->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class BottomWidget: public Ui_BottomWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BOTTOMFORMWIDGET_H
