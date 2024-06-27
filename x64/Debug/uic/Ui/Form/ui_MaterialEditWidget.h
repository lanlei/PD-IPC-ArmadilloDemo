/********************************************************************************
** Form generated from reading UI file 'MaterialEditWidget.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MATERIALEDITWIDGET_H
#define UI_MATERIALEDITWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MaterialEditWidget
{
public:
    QWidget *layoutWidget;
    QGridLayout *gridLayout_2;
    QVBoxLayout *texture_Layout;
    QLabel *texture_label;
    QHBoxLayout *horizontalLayout;
    QPushButton *ambient_texture_button;
    QLineEdit *ambient_texture_dir_lineEdit;
    QCheckBox *enable_amibent_texture_checkBox;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *diffuse_texture_button;
    QLineEdit *diffuse_texture_dir_lineEdit;
    QCheckBox *enable_diffuse_texture_checkBox;
    QHBoxLayout *horizontalLayout_6;
    QPushButton *specular_texture_button;
    QLineEdit *specular_texture_dir_lineEdit;
    QCheckBox *enable_specular_texture_checkBox;
    QHBoxLayout *horizontalLayout_7;
    QPushButton *bump_texture_button;
    QLineEdit *bump_texture_dir_lineEdit;
    QCheckBox *enable_bump_texture_checkBox;
    QVBoxLayout *listWidget_Layout;
    QLabel *material_list_label;
    QListWidget *material_listWidget;
    QVBoxLayout *uniform_color_Layout;
    QLabel *uniform_color_label;
    QGridLayout *gridLayout;
    QPushButton *ambient_color_button;
    QHBoxLayout *horizontalLayout_3;
    QLineEdit *amibent_color_r;
    QLineEdit *amibent_color_g;
    QLineEdit *amibent_color_b;
    QPushButton *diffuse_color_button;
    QHBoxLayout *horizontalLayout_4;
    QLineEdit *diffuse_color_r;
    QLineEdit *diffuse_color_g;
    QLineEdit *diffuse_color_b;
    QPushButton *specular_color_button;
    QHBoxLayout *horizontalLayout_5;
    QLineEdit *specular_color_r;
    QLineEdit *specular_color_g;
    QLineEdit *specular_color_b;
    QCheckBox *transparent_checkBox;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout_8;
    QLineEdit *shinness_lineEdit;
    QLabel *Shinness_label;
    QSpacerItem *horizontalSpacer;
    QSlider *ShinnessSlider;

    void setupUi(QWidget *MaterialEditWidget)
    {
        if (MaterialEditWidget->objectName().isEmpty())
            MaterialEditWidget->setObjectName(QString::fromUtf8("MaterialEditWidget"));
        MaterialEditWidget->resize(505, 347);
        MaterialEditWidget->setMinimumSize(QSize(0, 0));
        MaterialEditWidget->setMaximumSize(QSize(16777215, 16777215));
        layoutWidget = new QWidget(MaterialEditWidget);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 0, 479, 327));
        gridLayout_2 = new QGridLayout(layoutWidget);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        texture_Layout = new QVBoxLayout();
        texture_Layout->setObjectName(QString::fromUtf8("texture_Layout"));
        texture_label = new QLabel(layoutWidget);
        texture_label->setObjectName(QString::fromUtf8("texture_label"));
        QFont font;
        font.setFamily(QString::fromUtf8("Times New Roman"));
        font.setPointSize(12);
        texture_label->setFont(font);

        texture_Layout->addWidget(texture_label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        ambient_texture_button = new QPushButton(layoutWidget);
        ambient_texture_button->setObjectName(QString::fromUtf8("ambient_texture_button"));
        ambient_texture_button->setMinimumSize(QSize(75, 23));
        ambient_texture_button->setMaximumSize(QSize(75, 23));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Times New Roman"));
        font1.setPointSize(10);
        ambient_texture_button->setFont(font1);

        horizontalLayout->addWidget(ambient_texture_button);

        ambient_texture_dir_lineEdit = new QLineEdit(layoutWidget);
        ambient_texture_dir_lineEdit->setObjectName(QString::fromUtf8("ambient_texture_dir_lineEdit"));
        ambient_texture_dir_lineEdit->setMinimumSize(QSize(300, 25));
        ambient_texture_dir_lineEdit->setMaximumSize(QSize(1677215, 25));

        horizontalLayout->addWidget(ambient_texture_dir_lineEdit);

        enable_amibent_texture_checkBox = new QCheckBox(layoutWidget);
        enable_amibent_texture_checkBox->setObjectName(QString::fromUtf8("enable_amibent_texture_checkBox"));

        horizontalLayout->addWidget(enable_amibent_texture_checkBox);


        texture_Layout->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        diffuse_texture_button = new QPushButton(layoutWidget);
        diffuse_texture_button->setObjectName(QString::fromUtf8("diffuse_texture_button"));
        diffuse_texture_button->setMinimumSize(QSize(75, 23));
        diffuse_texture_button->setMaximumSize(QSize(75, 23));
        diffuse_texture_button->setFont(font1);

        horizontalLayout_2->addWidget(diffuse_texture_button);

        diffuse_texture_dir_lineEdit = new QLineEdit(layoutWidget);
        diffuse_texture_dir_lineEdit->setObjectName(QString::fromUtf8("diffuse_texture_dir_lineEdit"));
        diffuse_texture_dir_lineEdit->setMinimumSize(QSize(300, 25));
        diffuse_texture_dir_lineEdit->setMaximumSize(QSize(16777215, 25));

        horizontalLayout_2->addWidget(diffuse_texture_dir_lineEdit);

        enable_diffuse_texture_checkBox = new QCheckBox(layoutWidget);
        enable_diffuse_texture_checkBox->setObjectName(QString::fromUtf8("enable_diffuse_texture_checkBox"));

        horizontalLayout_2->addWidget(enable_diffuse_texture_checkBox);


        texture_Layout->addLayout(horizontalLayout_2);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        specular_texture_button = new QPushButton(layoutWidget);
        specular_texture_button->setObjectName(QString::fromUtf8("specular_texture_button"));
        specular_texture_button->setMinimumSize(QSize(75, 23));
        specular_texture_button->setMaximumSize(QSize(16777215, 16777215));
        specular_texture_button->setFont(font1);

        horizontalLayout_6->addWidget(specular_texture_button);

        specular_texture_dir_lineEdit = new QLineEdit(layoutWidget);
        specular_texture_dir_lineEdit->setObjectName(QString::fromUtf8("specular_texture_dir_lineEdit"));
        specular_texture_dir_lineEdit->setMinimumSize(QSize(300, 25));
        specular_texture_dir_lineEdit->setMaximumSize(QSize(16777215, 25));

        horizontalLayout_6->addWidget(specular_texture_dir_lineEdit);

        enable_specular_texture_checkBox = new QCheckBox(layoutWidget);
        enable_specular_texture_checkBox->setObjectName(QString::fromUtf8("enable_specular_texture_checkBox"));

        horizontalLayout_6->addWidget(enable_specular_texture_checkBox);


        texture_Layout->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        bump_texture_button = new QPushButton(layoutWidget);
        bump_texture_button->setObjectName(QString::fromUtf8("bump_texture_button"));
        bump_texture_button->setMinimumSize(QSize(75, 23));
        bump_texture_button->setMaximumSize(QSize(16777215, 16777215));
        bump_texture_button->setFont(font1);

        horizontalLayout_7->addWidget(bump_texture_button);

        bump_texture_dir_lineEdit = new QLineEdit(layoutWidget);
        bump_texture_dir_lineEdit->setObjectName(QString::fromUtf8("bump_texture_dir_lineEdit"));
        bump_texture_dir_lineEdit->setMinimumSize(QSize(300, 25));
        bump_texture_dir_lineEdit->setMaximumSize(QSize(16777215, 25));

        horizontalLayout_7->addWidget(bump_texture_dir_lineEdit);

        enable_bump_texture_checkBox = new QCheckBox(layoutWidget);
        enable_bump_texture_checkBox->setObjectName(QString::fromUtf8("enable_bump_texture_checkBox"));

        horizontalLayout_7->addWidget(enable_bump_texture_checkBox);


        texture_Layout->addLayout(horizontalLayout_7);


        gridLayout_2->addLayout(texture_Layout, 2, 0, 1, 2);

        listWidget_Layout = new QVBoxLayout();
        listWidget_Layout->setObjectName(QString::fromUtf8("listWidget_Layout"));
        material_list_label = new QLabel(layoutWidget);
        material_list_label->setObjectName(QString::fromUtf8("material_list_label"));
        material_list_label->setFont(font);

        listWidget_Layout->addWidget(material_list_label);

        material_listWidget = new QListWidget(layoutWidget);
        material_listWidget->setObjectName(QString::fromUtf8("material_listWidget"));

        listWidget_Layout->addWidget(material_listWidget);


        gridLayout_2->addLayout(listWidget_Layout, 0, 0, 1, 1);

        uniform_color_Layout = new QVBoxLayout();
        uniform_color_Layout->setObjectName(QString::fromUtf8("uniform_color_Layout"));
        uniform_color_label = new QLabel(layoutWidget);
        uniform_color_label->setObjectName(QString::fromUtf8("uniform_color_label"));
        uniform_color_label->setFont(font);

        uniform_color_Layout->addWidget(uniform_color_label);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        ambient_color_button = new QPushButton(layoutWidget);
        ambient_color_button->setObjectName(QString::fromUtf8("ambient_color_button"));
        ambient_color_button->setMinimumSize(QSize(75, 23));
        ambient_color_button->setMaximumSize(QSize(75, 23));
        ambient_color_button->setFont(font1);

        gridLayout->addWidget(ambient_color_button, 0, 0, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(3);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        amibent_color_r = new QLineEdit(layoutWidget);
        amibent_color_r->setObjectName(QString::fromUtf8("amibent_color_r"));
        amibent_color_r->setMinimumSize(QSize(40, 0));
        amibent_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(amibent_color_r);

        amibent_color_g = new QLineEdit(layoutWidget);
        amibent_color_g->setObjectName(QString::fromUtf8("amibent_color_g"));
        amibent_color_g->setMinimumSize(QSize(40, 0));
        amibent_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(amibent_color_g);

        amibent_color_b = new QLineEdit(layoutWidget);
        amibent_color_b->setObjectName(QString::fromUtf8("amibent_color_b"));
        amibent_color_b->setMinimumSize(QSize(40, 0));
        amibent_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(amibent_color_b);


        gridLayout->addLayout(horizontalLayout_3, 0, 1, 1, 1);

        diffuse_color_button = new QPushButton(layoutWidget);
        diffuse_color_button->setObjectName(QString::fromUtf8("diffuse_color_button"));
        diffuse_color_button->setMaximumSize(QSize(16777215, 16777215));
        diffuse_color_button->setFont(font1);

        gridLayout->addWidget(diffuse_color_button, 1, 0, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(3);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        diffuse_color_r = new QLineEdit(layoutWidget);
        diffuse_color_r->setObjectName(QString::fromUtf8("diffuse_color_r"));
        diffuse_color_r->setMinimumSize(QSize(40, 0));
        diffuse_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(diffuse_color_r);

        diffuse_color_g = new QLineEdit(layoutWidget);
        diffuse_color_g->setObjectName(QString::fromUtf8("diffuse_color_g"));
        diffuse_color_g->setMinimumSize(QSize(40, 0));
        diffuse_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(diffuse_color_g);

        diffuse_color_b = new QLineEdit(layoutWidget);
        diffuse_color_b->setObjectName(QString::fromUtf8("diffuse_color_b"));
        diffuse_color_b->setMinimumSize(QSize(40, 0));
        diffuse_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(diffuse_color_b);


        gridLayout->addLayout(horizontalLayout_4, 1, 1, 1, 1);

        specular_color_button = new QPushButton(layoutWidget);
        specular_color_button->setObjectName(QString::fromUtf8("specular_color_button"));
        specular_color_button->setMaximumSize(QSize(16777215, 16777215));
        specular_color_button->setFont(font1);

        gridLayout->addWidget(specular_color_button, 2, 0, 1, 1);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(3);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        specular_color_r = new QLineEdit(layoutWidget);
        specular_color_r->setObjectName(QString::fromUtf8("specular_color_r"));
        specular_color_r->setMinimumSize(QSize(40, 0));
        specular_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(specular_color_r);

        specular_color_g = new QLineEdit(layoutWidget);
        specular_color_g->setObjectName(QString::fromUtf8("specular_color_g"));
        specular_color_g->setMinimumSize(QSize(40, 0));
        specular_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(specular_color_g);

        specular_color_b = new QLineEdit(layoutWidget);
        specular_color_b->setObjectName(QString::fromUtf8("specular_color_b"));
        specular_color_b->setMinimumSize(QSize(40, 0));
        specular_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(specular_color_b);


        gridLayout->addLayout(horizontalLayout_5, 2, 1, 1, 1);


        uniform_color_Layout->addLayout(gridLayout);

        transparent_checkBox = new QCheckBox(layoutWidget);
        transparent_checkBox->setObjectName(QString::fromUtf8("transparent_checkBox"));

        uniform_color_Layout->addWidget(transparent_checkBox);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        uniform_color_Layout->addItem(verticalSpacer);


        gridLayout_2->addLayout(uniform_color_Layout, 0, 1, 1, 1);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        shinness_lineEdit = new QLineEdit(layoutWidget);
        shinness_lineEdit->setObjectName(QString::fromUtf8("shinness_lineEdit"));
        shinness_lineEdit->setMinimumSize(QSize(40, 0));
        shinness_lineEdit->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_8->addWidget(shinness_lineEdit);

        Shinness_label = new QLabel(layoutWidget);
        Shinness_label->setObjectName(QString::fromUtf8("Shinness_label"));
        Shinness_label->setMinimumSize(QSize(75, 0));
        Shinness_label->setMaximumSize(QSize(75, 16777215));
        Shinness_label->setFont(font);
        Shinness_label->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        horizontalLayout_8->addWidget(Shinness_label);

        horizontalSpacer = new QSpacerItem(80, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer);


        gridLayout_2->addLayout(horizontalLayout_8, 1, 1, 1, 1);

        ShinnessSlider = new QSlider(layoutWidget);
        ShinnessSlider->setObjectName(QString::fromUtf8("ShinnessSlider"));
        ShinnessSlider->setOrientation(Qt::Horizontal);

        gridLayout_2->addWidget(ShinnessSlider, 1, 0, 1, 1);


        retranslateUi(MaterialEditWidget);

        QMetaObject::connectSlotsByName(MaterialEditWidget);
    } // setupUi

    void retranslateUi(QWidget *MaterialEditWidget)
    {
        MaterialEditWidget->setWindowTitle(QCoreApplication::translate("MaterialEditWidget", "Form", nullptr));
        texture_label->setText(QCoreApplication::translate("MaterialEditWidget", "Texture", nullptr));
        ambient_texture_button->setText(QCoreApplication::translate("MaterialEditWidget", "Ambient", nullptr));
        enable_amibent_texture_checkBox->setText(QString());
        diffuse_texture_button->setText(QCoreApplication::translate("MaterialEditWidget", "Diffuse", nullptr));
        enable_diffuse_texture_checkBox->setText(QString());
        specular_texture_button->setText(QCoreApplication::translate("MaterialEditWidget", "Specular", nullptr));
        enable_specular_texture_checkBox->setText(QString());
        bump_texture_button->setText(QCoreApplication::translate("MaterialEditWidget", "Bump", nullptr));
        enable_bump_texture_checkBox->setText(QString());
        material_list_label->setText(QCoreApplication::translate("MaterialEditWidget", "Materials List", nullptr));
        uniform_color_label->setText(QCoreApplication::translate("MaterialEditWidget", "Uniform Color", nullptr));
        ambient_color_button->setText(QCoreApplication::translate("MaterialEditWidget", "Ambient", nullptr));
        diffuse_color_button->setText(QCoreApplication::translate("MaterialEditWidget", "Diffuse", nullptr));
        specular_color_button->setText(QCoreApplication::translate("MaterialEditWidget", "Specular", nullptr));
        transparent_checkBox->setText(QCoreApplication::translate("MaterialEditWidget", "Transparent", nullptr));
        Shinness_label->setText(QCoreApplication::translate("MaterialEditWidget", "Shiness", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MaterialEditWidget: public Ui_MaterialEditWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MATERIALEDITWIDGET_H
