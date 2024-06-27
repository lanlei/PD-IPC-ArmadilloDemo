/********************************************************************************
** Form generated from reading UI file 'RightFormWidget.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RIGHTFORMWIDGET_H
#define UI_RIGHTFORMWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_RightWidget
{
public:
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *Scene;
    QGroupBox *camera_groupBox;
    QVBoxLayout *verticalLayout;
    QCheckBox *ortho_mode_checkBox;
    QCheckBox *lock_camera_checkBox;
    QCheckBox *show_FOV_checkBox;
    QGroupBox *lightting_groupBox;
    QGridLayout *gridLayout_3;
    QPushButton *light_specular_color_button;
    QPushButton *light_ambient_color_button;
    QPushButton *light_diffuse_color_button;
    QLabel *light_pos_lable;
    QCheckBox *shadow_checkBox;
    QPushButton *add_light_button;
    QHBoxLayout *horizontalLayout_4;
    QLineEdit *light_diffuse_color_r;
    QLineEdit *light_diffuse_color_g;
    QLineEdit *light_diffuse_color_b;
    QHBoxLayout *horizontalLayout_5;
    QLineEdit *light_specular_color_r;
    QLineEdit *light_specular_color_g;
    QLineEdit *light_specular_color_b;
    QHBoxLayout *horizontalLayout;
    QLineEdit *light_pos_x;
    QLineEdit *light_pos_y;
    QLineEdit *light_pos_z;
    QHBoxLayout *horizontalLayout_3;
    QLineEdit *light_amibent_color_r;
    QLineEdit *light_amibent_color_g;
    QLineEdit *light_amibent_color_b;
    QLabel *Attenuation_pos_lable;
    QHBoxLayout *horizontalLayout_7;
    QLineEdit *light_attenuation_constant;
    QLineEdit *light_attenuation_linear;
    QLineEdit *light_attenuation_quadratic;
    QGroupBox *floor_groupBox;
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout_2;
    QLabel *floor_XZ_scale_label;
    QPushButton *floor_metrial_button;
    QHBoxLayout *horizontalLayout_6;
    QLineEdit *floor_xz_scale;
    QLabel *floor_plane_label;
    QLineEdit *floor_y_plane;
    QCheckBox *floor_hide_checkBox;
    QGroupBox *physical_environment_groupBox;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout_3;
    QGridLayout *gridLayout_8;
    QLabel *gravity_lable;
    QHBoxLayout *horizontalLayout_19;
    QLineEdit *gravity_lineEdit;
    QCheckBox *gravity_checkBox;
    QPushButton *add_extra_force_button;
    QLabel *time_step_lable;
    QLineEdit *time_step_lineEdit;
    QPushButton *setting_list_pushButton;
    QWidget *Model;
    QGroupBox *floor_groupBox_2;
    QGroupBox *floor_groupBox_3;
    QWidget *layoutWidget1;
    QGridLayout *gridLayout_5;
    QHBoxLayout *horizontalLayout_8;
    QPushButton *mesh_metrial_button_2;
    QCheckBox *mesh_hide_checkBox_2;
    QLabel *mesh_scale_label_2;
    QLineEdit *mesh_scale_lineEdit_2;
    QLabel *mesh_translate_label_2;
    QHBoxLayout *horizontalLayout_11;
    QLineEdit *mesh_translation_x_2;
    QLineEdit *mesh_translation_y_2;
    QLineEdit *mesh_translation_z_2;
    QLabel *mesh_rotate_label_2;
    QHBoxLayout *horizontalLayout_14;
    QSpacerItem *horizontalSpacer_2;
    QLineEdit *mesh_rotation_x_2;
    QLineEdit *mesh_rotation_y_2;
    QLineEdit *mesh_rotation_z_2;
    QWidget *layoutWidget2;
    QGridLayout *gridLayout_4;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *mesh_material_button;
    QCheckBox *mesh_hide_checkBox;
    QLabel *mesh_scale_label;
    QHBoxLayout *horizontalLayout_12;
    QLineEdit *mesh_scale_x;
    QLineEdit *mesh_scale_y;
    QLineEdit *mesh_scale_z;
    QSpacerItem *horizontalSpacer_6;
    QLabel *mesh_translate_label;
    QHBoxLayout *horizontalLayout_10;
    QLineEdit *mesh_translation_x;
    QLineEdit *mesh_translation_y;
    QLineEdit *mesh_translation_z;
    QSpacerItem *horizontalSpacer_5;
    QLabel *mesh_rotate_label;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_13;
    QLineEdit *mesh_rotation_x;
    QLineEdit *mesh_rotation_y;
    QLineEdit *mesh_rotation_z;
    QLineEdit *mesh_rotation_sita;
    QListWidget *mesh_listWidget;
    QGroupBox *floor_groupBox_4;
    QWidget *layoutWidget3;
    QHBoxLayout *horizontalLayout_20;
    QCheckBox *medial_mesh_hide_checkBox;
    QCheckBox *medial_mesh_transparent_checkBox;
    QListWidget *medial_mesh_set_listWidget;
    QGroupBox *floor_groupBox_6;
    QWidget *layoutWidget_4;
    QHBoxLayout *horizontalLayout_21;
    QCheckBox *tet_mesh_hide_checkBox;
    QCheckBox *tet_mesh_transparent_checkBox;
    QListWidget *element_set_listWidget;

    void setupUi(QWidget *RightWidget)
    {
        if (RightWidget->objectName().isEmpty())
            RightWidget->setObjectName(QString::fromUtf8("RightWidget"));
        RightWidget->resize(280, 904);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(RightWidget->sizePolicy().hasHeightForWidth());
        RightWidget->setSizePolicy(sizePolicy);
        RightWidget->setMinimumSize(QSize(265, 0));
        RightWidget->setMaximumSize(QSize(280, 16777215));
        gridLayout = new QGridLayout(RightWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tabWidget = new QTabWidget(RightWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(tabWidget->sizePolicy().hasHeightForWidth());
        tabWidget->setSizePolicy(sizePolicy1);
        tabWidget->setMinimumSize(QSize(0, 0));
        tabWidget->setMaximumSize(QSize(1677725, 1677725));
        QFont font;
        font.setFamily(QString::fromUtf8("Times New Roman"));
        font.setPointSize(12);
        tabWidget->setFont(font);
        tabWidget->setFocusPolicy(Qt::NoFocus);
        Scene = new QWidget();
        Scene->setObjectName(QString::fromUtf8("Scene"));
        camera_groupBox = new QGroupBox(Scene);
        camera_groupBox->setObjectName(QString::fromUtf8("camera_groupBox"));
        camera_groupBox->setGeometry(QRect(0, 0, 255, 108));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(camera_groupBox->sizePolicy().hasHeightForWidth());
        camera_groupBox->setSizePolicy(sizePolicy2);
        camera_groupBox->setMinimumSize(QSize(255, 0));
        camera_groupBox->setMaximumSize(QSize(255, 16777215));
        QFont font1;
        font1.setFamily(QString::fromUtf8("Times New Roman"));
        font1.setPointSize(12);
        font1.setStyleStrategy(QFont::PreferAntialias);
        camera_groupBox->setFont(font1);
        verticalLayout = new QVBoxLayout(camera_groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        ortho_mode_checkBox = new QCheckBox(camera_groupBox);
        ortho_mode_checkBox->setObjectName(QString::fromUtf8("ortho_mode_checkBox"));
        ortho_mode_checkBox->setEnabled(false);
        QFont font2;
        font2.setFamily(QString::fromUtf8("Times New Roman"));
        font2.setPointSize(10);
        font2.setStyleStrategy(QFont::PreferAntialias);
        ortho_mode_checkBox->setFont(font2);
        ortho_mode_checkBox->setLayoutDirection(Qt::LeftToRight);

        verticalLayout->addWidget(ortho_mode_checkBox);

        lock_camera_checkBox = new QCheckBox(camera_groupBox);
        lock_camera_checkBox->setObjectName(QString::fromUtf8("lock_camera_checkBox"));
        lock_camera_checkBox->setFont(font2);
        lock_camera_checkBox->setLayoutDirection(Qt::LeftToRight);

        verticalLayout->addWidget(lock_camera_checkBox);

        show_FOV_checkBox = new QCheckBox(camera_groupBox);
        show_FOV_checkBox->setObjectName(QString::fromUtf8("show_FOV_checkBox"));
        show_FOV_checkBox->setEnabled(false);
        show_FOV_checkBox->setFont(font2);
        show_FOV_checkBox->setLayoutDirection(Qt::LeftToRight);

        verticalLayout->addWidget(show_FOV_checkBox);

        lightting_groupBox = new QGroupBox(Scene);
        lightting_groupBox->setObjectName(QString::fromUtf8("lightting_groupBox"));
        lightting_groupBox->setGeometry(QRect(0, 120, 255, 231));
        sizePolicy2.setHeightForWidth(lightting_groupBox->sizePolicy().hasHeightForWidth());
        lightting_groupBox->setSizePolicy(sizePolicy2);
        lightting_groupBox->setMinimumSize(QSize(255, 0));
        lightting_groupBox->setMaximumSize(QSize(255, 16777215));
        lightting_groupBox->setFont(font1);
        gridLayout_3 = new QGridLayout(lightting_groupBox);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        light_specular_color_button = new QPushButton(lightting_groupBox);
        light_specular_color_button->setObjectName(QString::fromUtf8("light_specular_color_button"));
        light_specular_color_button->setMaximumSize(QSize(16777215, 16777215));
        light_specular_color_button->setFont(font2);

        gridLayout_3->addWidget(light_specular_color_button, 4, 0, 1, 1);

        light_ambient_color_button = new QPushButton(lightting_groupBox);
        light_ambient_color_button->setObjectName(QString::fromUtf8("light_ambient_color_button"));
        light_ambient_color_button->setMaximumSize(QSize(16777215, 16777215));
        light_ambient_color_button->setFont(font2);

        gridLayout_3->addWidget(light_ambient_color_button, 2, 0, 1, 1);

        light_diffuse_color_button = new QPushButton(lightting_groupBox);
        light_diffuse_color_button->setObjectName(QString::fromUtf8("light_diffuse_color_button"));
        light_diffuse_color_button->setMaximumSize(QSize(16777215, 16777215));
        light_diffuse_color_button->setFont(font2);

        gridLayout_3->addWidget(light_diffuse_color_button, 3, 0, 1, 1);

        light_pos_lable = new QLabel(lightting_groupBox);
        light_pos_lable->setObjectName(QString::fromUtf8("light_pos_lable"));
        light_pos_lable->setFont(font2);
        light_pos_lable->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(light_pos_lable, 0, 0, 1, 1);

        shadow_checkBox = new QCheckBox(lightting_groupBox);
        shadow_checkBox->setObjectName(QString::fromUtf8("shadow_checkBox"));
        shadow_checkBox->setFont(font2);
        shadow_checkBox->setLayoutDirection(Qt::LeftToRight);

        gridLayout_3->addWidget(shadow_checkBox, 5, 0, 1, 1);

        add_light_button = new QPushButton(lightting_groupBox);
        add_light_button->setObjectName(QString::fromUtf8("add_light_button"));
        add_light_button->setEnabled(false);
        add_light_button->setFont(font2);

        gridLayout_3->addWidget(add_light_button, 5, 1, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(3);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        light_diffuse_color_r = new QLineEdit(lightting_groupBox);
        light_diffuse_color_r->setObjectName(QString::fromUtf8("light_diffuse_color_r"));
        light_diffuse_color_r->setMinimumSize(QSize(40, 0));
        light_diffuse_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(light_diffuse_color_r);

        light_diffuse_color_g = new QLineEdit(lightting_groupBox);
        light_diffuse_color_g->setObjectName(QString::fromUtf8("light_diffuse_color_g"));
        light_diffuse_color_g->setMinimumSize(QSize(40, 0));
        light_diffuse_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(light_diffuse_color_g);

        light_diffuse_color_b = new QLineEdit(lightting_groupBox);
        light_diffuse_color_b->setObjectName(QString::fromUtf8("light_diffuse_color_b"));
        light_diffuse_color_b->setMinimumSize(QSize(40, 0));
        light_diffuse_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_4->addWidget(light_diffuse_color_b);


        gridLayout_3->addLayout(horizontalLayout_4, 3, 1, 1, 1);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(3);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        light_specular_color_r = new QLineEdit(lightting_groupBox);
        light_specular_color_r->setObjectName(QString::fromUtf8("light_specular_color_r"));
        light_specular_color_r->setMinimumSize(QSize(40, 0));
        light_specular_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(light_specular_color_r);

        light_specular_color_g = new QLineEdit(lightting_groupBox);
        light_specular_color_g->setObjectName(QString::fromUtf8("light_specular_color_g"));
        light_specular_color_g->setMinimumSize(QSize(40, 0));
        light_specular_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(light_specular_color_g);

        light_specular_color_b = new QLineEdit(lightting_groupBox);
        light_specular_color_b->setObjectName(QString::fromUtf8("light_specular_color_b"));
        light_specular_color_b->setMinimumSize(QSize(40, 0));
        light_specular_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_5->addWidget(light_specular_color_b);


        gridLayout_3->addLayout(horizontalLayout_5, 4, 1, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(3);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        light_pos_x = new QLineEdit(lightting_groupBox);
        light_pos_x->setObjectName(QString::fromUtf8("light_pos_x"));
        light_pos_x->setMinimumSize(QSize(40, 0));
        light_pos_x->setMaximumSize(QSize(40, 16777215));

        horizontalLayout->addWidget(light_pos_x);

        light_pos_y = new QLineEdit(lightting_groupBox);
        light_pos_y->setObjectName(QString::fromUtf8("light_pos_y"));
        light_pos_y->setMinimumSize(QSize(40, 0));
        light_pos_y->setMaximumSize(QSize(40, 16777215));

        horizontalLayout->addWidget(light_pos_y);

        light_pos_z = new QLineEdit(lightting_groupBox);
        light_pos_z->setObjectName(QString::fromUtf8("light_pos_z"));
        light_pos_z->setMinimumSize(QSize(40, 0));
        light_pos_z->setMaximumSize(QSize(40, 16777215));

        horizontalLayout->addWidget(light_pos_z);


        gridLayout_3->addLayout(horizontalLayout, 0, 1, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(3);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        light_amibent_color_r = new QLineEdit(lightting_groupBox);
        light_amibent_color_r->setObjectName(QString::fromUtf8("light_amibent_color_r"));
        light_amibent_color_r->setMinimumSize(QSize(40, 0));
        light_amibent_color_r->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(light_amibent_color_r);

        light_amibent_color_g = new QLineEdit(lightting_groupBox);
        light_amibent_color_g->setObjectName(QString::fromUtf8("light_amibent_color_g"));
        light_amibent_color_g->setMinimumSize(QSize(40, 0));
        light_amibent_color_g->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(light_amibent_color_g);

        light_amibent_color_b = new QLineEdit(lightting_groupBox);
        light_amibent_color_b->setObjectName(QString::fromUtf8("light_amibent_color_b"));
        light_amibent_color_b->setMinimumSize(QSize(40, 0));
        light_amibent_color_b->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_3->addWidget(light_amibent_color_b);


        gridLayout_3->addLayout(horizontalLayout_3, 2, 1, 1, 1);

        Attenuation_pos_lable = new QLabel(lightting_groupBox);
        Attenuation_pos_lable->setObjectName(QString::fromUtf8("Attenuation_pos_lable"));
        Attenuation_pos_lable->setFont(font2);
        Attenuation_pos_lable->setAlignment(Qt::AlignCenter);

        gridLayout_3->addWidget(Attenuation_pos_lable, 1, 0, 1, 1);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(3);
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        light_attenuation_constant = new QLineEdit(lightting_groupBox);
        light_attenuation_constant->setObjectName(QString::fromUtf8("light_attenuation_constant"));
        light_attenuation_constant->setMinimumSize(QSize(40, 0));
        light_attenuation_constant->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_7->addWidget(light_attenuation_constant);

        light_attenuation_linear = new QLineEdit(lightting_groupBox);
        light_attenuation_linear->setObjectName(QString::fromUtf8("light_attenuation_linear"));
        light_attenuation_linear->setMinimumSize(QSize(40, 0));
        light_attenuation_linear->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_7->addWidget(light_attenuation_linear);

        light_attenuation_quadratic = new QLineEdit(lightting_groupBox);
        light_attenuation_quadratic->setObjectName(QString::fromUtf8("light_attenuation_quadratic"));
        light_attenuation_quadratic->setMinimumSize(QSize(40, 0));
        light_attenuation_quadratic->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_7->addWidget(light_attenuation_quadratic);


        gridLayout_3->addLayout(horizontalLayout_7, 1, 1, 1, 1);

        floor_groupBox = new QGroupBox(Scene);
        floor_groupBox->setObjectName(QString::fromUtf8("floor_groupBox"));
        floor_groupBox->setGeometry(QRect(0, 370, 255, 101));
        sizePolicy2.setHeightForWidth(floor_groupBox->sizePolicy().hasHeightForWidth());
        floor_groupBox->setSizePolicy(sizePolicy2);
        floor_groupBox->setMinimumSize(QSize(255, 0));
        floor_groupBox->setMaximumSize(QSize(255, 16777215));
        QFont font3;
        font3.setPointSize(12);
        floor_groupBox->setFont(font3);
        verticalLayout_2 = new QVBoxLayout(floor_groupBox);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        floor_XZ_scale_label = new QLabel(floor_groupBox);
        floor_XZ_scale_label->setObjectName(QString::fromUtf8("floor_XZ_scale_label"));
        QFont font4;
        font4.setPointSize(10);
        floor_XZ_scale_label->setFont(font4);

        gridLayout_2->addWidget(floor_XZ_scale_label, 1, 0, 1, 1);

        floor_metrial_button = new QPushButton(floor_groupBox);
        floor_metrial_button->setObjectName(QString::fromUtf8("floor_metrial_button"));
        QSizePolicy sizePolicy3(QSizePolicy::Minimum, QSizePolicy::Expanding);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(floor_metrial_button->sizePolicy().hasHeightForWidth());
        floor_metrial_button->setSizePolicy(sizePolicy3);
        floor_metrial_button->setMinimumSize(QSize(75, 23));
        floor_metrial_button->setMaximumSize(QSize(75, 23));
        floor_metrial_button->setFont(font4);

        gridLayout_2->addWidget(floor_metrial_button, 0, 0, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(1);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        floor_xz_scale = new QLineEdit(floor_groupBox);
        floor_xz_scale->setObjectName(QString::fromUtf8("floor_xz_scale"));
        floor_xz_scale->setMinimumSize(QSize(40, 0));
        floor_xz_scale->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_6->addWidget(floor_xz_scale);

        floor_plane_label = new QLabel(floor_groupBox);
        floor_plane_label->setObjectName(QString::fromUtf8("floor_plane_label"));
        floor_plane_label->setFont(font4);

        horizontalLayout_6->addWidget(floor_plane_label);

        floor_y_plane = new QLineEdit(floor_groupBox);
        floor_y_plane->setObjectName(QString::fromUtf8("floor_y_plane"));
        floor_y_plane->setMinimumSize(QSize(40, 0));
        floor_y_plane->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_6->addWidget(floor_y_plane);


        gridLayout_2->addLayout(horizontalLayout_6, 1, 1, 1, 1);

        floor_hide_checkBox = new QCheckBox(floor_groupBox);
        floor_hide_checkBox->setObjectName(QString::fromUtf8("floor_hide_checkBox"));
        floor_hide_checkBox->setFont(font4);
        floor_hide_checkBox->setLayoutDirection(Qt::LeftToRight);

        gridLayout_2->addWidget(floor_hide_checkBox, 0, 1, 1, 1);


        verticalLayout_2->addLayout(gridLayout_2);

        physical_environment_groupBox = new QGroupBox(Scene);
        physical_environment_groupBox->setObjectName(QString::fromUtf8("physical_environment_groupBox"));
        physical_environment_groupBox->setGeometry(QRect(0, 490, 255, 141));
        sizePolicy2.setHeightForWidth(physical_environment_groupBox->sizePolicy().hasHeightForWidth());
        physical_environment_groupBox->setSizePolicy(sizePolicy2);
        physical_environment_groupBox->setMinimumSize(QSize(255, 0));
        physical_environment_groupBox->setMaximumSize(QSize(255, 16777215));
        physical_environment_groupBox->setFont(font3);
        layoutWidget = new QWidget(physical_environment_groupBox);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(12, 30, 224, 95));
        verticalLayout_3 = new QVBoxLayout(layoutWidget);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        gridLayout_8 = new QGridLayout();
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        gravity_lable = new QLabel(layoutWidget);
        gravity_lable->setObjectName(QString::fromUtf8("gravity_lable"));
        gravity_lable->setFont(font4);
        gravity_lable->setAlignment(Qt::AlignCenter);

        gridLayout_8->addWidget(gravity_lable, 0, 0, 1, 1);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setObjectName(QString::fromUtf8("horizontalLayout_19"));
        gravity_lineEdit = new QLineEdit(layoutWidget);
        gravity_lineEdit->setObjectName(QString::fromUtf8("gravity_lineEdit"));
        gravity_lineEdit->setMinimumSize(QSize(40, 0));
        gravity_lineEdit->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_19->addWidget(gravity_lineEdit);

        gravity_checkBox = new QCheckBox(layoutWidget);
        gravity_checkBox->setObjectName(QString::fromUtf8("gravity_checkBox"));
        gravity_checkBox->setFont(font4);
        gravity_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_19->addWidget(gravity_checkBox);


        gridLayout_8->addLayout(horizontalLayout_19, 0, 1, 1, 1);

        add_extra_force_button = new QPushButton(layoutWidget);
        add_extra_force_button->setObjectName(QString::fromUtf8("add_extra_force_button"));
        add_extra_force_button->setEnabled(true);
        add_extra_force_button->setFont(font4);

        gridLayout_8->addWidget(add_extra_force_button, 0, 2, 1, 1);

        time_step_lable = new QLabel(layoutWidget);
        time_step_lable->setObjectName(QString::fromUtf8("time_step_lable"));
        time_step_lable->setFont(font4);
        time_step_lable->setAlignment(Qt::AlignCenter);

        gridLayout_8->addWidget(time_step_lable, 1, 0, 1, 1);

        time_step_lineEdit = new QLineEdit(layoutWidget);
        time_step_lineEdit->setObjectName(QString::fromUtf8("time_step_lineEdit"));
        time_step_lineEdit->setMinimumSize(QSize(40, 0));
        time_step_lineEdit->setMaximumSize(QSize(40, 16777215));

        gridLayout_8->addWidget(time_step_lineEdit, 1, 1, 1, 1);


        verticalLayout_3->addLayout(gridLayout_8);

        setting_list_pushButton = new QPushButton(layoutWidget);
        setting_list_pushButton->setObjectName(QString::fromUtf8("setting_list_pushButton"));

        verticalLayout_3->addWidget(setting_list_pushButton);

        tabWidget->addTab(Scene, QString());
        Model = new QWidget();
        Model->setObjectName(QString::fromUtf8("Model"));
        Model->setContextMenuPolicy(Qt::NoContextMenu);
        floor_groupBox_2 = new QGroupBox(Model);
        floor_groupBox_2->setObjectName(QString::fromUtf8("floor_groupBox_2"));
        floor_groupBox_2->setGeometry(QRect(0, 140, 255, 161));
        sizePolicy2.setHeightForWidth(floor_groupBox_2->sizePolicy().hasHeightForWidth());
        floor_groupBox_2->setSizePolicy(sizePolicy2);
        floor_groupBox_2->setMinimumSize(QSize(255, 0));
        floor_groupBox_2->setMaximumSize(QSize(255, 16777215));
        floor_groupBox_2->setFont(font3);
        floor_groupBox_3 = new QGroupBox(floor_groupBox_2);
        floor_groupBox_3->setObjectName(QString::fromUtf8("floor_groupBox_3"));
        floor_groupBox_3->setGeometry(QRect(220, 160, 255, 161));
        sizePolicy2.setHeightForWidth(floor_groupBox_3->sizePolicy().hasHeightForWidth());
        floor_groupBox_3->setSizePolicy(sizePolicy2);
        floor_groupBox_3->setMinimumSize(QSize(255, 0));
        floor_groupBox_3->setMaximumSize(QSize(255, 16777215));
        floor_groupBox_3->setFont(font3);
        layoutWidget1 = new QWidget(floor_groupBox_3);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(12, 23, 186, 124));
        gridLayout_5 = new QGridLayout(layoutWidget1);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        gridLayout_5->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        mesh_metrial_button_2 = new QPushButton(layoutWidget1);
        mesh_metrial_button_2->setObjectName(QString::fromUtf8("mesh_metrial_button_2"));
        sizePolicy3.setHeightForWidth(mesh_metrial_button_2->sizePolicy().hasHeightForWidth());
        mesh_metrial_button_2->setSizePolicy(sizePolicy3);
        mesh_metrial_button_2->setMinimumSize(QSize(75, 23));
        mesh_metrial_button_2->setMaximumSize(QSize(75, 23));
        mesh_metrial_button_2->setFont(font4);

        horizontalLayout_8->addWidget(mesh_metrial_button_2);

        mesh_hide_checkBox_2 = new QCheckBox(layoutWidget1);
        mesh_hide_checkBox_2->setObjectName(QString::fromUtf8("mesh_hide_checkBox_2"));
        mesh_hide_checkBox_2->setFont(font4);
        mesh_hide_checkBox_2->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_8->addWidget(mesh_hide_checkBox_2);


        gridLayout_5->addLayout(horizontalLayout_8, 0, 0, 1, 3);

        mesh_scale_label_2 = new QLabel(layoutWidget1);
        mesh_scale_label_2->setObjectName(QString::fromUtf8("mesh_scale_label_2"));
        mesh_scale_label_2->setFont(font4);
        mesh_scale_label_2->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(mesh_scale_label_2, 1, 0, 1, 2);

        mesh_scale_lineEdit_2 = new QLineEdit(layoutWidget1);
        mesh_scale_lineEdit_2->setObjectName(QString::fromUtf8("mesh_scale_lineEdit_2"));
        mesh_scale_lineEdit_2->setMinimumSize(QSize(40, 0));
        mesh_scale_lineEdit_2->setMaximumSize(QSize(40, 16777215));

        gridLayout_5->addWidget(mesh_scale_lineEdit_2, 1, 2, 1, 1);

        mesh_translate_label_2 = new QLabel(layoutWidget1);
        mesh_translate_label_2->setObjectName(QString::fromUtf8("mesh_translate_label_2"));
        mesh_translate_label_2->setFont(font4);
        mesh_translate_label_2->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(mesh_translate_label_2, 2, 0, 1, 2);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setSpacing(3);
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        mesh_translation_x_2 = new QLineEdit(layoutWidget1);
        mesh_translation_x_2->setObjectName(QString::fromUtf8("mesh_translation_x_2"));
        mesh_translation_x_2->setMinimumSize(QSize(40, 0));
        mesh_translation_x_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_11->addWidget(mesh_translation_x_2);

        mesh_translation_y_2 = new QLineEdit(layoutWidget1);
        mesh_translation_y_2->setObjectName(QString::fromUtf8("mesh_translation_y_2"));
        mesh_translation_y_2->setMinimumSize(QSize(40, 0));
        mesh_translation_y_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_11->addWidget(mesh_translation_y_2);

        mesh_translation_z_2 = new QLineEdit(layoutWidget1);
        mesh_translation_z_2->setObjectName(QString::fromUtf8("mesh_translation_z_2"));
        mesh_translation_z_2->setMinimumSize(QSize(40, 0));
        mesh_translation_z_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_11->addWidget(mesh_translation_z_2);


        gridLayout_5->addLayout(horizontalLayout_11, 2, 2, 1, 1);

        mesh_rotate_label_2 = new QLabel(layoutWidget1);
        mesh_rotate_label_2->setObjectName(QString::fromUtf8("mesh_rotate_label_2"));
        mesh_rotate_label_2->setFont(font4);
        mesh_rotate_label_2->setAlignment(Qt::AlignCenter);

        gridLayout_5->addWidget(mesh_rotate_label_2, 3, 0, 1, 1);

        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setSpacing(3);
        horizontalLayout_14->setObjectName(QString::fromUtf8("horizontalLayout_14"));
        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_14->addItem(horizontalSpacer_2);

        mesh_rotation_x_2 = new QLineEdit(layoutWidget1);
        mesh_rotation_x_2->setObjectName(QString::fromUtf8("mesh_rotation_x_2"));
        mesh_rotation_x_2->setMinimumSize(QSize(40, 0));
        mesh_rotation_x_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_14->addWidget(mesh_rotation_x_2);

        mesh_rotation_y_2 = new QLineEdit(layoutWidget1);
        mesh_rotation_y_2->setObjectName(QString::fromUtf8("mesh_rotation_y_2"));
        mesh_rotation_y_2->setMinimumSize(QSize(40, 0));
        mesh_rotation_y_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_14->addWidget(mesh_rotation_y_2);

        mesh_rotation_z_2 = new QLineEdit(layoutWidget1);
        mesh_rotation_z_2->setObjectName(QString::fromUtf8("mesh_rotation_z_2"));
        mesh_rotation_z_2->setMinimumSize(QSize(40, 0));
        mesh_rotation_z_2->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_14->addWidget(mesh_rotation_z_2);


        gridLayout_5->addLayout(horizontalLayout_14, 3, 1, 1, 2);

        layoutWidget2 = new QWidget(floor_groupBox_2);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(11, 24, 241, 126));
        gridLayout_4 = new QGridLayout(layoutWidget2);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        gridLayout_4->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        mesh_material_button = new QPushButton(layoutWidget2);
        mesh_material_button->setObjectName(QString::fromUtf8("mesh_material_button"));
        sizePolicy3.setHeightForWidth(mesh_material_button->sizePolicy().hasHeightForWidth());
        mesh_material_button->setSizePolicy(sizePolicy3);
        mesh_material_button->setMinimumSize(QSize(75, 23));
        mesh_material_button->setMaximumSize(QSize(75, 23));
        mesh_material_button->setFont(font4);

        horizontalLayout_2->addWidget(mesh_material_button);

        mesh_hide_checkBox = new QCheckBox(layoutWidget2);
        mesh_hide_checkBox->setObjectName(QString::fromUtf8("mesh_hide_checkBox"));
        mesh_hide_checkBox->setFont(font4);
        mesh_hide_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_2->addWidget(mesh_hide_checkBox);


        gridLayout_4->addLayout(horizontalLayout_2, 0, 0, 1, 3);

        mesh_scale_label = new QLabel(layoutWidget2);
        mesh_scale_label->setObjectName(QString::fromUtf8("mesh_scale_label"));
        mesh_scale_label->setFont(font4);
        mesh_scale_label->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(mesh_scale_label, 1, 0, 1, 1);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setSpacing(3);
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        mesh_scale_x = new QLineEdit(layoutWidget2);
        mesh_scale_x->setObjectName(QString::fromUtf8("mesh_scale_x"));
        mesh_scale_x->setMinimumSize(QSize(40, 0));
        mesh_scale_x->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_12->addWidget(mesh_scale_x);

        mesh_scale_y = new QLineEdit(layoutWidget2);
        mesh_scale_y->setObjectName(QString::fromUtf8("mesh_scale_y"));
        mesh_scale_y->setMinimumSize(QSize(40, 0));
        mesh_scale_y->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_12->addWidget(mesh_scale_y);

        mesh_scale_z = new QLineEdit(layoutWidget2);
        mesh_scale_z->setObjectName(QString::fromUtf8("mesh_scale_z"));
        mesh_scale_z->setMinimumSize(QSize(40, 0));
        mesh_scale_z->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_12->addWidget(mesh_scale_z);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_12->addItem(horizontalSpacer_6);


        gridLayout_4->addLayout(horizontalLayout_12, 1, 2, 1, 1);

        mesh_translate_label = new QLabel(layoutWidget2);
        mesh_translate_label->setObjectName(QString::fromUtf8("mesh_translate_label"));
        mesh_translate_label->setFont(font4);
        mesh_translate_label->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(mesh_translate_label, 2, 0, 1, 1);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setSpacing(3);
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        mesh_translation_x = new QLineEdit(layoutWidget2);
        mesh_translation_x->setObjectName(QString::fromUtf8("mesh_translation_x"));
        mesh_translation_x->setMinimumSize(QSize(40, 0));
        mesh_translation_x->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_10->addWidget(mesh_translation_x);

        mesh_translation_y = new QLineEdit(layoutWidget2);
        mesh_translation_y->setObjectName(QString::fromUtf8("mesh_translation_y"));
        mesh_translation_y->setMinimumSize(QSize(40, 0));
        mesh_translation_y->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_10->addWidget(mesh_translation_y);

        mesh_translation_z = new QLineEdit(layoutWidget2);
        mesh_translation_z->setObjectName(QString::fromUtf8("mesh_translation_z"));
        mesh_translation_z->setMinimumSize(QSize(40, 0));
        mesh_translation_z->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_10->addWidget(mesh_translation_z);

        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_10->addItem(horizontalSpacer_5);


        gridLayout_4->addLayout(horizontalLayout_10, 2, 2, 1, 1);

        mesh_rotate_label = new QLabel(layoutWidget2);
        mesh_rotate_label->setObjectName(QString::fromUtf8("mesh_rotate_label"));
        mesh_rotate_label->setFont(font4);
        mesh_rotate_label->setAlignment(Qt::AlignCenter);

        gridLayout_4->addWidget(mesh_rotate_label, 3, 0, 1, 1);

        horizontalSpacer = new QSpacerItem(13, 24, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer, 3, 1, 1, 1);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setSpacing(3);
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        mesh_rotation_x = new QLineEdit(layoutWidget2);
        mesh_rotation_x->setObjectName(QString::fromUtf8("mesh_rotation_x"));
        mesh_rotation_x->setMinimumSize(QSize(40, 0));
        mesh_rotation_x->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_13->addWidget(mesh_rotation_x);

        mesh_rotation_y = new QLineEdit(layoutWidget2);
        mesh_rotation_y->setObjectName(QString::fromUtf8("mesh_rotation_y"));
        mesh_rotation_y->setMinimumSize(QSize(40, 0));
        mesh_rotation_y->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_13->addWidget(mesh_rotation_y);

        mesh_rotation_z = new QLineEdit(layoutWidget2);
        mesh_rotation_z->setObjectName(QString::fromUtf8("mesh_rotation_z"));
        mesh_rotation_z->setMinimumSize(QSize(40, 0));
        mesh_rotation_z->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_13->addWidget(mesh_rotation_z);

        mesh_rotation_sita = new QLineEdit(layoutWidget2);
        mesh_rotation_sita->setObjectName(QString::fromUtf8("mesh_rotation_sita"));
        mesh_rotation_sita->setMinimumSize(QSize(40, 0));
        mesh_rotation_sita->setMaximumSize(QSize(40, 16777215));

        horizontalLayout_13->addWidget(mesh_rotation_sita);


        gridLayout_4->addLayout(horizontalLayout_13, 3, 2, 1, 1);

        mesh_listWidget = new QListWidget(Model);
        mesh_listWidget->setObjectName(QString::fromUtf8("mesh_listWidget"));
        mesh_listWidget->setGeometry(QRect(10, 10, 241, 121));
        floor_groupBox_4 = new QGroupBox(Model);
        floor_groupBox_4->setObjectName(QString::fromUtf8("floor_groupBox_4"));
        floor_groupBox_4->setGeometry(QRect(0, 500, 255, 181));
        sizePolicy2.setHeightForWidth(floor_groupBox_4->sizePolicy().hasHeightForWidth());
        floor_groupBox_4->setSizePolicy(sizePolicy2);
        floor_groupBox_4->setMinimumSize(QSize(255, 0));
        floor_groupBox_4->setMaximumSize(QSize(255, 16777215));
        floor_groupBox_4->setFont(font3);
        layoutWidget3 = new QWidget(floor_groupBox_4);
        layoutWidget3->setObjectName(QString::fromUtf8("layoutWidget3"));
        layoutWidget3->setGeometry(QRect(10, 20, 144, 21));
        horizontalLayout_20 = new QHBoxLayout(layoutWidget3);
        horizontalLayout_20->setObjectName(QString::fromUtf8("horizontalLayout_20"));
        horizontalLayout_20->setContentsMargins(0, 0, 0, 0);
        medial_mesh_hide_checkBox = new QCheckBox(layoutWidget3);
        medial_mesh_hide_checkBox->setObjectName(QString::fromUtf8("medial_mesh_hide_checkBox"));
        medial_mesh_hide_checkBox->setEnabled(false);
        medial_mesh_hide_checkBox->setFont(font4);
        medial_mesh_hide_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_20->addWidget(medial_mesh_hide_checkBox);

        medial_mesh_transparent_checkBox = new QCheckBox(layoutWidget3);
        medial_mesh_transparent_checkBox->setObjectName(QString::fromUtf8("medial_mesh_transparent_checkBox"));
        medial_mesh_transparent_checkBox->setEnabled(false);
        medial_mesh_transparent_checkBox->setFont(font4);
        medial_mesh_transparent_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_20->addWidget(medial_mesh_transparent_checkBox);

        medial_mesh_set_listWidget = new QListWidget(floor_groupBox_4);
        medial_mesh_set_listWidget->setObjectName(QString::fromUtf8("medial_mesh_set_listWidget"));
        medial_mesh_set_listWidget->setGeometry(QRect(10, 50, 241, 121));
        floor_groupBox_6 = new QGroupBox(Model);
        floor_groupBox_6->setObjectName(QString::fromUtf8("floor_groupBox_6"));
        floor_groupBox_6->setGeometry(QRect(0, 310, 255, 181));
        sizePolicy2.setHeightForWidth(floor_groupBox_6->sizePolicy().hasHeightForWidth());
        floor_groupBox_6->setSizePolicy(sizePolicy2);
        floor_groupBox_6->setMinimumSize(QSize(255, 0));
        floor_groupBox_6->setMaximumSize(QSize(255, 16777215));
        floor_groupBox_6->setFont(font3);
        layoutWidget_4 = new QWidget(floor_groupBox_6);
        layoutWidget_4->setObjectName(QString::fromUtf8("layoutWidget_4"));
        layoutWidget_4->setGeometry(QRect(10, 20, 144, 21));
        horizontalLayout_21 = new QHBoxLayout(layoutWidget_4);
        horizontalLayout_21->setObjectName(QString::fromUtf8("horizontalLayout_21"));
        horizontalLayout_21->setContentsMargins(0, 0, 0, 0);
        tet_mesh_hide_checkBox = new QCheckBox(layoutWidget_4);
        tet_mesh_hide_checkBox->setObjectName(QString::fromUtf8("tet_mesh_hide_checkBox"));
        tet_mesh_hide_checkBox->setEnabled(false);
        tet_mesh_hide_checkBox->setFont(font4);
        tet_mesh_hide_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_21->addWidget(tet_mesh_hide_checkBox);

        tet_mesh_transparent_checkBox = new QCheckBox(layoutWidget_4);
        tet_mesh_transparent_checkBox->setObjectName(QString::fromUtf8("tet_mesh_transparent_checkBox"));
        tet_mesh_transparent_checkBox->setEnabled(false);
        tet_mesh_transparent_checkBox->setFont(font4);
        tet_mesh_transparent_checkBox->setLayoutDirection(Qt::LeftToRight);

        horizontalLayout_21->addWidget(tet_mesh_transparent_checkBox);

        element_set_listWidget = new QListWidget(floor_groupBox_6);
        element_set_listWidget->setObjectName(QString::fromUtf8("element_set_listWidget"));
        element_set_listWidget->setGeometry(QRect(10, 50, 241, 121));
        tabWidget->addTab(Model, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(RightWidget);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(RightWidget);
    } // setupUi

    void retranslateUi(QWidget *RightWidget)
    {
        RightWidget->setWindowTitle(QCoreApplication::translate("RightWidget", "Form", nullptr));
        camera_groupBox->setTitle(QCoreApplication::translate("RightWidget", "Camera", nullptr));
        ortho_mode_checkBox->setText(QCoreApplication::translate("RightWidget", "Ortho Mode", nullptr));
        lock_camera_checkBox->setText(QCoreApplication::translate("RightWidget", "Lock Camera", nullptr));
        show_FOV_checkBox->setText(QCoreApplication::translate("RightWidget", "Show FOV", nullptr));
        lightting_groupBox->setTitle(QCoreApplication::translate("RightWidget", "Lighting", nullptr));
        light_specular_color_button->setText(QCoreApplication::translate("RightWidget", "Specular", nullptr));
        light_ambient_color_button->setText(QCoreApplication::translate("RightWidget", "Ambient", nullptr));
        light_diffuse_color_button->setText(QCoreApplication::translate("RightWidget", "Diffuse", nullptr));
        light_pos_lable->setText(QCoreApplication::translate("RightWidget", "          Pos", nullptr));
        shadow_checkBox->setText(QCoreApplication::translate("RightWidget", "Shadow", nullptr));
        add_light_button->setText(QCoreApplication::translate("RightWidget", "Add Light", nullptr));
        Attenuation_pos_lable->setText(QCoreApplication::translate("RightWidget", "Attenuation", nullptr));
        floor_groupBox->setTitle(QCoreApplication::translate("RightWidget", "Floor", nullptr));
        floor_XZ_scale_label->setText(QCoreApplication::translate("RightWidget", "      XZ-Scale", nullptr));
        floor_metrial_button->setText(QCoreApplication::translate("RightWidget", "Render", nullptr));
        floor_plane_label->setText(QCoreApplication::translate("RightWidget", "      Y-Plane", nullptr));
        floor_hide_checkBox->setText(QCoreApplication::translate("RightWidget", "Hide", nullptr));
        physical_environment_groupBox->setTitle(QCoreApplication::translate("RightWidget", "Environment", nullptr));
        gravity_lable->setText(QCoreApplication::translate("RightWidget", "Gravity", nullptr));
        gravity_checkBox->setText(QString());
        add_extra_force_button->setText(QCoreApplication::translate("RightWidget", "Add Extra Force", nullptr));
        time_step_lable->setText(QCoreApplication::translate("RightWidget", "Time Step", nullptr));
        setting_list_pushButton->setText(QCoreApplication::translate("RightWidget", "Setting List", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Scene), QCoreApplication::translate("RightWidget", "Scene", nullptr));
        floor_groupBox_2->setTitle(QCoreApplication::translate("RightWidget", "Surface Mesh", nullptr));
        floor_groupBox_3->setTitle(QCoreApplication::translate("RightWidget", "Mesh", nullptr));
        mesh_metrial_button_2->setText(QCoreApplication::translate("RightWidget", "Render", nullptr));
        mesh_hide_checkBox_2->setText(QCoreApplication::translate("RightWidget", "Hide", nullptr));
        mesh_scale_label_2->setText(QCoreApplication::translate("RightWidget", "      Scale", nullptr));
        mesh_translate_label_2->setText(QCoreApplication::translate("RightWidget", "Translate", nullptr));
        mesh_rotate_label_2->setText(QCoreApplication::translate("RightWidget", "Rotate", nullptr));
        mesh_material_button->setText(QCoreApplication::translate("RightWidget", "Render", nullptr));
        mesh_hide_checkBox->setText(QCoreApplication::translate("RightWidget", "Hide", nullptr));
        mesh_scale_label->setText(QCoreApplication::translate("RightWidget", "      Scale", nullptr));
        mesh_translate_label->setText(QCoreApplication::translate("RightWidget", "Translate", nullptr));
        mesh_rotate_label->setText(QCoreApplication::translate("RightWidget", "    Rotate", nullptr));
        floor_groupBox_4->setTitle(QCoreApplication::translate("RightWidget", "Medial Mesh", nullptr));
        medial_mesh_hide_checkBox->setText(QCoreApplication::translate("RightWidget", "Hide", nullptr));
        medial_mesh_transparent_checkBox->setText(QCoreApplication::translate("RightWidget", "Transparent", nullptr));
        floor_groupBox_6->setTitle(QCoreApplication::translate("RightWidget", "Volumetric Mesh", nullptr));
        tet_mesh_hide_checkBox->setText(QCoreApplication::translate("RightWidget", "Hide", nullptr));
        tet_mesh_transparent_checkBox->setText(QCoreApplication::translate("RightWidget", "Transparent", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Model), QCoreApplication::translate("RightWidget", "Model", nullptr));
    } // retranslateUi

};

namespace Ui {
    class RightWidget: public Ui_RightWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RIGHTFORMWIDGET_H
