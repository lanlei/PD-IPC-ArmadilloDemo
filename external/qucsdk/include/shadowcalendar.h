﻿#ifndef SHADOWCALENDAR_H
#define SHADOWCALENDAR_H

/**
 * 光晕日历控件 作者:雨田哥(QQ:3246214072) 整理:feiyangqingyun(QQ:517216493) 2019-10-07
 * 1:可设置背景颜色
 * 2:可设置光晕颜色
 * 3:可设置文字颜色
 * 4:可设置选中日期背景
 * 5:光晕跟随鼠标移动
 */

#include <QWidget>
#include <QDate>

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT ShadowCalendar : public QWidget
#else
class ShadowCalendar : public QWidget
#endif

{
    Q_OBJECT
    Q_PROPERTY(QColor bgColor READ getBgColor WRITE setBgColor)
    Q_PROPERTY(QColor textColor READ getTextColor WRITE setTextColor)
    Q_PROPERTY(QColor shadowColor READ getShadowColor WRITE setShadowColor)
    Q_PROPERTY(QColor selectColor READ getSelectColor WRITE setSelectColor)

public:
    struct DateItem {
        int year;
        int month;
        int day;

        DateItem()
        {
            year = -1;
            month = -1;
            day = -1;
        }
    };

    explicit ShadowCalendar(QWidget *parent = 0);
    ~ShadowCalendar();

public:
    void updateCalendar(const QDate &selectDate);

protected:
    void leaveEvent(QEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void paintEvent(QPaintEvent *);

private:
    QColor bgColor;             //背景颜色
    QColor textColor;           //文字颜色
    QColor shadowColor;         //光晕颜色
    QColor selectColor;         //选中颜色

    QDate selectDate;           //今天日期
    DateItem dateItem[6][7];    //日期数组

public:
    QColor getBgColor()         const;
    QColor getTextColor()       const;
    QColor getShadowColor()     const;
    QColor getSelectColor()     const;

    QSize sizeHint()            const;
    QSize minimumSizeHint()     const;

public Q_SLOTS:
    //设置背景颜色+文字颜色+光晕颜色+选中颜色
    void setBgColor(const QColor &bgColor);
    void setTextColor(const QColor &textColor);
    void setShadowColor(const QColor &shadowColor);
    void setSelectColor(const QColor &selectColor);
};

#endif // SHADOWCALENDAR_H
