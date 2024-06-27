﻿#ifndef TIMEAXIS_H
#define TIMEAXIS_H

/**
 * 垂直时间轴控件 作者:雨田哥(QQ:3246214072) 整理:feiyangqingyun(QQ:517216493) 2019-10-07
 * 1:可设置节点边距
 * 2:可设置节点高度
 * 3:可设置信息边框边距
 * 4:可设置信息所占高度
 * 5:可设置基准颜色/线条颜色
 * 6:可设置标题/信息集合
 * 7:自动产生滚动条
 * 8:支持字符串形式设置数据
 */

#include <QScrollArea>
class TimeAxisWidget;

#ifdef quc
#if (QT_VERSION < QT_VERSION_CHECK(5,7,0))
#include <QtDesigner/QDesignerExportWidget>
#else
#include <QtUiPlugin/QDesignerExportWidget>
#endif

class QDESIGNER_WIDGET_EXPORT TimeAxis : public QScrollArea
#else
class TimeAxis : public QScrollArea
#endif

{
    Q_OBJECT
    Q_PROPERTY(int itemMargin READ getItemMargin WRITE setItemMargin)
    Q_PROPERTY(int itemHeight READ getItemHeight WRITE setItemHeight)
    Q_PROPERTY(int infoPadding READ getInfoPadding WRITE setInfoPadding)
    Q_PROPERTY(int infoHeight READ getInfoHeight WRITE setInfoHeight)

    Q_PROPERTY(QColor baseColor READ getBaseColor WRITE setBaseColor)
    Q_PROPERTY(QColor lineColor READ getLineColor WRITE setLineColor)

    Q_PROPERTY(QString title READ getTitle WRITE setTitle)
    Q_PROPERTY(QString infos READ getInfos WRITE setInfos)

public:
    explicit TimeAxis(QWidget *parent = 0);

private:
    int itemMargin;         //节点边距
    int itemHeight;         //节点高度
    int infoPadding;        //信息边距
    int infoHeight;         //信息高度

    QColor baseColor;       //基准颜色
    QColor lineColor;       //线条颜色

    QString title;          //标题
    QString infos;          //信息集合

    //时间轴主控件
    TimeAxisWidget *timeAxisWidget;

public:
    int getItemMargin()     const;
    int getItemHeight()     const;
    int getInfoPadding()    const;
    int getInfoHeight()     const;

    QColor getBaseColor()   const;
    QColor getLineColor()   const;

    QString getTitle()      const;
    QString getInfos()      const;

    QSize sizeHint()        const;
    QSize minimumSizeHint() const;

    TimeAxisWidget *getWidget();

public Q_SLOTS:
    //设置节点边距+节点高度
    void setItemMargin(int itemMargin);
    void setItemHeight(int itemHeight);

    //设置信息边距+信息高度
    void setInfoPadding(int infoPadding);
    void setInfoHeight(int infoHeight);

    //设置基准颜色+线条颜色
    void setBaseColor(const QColor &baseColor);
    void setLineColor(const QColor &lineColor);

    //设置标题+信息集合
    void setTitle(const QString &title);
    void setInfos(const QString &infos);
};

class TimeAxisWidget : public QWidget
{
    Q_OBJECT

public:
    //可以自行拓展其他信息
    struct TimeAxisInfo {
        QString time;   //时间
        QString info;   //信息
    };

    explicit TimeAxisWidget(QWidget *parent = 0);

protected:
    void paintEvent(QPaintEvent *);
    void drawTitle(QPainter *painter);
    void drawLine(QPainter *painter);
    void drawInfo(QPainter *painter);
    void drawInfoRight(QPainter *painter, const QRectF &infoRect, int infoHeight);
    void drawInfoLeft(QPainter *painter, const QRectF &infoRect, int infoHeight);

private:
    int itemMargin;         //节点边距
    int itemHeight;         //节点高度
    int infoPadding;        //信息边距
    int infoHeight;         //信息高度

    QColor baseColor;       //基准颜色
    QColor lineColor;       //线条颜色

    QString title;          //标题
    QString infos;          //信息集合

    //信息集合结构体
    QList<TimeAxisInfo> itemInfos;

public:
    int getItemMargin()     const;
    int getItemHeight()     const;
    int getInfoPadding()    const;
    int getInfoHeight()     const;

    QColor getBaseColor()   const;
    QColor getLineColor()   const;

    QString getTitle()      const;
    QString getInfos()      const;

    QSize sizeHint()        const;
    QSize minimumSizeHint() const;

public Q_SLOTS:
    //设置节点边距+节点高度
    void setItemMargin(int itemMargin);
    void setItemHeight(int itemHeight);

    //设置信息边距+信息高度
    void setInfoPadding(int infoPadding);
    void setInfoHeight(int infoHeight);

    //设置基准颜色+线条颜色
    void setBaseColor(const QColor &baseColor);
    void setLineColor(const QColor &lineColor);

    //设置标题+信息集合
    void setTitle(const QString &title);
    void setInfos(const QString &infos);

    //设置信息集合,结构体方式
    void setItemInfos(const QList<TimeAxisInfo> &itemInfos);
};

#endif // TIMEAXIS_H
