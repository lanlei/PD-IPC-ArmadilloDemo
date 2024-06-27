/****************************************************************************
** Meta object code from reading C++ file 'BaseMainWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.14.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../Ui/BaseMainWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'BaseMainWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.14.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_BaseMainWidget_t {
    QByteArrayData data[13];
    char stringdata0[279];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_BaseMainWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_BaseMainWidget_t qt_meta_stringdata_BaseMainWidget = {
    {
QT_MOC_LITERAL(0, 0, 14), // "BaseMainWidget"
QT_MOC_LITERAL(1, 15, 14), // "addFileSuccess"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 20), // "handleOpenFileAction"
QT_MOC_LITERAL(4, 52, 20), // "handleSaveFileAction"
QT_MOC_LITERAL(5, 73, 18), // "handleTetgenAction"
QT_MOC_LITERAL(6, 92, 26), // "handleRenderLineModeAction"
QT_MOC_LITERAL(7, 119, 31), // "handleSelectSurfacePointsAction"
QT_MOC_LITERAL(8, 151, 27), // "handleSelectTetPointsAction"
QT_MOC_LITERAL(9, 179, 22), // "handlePlayAnimationBtn"
QT_MOC_LITERAL(10, 202, 23), // "handleRecorAnimationBtn"
QT_MOC_LITERAL(11, 226, 23), // "handleResetAnimationBtn"
QT_MOC_LITERAL(12, 250, 28) // "handleMultiViewerButtonGroup"

    },
    "BaseMainWidget\0addFileSuccess\0\0"
    "handleOpenFileAction\0handleSaveFileAction\0"
    "handleTetgenAction\0handleRenderLineModeAction\0"
    "handleSelectSurfacePointsAction\0"
    "handleSelectTetPointsAction\0"
    "handlePlayAnimationBtn\0handleRecorAnimationBtn\0"
    "handleResetAnimationBtn\0"
    "handleMultiViewerButtonGroup"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_BaseMainWidget[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   69,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   72,    2, 0x0a /* Public */,
       4,    0,   73,    2, 0x0a /* Public */,
       5,    0,   74,    2, 0x0a /* Public */,
       6,    0,   75,    2, 0x0a /* Public */,
       7,    0,   76,    2, 0x0a /* Public */,
       8,    0,   77,    2, 0x0a /* Public */,
       9,    1,   78,    2, 0x0a /* Public */,
      10,    1,   81,    2, 0x0a /* Public */,
      11,    0,   84,    2, 0x0a /* Public */,
      12,    1,   85,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Bool,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,

       0        // eod
};

void BaseMainWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<BaseMainWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->addFileSuccess((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->handleOpenFileAction(); break;
        case 2: _t->handleSaveFileAction(); break;
        case 3: _t->handleTetgenAction(); break;
        case 4: _t->handleRenderLineModeAction(); break;
        case 5: _t->handleSelectSurfacePointsAction(); break;
        case 6: _t->handleSelectTetPointsAction(); break;
        case 7: _t->handlePlayAnimationBtn((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: _t->handleRecorAnimationBtn((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: _t->handleResetAnimationBtn(); break;
        case 10: _t->handleMultiViewerButtonGroup((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (BaseMainWidget::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseMainWidget::addFileSuccess)) {
                *result = 0;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject BaseMainWidget::staticMetaObject = { {
    QMetaObject::SuperData::link<QWidget::staticMetaObject>(),
    qt_meta_stringdata_BaseMainWidget.data,
    qt_meta_data_BaseMainWidget,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *BaseMainWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *BaseMainWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_BaseMainWidget.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int BaseMainWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}

// SIGNAL 0
void BaseMainWidget::addFileSuccess(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
