/****************************************************************************
** Meta object code from reading C++ file 'BaseScene.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.14.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../Scene/BaseScene.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'BaseScene.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.14.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_BaseScene_t {
    QByteArrayData data[6];
    char stringdata0[71];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_BaseScene_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_BaseScene_t qt_meta_stringdata_BaseScene = {
    {
QT_MOC_LITERAL(0, 0, 9), // "BaseScene"
QT_MOC_LITERAL(1, 10, 15), // "initSceneSignal"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 15), // "changeSimulator"
QT_MOC_LITERAL(4, 43, 15), // "recordAnimation"
QT_MOC_LITERAL(5, 59, 11) // "needInitSim"

    },
    "BaseScene\0initSceneSignal\0\0changeSimulator\0"
    "recordAnimation\0needInitSim"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_BaseScene[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   34,    2, 0x06 /* Public */,
       3,    0,   35,    2, 0x06 /* Public */,
       4,    1,   36,    2, 0x06 /* Public */,
       5,    0,   39,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,

       0        // eod
};

void BaseScene::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<BaseScene *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->initSceneSignal(); break;
        case 1: _t->changeSimulator(); break;
        case 2: _t->recordAnimation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: _t->needInitSim(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (BaseScene::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseScene::initSceneSignal)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (BaseScene::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseScene::changeSimulator)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (BaseScene::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseScene::recordAnimation)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (BaseScene::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseScene::needInitSim)) {
                *result = 3;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject BaseScene::staticMetaObject = { {
    QMetaObject::SuperData::link<QGLViewer::staticMetaObject>(),
    qt_meta_stringdata_BaseScene.data,
    qt_meta_data_BaseScene,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *BaseScene::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *BaseScene::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_BaseScene.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "QOpenGLFunctions"))
        return static_cast< QOpenGLFunctions*>(this);
    return QGLViewer::qt_metacast(_clname);
}

int BaseScene::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLViewer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void BaseScene::initSceneSignal()
{
    QMetaObject::activate(this, &staticMetaObject, 0, nullptr);
}

// SIGNAL 1
void BaseScene::changeSimulator()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}

// SIGNAL 2
void BaseScene::recordAnimation(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void BaseScene::needInitSim()
{
    QMetaObject::activate(this, &staticMetaObject, 3, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
