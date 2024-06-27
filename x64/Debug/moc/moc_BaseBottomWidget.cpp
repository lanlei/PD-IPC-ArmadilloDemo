/****************************************************************************
** Meta object code from reading C++ file 'BaseBottomWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.14.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../Ui/BaseBottomWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'BaseBottomWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.14.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_BaseBottomWidget_t {
    QByteArrayData data[17];
    char stringdata0[307];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_BaseBottomWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_BaseBottomWidget_t qt_meta_stringdata_BaseBottomWidget = {
    {
QT_MOC_LITERAL(0, 0, 16), // "BaseBottomWidget"
QT_MOC_LITERAL(1, 17, 15), // "changeViewerNum"
QT_MOC_LITERAL(2, 33, 0), // ""
QT_MOC_LITERAL(3, 34, 15), // "recordAnimation"
QT_MOC_LITERAL(4, 50, 13), // "saveAnimation"
QT_MOC_LITERAL(5, 64, 13), // "playAnimation"
QT_MOC_LITERAL(6, 78, 14), // "resetAnimation"
QT_MOC_LITERAL(7, 93, 20), // "bindMultiViewerGroup"
QT_MOC_LITERAL(8, 114, 18), // "bindAnimationGroup"
QT_MOC_LITERAL(9, 133, 28), // "handleMultiViewerButtonGroup"
QT_MOC_LITERAL(10, 162, 16), // "QAbstractButton*"
QT_MOC_LITERAL(11, 179, 19), // "handleInitSimulator"
QT_MOC_LITERAL(12, 199, 21), // "handleRecordAnimation"
QT_MOC_LITERAL(13, 221, 19), // "handleSaveAnimation"
QT_MOC_LITERAL(14, 241, 19), // "handlePlayAnimation"
QT_MOC_LITERAL(15, 261, 20), // "handleResetAnimation"
QT_MOC_LITERAL(16, 282, 24) // "resetInitSimulatorStatus"

    },
    "BaseBottomWidget\0changeViewerNum\0\0"
    "recordAnimation\0saveAnimation\0"
    "playAnimation\0resetAnimation\0"
    "bindMultiViewerGroup\0bindAnimationGroup\0"
    "handleMultiViewerButtonGroup\0"
    "QAbstractButton*\0handleInitSimulator\0"
    "handleRecordAnimation\0handleSaveAnimation\0"
    "handlePlayAnimation\0handleResetAnimation\0"
    "resetInitSimulatorStatus"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_BaseBottomWidget[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   84,    2, 0x06 /* Public */,
       3,    1,   87,    2, 0x06 /* Public */,
       4,    1,   90,    2, 0x06 /* Public */,
       5,    1,   93,    2, 0x06 /* Public */,
       6,    0,   96,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       7,    0,   97,    2, 0x0a /* Public */,
       8,    0,   98,    2, 0x0a /* Public */,
       9,    1,   99,    2, 0x0a /* Public */,
      11,    0,  102,    2, 0x0a /* Public */,
      12,    0,  103,    2, 0x0a /* Public */,
      13,    0,  104,    2, 0x0a /* Public */,
      14,    0,  105,    2, 0x0a /* Public */,
      15,    0,  106,    2, 0x0a /* Public */,
      16,    0,  107,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 10,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void BaseBottomWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<BaseBottomWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->changeViewerNum((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->recordAnimation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->saveAnimation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: _t->playAnimation((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: _t->resetAnimation(); break;
        case 5: _t->bindMultiViewerGroup(); break;
        case 6: _t->bindAnimationGroup(); break;
        case 7: _t->handleMultiViewerButtonGroup((*reinterpret_cast< QAbstractButton*(*)>(_a[1]))); break;
        case 8: _t->handleInitSimulator(); break;
        case 9: _t->handleRecordAnimation(); break;
        case 10: _t->handleSaveAnimation(); break;
        case 11: _t->handlePlayAnimation(); break;
        case 12: _t->handleResetAnimation(); break;
        case 13: _t->resetInitSimulatorStatus(); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 7:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAbstractButton* >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (BaseBottomWidget::*)(int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseBottomWidget::changeViewerNum)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (BaseBottomWidget::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseBottomWidget::recordAnimation)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (BaseBottomWidget::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseBottomWidget::saveAnimation)) {
                *result = 2;
                return;
            }
        }
        {
            using _t = void (BaseBottomWidget::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseBottomWidget::playAnimation)) {
                *result = 3;
                return;
            }
        }
        {
            using _t = void (BaseBottomWidget::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&BaseBottomWidget::resetAnimation)) {
                *result = 4;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject BaseBottomWidget::staticMetaObject = { {
    QMetaObject::SuperData::link<QWidget::staticMetaObject>(),
    qt_meta_stringdata_BaseBottomWidget.data,
    qt_meta_data_BaseBottomWidget,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *BaseBottomWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *BaseBottomWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_BaseBottomWidget.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "Ui::BottomWidget"))
        return static_cast< Ui::BottomWidget*>(this);
    return QWidget::qt_metacast(_clname);
}

int BaseBottomWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    }
    return _id;
}

// SIGNAL 0
void BaseBottomWidget::changeViewerNum(int _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void BaseBottomWidget::recordAnimation(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void BaseBottomWidget::saveAnimation(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void BaseBottomWidget::playAnimation(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void BaseBottomWidget::resetAnimation()
{
    QMetaObject::activate(this, &staticMetaObject, 4, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
