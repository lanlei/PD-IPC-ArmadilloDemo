/****************************************************************************
** Meta object code from reading C++ file 'ColorSelectionPanel.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.14.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../Ui/ColorSelectionPanel.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ColorSelectionPanel.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.14.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ColorSelectionPanel_t {
    QByteArrayData data[19];
    char stringdata0[294];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ColorSelectionPanel_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ColorSelectionPanel_t qt_meta_stringdata_ColorSelectionPanel = {
    {
QT_MOC_LITERAL(0, 0, 19), // "ColorSelectionPanel"
QT_MOC_LITERAL(1, 20, 11), // "colorChange"
QT_MOC_LITERAL(2, 32, 0), // ""
QT_MOC_LITERAL(3, 33, 7), // "QColor*"
QT_MOC_LITERAL(4, 41, 1), // "c"
QT_MOC_LITERAL(5, 43, 10), // "panelClose"
QT_MOC_LITERAL(6, 54, 23), // "updateFromColorPanelBtn"
QT_MOC_LITERAL(7, 78, 5), // "color"
QT_MOC_LITERAL(8, 84, 23), // "updateFromColorPanelHSB"
QT_MOC_LITERAL(9, 108, 3), // "hue"
QT_MOC_LITERAL(10, 112, 3), // "sat"
QT_MOC_LITERAL(11, 116, 21), // "updateFromRedEditLine"
QT_MOC_LITERAL(12, 138, 23), // "updateFromGreenEditLine"
QT_MOC_LITERAL(13, 162, 22), // "updateFromBlueEditLine"
QT_MOC_LITERAL(14, 185, 23), // "updateFromAlphaEditLine"
QT_MOC_LITERAL(15, 209, 19), // "updateFromRedSlider"
QT_MOC_LITERAL(16, 229, 21), // "updateFromGreenSlider"
QT_MOC_LITERAL(17, 251, 20), // "updateFromBlueSlider"
QT_MOC_LITERAL(18, 272, 21) // "updateFromAlphaSlider"

    },
    "ColorSelectionPanel\0colorChange\0\0"
    "QColor*\0c\0panelClose\0updateFromColorPanelBtn\0"
    "color\0updateFromColorPanelHSB\0hue\0sat\0"
    "updateFromRedEditLine\0updateFromGreenEditLine\0"
    "updateFromBlueEditLine\0updateFromAlphaEditLine\0"
    "updateFromRedSlider\0updateFromGreenSlider\0"
    "updateFromBlueSlider\0updateFromAlphaSlider"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ColorSelectionPanel[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   74,    2, 0x06 /* Public */,
       5,    0,   77,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       6,    1,   78,    2, 0x0a /* Public */,
       8,    3,   81,    2, 0x0a /* Public */,
      11,    0,   88,    2, 0x0a /* Public */,
      12,    0,   89,    2, 0x0a /* Public */,
      13,    0,   90,    2, 0x0a /* Public */,
      14,    0,   91,    2, 0x0a /* Public */,
      15,    1,   92,    2, 0x0a /* Public */,
      16,    1,   95,    2, 0x0a /* Public */,
      17,    1,   98,    2, 0x0a /* Public */,
      18,    1,  101,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QColor,    7,
    QMetaType::Void, QMetaType::QColor, QMetaType::Double, QMetaType::Double,    7,    9,   10,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,

       0        // eod
};

void ColorSelectionPanel::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<ColorSelectionPanel *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->colorChange((*reinterpret_cast< QColor*(*)>(_a[1]))); break;
        case 1: _t->panelClose(); break;
        case 2: _t->updateFromColorPanelBtn((*reinterpret_cast< const QColor(*)>(_a[1]))); break;
        case 3: _t->updateFromColorPanelHSB((*reinterpret_cast< const QColor(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        case 4: _t->updateFromRedEditLine(); break;
        case 5: _t->updateFromGreenEditLine(); break;
        case 6: _t->updateFromBlueEditLine(); break;
        case 7: _t->updateFromAlphaEditLine(); break;
        case 8: _t->updateFromRedSlider((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->updateFromGreenSlider((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->updateFromBlueSlider((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->updateFromAlphaSlider((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (ColorSelectionPanel::*)(QColor * );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ColorSelectionPanel::colorChange)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (ColorSelectionPanel::*)();
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&ColorSelectionPanel::panelClose)) {
                *result = 1;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject ColorSelectionPanel::staticMetaObject = { {
    QMetaObject::SuperData::link<QWidget::staticMetaObject>(),
    qt_meta_stringdata_ColorSelectionPanel.data,
    qt_meta_data_ColorSelectionPanel,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *ColorSelectionPanel::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ColorSelectionPanel::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ColorSelectionPanel.stringdata0))
        return static_cast<void*>(this);
    if (!strcmp(_clname, "Ui::ColorSelectionWidget"))
        return static_cast< Ui::ColorSelectionWidget*>(this);
    return QWidget::qt_metacast(_clname);
}

int ColorSelectionPanel::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void ColorSelectionPanel::colorChange(QColor * _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void ColorSelectionPanel::panelClose()
{
    QMetaObject::activate(this, &staticMetaObject, 1, nullptr);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
