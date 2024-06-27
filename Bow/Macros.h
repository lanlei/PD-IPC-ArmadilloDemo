#ifndef BOW_MACROS
#define BOW_MACROS

#include <iostream>
#define BOW_STATIC_LIBRARY

#ifndef BOW_STATIC_LIBRARY
#define BOW_INLINE inline
#else
#define BOW_INLINE
#endif



#define BOW_COMPILE_FLOAT
#define BOW_COMPILE_DOUBLE
#define BOW_COMPILE_2D
#define BOW_COMPILE_3D

#endif