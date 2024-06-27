#pragma once
#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H
#include "MatrixCore.h"
// Ax = b
void gaussElimination(MatrixX& A, VectorX& b, VectorX& x, MatrixX& U);

#endif