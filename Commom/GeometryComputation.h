#pragma once
#ifndef GEOMETRY_COMPUTATION_H
#define GEOMETRY_COMPUTATION_H
#include "MatrixCore.h"

qeal TriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2);

void TriangleNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& n);

bool DistanceToLine(const Vector3 & p, const Vector3 & v0, const Vector3 & v1, qeal& dist, qeal & t, Vector3 & fp);

bool SameSize(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& b);

bool InsideTriangle(const Vector3& p, const Vector3& v0, const Vector3& v1, const Vector3& v2);


void TriBarycentericCoordinate(const Vector3& p, const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& coordinate);

bool PointinTriangle1(Vector3 A, Vector3 B, Vector3 C, Vector3 P);

bool SameSide(Vector3 A, Vector3 B, Vector3 C, Vector3 P);


qeal VertexEdgeSqDistance(Vector3 v, Vector3 ea, Vector3 eb, qeal& r, Vector3& dir);

qeal EdgeEdgeSqDistance(Vector3 xi, Vector3 xj, Vector3 xa, Vector3 xb, qeal& r, qeal& s, Vector3& dir);

qeal VertexTriangleDistance(Vector3 xi, Vector3 xa, Vector3 xb, Vector3 xc, qeal& ba, qeal& bb, qeal& bc, Vector3& dir);



#endif