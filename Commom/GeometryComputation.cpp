#include "GeometryComputation.h"

qeal TriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
	return 0.5 * fabs(((v1 - v0).cross(v2 - v0)).norm());
}

void TriangleNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& n)
{
	Vector3 vec1 = v1 - v0;
	Vector3 vec2 = v2 - v0;
	n = vec1.cross(vec2);
	n.normalize();
}

bool DistanceToLine(const Vector3 & p, const Vector3 & v0, const Vector3 & v1, qeal& dist, qeal & t, Vector3 & fp)
{
	Vector3 v0v1 = v1 - v0;
	Vector3 pv0 = v0 - p;
	Vector3 pv1 = v1 - p;

	qeal area = abs((v0v1.cross(pv0)).norm());

	if (!IS_QEAL_ZERO(v0v1.norm()))
	{
		dist = area / v0v1.norm();
		t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
		fp = (1.0 - t) * v0 + t * v1;
		return true;
	}
	else return false;
}

bool SameSize(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& b)
{
	Vector3 cp1 = (b - a).cross(p1 - a);
	Vector3 cp2 = (b - a).cross(p2 - a);
	if (!IS_QEAL_ZERO(cp1.dot(cp2)))
		return true;
	else
		return false;
}

bool InsideTriangle(const Vector3& p, const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
	if (SameSize(p, v0, v1, v2) && SameSize(p, v1, v0, v2) && SameSize(p, v2, v1, v0))
		return true;
	return false;

}

void TriBarycentericCoordinate(const Vector3& p, const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& coordinate)
{
	qeal area = TriangleArea(v0, v1, v2);
	coordinate[0] = TriangleArea(p, v1, v2) / area;
	coordinate[1] = TriangleArea(p, v0, v2) / area;
	coordinate[2] = TriangleArea(p, v1, v0) / area;
}


bool SameSide(Vector3 A, Vector3 B, Vector3 C, Vector3 P)
{
	Vector3 AB = B - A;
	Vector3 AC = C - A;
	Vector3 AP = P - A;

	Vector3 v1 = AB.cross(AC);
	Vector3 v2 = AB.cross(AP);

	// v1 and v2 should point to the same direction
	return v1.dot(v2) >= 0;
}

// Same side method
// Determine whether point P in triangle ABC
bool PointinTriangle1(Vector3 A, Vector3 B, Vector3 C, Vector3 P)
{
	return SameSide(A, B, C, P) &&
		SameSide(B, C, A, P) &&
		SameSide(C, A, B, P);
}


qeal VertexEdgeSqDistance(Vector3 v, Vector3 ea, Vector3 eb, qeal& r, Vector3& dir)
{
	Vector3 va, ba;
	va = v - ea;
	ba = eb - ea;

	qeal va_ba = va.dot(ba);
	qeal ba_ba = ba.dot(ba);

	if (va_ba < 0)
		r = 0;
	else if (va_ba > ba_ba)
		r = 1;
	else r = va_ba / ba_ba;

	dir = v - ((1.0 - r) * ea + r * eb);

	return dir.dot(dir);
}

qeal EdgeEdgeSqDistance(Vector3 xi, Vector3 xj, Vector3 xa, Vector3 xb, qeal& r, qeal& s, Vector3& dir)
{
	Vector3 xba, xji, xai;
	xba = xb - xa;
	xji = xj - xi;
	xai = xa - xi;
	dir = xji.cross(xba);
	qeal nn = dir.dot(dir);
	Vector3 temp;
	temp = xai.cross(xji);

	qeal weight_aiji = dir.dot(temp);

	temp = xai.cross(xba);
	qeal weight_aiba = dir.dot(temp);

	if (nn > 1e-24f && weight_aiji >= 0 && weight_aiji <= nn && weight_aiba >= 0 && weight_aiba <= nn)
	{
		r = weight_aiba / nn;
		s = weight_aiji / nn;
	}
	else
	{
		qeal minDistance = 999999999;
		qeal distance, v;
		Vector3 Ndir;
		if (weight_aiba < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xb, v, Ndir)) < minDistance))
		{
			minDistance = distance;
			r = 0;
			s = v;
		}
		if (weight_aiba > nn && ((distance = VertexEdgeSqDistance(xj, xa, xb, v, Ndir)) < minDistance))
		{
			minDistance = distance;
			r = 1;
			s = v;
		}
		if (weight_aiji < 0 && ((distance = VertexEdgeSqDistance(xa, xi, xj, v, Ndir)) < minDistance))
		{
			minDistance = distance;
			r = v;
			s = 0;
		}
		if (weight_aiji > nn && ((distance = VertexEdgeSqDistance(xb, xi, xj, v, Ndir)) < minDistance))
		{
			minDistance = distance;
			r = v;
			s = 1;
		}
	}
	dir = xi * (1.0 - r) + xj * r - xa * (1.0 - s) - xb * s;
	return dir.dot(dir);
}

qeal VertexTriangleDistance(Vector3 xi, Vector3 xa, Vector3 xb, Vector3 xc, qeal& ba, qeal& bb, qeal& bc, Vector3& dir)
{
	Vector3 xba, xca, xia;
	xba = xb - xa;
	xca = xc - xa;
	xia = xi - xa;
	dir = xba.cross(xca);

	qeal nn = dir.dot(dir);

	Vector3 temp;
	temp = xia.cross(xca);
	qeal weight_iaca = dir.dot(temp);
	temp = xba.cross(xia);
	qeal weight_baia = dir.dot(temp);

	if (nn > 1e-24f && weight_iaca >= 0 && weight_baia >= 0 && nn - weight_iaca - weight_baia >= 0)
	{
		bb = weight_iaca / nn;
		bc = weight_baia / nn;
		ba = 1 - bb - bc;
	}
	else
	{
		qeal minDistance = 999999999;
		qeal r, distance;
		Vector3 N;
		if (nn - weight_iaca - weight_baia < 0 && ((distance = VertexEdgeSqDistance(xi, xb, xc, r, N)) < minDistance))
		{
			minDistance = distance;
			bb = 1 - r;
			bc = r;
			ba = 0;
		}
		if (weight_iaca < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xc, r, N)) < minDistance))
		{
			minDistance = distance;
			bb = 0;
			bc = r;
			ba = 1 - bb - bc;
		}
		if (weight_baia < 0 && ((distance = VertexEdgeSqDistance(xi, xa, xb, r, N)) < minDistance))
		{
			minDistance = distance;
			bb = r;
			bc = 0;
			ba = 1 - bb - bc;
		}
	}
	dir = xi - ba * xa - bb * xb - bc * xc;
	return dir.dot(dir);
}