#pragma once
#ifndef Utils_H
#define Utils_H
#include "MatrixCore.h"

bool v2IndComp(const Vector2i& lhs, const Vector2i& rhs);

bool v3IndComp(const Vector3i& lhs, const Vector3i& rhs);
struct Vec2IndComp {
	bool operator() (const Vector2i& lhs, const Vector2i& rhs) const {
		return v2IndComp(lhs, rhs);
	}
};

struct Vec3IndComp {
	bool operator() (const Vector3i& lhs, const Vector3i& rhs) const {
		return v3IndComp(lhs, rhs);
	}
};

struct VectorHash2i {
	typedef Vector2i IV;
	size_t operator()(const IV& a) const
	{
		std::size_t h = 0;
		for (int d = 0; d < 2; ++d) {
			h ^= std::hash<int>{}(a(d)) + 0x9e3779b9 + (h << 6) + (h >> 2);
		}
		return h;
	}
};

struct VectorHash3i {
	typedef Vector3i IV;
	size_t operator()(const IV& a) const
	{
		std::size_t h = 0;
		for (int d = 0; d < 3; ++d) {
			h ^= std::hash<int>{}(a(d)) + 0x9e3779b9 + (h << 6) + (h >> 2);
		}
		return h;
	}
};



#endif