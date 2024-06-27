#include "Utils.h"

bool v2IndComp(const Vector2i& lhs, const Vector2i& rhs) {
	if (lhs[0] == rhs[0])
		return lhs[1] < rhs[1];
	else return lhs[0] < rhs[0];
}

bool v3IndComp(const Vector3i& lhs, const Vector3i& rhs) {
	if (lhs[0] == rhs[0])
	{
		if (lhs[1] == rhs[1])
			return lhs[2] < rhs[2];
		else return lhs[1] < rhs[1];
	}
	else return lhs[0] < rhs[0];
}

