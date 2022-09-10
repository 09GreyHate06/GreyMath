#pragma once
#include "Types.h"
#include "Operators.h"

namespace GreyMath
{
	// =========================================== Vector =================================================

	// = sqrt(x^2 + y^2 ... n^2)
	float Magnitude(const Vector& v)
	{
		return std::sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	}

	// v / magnitude(v)
	Vector Normalized(const Vector& v)
	{
		return v / Magnitude(v);
	}

	// ||P|| * ||Q|| * cos(a)
	float Dot(const Vector& lhs, const Vector& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}


	// =========================================== Matrix ================================================= 
}