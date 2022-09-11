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

	// ||lhs|| * ||rhs|| * cos(a)
	float Dot(const Vector& lhs, const Vector& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}

	// (v0 . v1) / ||v1||^2) * v1
	Vector ProjectionOfV0OntoV1(const Vector& v0, const Vector& v1)
	{
		float v1MagSqrd = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z + v1.w * v1.w;
		return (Dot(v0, v1) / v1MagSqrd) * v1;
	}


	// =========================================== Matrix ================================================= 
}