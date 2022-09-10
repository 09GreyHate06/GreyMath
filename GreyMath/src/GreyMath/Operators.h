#pragma once
#include "Types.h"
#include <iostream>
#include <iomanip>
#include <string>

namespace GreyMath
{
	// =========================================== Vector =================================================

	Vector operator+(const Vector& lhs, const Vector& rhs)
	{
		return Vector(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
	}

	Vector operator-(const Vector& lhs, const Vector& rhs)
	{
		return Vector(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, lhs.w - rhs.w);
	}

	Vector operator*(const Vector& v, float s)
	{
		return Vector(v.x * s, v.y * s, v.z * s, v.w * s);
	}

	Vector operator*(float s, const Vector& v)
	{
		return Vector(v.x * s, v.y * s, v.z * s, v.w * s);
	}

	Vector operator/(const Vector& v, float s)
	{
		return Vector(v.x / s, v.y / s, v.z / s, v.w / s);
	}


	Vector& operator+=(Vector& lhs, const Vector& rhs)
	{
		lhs.x += rhs.x;
		lhs.y += rhs.y;
		lhs.z += rhs.z;
		lhs.w += rhs.w;

		return lhs;
	}

	Vector& operator-=(Vector& lhs, const Vector& rhs)
	{
		lhs.x -= rhs.x;
		lhs.y -= rhs.y;
		lhs.z -= rhs.z;
		lhs.w -= rhs.w;

		return lhs;
	}

	Vector& operator*=(Vector& lhs, float s)
	{
		lhs.x *= s;
		lhs.y *= s;
		lhs.z *= s;
		lhs.w *= s;

		return lhs;
	}

	Vector& operator/=(Vector& lhs, float s)
	{
		lhs.x /= s;
		lhs.y /= s;
		lhs.z /= s;
		lhs.w /= s;

		return lhs;
	}


	Vector operator-(const Vector& v)
	{
		return Vector(-v.x, -v.y, -v.z, -v.w);
	}

	std::ostream& operator<<(std::ostream& os, const Vector& v)
	{
		os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
		return os;
	}





	// =========================================== Matrix ================================================= 

	std::ostream& operator<<(std::ostream& os, const Matrix& m)
	{
		size_t maxLenPerCol[4] = { 0, 0, 0, 0 };

		for (int y = 0; y < 4; y++)
		{
			for (int x = 0; x < 4; x++)
			{
				size_t curLength = std::to_string(m[x][y]).size();
				if (curLength > maxLenPerCol[y])
					maxLenPerCol[y] = curLength;
			}
		}

		for (int x = 0; x < 4; x++)
		{
			os << "|";
			for (int y = 0; y < 4; y++)
			{
				os << std::setw(maxLenPerCol[y]) << m[x][y];
				if (y == 3)
					os << "|";
			}

			os << '\n';
		}

		//os << "****** Matrix ******\n";
		//os << m[0] << '\n' << m[1] << '\n' << m[2] << '\n' << m[3] << '\n';
		//os << "********************";

		return os;
	}
}