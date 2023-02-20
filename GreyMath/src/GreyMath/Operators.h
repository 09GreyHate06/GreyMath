#pragma once
#include "Types.h"
#include <iostream>
#include <iomanip>
#include <string>

namespace GreyMath
{
	struct Matrix;
	constexpr float Dot(const Vector& lhs, const Vector& rhs);

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

	Vector operator*(const Vector& v, const Matrix& m)
	{
		Vector result = {};
		for (int i = 0; i < 4; i++)
		{
			result[i] = Dot(v, m.GetColumnVector(i));
		}

		return result;
	}

	Vector operator*(const Matrix& m, const Vector& v)
	{
		Vector result = {};
		for (int i = 0; i < 4; i++)
		{
			result[i] = Dot(m[i], v);
		}

		return result;
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


	Matrix operator+(const Matrix& a, const Matrix& b)
	{
		return Matrix(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]);
	}

	Matrix operator-(const Matrix& a, const Matrix& b)
	{
		return Matrix(a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]);
	}

	Matrix operator-(const Matrix& m)
	{
		return Matrix(-m[0], -m[1], -m[2], -m[3]);
	}

	Matrix operator*(const Matrix& m, float s)
	{
		return Matrix(m[0] * s, m[1] * s, m[2] * s, m[3] * s);
	}

	Matrix operator*(float s, const Matrix& m)
	{
		return Matrix(m[0] * s, m[1] * s, m[2] * s, m[3] * s);
	}

	Matrix operator*(const Matrix& a, const Matrix& b)
	{
		Matrix ab;
		
		// algorithm
		//for (int i = 0; i < 4; i++)
		//{
		//	for (int j = 0; j < 4; j++)
		//	{
		//		float ij = 0.0f;
		//		for (int k = 0; k < 4; k++)
		//		{
		//			ij += a[i][k] * b[k][j];
		//		}
		//
		//		ab[i][j] = ij;
		//	}
		//}

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ab[i][j] = Dot(a[i], b.GetColumnVector(j));
			}
		}

		return ab;
	}

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