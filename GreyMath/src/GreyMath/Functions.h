#pragma once
#include "Types.h"
#include "Operators.h"
#include <array>
#include <vector>
#include <assert.h>

namespace GreyMath
{
	Quaternion QuatConjugate(const Quaternion&);
	Quaternion QuatMultiply(const Quaternion&, const Quaternion&);

	// =========================================== Vector =================================================
	// = sqrt(x^2 + y^2 ... n^2)
	float Vec4Magnitude(const Vector& v)
	{
		return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	}

	// = sqrt(x^2 + y^2 ... n^2)
	float Vec3Magnitude(const Vector& v)
	{
		return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	}

	// = sqrt(x^2 + y^2 ... n^2)
	float Vec2Magnitude(const Vector& v)
	{
		return sqrtf(v.x * v.x + v.y * v.y);
	}

	// v / |v|
	Vector Vec4Normalized(const Vector& v)
	{
		return v / Vec4Magnitude(v);
	}

	// v / |v|
	Vector Vec3Normalized(const Vector& v)
	{
		return v / Vec3Magnitude(v);
	}

	// v / |v|
	Vector Vec2Normalized(const Vector& v)
	{
		return v / Vec2Magnitude(v);
	}


	// v0.x^2 + v0.y^2 + v0.z^2 + v0.w^2 == ||v0|| * ||v1|| * cos(a)
	float Vec4Dot(const Vector& v0, const Vector& v1)
	{
		return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
	}

	// v0.x^2 + v0.y^2 + v0.z^2== ||v0|| * ||v1|| * cos(a)
	float Vec3Dot(const Vector& v0, const Vector& v1)
	{
		return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
	}

	// v0.x^2 + v0.y^2 == ||v0|| * ||v1|| * cos(a)
	float Vec2Dot(const Vector& v0, const Vector& v1)
	{
		return v0.x * v1.x + v0.y * v1.y;
	}


	// Pseudodeterminant
	// |   i      j       k   |
	// | v0.x  v0.y   v0.z | 
	// | v1.x  v1.y   v1.z |
	// = i(v0.y * v1.z - v0.z * v1.y) - j(v0.x * v1.z - v0.z * v1.x) + k(v0.x * v1.y - v0.y * v1.x)
	// ||result|| = ||v0|| * ||v1|| * sin(a)
	Vector Vec3Cross(const Vector& v0, const Vector& v1)
	{
		//Vector i(1.0f, 0.0f, 0.0f, 0.0f);
		//Vector j(0.0f, 1.0f, 0.0f, 0.0f);
		//Vector k(0.0f, 0.0f, 1.0f, 0.0f);
		//
		//return i * (v0.y * v1.z - v0.z * v1.y) - j * (v0.x * v1.z - v0.z * v1.x) + k * (v0.x * v1.y - v0.y * v1.x);

		return Vector(v0.y * v1.z - v0.z * v1.y, v0.z * v1.x - v0.x * v1.z, v0.x * v1.y - v0.y * v1.x, 1.0f);
	}

	// (v0 . v1) / ||v1||^2) * v1
	Vector Vec3ProjectionOfV0OntoV1(const Vector& v0, const Vector& v1)
	{
		//return (Dot(v0, v1) / Magnitude(v1)) * Normalized(v1);

		float v1MagSqrd = Vec3Dot(v1, v1);
		return (Vec3Dot(v0, v1) / v1MagSqrd) * v1;
	}

	Vector Vec4Transform(const Vector& v, const Matrix& m)
	{
		Vector v_(0.0f, 0.0f, 0.0f, 1.0f);

		for (int i = 0; i < 4; i++)
		{
			v_[i] = Vec4Dot(v, m.GetColumnVector(i));
		}

		return v_;
	}

	Vector Vec3Transform(const Vector& v, const Matrix& m)
	{
		Vector v_(0.0f, 0.0f, 0.0f, 1.0f);

		for (int i = 0; i < 3; i++)
		{
			v_[i] = Vec3Dot(v, m.GetColumnVector(i));
		}

		return v_;
	}

	Vector Vec2Transform(const Vector& v, const Matrix& m)
	{
		Vector v_(0.0f, 0.0f, 0.0f, 1.0f);

		for (int i = 0; i < 2; i++)
		{
			v_[i] = Vec2Dot(v, m.GetColumnVector(i));
		}

		return v_;
	}

	/// <summary>
	/// Rotate a vector using quaternion
	/// </summary>
	/// <param name="v">vector to ratate</param>
	/// <param name="q">unit quaternion</param>
	/// <returns></returns>
	Vector Vec3Rotate(const Vector& v, const Quaternion& q)
	{
		return QuatMultiply(QuatMultiply(q, v), QuatConjugate(q));
	}


















	





	// =========================================== Quaternion =============================================

	Quaternion QuatIdentity()
	{
		return Quaternion(0.0f, 0.0f, 0.0f, 1.0f);
	}

	/// <summary>
	/// lets say that q0 = (s0 + v0) and q1 = (s1 + v1) s0 and s1 is the scalar part(the w component) and v0 and v1 is the vector part(the x, y, z component) then 
	/// q0*q1 = s0 * s1 - Dot(v0, v1) + s0 * v1 + s1 * v0 + Cross(v0, v1) 
	/// = Vector((s0 * v1 + s1 * v0 + Cross(v0, v1)).xyz, s0 * s1 - Dot(v0, v1))
	/// </summary>
	/// <param name="q0"></param>
	/// <param name="q1"></param>
	/// <returns></returns>
	Quaternion QuatMultiply(const Quaternion& q0, const Quaternion& q1)
	{
		// i^2 = j^2 = k^2 = -1
		// ij = -ji = k; jk = -kj = i; ki = -ik = j
		// 
		// (q0.w + q0.x*i + q0.y*j + q0.z*k) * (q1.w + q1.x*i + q1.y*j + q1.z*k)
		// =
		// a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z +
		// (a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y)*i +
		// (a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x)*j +
		// (a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w)*k
		// 
		// compact way:
		// lets say that q0 = (s0 + v0) and q1 = (s1 + v1) s0 and s1 is the scalar part(the w component) and v0 and v1 is the vector part(the x, y, z component) then
		// q0*q1 = s0*s1 - Dot(v0, v1) + s0*v1 + s1*v0 + Cross(v0, v1)
		// = Vector((s0*v1 + s1*v0 + Cross(v0, v1)).xyz, s0*s1 - Dot(v0, v1))

		float s = q0.w * q1.w - Vec3Dot(q0, q1);
		Vector v = q0.w * q1 + q1.w * q0 + Vec3Cross(q0, q1);
		v.w = s;
		return v;
	}

	Quaternion QuatNormalized(const Quaternion& q)
	{
		return Vec4Normalized(q);
	}

	/// <summary>
	/// = (-q.x, -q.y, - q.z, q.w)
	/// </summary>
	/// <param name="q"></param>
	/// <returns></returns>
	Quaternion QuatConjugate(const Quaternion& q)
	{
		return Quaternion(-q.x, -q.y, -q.z, q.w);
	}

	/// <summary>
	/// Inverse(q) = Conjugate(q)/|q|^2
	/// </summary>
	/// <param name="q"></param>
	/// <returns></returns>
	Quaternion QuatInverse(const Quaternion& q)
	{
		//q*Conjugate(q) = Conjugate(q)*q = Dot(q, q) = |q|^2
		// given that
		// Inverse(q) = Conjugate(q)/|q|^2
		// q*Inverse(q) = q*Conjugate(q) / |q|^2 = |q|^2 / |q|^2 = 1
		return QuatConjugate(q) / Vec4Dot(q, q);
	}

	Quaternion QuatRotationAxis(const Vector& axis, float angle)
	{
		float s = sinf(0.5f * angle);
		float c = cosf(0.5f * angle);

		return Quaternion(s * axis.x, s * axis.y, s * axis.z, c);
	}

	Quaternion QuatRotationRollPitchYaw(float pitch, float yaw, float roll)
	{
		float sP = sinf(0.5f * pitch);
		float sY = sinf(0.5f * yaw);
		float sR = sinf(0.5f * roll);

		float cP = cosf(0.5f * pitch);
		float cY = cosf(0.5f * yaw);
		float cR = cosf(0.5f * roll);


		Quaternion p(sP, 0.0f, 0.0f, cP);
		Quaternion y(0.0f, sY, 0.0f, cY);
		Quaternion r(0.0f, 0.0f, sR, cR);

		return QuatMultiply(QuatMultiply(y, p), r);
	}








































	// =========================================== Matrix =================================================

	Matrix MatIdentity()
	{
		return Matrix
		(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	Matrix MatTranspose(const Matrix& m)
	{
		return Matrix
		(
			m.GetColumnVector(0),
			m.GetColumnVector(1),
			m.GetColumnVector(2),
			m.GetColumnVector(3)
		);
	}




	/// <summary>
	/// Calculate the determinant using elementary row operations
	/// </summary>
	/// <param name="m"> Linearly independent matrix</param>
	/// <param name="n"> Size of row and column</param>
	/// <returns>Determinant of Matrix m</returns>
	float MatDeterminant(const Matrix& m, int n)
	{
		assert(0 < n && n < 5 && "size:n not supported");

		if (n == 1)
			return m[0][0];

		// algorithm: make the matrix triangular using row operations then multiply the diagonal
		Matrix M = m;
		int sign = 1;

		int i = 0;
		for (int j = 0; j < n; j++)
		{
			int rowOfLargestLeadCoef = i;
			float curLargestCoef = abs(M[i][j]);
			for (int k = i; k < n; k++)
			{
				float coef = abs(M[k][j]);
				if (coef > curLargestCoef)
				{
					curLargestCoef = coef;
					rowOfLargestLeadCoef = k;
				}
			}

			// swap then change sign
			if (rowOfLargestLeadCoef != i)
			{
				if (i == 0)
				{
					switch (rowOfLargestLeadCoef)
					{
					case 1:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[i],
							M[2],
							M[3]
						);

						break;
					case 2:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[i],
							M[3]
						);

						break;
					case 3:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[2],
							M[i]
						);

						break;
					}
				}
				else if (i == 1)
				{
					switch (rowOfLargestLeadCoef)
					{
					case 2:
						M = Matrix(
							M[0],
							M[rowOfLargestLeadCoef],
							M[i],
							M[3]
						);

						break;
					case 3:
						M = Matrix(
							M[0],
							M[rowOfLargestLeadCoef],
							M[2],
							M[i]
						);

						break;
					}
				}
				else
				{
					// largest coef row = 3
					// i == 2
					M = Matrix(
						M[0],
						M[1],
						M[rowOfLargestLeadCoef],
						M[i]
					);
				}

				sign = -sign;
			}

			for (int r = i + 1; r < n; r++)
			{
				M[r] = (-M[r][j] * (1.0f / M[i][j] * M[i])) + M[r];
			}

			i++;
		}

		float result = 1.0f;
		for (i = 0; i < n; i++)
		{
			result *= M[i][i];
		}

		return sign * result;
	}

	/// <summary>
	/// Algorithm:
	/// Mathematics for 3D Game Programming  and Computer Graphics; Page 42: Algorthm 3.11 Gauss Jordan Elimination
	/// https://canvas.projekti.info/ebooks/Mathematics%20for%203D%20Game%20Programming%20and%20Computer%20Graphics,%20Third%20Edition.pdf
	/// </summary>
	/// <param name="m"> Invertible 4x4 Matrix</param>
	/// <param name="n"> size row and column</param>
	/// <param name="elementaryRows"> nullptr if you dont want to get 'm' elementary rows</param>
	/// <returns></returns>
	Matrix MatInverse(const Matrix& m, int n, std::vector<Matrix>* elementaryRows)
	{
		assert(1 < n && n < 5 && "size:n not supported");

		// alternative algorithm
		// Matrix CF = CofactorMatrix(Transpose(m), 4)
		// (1/Determinant(m, n)) * CF;
		// since the determinant can be evaluated by choosing any row of k of the matrix m
		// CofactorMatrix(Transpose(m), 4) and summing the products(dot product) 
		// with the entries of the k-th column(because its been transpose) of the matrix CF we can:
		// Matrix CF = CofactorMatrix(Transpose(m), 4);
		// (1/Dot(m[k], CF.GetColumnVector(k))) * CF;

		Matrix M = m;
		Matrix M_ = MatIdentity();

		for (int j = 0; j < n; j++)
		{
			int rowOfLargestLeadCoef = j;
			float curLargestCoef = abs(M[j][j]);
			for (int i = j + 1; i < n; i++)
			{
				float coef = abs(M[i][j]);
				if (coef > curLargestCoef)
				{
					curLargestCoef = coef;
					rowOfLargestLeadCoef = i;
				}
			}

			if (rowOfLargestLeadCoef != j)
			{
				if (j == 0)
				{
					switch (rowOfLargestLeadCoef)
					{
					case 1:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[j],
							M[2],
							M[3]
						);

						M_ = Matrix(
							M_[rowOfLargestLeadCoef],
							M_[j],
							M_[2],
							M_[3]
						);

						if (elementaryRows)
						{
							elementaryRows->push_back(
								Matrix(
								0, 1, 0, 0,
								1, 0, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1
							));
						}
						break;
					case 2:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[j],
							M[3]
						);

						M_ = Matrix(
							M_[rowOfLargestLeadCoef],
							M_[1],
							M_[j],
							M_[3]
						);

						if (elementaryRows)
						{
							elementaryRows->push_back(
								Matrix(
								0, 0, 1, 0,
								0, 1, 0, 0,
								1, 0, 0, 0,
								0, 0, 0, 1
							));
						}

						break;
					case 3:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[2],
							M[j]
						);

						M_ = Matrix(
							M_[rowOfLargestLeadCoef],
							M_[1],
							M_[2],
							M_[j]
						);

						if (elementaryRows)
						{
							elementaryRows->push_back(
								Matrix(
								0, 0, 0, 1,
								0, 1, 0, 0,
								0, 0, 1, 0,
								1, 0, 0, 0
							));
						}

						break;
					}
				}
				else if (j == 1)
				{
					switch (rowOfLargestLeadCoef)
					{
					case 2:
						M = Matrix(
							M[0],
							M[rowOfLargestLeadCoef],
							M[j],
							M[3]
						);

						M_ = Matrix(
							M_[0],
							M_[rowOfLargestLeadCoef],
							M_[j],
							M_[3]
						);

						if (elementaryRows)
						{
							elementaryRows->push_back(
								Matrix(
								1, 0, 0, 0,
								0, 0, 1, 0,
								0, 1, 0, 0,
								0, 0, 0, 1
							));
						}

						break;
					case 3:
						M = Matrix(
							M[0],
							M[rowOfLargestLeadCoef],
							M[2],
							M[j]
						);

						M_ = Matrix(
							M_[0],
							M_[rowOfLargestLeadCoef],
							M_[2],
							M_[j]
						);

						if (elementaryRows)
						{
							elementaryRows->push_back(
								Matrix(
								1, 0, 0, 0,
								0, 0, 0, 1,
								0, 0, 1, 0,
								0, 1, 0, 0
							));
						}

						break;
					}
				}
				else
				{
					// largest coef row = 3
					// j == 2
					M = Matrix(
						M[0],
						M[1],
						M[rowOfLargestLeadCoef],
						M[j]
					);

					M_ = Matrix(
						M_[0],
						M_[1],
						M_[rowOfLargestLeadCoef],
						M_[j]
					);

					if (elementaryRows)
					{
						elementaryRows->push_back(
							Matrix(
							1, 0, 0, 0,
							0, 1, 0, 0,
							0, 0, 0, 1,
							0, 0, 1, 0
						));
					}
				}
			}


			float temp = M[j][j];
			M[j] = 1.0f / temp * M[j];
			M_[j] = 1.0f / temp * M_[j];

			if (elementaryRows)
			{
				Matrix I = MatIdentity();
				I[j] = 1.0f / temp * I[j];
				elementaryRows->push_back(I);
			}

			for (int r = 0; r < n; r++)
			{
				if (r == j)
					continue;

				float temp = -M[r][j];
				M[r] = temp * M[j] + M[r];
				M_[r] = temp * M_[j] + M_[r];

				if (elementaryRows)
				{
					Matrix I = MatIdentity();
					I[r] = temp * I[j] + I[r];
					elementaryRows->push_back(I);
				}
			}
		}

		return M_;
	}

	/// <param name="m">Matrix</param>
	/// <param name="n">row and column of matrix</param>
	/// <param name="i">row of number to be calculated</param>
	/// <param name="j">column of number to be calculated</param>
	/// <returns></returns>
	float MatCofactor(const Matrix& m, int n, int i, int j)
	{
		assert(1 < n && n < 5 && "size:n not supported");

		float exp = (i + 1.0f) + (j + 1.0f);
		float sign = pow(-1.0f, exp);

		Matrix M = {};

		int r = 0;
		for (int x = 0; x < n; x++)
		{
			if (x == i)
				continue;

			int s = 0;
			for (int y = 0; y < n; y++)
			{
				if (y == j)
					continue;

				M[r][s] = m[x][y];

				s++;
			}

			r++;
		}


		return sign * MatDeterminant(M, n - 1);
	}

	/// <param name="m">matrix</param>
	/// <param name="n">row and column size</param>
	/// <returns></returns>
	Matrix MatCofactorMatrix(const Matrix& m, int n)
	{
		assert(1 < n && n < 5 && "size:n not supported");

		Matrix M = {};

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = MatCofactor(m, n, i, j);
			}
		}

		return M;
	}

	// Matrix Transformation

	Matrix MatScale(float x, float y, float z)
	{
		Matrix res = MatIdentity();
		res[0][0] = x;
		res[1][1] = y;
		res[2][2] = z;

		return res;
	}

	Matrix MatScale(const Vector& v)
	{
		return MatScale(v.x, v.y, v.z);
	}

	Matrix MatTranslate(float x, float y, float z)
	{
		Matrix res = MatIdentity();
		res[3][0] = x;
		res[3][1] = y;
		res[3][2] = z;

		return res;
	}

	Matrix MatTranslate(const Vector& v)
	{
		return MatTranslate(v.x, v.y, v.z);
	}

	Matrix MatRotationX(float angle)
	{
		return Matrix(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, cosf(angle), sinf(angle), 0.0f,
			0.0f, -sinf(angle), cosf(angle), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	Matrix MatRotationY(float angle)
	{
		return Matrix(
			cosf(angle), 0.0f, -sinf(angle), 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			sinf(angle), 0.0f, cosf(angle), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	Matrix MatRotationZ(float angle)
	{
		return Matrix(
			cosf(angle), sinf(angle), 0.0f, 0.0f,
			-sinf(angle), cosf(angle), 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	/// <summary>
	/// Rodrigues' Rotation Formula 
	/// P_rot = P * cos(angle) + Cross(axis, P) * sin(angle) + axis * Dot(axis, P) * (1 - cos(angle))
	/// </summary>
	/// <param name="angle"></param>
	/// <param name="axis"></param>
	/// <returns></returns>
	Matrix MatRotationAxis(float angle, const Vector& axis)
	{
		// P_rot = P_para + cos(angle)P_perp + sin(angle)axis x P_perp
		// P_rot = Pcos(angle) + (axis x P)sin(angle) + axis(axis . P)(1 - cos(angle))
		// =
		//             | 1  0  0 |	              |  0		  axis.z  -axis.y |						|   axis.x^2     axis.x*axis.y  axis.x*axis.z  |
		// Pcos(angle) | 0  1  0 |  + Psin(angle) | -axis.z     0	   axis.x | + P(1 - cos(angle)) | axis.x*axis.y    axis.y^2     axis.y*axis.z  |
		//             | 0  0  1 |	              |  axis.y	 -axis.x	 0    |                     | axis.x*axis.z  axis.y*axis.z    axis.z^2	   |
		// =
		// note: c = cos(angle) s = sin(angle)
		//   | c + (1-c)axis.x^2              (1-c)axis.x*axis.y + s*axis.z  (1-c)axis.x*axis.z - s*axis.y  |
		// P | (1-c)axis.x*axis.y - s*axis.z  c + (1-c)axis.y^2              (1-c)axis.y*axis.z + s*axis.x	|
		//	 | (1-c)axis.x*axis.z + s*axis.y  (1-c)axis.y*axis.z - s*axis.x  c + (1-c)axis.z^2				|

		Vector a = Vec3Normalized(axis);
		float s = sinf(angle);
		float c = cosf(angle);
		Vector temp = (1 - c) * a;

		return Matrix(
			c + temp.x * a.x, temp.x * a.y + s * a.z, temp.x * a.z - s * a.y, 0,
			temp.x * a.y - s * a.z, c + temp.y * a.y, temp.y * a.z + s * a.x, 0,
			temp.x * a.z + s * a.y, temp.y * a.z - s * a.x, c + temp.z * a.z, 0,
			0, 0, 0, 1
		);
	}

	/// <summary>
	/// Matrix representation of quaternion rotation. Read the "Quaternion section of Mathematics for 3D Game Programming and Computer Graphics"
	/// </summary>
	/// <param name="q"></param>
	/// <returns></returns>
	Matrix MatRotationQuaternion(const Quaternion& q)
	{
		Quaternion temp = 2 * q;
		return Matrix
		(
			1.0f - temp.y * q.y - temp.z*q.z, temp.x * q.y + temp.w * q.z,        temp.x * q.z - temp.w * q.y,        0.0f,
			temp.x * q.y - temp.w * q.z,      1.0f - temp.x * q.x - temp.z * q.z, temp.y * q.z + temp.w * q.x,        0.0f,
			temp.x * q.z + temp.w * q.y,      temp.y * q.z - temp.w * q.x,        1.0f - temp.x * q.x - temp.y * q.y, 0.0f,
			0.0f,                             0.0f,                               0.0f,                               1.0f
		);
	}
}