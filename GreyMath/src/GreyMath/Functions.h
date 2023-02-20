#pragma once
#include "Types.h"
#include "Operators.h"
#include <array>
#include <vector>
#include <assert.h>

namespace GreyMath
{
	float Determinant(const Matrix&, int);
	Matrix Inverse(const Matrix&, int, std::vector<Matrix>*);
	Matrix Transpose(const Matrix&);
	Matrix CofactorMatrix(const Matrix&, int);

	// =========================================== Vector =================================================

	// = sqrt(x^2 + y^2 ... n^2)
	float Magnitude(const Vector& v)
	{
		return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	}

	// v / magnitude(v)
	Vector Normalized(const Vector& v)
	{
		return v / Magnitude(v);
	}

	// lhs.x^2 + lhs.y^2 + lhs.z^2 + lhs.w^2 == ||lhs|| * ||rhs|| * cos(a)
	constexpr float Dot(const Vector& lhs, const Vector& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
	}

	// Pseudodeterminant
	// |   i      j       k   |
	// | lhs.x  lhs.y   lhs.z | 
	// | rhs.x  rhs.y   rhs.z |
	// = i(lhs.y * rhs.z - lhs.z * rhs.y) - j(lhs.x * rhs.z - lhs.z * rhs.x) + k(lhs.x * rhs.y - lhs.y * rhs.x)
	// ||result|| = ||lhs|| * ||rhs|| * sin(a)
	Vector Cross(const Vector& lhs, const Vector& rhs)
	{
		//Vector i(1.0f, 0.0f, 0.0f, 0.0f);
		//Vector j(0.0f, 1.0f, 0.0f, 0.0f);
		//Vector k(0.0f, 0.0f, 1.0f, 0.0f);
		//
		//return i * (lhs.y * rhs.z - lhs.z * rhs.y) - j * (lhs.x * rhs.z - lhs.z * rhs.x) + k * (lhs.x * rhs.y - lhs.y * rhs.x);

		return Vector(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x, 0.0f);
	}

	// (v0 . v1) / ||v1||^2) * v1
	Vector ProjectionOfV0OntoV1(const Vector& v0, const Vector& v1)
	{
		//return (Dot(v0, v1) / Magnitude(v1)) * Normalized(v1);

		float v1MagSqrd = Dot(v1, v1);
		return (Dot(v0, v1) / v1MagSqrd) * v1;
	}

	// Algorithm: https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
	std::array<Vector, 3> GramSchmidtOrthogonalization(const std::array<Vector, 3>& basis)
	{
		std::array<Vector, 3> result;
		result[0] = basis[0];

		for (int i = 1; i < 3; i++)
		{
			Vector projectionSum(0.0f, 0.0f, 0.0f, 0.0f);
			for (int k = 0; k < i; k++)
			{
				projectionSum = projectionSum + ProjectionOfV0OntoV1(basis[i], result[k]);
			}

			result[i] = basis[i] - projectionSum;
		}

		return result;
	}

	/// <summary>
	/// Algorithm:
	/// Mathematics for 3D Game Programming  and Computer Graphics; Page 36: Algorthm 3.6
	/// https://canvas.projekti.info/ebooks/Mathematics%20for%203D%20Game%20Programming%20and%20Computer%20Graphics,%20Third%20Edition.pdf
	/// </summary>
	/// <param name="coefs">Coefficient Matrix</param>
	/// <param name="constants">Constants or the right hand side of the system</param>
	/// <param name="n">number of variables or size of Matrix_nxn</param>
	/// <returns></returns>
	Vector GussianElimination(const Matrix& coefs, const Vector& constants, int n)
	{
		Matrix M = coefs;
		Vector C = constants;

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

						C = Vector(C[rowOfLargestLeadCoef], C[i], C[2], C[3]);
						break;
					case 2:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[i],
							M[3]
						);

						C = Vector(C[rowOfLargestLeadCoef], C[1], C[i], C[3]);
						break;
					case 3:
						M = Matrix(
							M[rowOfLargestLeadCoef],
							M[1],
							M[2],
							M[i]
						);

						C = Vector(C[rowOfLargestLeadCoef], C[1], C[2], C[i]);
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

						C = Vector(C[0], C[rowOfLargestLeadCoef], C[i], C[3]);
						break;
					case 3:
						M = Matrix(
							M[0],
							M[rowOfLargestLeadCoef],
							M[2],
							M[i]
						);

						C = Vector(C[0], C[rowOfLargestLeadCoef], C[2], C[i]);
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

					C = Vector(C[0], C[1], C[rowOfLargestLeadCoef], C[i]);
				}
			}

			// set leading coef to 1
			{
				float temp = 1.0f / M[i][j];
				M[i] = temp * M[i];
				C[i] = temp * C[i];
			}

			// set all under the leading coef to 0
			for (int r = 0; r < n; r++)
			{
				if (r == i)
					continue;

				float temp = -M[r][j];
				M[r] = temp * M[i] + M[r];
				C[r] = temp * C[i] + C[r];
			}

			i++;
		}

		return C;
	}

	/// <summary>
	/// Algorithm:
	/// Mathematics for 3D Game Programming  and Computer Graphics; Page 53: Equation 3.69:
	/// https://canvas.projekti.info/ebooks/Mathematics%20for%203D%20Game%20Programming%20and%20Computer%20Graphics,%20Third%20Edition.pdf
	/// </summary>
	/// <param name="coefs">Coefficient Matrix</param>
	/// <param name="constants">Constants or the right hand side of the system</param>
	/// <param name="n">number of variables or size of Matrix_nxn</param>
	/// <returns></returns>
	Vector CramersRule(const Matrix& coefM, const Vector& constants, int n)
	{
		// alt algorithm
		// ax + by + cz = d
		// D = Determinant of 3x3 matrix of coef [a_i, b_i, c_i]
		// Dx = Determinant of 3x3 matrix of coef and constant [d_i, b_i, c_i]
		// Dy = Determinant of 3x3 matrix of coef and constant [a_i, d_i, c_i]
		// Dz = Determinant of 3x3 matrix of coef and constant [a_i, b_i, d_i]
		// x = Dx/D; y = Dy/D; z= Dz/D;
		//
		//Matrix CM = CofactorMatrix(coefM, 3);
		//float detM = Determinant(coefM, 3);
		//Vector res;
		//for (int i = 0; i < n; i++)
		//{
		//	for (int j = 0; j < n; j++)
		//	{
		//		res[i] += CM[j][i] * constants[j];
		//	}
		//
		//	res[i] *= 1 / detM;
		//}
		//
		//return res;

		Matrix coefMI = Inverse(coefM, n, nullptr);
		return constants * Transpose(coefMI);
	}

	// =========================================== Matrix =================================================

	Matrix Identity()
	{
		return Matrix
		(
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	Matrix Transpose(const Matrix& m)
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
	float Determinant(const Matrix& m, int n)
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
	Matrix Inverse(const Matrix& m, int n, std::vector<Matrix>* elementaryRows)
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
		Matrix M_ = Identity();

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
				Matrix I = Identity();
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
					Matrix I = Identity();
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
	float Cofactor(const Matrix& m, int n, int i, int j)
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

		
		return sign * Determinant(M, n - 1);
	}

	/// <param name="m">matrix</param>
	/// <param name="n">row and column size</param>
	/// <returns></returns>
	Matrix CofactorMatrix(const Matrix& m, int n)
	{
		assert(1 < n && n < 5 && "size:n not supported");

		Matrix M = {};

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = Cofactor(m, n, i, j);
			}
		}

		return M;
	}





	// Matrix Transformation

	Matrix Scale(float x, float y, float z)
	{
		Matrix res = Identity();
		res[0][0] = x;
		res[1][1] = y;
		res[2][2] = z;

		return res;
	}

	Matrix Scale(const Vector& v)
	{
		return Scale(v.x, v.y, v.z);
	}

	Matrix Translate(float x, float y, float z)
	{
		Matrix res = Identity();
		res[3][0] = x;
		res[3][1] = y;
		res[3][2] = z;

		return res;
	}

	Matrix Translate(const Vector& v)
	{
		return Translate(v.x, v.y, v.z);
	}

	Matrix RotateY(float angle)
	{
		return Matrix(
			cosf(angle), 0.0f, -sinf(angle), 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			sinf(angle), 0.0f, cosf(angle), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);
	}

	Matrix RotateZ(float angle)
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
	Matrix RotateAxis(float angle, const Vector& axis)
	{
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

		Vector a = Normalized(axis);
		float s = sinf(angle);
		float c = cosf(angle);
		Vector temp = (1 - c) * a;

		return Matrix(
			c + temp.x * a.x,       temp.x * a.y + s * a.z, temp.x * a.z - s * a.y, 0,
			temp.x * a.y - s * a.z, c + temp.y * a.y,       temp.y * a.z + s * a.x, 0,
			temp.x * a.z + s * a.y, temp.y * a.z - s * a.x, c + temp.z * a.z,       0,
			0,                      0,                      0,                      1
		);
	}

	Vector CalcTriangleNormal(const Vector& a, const Vector& b, const Vector& c)
	{
		return Cross(b - a, c - a);
	}

	Vector RotateEulerFormula(float angle, const Vector& axis, const Vector& v)
	{
		return cosf(angle) * v + sinf(angle) * Cross(axis, v);
	}
}