#pragma once
#include <array>
#include "Types.h"
#include "Functions.h"

namespace GreyMath
{
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
				projectionSum = projectionSum + Vec3ProjectionOfV0OntoV1(basis[i], result[k]);
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
		assert(1 < n && n < 5 && "size:n not supported");

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
		assert(1 < n && n < 5 && "size:n not supported");

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

		Matrix coefMI = MatInverse(coefM, n, nullptr);
		return Vec4Transform(constants, MatTranspose(coefMI));
	}
}