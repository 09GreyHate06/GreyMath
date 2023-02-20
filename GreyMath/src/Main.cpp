#include "GreyMath.h"
#include <DirectXMath.h>
#include <array>
using namespace GreyMath;
using namespace DirectX;

void DotProductAndProjections()
{
	Vector P(6.0f, -75.64f, 32.75f, 0.0f);
	Vector Q(-123.2f, 23.12f, 90.5f, 0.0f);
	Vector R(6.0f, -3.0f, 7.0f, 0.0f);

	std::cout << "=========== Dot Product Properties ===========\n";
	std::cout << "P" << P << '\n';
	std::cout << "Q" << Q << '\n';
	std::cout << "R" << R << '\n';
	std::cout << "a = 5\n\n";

	std::cout << "P . Q = Q . P" << '\n';
	std::cout << Dot(P, Q) << " = " << Dot(Q, P) << "\n\n";

	std::cout << "(a P) . Q = a (P . Q)\n";
	std::cout << Dot(5.0f * P, Q) << " = " << 5.0f * Dot(P, Q) << "\n\n";

	std::cout << "P . (Q + R) = P . Q + P . R\n";
	std::cout << Dot(P, Q + R) << " = " << Dot(P, Q) + Dot(P, R) << "\n\n";

	std::cout << "P . P = ||P||^2\n";
	std::cout << Dot(P, P) << " = " << powf(Magnitude(P), 2.0f) << "\n\n";

	std::cout << "P . Q <= ||P|| ||Q||\n";
	std::cout << Dot(P, Q) << " <= " << Magnitude(P) * Magnitude(Q) << '\n';


	Vector U(6.0f, -3.0f, 9.0f, 0.0f);
	Vector V(4.0f, -1.0f, 8.0f, 0.0f);
	std::cout << "\n================= Projection =================\n";
	std::cout << "U" << U << '\n';
	std::cout << "V" << V << "\n\n";

	Vector projUontoV = ProjectionOfV0OntoV1(U, V);

	std::cout << "Projection of U onto V\n";
	std::cout << projUontoV << "\n\n";

	std::cout << "The vector component of U orthogonal to V\n";
	std::cout << U - projUontoV << '\n';
	std::cout << "Check: V . projUontoV = " << Dot(Normalized(V), Normalized(U - projUontoV)) << "\n\n"; // floating point error
}

template<typename ...ArgsT>
float LCM(ArgsT... fs)
{
	std::array<float, sizeof...(ArgsT)> const nums{ fs... };
	return nums[0];
}

Matrix RotateX(float angle)
{
	return Matrix(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, cosf(angle), sinf(angle), 0.0f,
		0.0f, -sinf(angle), cosf(angle), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}



#define G_PI 3.14159265359

float f(float x)
{
	return powf(2.0f, x);
}

int main()
{
	//Quaternion a(2, 3, 4, 1);
	//Quaternion b(-2, -3, -4, 1);

	//float real = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
	//float i = (a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y);
	//float j = (a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x);
	//float k = (a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w);

	//std::cout << "real: " << real << '\n';
	//std::cout << "i: " << i << '\n';
	//std::cout << "j: " << j << '\n';
	//std::cout << "k: " << k << '\n';
	//std::cout << "result: " << real + i + j + k << '\n';

	//Quaternion uq = (Quaternion(2, 3, 4, 1));
	//std::cout << uq << '\n';
	//std::cout << Vector(-uq.x, -uq.y, -uq.z, uq.w) / Dot(uq, uq) << '\n';
	//std::cout << Vector(-uq.x, -uq.y, -uq.z, uq.w) << '\n';

	//Vector a(1, 2, 3, 4);
	//Vector b(5, 6, 7, 8);
	//std::cout << Dot(Cross(a, b), a) << '\n';
	//std::cout << Dot(4*b, a) << '\n';

	Vector P(2, 1, 0, 0);
	std::cout << P * cosf(45.0f * G_PI / 180.0f)  + Vector(-P.y, P.x, 0.0f, 0.0f) * sinf(45.0f * G_PI / 180.0f) << '\n';
	std::cout << P * cosf(90.0f * G_PI / 180.0f)  + Vector(-P.y, P.x, 0.0f, 0.0f) * sinf(90.0f * G_PI / 180.0f) << '\n';
	return 0;
}