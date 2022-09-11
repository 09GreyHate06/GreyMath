#include "GreyMath.h"
#include <DirectXMath.h>

using namespace GreyMath;
int main()
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

	return 0;
}