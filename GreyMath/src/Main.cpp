#include "GreyMath.h"
#include <DirectXMath.h>

using namespace GreyMath;
int main()
{
	Vector v0(0.5, 0.0f, 23.0f, 0.0f);
	Vector v1(0.5, 0.0f, 23.0f, 0.0f);

	//v0 = Normalized(v0);
	//v1 = Normalized(v1);
	std::cout << Dot(5 * v0, v1) << '\n';
	std::cout << 5 * Dot(v1, v0) << '\n';

	return 0;
}