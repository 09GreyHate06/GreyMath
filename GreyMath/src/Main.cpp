#include "GreyMath.h"
#include <DirectXMath.h>
#include <array>

using namespace GreyMath;
using namespace DirectX;


#define G_PI 3.14159265359f

int main()
{
	//Vector N(0, 1, 0, 0);
	//Vector P(0, 6, 0, 0);
	//Vector Q(1, 2, 23, 0);

	//float D = Vec3Dot(-N, P);

	//Matrix planeTransform = MatRotationQuaternion(QuatRotationRollPitchYaw(0.0f, 0.0f, -90.0f * G_PI / 180.0f)) *
	//	MatTranslate(6.0f, 0.0f, 0.0f);

	//Vector L = N;
	//L.w = D;

	//std::cout << L << '\n';
	////std::cout << MatInverse(planeTransform, 4, nullptr) << '\n';
	//std::cout << Vec4Transform(L, MatTranspose(MatInverse(planeTransform, 4, nullptr))) << '\n';

	Vector v(255, 0, 0, 0);
	Vector n(255, 255, 255, 0);
	float w = 1.0f / 4.0f;

	Vector res;
	for (int i = 0; i < 4; i++)
	{
		if (i < 2)
			res += w * v;
		else
			res += w * n;
	}
	std::cout << res << '\n';
	
	return 0;
}