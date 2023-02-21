#include "GreyMath.h"
#include <DirectXMath.h>
#include <array>

using namespace GreyMath;
using namespace DirectX;


#define G_PI 3.14159265359f

Quaternion Test(const Vector& v)
{
	float theta = 360.0f * G_PI / 180.0f;

	Vector n = Vec3Normalized(Vector(1.0f, 1.0f, 1.0f, 0.0f));
	Quaternion q(sinf(theta / 2.0f) * n.x, sinf(theta / 2.0f) * n.y, sinf(theta / 2.0f) * n.z, cosf(theta / 2.0f));

	return QuatMultiply(QuatMultiply(q, v), QuatInverse(q));
}

int main()
{
	Vector v(1.0f, 3.0f, 4.0f, 0.0f);
	Matrix rot = MatRotationQuaternion(QuatRotationRollPitchYaw(90.0f * G_PI / 180.0f, 90.0f * G_PI / 180.0f, 90.0f * G_PI / 180.0f));
	std::cout << Vec3Transform(v, rot) << '\n';

	XMVECTOR v_ = XMVectorSet(1.0f, 3.0f, 4.0f, 0.0f);
	XMMATRIX rot_ = XMMatrixRotationQuaternion(XMQuaternionRotationRollPitchYaw(90.0f * G_PI / 180.0f, 90.0f * G_PI / 180.0f, 90.0f * G_PI / 180.0f));


	XMFLOAT4 v__;
	XMFLOAT4X4 rot__;
	XMStoreFloat4x4(&rot__, rot_);
	XMStoreFloat4(&v__, XMVector3Transform(v_, rot_));
	std::cout << (*(Vector*)&v__) << '\n';

	return 0;
}