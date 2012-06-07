//Standard math class. I have not created this completly by myself. All SSE parts are done by me.


#ifndef miscmatte
#define miscmatte

#include "math.h"
#include "stdlib.h"
#include "xmmintrin.h"

inline float Rand( float a_Range ) { return ((float)rand() / RAND_MAX) * a_Range; }

#define DOT(A,B)		(A.x*B.x+A.y*B.y+A.z*B.z)
#define NORMALIZE(A)	{float l=InvSqrt(A.x*A.x+A.y*A.y+A.z*A.z);A.x*=l;A.y*=l;A.z*=l;}
#define LENGTH(A)		(sqrtf(A.x*A.x+A.y*A.y+A.z*A.z))
#define SQRLENGTH(A)	(A.x*A.x+A.y*A.y+A.z*A.z)
#define SQRDISTANCE(A,B) ((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z))

#define EPSILON			0.0001f
//#define MAXDEPTH		2

#define PI				3.141592653589793238462f

float MAX(float,float);
float Q_rsqrt(float);
float InvSqrt(float x);

typedef struct vec4{
    float x,y,z,w;
} Vec4;

__m128 dot_sse(__m128 &x1, __m128 &y1, __m128 &z1,__m128 &x2, __m128 &y2, __m128 &z2);
__m128 length_sse(__m128 &x, __m128 &y, __m128 &z);
void normalize_sse(__m128 &x, __m128 &y,__m128 &z);

class Vec3
{
public:
	Vec3() : x( 0.0f ), y( 0.0f ), z( 0.0f ) {};
	Vec3( float a_X, float a_Y, float a_Z ) : x( a_X ), y( a_Y ), z( a_Z ) {};
	void Set( float a_X, float a_Y, float a_Z ) { x = a_X; y = a_Y; z = a_Z; }
	void Normalize() { float l = 1.0f / Length(); x *= l; y *= l; z *= l; }
	float Length() { return (float)sqrt( x * x + y * y + z * z ); }
	float SqrLength() { return x * x + y * y + z * z; }
	float Dot( Vec3 a_V ) { return x * a_V.x + y * a_V.y + z * a_V.z; }
	Vec3 Cross( Vec3 b ) { return Vec3( y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x ); }
	void operator += ( Vec3& a_V ) { x += a_V.x; y += a_V.y; z += a_V.z; }
	void operator += ( Vec3* a_V ) { x += a_V->x; y += a_V->y; z += a_V->z; }
	void operator -= ( Vec3& a_V ) { x -= a_V.x; y -= a_V.y; z -= a_V.z; }
	void operator -= ( Vec3* a_V ) { x -= a_V->x; y -= a_V->y; z -= a_V->z; }
	void operator *= ( float f ) { x *= f; y *= f; z *= f; }
	void operator *= ( Vec3& a_V ) { x *= a_V.x; y *= a_V.y; z *= a_V.z; }
	void operator *= ( Vec3* a_V ) { x *= a_V->x; y *= a_V->y; z *= a_V->z; }
	Vec3 operator- () const { return Vec3( -x, -y, -z ); }
	friend Vec3 operator + ( const Vec3& v1, const Vec3& v2 ) { return Vec3( v1.x + v2.x, v1.y + v2.y, v1.z + v2.z ); }
	friend Vec3 operator - ( const Vec3& v1, const Vec3& v2 ) { return Vec3( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z ); }
	friend Vec3 operator + ( const Vec3& v1, Vec3* v2 ) { return Vec3( v1.x + v2->x, v1.y + v2->y, v1.z + v2->z ); }
	friend Vec3 operator - ( const Vec3& v1, Vec3* v2 ) { return Vec3( v1.x - v2->x, v1.y - v2->y, v1.z - v2->z ); }
	friend Vec3 operator * ( const Vec3& v, float f ) { return Vec3( v.x * f, v.y * f, v.z * f ); }
	friend Vec3 operator * ( const Vec3& v1, Vec3& v2 ) { return Vec3( v1.x * v2.x, v1.y * v2.y, v1.z * v2.z ); }
	friend Vec3 operator * ( float f, const Vec3& v ) { return Vec3( v.x * f, v.y * f, v.z * f ); }
	union
	{
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float cell[3]; };
	};
};

typedef Vec3 Color;

struct Plane{
    Vec3 N;
    float D;
    float cell[4];
    
	Plane() : N( 0, 0, 0 ), D( 0 ) {};
	Plane( Vec3 a_Normal, float a_D ) : N( a_Normal ), D( a_D ) {};
};

#endif