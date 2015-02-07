// Standard include files
#ifndef _VECTOR_3_4
#define _VECTOR_3_4 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#include "global.h"
class Vector3;

class Vector4 {
public: 
	float x, y, z,w;
	Vector4();
	Vector4(Vector3 &p);
};

class Vector3 {
public: 
	float x, y, z;
	Vector3(void);
	Vector3(float xx,float yy,float zz);
	Vector3(Vector4 &p);
	float amplitude(Vector3 vector);
	float Distance(Vector3 vector1,Vector3 vector2);
	Vector3 Direction(Vector3 vector1,Vector3 vector2);
	Vector3 Normalize(Vector3 vector1);
	float Max(Vector3 &v);
	Vector3 mean(Vector3 *vertex,long nverts);
	Vector3 CrossOver(Vector3 vector1,Vector3 vector2);//vector1 cross vector2
	Vector4 Mat4MulV4(Matrix4 &m4,Vector4 &v4);//multiply matrix4 with vector4
	Vector3 Mat3MulV3(Matrix3 &m3,Vector3 &v3);//multiply matrix3 with vector3
	Vector3 Rotate(Vector3 v,float a,float b,float r,float dx,float dy,float dz);
	Vector3 operator + (Vector3& vector1);	// 'Vector3 + Vector3'
	Vector3 operator - (Vector3& vector1);	// 'Vector3 - Vector3'
	Matrix3 operator * (Vector3& vector1);	// 'Vector3 * Vector3'
	Vector3 operator * (float dValue);	// 'Vector3 * double'	
	//Vector3 operator *= (float dValue);	// 'Vector3 *= double'	
	Vector3& operator = (Vector3& vector1);	// 'Vector3 = Vector3'	
	Vector3& operator += (Vector3& vector1);	// 'Vector3 += Vector3'	
	float operator / (Vector3& vector1);	// 'dot metrixs of Vector3 and Vector3'	
	Vector3 operator / (float dValue);// Vector3 / double	
	//Vector3 operator /= (float dValue);// Vector3 /= double
} ;




#endif
