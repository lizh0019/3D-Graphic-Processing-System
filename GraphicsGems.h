/* 
 * GraphicsGems.h  
 * Version 1.0 - Andrew Glassner
 * from "Graphics Gems", Academic Press, 1990
 */
#include <MATH.H>
#ifndef _GRAPHICS_GEMS

#define _GRAPHICS_GEMS 1

/*********************/
/* 2d geometry types */
/*********************/

typedef struct Point2Struct {	/* 2d point */
	double x, y;
	} Point2;
typedef Point2 Vector2;

typedef struct IntPoint2Struct {	/* 2d integer point */
	int x, y;
	} IntPoint2;

typedef struct Box2dStruct {		/* 2d box */
	Point2 min, max;
	} Box2;
	

/*********************/
/* 3d geometry types */
/*********************/


typedef struct IntPoint3Struct {	/* 3d integer point */
	int x, y, z;
	} IntPoint3;

class Matrix4{
public:
    Matrix4(Matrix4 &M){/* 4-by-4 matrix */
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
			element[i][j]=M.element[i][j];
			}
		}
	}
	Matrix4(){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
			element[i][j]=0;
			}
		}
	}
	Matrix4(float **m){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
			element[i][j]=m[i][j];
			}
		}
	}
	Matrix4(float a,float b,float r,float dx,float dy,float dz){
		element[0][0]=cos(a)*cos(b);
		element[0][1]=(cos(a)*sin(b)*sin(r)-sin(a)*cos(r));
		element[0][2]=(cos(a)*sin(b)*cos(r)+sin(a)*sin(r));
		element[0][3]=dx;
		element[1][0]=sin(a)*cos(b);
		element[1][1]=(sin(a)*sin(b)*sin(r)+cos(a)*cos(r));
		element[1][2]=(sin(a)*sin(b)*cos(r)-cos(a)*sin(r));
		element[1][3]=dy;
		element[2][0]=(-1*sin(b));
		element[2][1]=cos(b)*sin(r);
		element[2][2]=cos(b)*cos(r);
		element[2][3]=dz;
		element[3][0]=element[3][1]=element[3][2]=0;
		element[3][3]=1;
		}
	double element[4][4];
};
class Matrix3{
public:
    Matrix3(Matrix3 &M){/* 3-by-3 matrix */
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			element[i][j]=M.element[i][j];
			}
		}
	}
	Matrix3(){
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			element[i][j]=0;
			}
		}
	}
	Matrix3(float **m){
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			element[i][j]=m[i][j];
			}
		}
	}
	double element[3][3];
};



extern double V2SquaredLength(), V2Length();
extern double V2Dot(), V2DistanceBetween2Points(); 
extern Vector2 *V2Negate(), *V2Normalize(), *V2Scale(), *V2Add(), *V2Sub();
extern Vector2 *V2Lerp(), *V2Combine(), *V2Mul(), *V2MakePerpendicular();
extern Vector2 *V2New(), *V2Duplicate();
extern Point2 *V2MulPointByMatrix();
extern Matrix3 *V2MatMul();

extern double V3SquaredLength(), V3Length();
extern double V3Dot(), V3DistanceBetween2Points();

extern Matrix4 *V3MatMul();

extern double RegulaFalsi(), NewtonRaphson(), findroot();

#endif