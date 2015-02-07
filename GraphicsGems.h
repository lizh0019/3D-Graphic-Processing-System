/* 
 * GraphicsGems.h  
 * Version 1.0 - Andrew Glassner
 * from "Graphics Gems", Academic Press, 1990
 */
#include <MATH.H>
#ifndef GG_H

#define GG_H 1
#define   MAXDATA 1.0e+8
#define   TINY 1.0e-20
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

//typedef struct Matrix3Struct {	/* 3-by-3 matrix */
//	double element[3][3];
//	} Matrix3;

typedef struct Box2dStruct {		/* 2d box */
	Point2 min, max;
	} Box2;
	

/*********************/
/* 3d geometry types */
/*********************/

/*
class Vector3 {	// 3d point 
public: 
	double x, y, z;
	Vector3(void):x(0),y(0),z(0){};
	Vector3 Scale(Vector3 vertex,float scale){vertex.x=vertex.x*scale;vertex.y=vertex.y*scale;vertex.z=vertex.z*scale;return vertex;}
	Vector3 AddVector3(Vector3 vector1,Vector3 vector2){Vector3 vector;vector.x=vector1.x+vector2.x;vector.y=vector1.y+vector2.y;vector.z=vector1.z+vector2.z;return vector;}
	Vector3 SubVector3(Vector3 vector1,Vector3 vector2){Vector3 vector;vector.x=vector1.x-vector2.x;vector.y=vector1.y-vector2.y;vector.z=vector1.z-vector2.z;return vector;}
	float Magnitude(Vector3 vector){return sqrt((vector.x)*(vector.x)+(vector.y)*(vector.y)+(vector.z)*(vector.z));}
	float Distance(Vector3 vector1,Vector3 vector2){return Magnitude(SubVector3(vector1,vector2));}
	Vector3 Direction(Vector3 vector1,Vector3 vector2){return Scale(SubVector3(vector1,vector2),1.0/Distance(vector1,vector2));}
	Vector3 Normalize(Vector3 vector1){
		if (vector.Magnitude(vector1)<1e-5) return vector.Scale(vector1,1e+5);
		return vector.Scale(vector1,1.0/vector.Magnitude(vector1));
	}
} vector,Point3;*/

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


/***********************/
/* one-argument macros */
/***********************/

/* absolute value of a */
#define ABS(a)		(((a)<0) ? -(a) : (a))

/* round a to nearest int */
#define ROUND(a)	((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a)		(((a)<0) ? -1 : (a)>0 ? 1 : 0)	

/* take binary sign of a, either -1, or 1 if >= 0 */
#define SGN(a)		(((a)<0) ? -1 : 1)

/* shout if something that should be true isn't */
#define ASSERT(x) \
if (!(x)) fprintf(stderr," Assert failed: x\n");

/* square a */
#define SQR(a)		((a)*(a))	


/***********************/
/* two-argument macros */
/***********************/

/* find minimum of a and b */
#define MIN(a,b)	(((a)<(b))?(a):(b))	

/* find maximum of a and b */
#define MAX(a,b)	(((a)>(b))?(a):(b))	

/* swap a and b (see Gem by Wyvill) */
#define SWAP(a,b)	{ (a)=a^(b); (b)=b^(a); (a)^=a^(b); }

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) */
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

/* clamp the input to the specified range */
#define CLAMP(v,l,h)	((v)<(l) ? (l) : (v) > (h) ? (h) : v)


/****************************/
/* memory allocation macros */
/****************************/

/* create a new instance of a structure (see Gem by Hultquist) */
#define NEWSTRUCT(x)	(struct x *)(malloc((unsigned)sizeof(struct x)))

/* create a new instance of a type */
#define NEWTYPE(x)	(x *)(malloc((unsigned)sizeof(x)))


/********************/
/* useful constants */
/********************/

#define PI		3.141592653589793	/* the venerable pi */
#define PITIMES2	6.283185	/* 2 * pi */
#define PIOVER2		1.570796	/* pi / 2 */
#define E		2.718282	/* the venerable e */
#define SQRT2		1.414214	/* sqrt(2) */
#define SQRT3		1.732051	/* sqrt(3) */
#define GOLDEN		1.618034	/* the golden ratio */
#define DTOR		0.017453	/* convert degrees to radians */
#define RTOD		57.29578	/* convert radians to degrees */


/************/
/* booleans */
/************/

#define TRUE		1
#define FALSE		0
#define ON		1
#define OFF 		0
//typedef int boolean;			/* boolean data type */
//typedef boolean flag;			/* flag data type */

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