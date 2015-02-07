#ifndef _MATRIX_3_4
#define _MATRIX_3_4 3

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#include "global.h"
#include "vector.h"
#include "mesh.h"


class Matrix{
public:
	double **element;
	Matrix();
	Matrix(Matrix3 &m3);
	Matrix(Matrix4 &m4);
	void nrerror(char error_text[]);
	double *vectorm(int nl,int nh);
	double **matrix(int nrl=0,int nrh=3,int ncl=0,int nch=3) ;
	void free_vector(double *v,int nl,int nh);
	void free_matrix(double **m,int nrl=0,int nrh=3,int ncl=0,int nch=3);
	void inverse(double **mat, int dim=3) ;
	void ludcmp(double **a, int n, int *indx, double *d) ;
	void lubksb(double **a, int n, int *indx, double *b) ;
	void print(double **m,int nrl=0,int nrh=3,int ncl=0,int nch=3);
	Matrix3 AddMatrix3(Matrix3 &M1,Matrix3 &M2);
	Matrix3 Add(Matrix3 &M1,float fValue);
	Matrix3 Mat3Div(Matrix3 &M,float n);
	Matrix3 Sqrt(Matrix3 &M);
	float Mat3lemda(Matrix3 &M);
	float ComMat3lemda(Matrix3 &M,Vector3 *v,int iteration);
};




#endif
