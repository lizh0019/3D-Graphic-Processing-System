#ifndef _ANNEAL_REGISTRATION
#define _ANNEAL_REGISTRATION 5

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#include "global.h"
#include "mesh.h"
#include "globalfunc.h"
#include "vector.h"
#include "watermark.h"
class Anneal{
public:
#define TFACTR  0.5
#define MBIG 1000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0f/MBIG)
#define MMAX 2
float y,T;
long idum;
Anneal();
int bolziman(float de,float t);
float ran(long *idum);
int irbit(unsigned long *iseed);
float move(Mesh *mesh1,Mesh *mesh2,Mesh *tempmesh,int n,int s,float x[],float t[],float c[],float min,float *up,float *low);
Mesh *registration(Mesh *mesh1,Mesh *mesh2);
};





#endif
