#ifndef _WATERMARK_OPERATIONS
#define _WATERMARK_OPERATIONS 1
// Standard include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#ifndef MaxNeighb
#define MaxNeighb 250
#endif
#include "global.h"
#include "vector.h"
#include "globalfunc.h"
////////////////////////////////////////////////////////////
// Classes
////////////////////////////////////////////////////////////

class Recover;
class Embed;
class extract;
class Permute;
class Attack;
class Embed{
public:
Mesh *EmbedWatermark0(Mesh *mesh,unsigned char *Watermark,float alpha);//embedding method of literature [33]
Mesh *EmbedWatermark1(Mesh *mesh,unsigned char *Watermark,float alpha);//my embedding method 1
Mesh *EmbedWatermark2(Mesh *mesh,unsigned char *Watermark,float alpha);//my embedding method 2
};

class Extract{
public:
	void ExtractWatermark0(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark);
	void ExtractWatermark1(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark);
};
class Attack{
public:
	Mesh *ChangeView(Mesh *mesh,Mesh *tempmesh,float a=0,float b=0,float r=0,float dx=0,float dy=0,float dz=0);
	Mesh *Cut(Mesh *mesh);
	Mesh *AddPulseNoise(Mesh *mesh,float AlphaNoise);
	Mesh *ChangeOrder(Mesh *mesh,int RotateNum);
	Mesh *Disturb(Mesh *mesh);
	Mesh *Move(Mesh *mesh);
	Mesh *Scale(Mesh *mesh);
	Mesh *DownSample(Mesh *mesh);
	Mesh *DownSample1(Mesh *mesh,float prop);
	Mesh *DownSample2(Mesh *mesh);
	Mesh *AddNoise(Mesh *mesh,float *noise);
};
class Recover{
public:
	Mesh *resort(Mesh *mesh1,Mesh *mesh2);//mesh1 is orignal,mesh2 is disturbed
	Mesh *resample(Mesh *mesh1,Mesh *mesh2);//mesh1 is orignal,mesh2 is downsampled
};

class Permute{
public:
	void PSEUDO(unsigned char *origwmdata,unsigned char *pmwmdata,int N1,int N2,int key);
	void DPSEUDO(unsigned char *pmwmdata,unsigned char *rewmdata,int N1,int N2,int key);
};


#endif





