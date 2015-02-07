    // Source file for the OFF file viewer



////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////

// Standard include files
#include "GraphicsGems.h"
//#include "Matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#define WatermarkLength 256
#define WatermarkColumn 16
#define MaxNeighb 250
static int key=88;
static int iMethod;
static float alpha;
// Windows include files 

#ifdef _WIN32
#include <windows.h>
#endif



// OpenGL include files 

#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"




////////////////////////////////////////////////////////////
// CLASSES
////////////////////////////////////////////////////////////
class Attack;
class Face;
class Matrix;
class Mesh;
class Recover;
class Tools;
class Vector3;
class Vector4;
class Anneal;
class Permute;
class Embed;
class Attract;
class Embed{
public:
Mesh *EmbedWatermark0(Mesh *mesh,unsigned char *Watermark,float alpha);//embedding method of literature [33]
Mesh *EmbedWatermark1(Mesh *mesh,unsigned char *Watermark,float alpha);//my embedding method 1
Mesh *EmbedWatermark2(Mesh *mesh,unsigned char *Watermark,float alpha);//my embedding method 2
}embed;
class Attract{
public:
	void AttractWatermark0(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark);
	void AttractWatermark1(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark);
}attract;
class Permute{
public:
	void PSEUDO(unsigned char *origwmdata,unsigned char *pmwmdata,int N1,int N2,int seed);
	void DPSEUDO(unsigned char *pmwmdata,unsigned char *rewmdata,int N1,int N2,int seed);
}permute;
class Anneal{
public:
#define TFACTR  0.5
#define MBIG 1000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define MMAX 2
float y,T;
long idum;
Anneal();
int bolziman(float de,float t);
float ran(long *idum);
int irbit(unsigned long *iseed);
float move(Mesh *mesh1,Mesh *mesh2,Mesh *tempmesh,int n,int s,float x[],float t[],float c[],float min,float *up,float *low);
Mesh *registration(Mesh *mesh1,Mesh *mesh2);
}anneal;


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
	float Magnitude(Vector3 vector);
	float Distance(Vector3 vector1,Vector3 vector2);
	Vector3 Direction(Vector3 vector1,Vector3 vector2);
	Vector3 Normalize(Vector3 vector1);
	float Max(Vector3 &v);
	Vector3 operator + (Vector3& vector1);	// 'CVector + CVector'
	Vector3 operator - (Vector3& vector1);	// 'CVector - CVector'
	Matrix3 operator * (Vector3& vector1);	// 'CVector * CVector'
	Vector3 operator * (double dValue);	// 'CVector * double'	
	//Vector3 operator *= (double dValue);	// 'CVector *= double'	
	Vector3& operator = (Vector3& vector1);	// 'CVector = CVector'	
	Vector3& operator += (Vector3& vector1);	// 'CVector += CVector'	
	float operator / (Vector3& vector1);	// 'dot metrixs of CVector and CVector'	
	Vector3 operator / (double dValue);// CVector / double	
	//Vector3 operator /= (double dValue);// CVector /= double
} vector;
class Tools{
public:
	double SubEnergy(Mesh *mesh1,Mesh *mesh2);
	Vector3 mean(Vector3 *vertex,long nverts);
	Vector3 CrossOver(Vector3 vector1,Vector3 vector2);
	Vector3 FaceNormal(Mesh *mesh,Face &face);
	Vector4 Mat4MulV4(Matrix4 &m4,Vector4 &v4);//multiply matrix4 with vector4
	Vector3 Mat3MulV3(Matrix3 &m3,Vector3 &v3);//multiply matrix3 with vector3
	Matrix3 Variance(Mesh *mesh,int vertex);//compute the second order moment of the vertex's neighbourhoods
	void DeleteVertex(Mesh *mesh,int iVertex);//delete a vertex of the model and compute its aftermath 
	void DeleteFace(Mesh *mesh,int iFace);//delete a face of the model
	void ReCompNeighb(Mesh *mesh);//recompute neighbourhoods of the model
	bool bIsNeighb(Mesh *mesh,int a,int b);//found out if vertex cut's neighbourhoods a and b is neighbourhood 
	void AppendFaces(Mesh *mesh,int cut,int (*a)[2],int countface);//append faces when faces including the vertex is cut
	void NeighbFaces(Mesh *mesh,int cut,int (*a)[2],int *countface);//compute neighbourhood face vertices
	bool IsOdd(int vertex,int (*a)[2],int countface);//find out if vertex is single in the cut vertex's neighbourhood face vertices
	float FindND(Mesh *mesh1,Mesh *mesh2,int ivertex,int *matchi);//find the nearest distance between mesh1's vertex ivertex and mesh2's vertices 
	void DeleteOverlapFaces(Mesh *mesh);//delete overlap faces in mesh
	void DeleteAbnormal(Mesh *mesh);//delete all the abnormal faces in mesh
}tools;
class Recover{
public:
	Mesh *resort(Mesh *mesh1,Mesh *mesh2);//mesh1 is orignal,mesh2 is disturbed
	Mesh *resample(Mesh *mesh1,Mesh *mesh2);//mesh1 is orignal,mesh2 is downsampled
}recover;
class Mesh {
public:  
	Mesh(void);
	Mesh(Mesh *mesh);
	~Mesh();
  Vector3 centriod;
  int nverts;
  Vector3 *verts;
  int nfaces;
  Face *faces;
  Vector3 *VertexNormal;
  int (*neighbourhood)[MaxNeighb];
  bool *VertexExcluded;
  float *weight;
  bool *VNDirection;
  float magnitude;
} ;
class Matrix{
public:
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
	Matrix3 Mat3Div(Matrix3 &M,float n);
	float Mat3lemda(Matrix3 &M);
	float ComMat3lemda(Matrix3 &M,Vector3 *v,int iteration);
}mat;
class Face {
public:
  Face(void);
  Face(Face &face);
  Face(Face *face);
  ~Face();
  int nverts;
  int *verts;
  Vector3 normal;
}face;
class Attack{
public:
	Mesh *ChangeView(Mesh *mesh,Mesh *tempmesh,float a=0,float b=0,float r=0,float dx=0,float dy=0,float dz=0);
	Vector3 Rotate(Vector3 v,float a,float b,float r,float dx,float dy,float dz);
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
}attack;

Anneal::Anneal():T(0),idum(-1){}
int Anneal::bolziman(float de,float t)
 {
      //float ran(long *idum);
      static long gljdum=-1;
      return de<0.0 ||ran(&gljdum)<exp(-de/(t));
  }
float Anneal::ran(long *idum)
{static int inext,inextp;
 static long ma[56];
 static int iff=0;
 long mj,mk;
 int i,ii,k;
 if(*idum<0||iff==0) {
		      iff=1;
		      mj=MSEED-(*idum<0?-*idum:*idum);
		      mj%=MBIG;
		      ma[55]=mj;
		      mk=1;
		      for(i=1;i<54;i++){
					ii=(21*i)%55;
					ma[ii]=mk;
					mk=mj-mk;
					if(mk<MZ)mk+=MBIG;
					mj=ma[ii];
			}
			for(k=1;k<=4;k++)
			     for(i=1;i<=55;i++){
				   ma[i]-=ma[1+(i+30)%55];
				   if(ma[i]<MZ) ma[i]+=MBIG;
				   }
			 inext=0;
			 inextp=31;
			 *idum=1;
		     }
    if(++inext==56)inext=1;
    if(++inextp==56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if(mj<MZ) mj+=MBIG;
    ma[inext]=mj;
    return mj*FAC;
    }
 int Anneal::irbit(unsigned long *iseed)
 {
     unsigned long newbit;
      newbit=(*iseed &131072)>>17
	    ^(*iseed&16)>>4
	      ^(*iseed &2)>>1
	      ^(*iseed &1);
	      *iseed=(*iseed<<1)|newbit;
	      return newbit;
}
float Anneal::move(Mesh *mesh1,Mesh *mesh2,Mesh *tempmesh,int n,int s,float x[],\
				   float t[],float c[],float min,float *up,float *low)
{
 int i[6],j,k,flag;
 float sig(0),temp1(0),temp2(0),p1(0),p2(0),cost,mincost;

 float xtemp[6],xmin[6];
 for(j=0;j<6;j++)xtemp[j]=x[j];
   mincost=min;  flag=0;
for(i[5]=1;i[5]<=2;i[5]++)
    {  for(i[4]=1;i[4]<=2;i[4]++)
     {  for(i[3]=1;i[3]<=2;i[3]++)
       {for(i[2]=1;i[2]<=2;i[2]++)
	{for(i[1]=1;i[1]<=2;i[1]++)
	 {for(i[0]=1;i[0]<=2;i[0]++)
	  { for(j=0;j<n;j++)
	   { k=(j+s)%6;
	   if (i[j]==1) {
	   xtemp[k]=x[k]+(int)((up[k]-x[k])*ran(&idum)*MBIG*T)*FAC;
	   if(xtemp[k]>up[k]) xtemp[k]=up[k];
	   if(xtemp[k]<low[k])xtemp[k]=low[k];}
	    if(i[j]==2) {
	    xtemp[k]=x[k]-(int)((x[k]-low[k])*ran(&idum)*MBIG*T)*FAC;
            	   if(xtemp[k]>up[k]) xtemp[k]=up[k];
	   if(xtemp[k]<low[k])xtemp[k]=low[k];}

	    }

      float a=xtemp[0];float b=xtemp[1];float r=xtemp[2];
	  float dx=xtemp[3];float dy=xtemp[4];float dz=xtemp[5];
	  tempmesh=attack.ChangeView(mesh1,tempmesh,a,b,r,dx,dy,dz);
	  //fitness function
	 cost=tools.SubEnergy(mesh2,tempmesh);
      if ((mincost-cost)>FAC)
	   {flag=1;
	    mincost=cost;
	    for(j=0;j<6;j++)xmin[j]=xtemp[j];
	    }
	 } if(n==1) break;
	 }if(n==1||n==2) break;
	 }if(n==1||n==2||n==3) break;
	 }if(n==1||n==2||n==3||n==4) break;
	 }if(n==1||n==2||n==3||n==4||n==5)break;
	 }

      if(flag)
      for(j=0;j<6;j++)
	x[j]=xmin[j];

	return mincost;
     }

Mesh *Anneal::registration(Mesh *mesh1,Mesh *mesh2){
		printf("registration the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh2;
	int i;
	Vector3 dCentriod=mesh2->centriod-mesh1->centriod;
	for(i=0;i<mesh1->nverts;i++){
		mesh2->verts[i]=mesh2->verts[i]-dCentriod;
	}
	mesh2->centriod=mesh1->centriod;
	float num;
	printf("please input the anneal number:");
	 scanf("%f",&num);printf("\n");
	 printf("registrating using anneal....\n");
	Mesh *tempmesh=new Mesh(mesh1);
	float x[6],t[6];
	float a,b,r,dx,dy,dz;
  printf("please input the approximate rotate parameters(a,b,r is degree and dx,dy,dz is percent,relative error less than 20%s):a,b,r,dx,dy,dz\n","%");
  scanf("%f,%f,%f,%f,%f,%f",&a,&b,&r,&dx,&dy,&dz);printf("\n");
  float mag=mesh2->magnitude;
  float mx[6]={a*DTOR,b*DTOR,r*DTOR,0,0,0};//initial condition
    float c[6],mincost=MBIG,cost,dc;
    float up[6]={mx[0]*1.2,mx[1]*1.2,mx[2]*1.2,mx[3]*1.2,mx[4]*1.2,mx[5]*1.2};//up bound
	float low[6]={mx[0]*0.8,mx[1]*0.8,mx[2]*0.8,mx[3]*0.8,mx[4]*0.8,mx[5]*0.8};//low bound
    int m=0,j,n,s,k;
	for(i=0;i<6;i++){
		t[i]=0.005;c[i]=25;
	}
		for(k=0;k<6;k++){
			mx[k]=low[k]+(int)((up[k]-low[k])*anneal.ran(&anneal.idum)/FAC)*FAC;
			    if (mx[k]>up[k]) mx[k]=up[k];
			    if(mx[k]<low[k])mx[k]=low[k];
		}
		anneal.T=1.0;mincost=MBIG;
		for(i=1;i<=num;i++){
			for(j=1;j<=200;j++){
				n=(int)(anneal.ran(&anneal.idum)*4)%4+1;
				s=(1+(int)(anneal.ran(&anneal.idum)*6))%6;
					for(k=0;k<6;k++)x[k]=mx[k];
						cost=anneal.move(mesh1,mesh2,tempmesh,n,s,x,t,c,mincost,up,low);

					dc=cost-mincost;
					if(fabs(dc)>=FAC && anneal.bolziman(dc,anneal.T)){// && ??????????????
					for(k=0;k<6;k++)
						mx[k]=x[k];
					mincost=cost;
					}
			}
			printf("\nnumber:%d,T=%f,mincost=%f",i,anneal.T,mincost);
			printf("\n %f,%f,%f,%f,%f,%f",mx[0]*RTOD,mx[1]*RTOD,mx[2]*RTOD,dx,dy,dz);
			anneal.T*=TFACTR;
		}/* end of i */
		double **m1 = mat.matrix(0,4,0,4); //用matrix函数，分配一个4*4的矩阵,通过m[i][j]赋值，对矩阵初始化 
		//Matrix4 RotMat(PI/6,PI/3,PI/4,0.1,0.2,0.3);//generate true Rotate Matrix according to a,b,r,dx,dy,dz		
		Matrix4 RotMat(mx[0],mx[1],mx[2],mx[3],mx[4],mx[5]);//generate computed Rotate Matrix according to a,b,r,dx,dy,dz
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				m1[i][j]=RotMat.element[i][j];//transfer matrix to array 4*4
			}
		}
		printf("\n");
		mat.print(m1,0,4,0,4);
		mat.inverse(m1,4); //矩阵求逆，结果也在m中 
		//Matrix4 iRotMat(m1);???????????
		Matrix4 iRotMat;
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
			iRotMat.element[i][j]=m1[i][j];
			}
		}
		mat.print(m1,0,4,0,4);
		for(i=0;i<mesh2->nverts;i++){
			Vector4 inv4(mesh2->verts[i]);
			Vector3 outv3(tools.Mat4MulV4(iRotMat,inv4));
			tempmesh->verts[i]=outv3;//实现重定位
		}
		mat.free_matrix(m1,0,4,0,4); //释放矩阵。 
		return tempmesh;
}



    Vector3::Vector3(void):x(0),y(0),z(0){}
	Vector3::Vector3(float xx,float yy,float zz){x=xx;y=yy;z=zz;}
	Vector3& Vector3::operator =(Vector3& vector1){
		 x=vector1.x;y=vector1.y;z=vector1.z;
		 return *this;
	}
	Vector3::Vector3(Vector4 &p){x=p.x;y=p.y;z=p.z;}
	Vector3 Vector3::operator + (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x+vector1.x;
		vector.y=y+vector1.y;
		vector.z=z+vector1.z;
		return	vector;
	}
	Vector3 Vector3::operator - (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x-vector1.x;
		vector.y=y-vector1.y;
		vector.z=z-vector1.z;
		return	vector;
	}
	float Vector3::operator / (Vector3& vector1){
		Vector3	vector = *this;
		vector.x=x*vector1.x;
		vector.y=y*vector1.y;
		vector.z=z*vector1.z;
		float sum=vector.x+vector.y+vector.z;
		return	sum;
	}
	Vector3 Vector3::operator / (double dValue){
		Vector3 vector = *this;
		vector.x=x/dValue;
		vector.y=y/dValue;
		vector.z=z/dValue;
		return vector;
	}
/*
		Vector3& Vector3::operator /= (double dValue){
			x/=dValue;y/=dValue;z/=dValue;
			return *this;
		}*/
	
	Matrix3 Vector3::operator * (Vector3& vector1){
		Matrix3 matrix3;
		matrix3.element[0][0]=x*vector1.x;
 		matrix3.element[0][1]=x*vector1.y;
 		matrix3.element[0][2]=x*vector1.z;
 		matrix3.element[1][0]=y*vector1.x;
 		matrix3.element[1][1]=y*vector1.y;
 		matrix3.element[1][2]=y*vector1.z;
 		matrix3.element[2][0]=z*vector1.x;
 		matrix3.element[2][1]=z*vector1.y;
 		matrix3.element[2][2]=z*vector1.z;
		return matrix3;
	}
	Vector3 Vector3::operator * (double dValue){
		Vector3	vector = *this;
		vector.x=x*dValue;
		vector.y=y*dValue;
		vector.z=z*dValue;
		return	vector;
	}
/*
		Vector3& Vector3::operator *= (double dValue){
			x*=dValue;
			y*=dValue;
			z*=dValue;
			return	*this;
		}*/
	
	Vector3& Vector3::operator += (Vector3& vector1){
		x+=vector1.x;
		y+=vector1.y;
		z+=vector1.z;
		return	*this;
	}
	float Vector3::Magnitude(Vector3 vector){
		float result=sqrt((vector.x)*(vector.x)+(vector.y)*(vector.y)+(vector.z)*(vector.z));
		return (result<1e-8)?1e-8:result;
	}
	float Vector3::Distance(Vector3 vector1,Vector3 vector2){return Magnitude(vector1-vector2);}
	Vector3 Vector3::Direction(Vector3 vector1,Vector3 vector2){return (vector1-vector2)/Distance(vector1,vector2);}
	Vector3 Vector3::Normalize(Vector3 vector1){return vector1/vector.Magnitude(vector1);}
	float Vector3::Max(Vector3 &vector){
		return max(max(vector.x,vector.y),vector.z);
	}
Vector4::Vector4(){x=y=z=0;w=1;}
Vector4::Vector4(Vector3 &p){x=p.x;y=p.y;z=p.z;w=1;}

Face::Face(void) : nverts(0), verts(0) {normal.x=normal.y=normal.z=0;}
Face::Face(Face &face){
	 nverts=face.nverts;
	 normal=face.normal;
	 verts=new int[face.nverts];assert(verts);
	 for (int i=0;i<face.nverts;i++){verts[i]=face.verts[i];}
}
Face::Face(Face *face){
	 nverts=face->nverts;
	 normal=face->normal;
	 verts=new int[face->nverts];assert(verts);
	 for (int i=0;i<face->nverts;i++){verts[i]=face->verts[i];}
}
Face::~Face(){
	delete []verts;
}
Mesh::Mesh(void) : nverts(0), verts(0), nfaces(0), faces(0) ,magnitude(0),\
		VertexNormal(0),neighbourhood(0),VertexExcluded(0),weight(0),VNDirection(0){};
Mesh::Mesh(Mesh *mesh):nverts(mesh->nverts),nfaces(mesh->nfaces),\
			centriod(mesh->centriod),magnitude(mesh->magnitude){
		verts=new Vector3[mesh->nverts];assert(verts);
		VertexNormal=new Vector3[mesh->nverts];assert(VertexNormal);
		VertexExcluded=new bool[mesh->nverts];assert(VertexExcluded);
		weight=new float[mesh->nverts];assert(weight);
		neighbourhood=new int[mesh->nverts][MaxNeighb];assert(neighbourhood);
		VNDirection=new bool[mesh->nverts];assert(VNDirection);
		for (int i=0;i<nverts;i++) {
			verts[i]=mesh->verts[i];
			VertexNormal[i]=mesh->VertexNormal[i];
			VertexExcluded[i]=mesh->VertexExcluded[i];
			weight[i]=mesh->weight[i];
			VNDirection[i]=mesh->VNDirection[i];
			for(int j=0;j<MaxNeighb;j++){neighbourhood[i][j]=mesh->neighbourhood[i][j];}
		}
		faces=new Face[nfaces];assert(faces);
		for (i=0;i<nfaces;i++){
			//faces[i]=mesh->faces[i];
			 faces[i].nverts=mesh->faces[i].nverts;
			 faces[i].normal=mesh->faces[i].normal;
			 faces[i].verts=new int[mesh->faces[i].nverts];
			 for (int j=0;j<faces[i].nverts;j++){faces[i].verts[j]=mesh->faces[i].verts[j];}
		}
	}
Mesh::~Mesh(){
  for (int i=0;i<nfaces;i++){
	delete []faces[i].verts;
  }
  for (i=0;i<MaxNeighb;i++){
	delete []&neighbourhood[i];
  }
  delete []verts;
  delete []VertexExcluded;
  delete []weight;
  delete []VNDirection;
  delete []neighbourhood;
}
	double Tools::SubEnergy(Mesh *mesh2,Mesh *mesh3){
		double energy=0;
		for (int i=0;i<mesh2->nverts;i++){
			energy+=vector.Distance(mesh2->verts[i],mesh3->verts[i])*vector.Distance(mesh2->verts[i],mesh3->verts[i]);
		}
		return energy;
	}
	Vector3 Tools::mean(Vector3 *vertex,long nverts){
		Vector3 sum;
		for (int i=0;i<nverts;i++){
			sum+=vertex[i];
		}
		return sum/nverts;
	};
	Vector3 Tools::CrossOver(Vector3 vector1,Vector3 vector2){
		Vector3 result;
		result.x=vector1.y*vector2.z-vector2.y*vector1.z;
		result.y=vector1.x*vector2.z-vector2.x*vector1.z;
		result.z=vector1.x*vector2.y-vector2.x*vector1.y;
		return result;
	}
	Vector3 Tools::FaceNormal(Mesh *mesh,Face &face){
		Vector3 *v1 = &(mesh->verts[face.verts[face.nverts-1]]);
		for (int i = 0; i < face.nverts; i++) {
			Vector3 *v2 = &(mesh->verts[face.verts[i]]);
			face.normal.x += (v1->y - v2->y) * (v1->z + v2->z);
			face.normal.y += (v1->z - v2->z) * (v1->x + v2->x);
			face.normal.z += (v1->x - v2->x) * (v1->y + v2->y);
			v1 = v2;
		}
		return face.normal;
	}
void Tools::DeleteVertex(Mesh *mesh,int iVertex){
		int i,j,k,k1,temp;bool endsearch=FALSE;
/*if the comments below is programmed rightly ,the neighbourhoods computing will be much more time saveing*/

		bool bFound=FALSE;int iNeighb=-1;
		for (i=0;i<mesh->nverts;i++){//change neighbourhoods concerning iVertex
			iNeighb=-1;bFound=FALSE;
			while (mesh->neighbourhood[i][++iNeighb]!=-1 && !bFound) {
				if(mesh->neighbourhood[i][iNeighb]==iVertex){
					bFound=TRUE;
					for(j=iNeighb;j<MaxNeighb-1;j++){
						temp=mesh->neighbourhood[i][j+1];
						mesh->neighbourhood[i][j]=(temp>iVertex)?(temp-1):temp;
					}
					for(j=0;j<iNeighb;j++){
						if(mesh->neighbourhood[i][j]>iVertex)
							mesh->neighbourhood[i][j]--;
					}
				}
			}
			if(!bFound){
				iNeighb=-1;
				while (mesh->neighbourhood[i][++iNeighb]!=-1) {
					if(mesh->neighbourhood[i][iNeighb]>iVertex)
						mesh->neighbourhood[i][iNeighb]--;
				}
			}
		}
		for(i=iVertex;i<mesh->nverts-1;i++)
			for(j=0;j<MaxNeighb;j++)
				mesh->neighbourhood[i][j]=mesh->neighbourhood[i+1][j];


		
		int nfaces=mesh->nfaces;bool endsearch1=FALSE;
		for (i=0;i<nfaces;i++){//change faces concerning iVertex
			if(endsearch1) i--;
			for (j=0;j<mesh->faces[i].nverts;j++){
				endsearch1=FALSE;
				if(mesh->faces[i].verts[j]==iVertex && !endsearch1){//delete the faces connected to iVertex
					mesh->nfaces--;
					if(i>=mesh->nfaces) endsearch=TRUE;
					if(!endsearch) {
						for(k=i;k<nfaces-1;k++) {
/*
								Face &face=mesh->faces[k+1];
							mesh->faces[k]=Face(face);
							delete []mesh->faces[k].verts;
							mesh->faces[k]=new Face(&mesh->faces[k+1]);
*/

							 mesh->faces[k].nverts=mesh->faces[k+1].nverts;
							 mesh->faces[k].normal=mesh->faces[k+1].normal;
							 for (k1=0;k1<mesh->faces[k+1].nverts;k1++){mesh->faces[k].verts[k1]=mesh->faces[k+1].verts[k1];}


						}
						j=mesh->faces[i].nverts;//end searching j
						i--;//researching iVertex
						if(i==-1) {i=0;endsearch1=TRUE;}
					}
				}
			}
			if(endsearch) i=nfaces;//end searching i
		}
		for (i=0;i<mesh->nfaces;i++){//change faces indexes concerning iVertex
			for(j=0;j<mesh->faces[i].nverts;j++){
				if(mesh->faces[i].verts[j]>iVertex)
					mesh->faces[i].verts[j]--;
			}
		}
								
		for (i=iVertex;i<mesh->nverts-1;i++){//delete the vertex and its properities
			mesh->weight[i]=mesh->weight[i+1];
			mesh->VertexNormal[i]=mesh->VertexNormal[i+1];
			mesh->VertexExcluded[i]=mesh->VertexExcluded[i+1];
			mesh->verts[i]=mesh->verts[i+1];
			//for(j=iNeighb;j<MaxNeighb-1;j++)
				//mesh->neighbourhood[i][j]=mesh->neighbourhood[i+1][j];
		}
		mesh->nverts--;


}
	void Tools::ReCompNeighb(Mesh *mesh){//recompute neighbourhoods
		int i,j,k,k1,iNeighb;bool endsearch=FALSE;
		for (i=0;i<mesh->nverts;i++)
			for(j=0;j<MaxNeighb;j++)
				mesh->neighbourhood[i][j]=-1;//initiate the neighbourhoods
		for(i=0;i<mesh->nverts;i++){
			bool m_bFound=FALSE;
			for (j=0;j<mesh->nfaces;j++){
				for(k=0;k<mesh->faces[j].nverts;k++){
					  if (i==mesh->faces[j].verts[k]){//match a vertex and its face
						  k=mesh->faces[j].nverts;//end searching
						  //find the neighbourhoods
						  m_bFound=FALSE;
						  for(k1=0;k1<mesh->faces[j].nverts;k1++){//neighbourhood face of vertex i is found
							iNeighb=-1;while(mesh->neighbourhood[i][++iNeighb]!=-1){
								if(mesh->neighbourhood[i][iNeighb]==mesh->faces[j].verts[k1])
									m_bFound=TRUE;//if mesh->faces[j].verts[kk] is found in vertex i's neighbourhood,do nothing
							}
							if(!m_bFound && i!=mesh->faces[j].verts[k1])
							{
								iNeighb=-1;while(mesh->neighbourhood[i][++iNeighb]!=-1);//trace the iNeighb
								mesh->neighbourhood[i][iNeighb]=mesh->faces[j].verts[k1];//append vertex i's neighbourhood
								iNeighb=-1;while(mesh->neighbourhood[mesh->faces[j].verts[k1]][++iNeighb]!=-1);//trace the tempneighb
								mesh->neighbourhood[mesh->faces[j].verts[k1]][iNeighb]=i;//append vertex mesh->faces[j].verts[kk]'s neighbourhood
							}
						  }
					  }
				  }
			  }
		  }	
	}
	void Tools::DeleteFace(Mesh *mesh,int iFace){
		int i;
		for (i=iFace;i<mesh->nfaces;i++){
			if(i==(mesh->nfaces-1)){//face to be cut is the last face
				mesh->nfaces--;
				return;
			}
			//delete face
			mesh->faces[i].normal =mesh->faces[i+1].normal;
			mesh->faces[i].nverts=mesh->faces[i+1].nverts;
			mesh->faces[i].verts[0]=mesh->faces[i+1].verts[0];
			mesh->faces[i].verts[1]=mesh->faces[i+1].verts[1];
			mesh->faces[i].verts[2]=mesh->faces[i+1].verts[2];
		}
	}
	void Tools::NeighbFaces(Mesh *mesh,int cut,int (*a)[2],int *countface){
		int t1,t2,t3;

		int icountface=0;
		for (t1=0;t1<mesh->nfaces;t1++){
			for(t2=0,t3=0;t2<mesh->faces[t1].nverts;t2++){//find all the neighbourhood faces vertices
				if (cut==mesh->faces[t1].verts[t2]){
						while(cut==mesh->faces[t1].verts[t3++]); 
						a[icountface][0]=mesh->faces[t1].verts[--t3];
						while(cut==mesh->faces[t1].verts[++t3]);
						a[icountface++][1]=mesh->faces[t1].verts[t3];
				}
			}
		}
	}
	void Tools::AppendFaces(Mesh *mesh,int cut,int (*a)[2],int countface){
	//	if(countface<1) return;
		int i,j;int first,mid,last,icountface=0;
		bool *NeighfaceUsed=new bool[countface];
		for(i=0;i<countface;i++){NeighfaceUsed[i]=FALSE;}
		first=a[0][0];mid=a[0][1];icountface=0;i=0;bool endflag=FALSE;
		NeighfaceUsed[0]=TRUE;
		do {
			bool tempflag=IsOdd(mid,a,countface);
			if(tempflag){
				tempflag=FALSE;
				for(int t1=0;t1<countface;t1++){
					for(int t2=0;t2<2;t2++){
						if(IsOdd(a[t1][t2],a,countface)&&!NeighfaceUsed[t1]&&(a[t1][t2]!=mid)){
							NeighfaceUsed[t1]=TRUE;
							mid=a[t1][t2];
							tempflag=TRUE;
						}
					}
				}
			if (!tempflag) return;
			}//process the odd condition

			if(i){
				if(mid==a[i][0]) {
					NeighfaceUsed[i]=TRUE;
					last=a[i][1];
					if(first!=last){
						mesh->faces[mesh->nfaces].verts[0]=first;
						mesh->faces[mesh->nfaces].verts[1]=mid;
						mesh->faces[mesh->nfaces].verts[2]=last;
						mesh->nfaces++;
					}
					mid=last;
				}
				else if(mid==a[i][1]){
					NeighfaceUsed[i]=TRUE;
					last=a[i][0];
					if(first!=last){
						mesh->faces[mesh->nfaces].verts[0]=first;
						mesh->faces[mesh->nfaces].verts[1]=mid;
						mesh->faces[mesh->nfaces].verts[2]=last;
						mesh->nfaces++;
					}
					mid=last;
				}
			}
			for(j=0,endflag=TRUE;j<countface;j++){endflag&=NeighfaceUsed[j];}
			i=++i%countface;
			if(first==last) endflag=TRUE;
		} while(!endflag);
	}
	void Tools::DeleteOverlapFaces(Mesh *mesh){
		int t1,t2,t3,t4;
				int equalnum;
				for(t1=0;t1<mesh->nfaces-1;t1++){
					for(t2=t1+1;t2<mesh->nfaces;t2++){
						equalnum=0;
						for(t3=0;t3<3;t3++){
							for(t4=0;t4<3;t4++){
								if(mesh->faces[t1].verts[t3]==mesh->faces[t2].verts[t4]){
									t4=3;
									equalnum++;
								}								
							}
						}
						if(equalnum==3){
							tools.DeleteFace(mesh,t2);
						}
					}
				}
	}
	void Tools::DeleteAbnormal(Mesh *mesh){
		printf("deleting abnormal faces and vertices....\n");
		int i,t1,t2;
		int countvert1=0,countface;int deleteface;bool abnormalend=TRUE;
		do{
			for (i=0;i<mesh->nverts;i++){
				countface=0;
				for (t1=0;t1<mesh->nfaces;t1++){
						for(t2=0;t2<3;t2++){
							if (i==mesh->faces[t1].verts[t2]){
								countface++;
								deleteface=t1;
							}
						}
					}
				if(countface==1) {
					tools.DeleteFace(mesh,deleteface);
					tools.DeleteVertex(mesh,i-countvert1++);
					abnormalend=FALSE;
				}
			}
		}while(!abnormalend);
	}
	bool Tools::IsOdd(int vertex,int (*a)[2],int countface){
		bool isodd=FALSE;
		for(int i=0;i<countface;i++){
			if((vertex==a[i][0])||(vertex==a[i][1]))
				isodd=!isodd;
		}
		return isodd;
	}
	bool Tools::bIsNeighb(Mesh *mesh,int a,int b){
		int iNeighb1=-1;
		while(mesh->neighbourhood[a][++iNeighb1]!=-1){
			if((b!=a)&&(mesh->neighbourhood[a][iNeighb1]==b))
				return TRUE;
			}
		return FALSE;
	}
	Vector4 Tools::Mat4MulV4(Matrix4 &m4,Vector4 &v4){
		Vector4 rv4;
		rv4.x=m4.element[0][0]*v4.x+m4.element[0][1]*v4.y+m4.element[0][2]*v4.z+m4.element[0][3]*v4.w;
		rv4.y=m4.element[1][0]*v4.x+m4.element[1][1]*v4.y+m4.element[1][2]*v4.z+m4.element[1][3]*v4.w;
		rv4.z=m4.element[2][0]*v4.x+m4.element[2][1]*v4.y+m4.element[2][2]*v4.z+m4.element[2][3]*v4.w;
		rv4.w=m4.element[3][0]*v4.x+m4.element[3][1]*v4.y+m4.element[3][2]*v4.z+m4.element[3][3]*v4.w;
		return rv4;
	}
	Vector3 Tools::Mat3MulV3(Matrix3 &m3,Vector3 &v3){
		Vector3 rv3;
		rv3.x=m3.element[0][0]*v3.x+m3.element[0][1]*v3.y+m3.element[0][2]*v3.z;
		rv3.y=m3.element[1][0]*v3.x+m3.element[1][1]*v3.y+m3.element[1][2]*v3.z;
		rv3.z=m3.element[2][0]*v3.x+m3.element[2][1]*v3.y+m3.element[2][2]*v3.z;
		return rv3;
	}
	Matrix3 Tools::Variance(Mesh *mesh,int vertex){
	  int iNeighb=0;Vector3 average,temp;
	  while(mesh->neighbourhood[vertex][iNeighb]!=-1){
		  temp+=mesh->verts[mesh->neighbourhood[vertex][iNeighb++]];
	  }
	  temp=temp+mesh->verts[vertex];
	  average=temp/(iNeighb+1);
	  Vector3 vector1,vector2;
	  Matrix3 matrix3;
	  iNeighb=0;temp.x=temp.y=temp.z=0;
	  while(mesh->neighbourhood[vertex][iNeighb]!=-1){
	    vector1=mesh->verts[mesh->neighbourhood[vertex][iNeighb]]-average;
		vector2=vector1;
		matrix3=mat.AddMatrix3(matrix3,(vector1*vector2));
		iNeighb++;
	  }
	  vector1=mesh->verts[vertex]-average;vector2=vector1;
	  matrix3=mat.AddMatrix3(matrix3,(vector1*vector2));
	  matrix3=mat.Mat3Div(matrix3,iNeighb+1);
	  return matrix3;
	}
	float Tools::FindND(Mesh *mesh1,Mesh *mesh2,int ivertex,int *matchi){
		long int i;float min=MBIG;float distance;
		for(i=0;i<mesh2->nverts;i++){
			distance=vector.Distance(mesh1->verts[ivertex],mesh2->verts[i]);
			if(min>distance){
				min=distance;
				*matchi=i;
			}
		}
		return min;
	}

	Vector3 Attack::Rotate(Vector3 v,float a,float b,float r,float dx,float dy,float dz){
		Vector4 inv4(v);Matrix4 RotMat(a,b,r,dx,dy,dz);
		Vector4 outv4=tools.Mat4MulV4(RotMat,inv4);
		Vector3 outv3(outv4);
		return outv3;
	}
	Mesh *Attack::Cut(Mesh *mesh){
		printf("cut the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		Mesh *tempmesh=new Mesh(mesh);
		float CutRate=0;//insection rate
	    printf("please input the cut rate:");scanf("%f",&CutRate);
		char axis;
again1:	printf("\nplease input the axis the model is to be cut along(x,y,z) or according to order(o):");scanf("%s",&axis);
		printf("\nCutting the model,please wait...");
		int nverts=tempmesh->nverts;Vector3 *verts=tempmesh->verts;int VertsCut=nverts;
		int i;
 		int n=int(nverts*CutRate);//compute the vertices number to be cut
		float bbox[2][3] = { { 1.0E30F, 1.0E30F, 1.0E30F }, { -1.0E30F, -1.0E30F, -1.0E30F } };
		for (i = 0; i < nverts; i++) {//compute the bounding box
			Vector3& vert = verts[i];
			if (vert.x < bbox[0][0]) bbox[0][0] = vert.x;
			else if (vert.x > bbox[1][0]) bbox[1][0] = vert.x;
			if (vert.y < bbox[0][1]) bbox[0][1] = vert.y;
			else if (vert.y > bbox[1][1]) bbox[1][1] = vert.y;
			if (vert.z < bbox[0][2]) bbox[0][2] = vert.z;
			else if (vert.z > bbox[1][2]) bbox[1][2] = vert.z;
		}

		float ratio=0;float step;long unsigned int nstep=0;
		float bound;
		switch(axis) {//find the bound and the actual vertices cut number:VertsCut
		case 'x':
			step=(bbox[1][0]-bbox[0][0])/(3*nverts);//resolution
			while (VertsCut>n) {
				VertsCut=0;
				bound=bbox[0][0]+nstep++*step;
				for (i=0;i<nverts;i++){
					if(verts[i].x>bound) VertsCut++;
				}
			}
			break;
		case 'y':			
			step=(bbox[1][1]-bbox[0][1])/(3*nverts);//resolution
			while (VertsCut>n) {
				VertsCut=0;
				bound=bbox[0][1]+nstep++*step;
				for (i=0;i<nverts;i++){
					if(verts[i].y>bound) VertsCut++;
				}
			}
			break;
		case 'z':
			step=(bbox[1][2]-bbox[0][2])/(3*nverts);//resolution
			while (VertsCut>n) {
				VertsCut=0;
				bound=bbox[0][2]+nstep++*step;
				for (i=0;i<nverts;i++){
					if(verts[i].z>bound) VertsCut++;
				}
			}
			break;
		case 'o':
			VertsCut=n;
			break;
		default:
			printf("axis input error\n");goto again1;
		}

		int *indexcut=new int[VertsCut];//restore the indexes cut
		switch(axis) {
		case 'x':
			bound=bbox[0][0]+(nstep-1)*step;
			VertsCut=0;
			for (i=0;i<nverts;i++){
				if(verts[i].x>bound) indexcut[VertsCut++]=i;
			}
			break;
		case 'y':
			bound=bbox[0][1]+(nstep-1)*step;
			VertsCut=0;
			for (i=0;i<nverts;i++){
				if(verts[i].y>bound) indexcut[VertsCut++]=i;
			}
			break;
		case 'z':
			bound=bbox[0][2]+(nstep-1)*step;
			VertsCut=0;
			for (i=0;i<nverts;i++){
				if(verts[i].z>bound) indexcut[VertsCut++]=i;
			}
			break;
		case 'o':
			for(i=0;i<VertsCut;i++){
				indexcut[i]=nverts-VertsCut+i;
			}
			break;
		default:printf("this is impossible!");
		}


		
		int displacement=0;//vertices cut indexes displacement
			for (i=0;i<VertsCut;i++){
				tools.DeleteVertex(tempmesh,indexcut[i]-displacement++);
			}
		printf("end\n");
		printf("%d vertices have been cut from model,%2.2f%s\n",VertsCut,(float)VertsCut/nverts*100,"%");//inform the cut number
		delete []indexcut;
		return tempmesh;
	}


	Mesh *Attack::AddPulseNoise(Mesh *mesh,float AlphaNoise){//add pulse noise
		printf("add pulse noise? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("Adding noise....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
	    srand((unsigned)time(NULL));
		Vector3 temp;float Noise;
		for (int i=0;i<mesh->nverts;i++){
			Noise=rand();
			if (Noise>=16383.5) Noise=AlphaNoise;
			else Noise =-AlphaNoise;
			temp=mesh->verts[i]-mesh->centriod;
			temp=temp*Noise;
			tempmesh->verts[i]=mesh->verts[i]+temp;
		}
		return tempmesh;
	};
	Mesh *Attack::AddNoise(Mesh *mesh,float *noise){
		printf("add gaussian noise? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("adding gaussian noise.....\n");
		float noisep;
		printf("please input the noise magnitude:");
	    scanf("%f",&noisep);printf("\n");
		Mesh *tempmesh=new Mesh(mesh);
		int i;float max=0;Vector3 temp;
		for(i=0;i<mesh->nverts;i++){
			temp=mesh->verts[i]-mesh->centriod;
			float t=vector.Magnitude(temp);
			if(max<t)
				max=t;
		}
		for(i=0;i<mesh->nverts;i++){
			noise[i]*=max*noisep/100;
		    //printf("noise[%3d]=%f\n",i,noise[i]);
			temp=mesh->verts[i]-mesh->centriod;
			temp=vector.Normalize(temp);
			temp=temp*noise[i];
			tempmesh->verts[i]=mesh->verts[i]+temp;
		}
		return tempmesh;
	}
	Mesh *Attack::ChangeOrder(Mesh *mesh,int RotateNum){
		printf("change vertices order? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("Rotating....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		RotateNum=RotateNum%mesh->nverts;
		for(int i=0;i<mesh->nverts;i++){
			tempmesh->verts[i]=mesh->verts[(i+RotateNum+mesh->nverts)%mesh->nverts];
		}
		for(i=0;i<mesh->nfaces;i++){
			for(int j=0;j<mesh->faces[i].nverts;j++){
				tempmesh->faces[i].verts[j]=(mesh->faces[i].verts[j]-RotateNum+mesh->nverts)%mesh->nverts;
			}
		}
		return tempmesh;
	}
	Mesh *Attack::ChangeView(Mesh *mesh,Mesh *tempmesh,float a,float b,float r,float dx,float dy,float dz){
		//printf("%f %f %f %f %f %f \n",a*RTOD,b*RTOD,r*RTOD,dx,dy,dz);
		for(int i=0;i<mesh->nverts;i++){
			tempmesh->verts[i]=attack.Rotate(mesh->verts[i],a,b,r,dx,dy,dz);
		}
		return tempmesh;
	}
	Mesh *Attack::Disturb(Mesh *mesh){
		printf("disturb the order of the model randomly? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("disturbing....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		srand((unsigned)time(NULL));
		for(int i=0;i<tempmesh->nverts;i++){
			int RotateNum=rand()%tempmesh->nverts;
			Vector3 tempvector;
			tempvector=tempmesh->verts[i];
			tempmesh->verts[i]=tempmesh->verts[(i+RotateNum+tempmesh->nverts)%tempmesh->nverts];
			tempmesh->verts[(i+RotateNum+tempmesh->nverts)%tempmesh->nverts]=tempvector;
			for(int j=0;j<tempmesh->nfaces;j++){
				for(int k=0;k<tempmesh->faces[j].nverts;k++){
					if (tempmesh->faces[j].verts[k]==i){
						tempmesh->faces[j].verts[k]=(tempmesh->faces[j].verts[k]+RotateNum+tempmesh->nverts)%tempmesh->nverts;
					}
					else if (tempmesh->faces[j].verts[k]==((i+RotateNum)%tempmesh->nverts)){
						tempmesh->faces[j].verts[k]=(tempmesh->faces[j].verts[k]-RotateNum+tempmesh->nverts)%tempmesh->nverts;
					}
				}
			}
		}
		return tempmesh;
	}
	Mesh *Attack::Move(Mesh *mesh){
		printf("move the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		Vector3 dCentriod(1,2,3);
		printf("moving whole model....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		for(int i=0;i<mesh->nverts;i++){
			tempmesh->verts[i]=mesh->verts[i]+dCentriod;
		}
		tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);
		return tempmesh;
	}
	Mesh *Attack::Scale(Mesh *mesh){
		printf("scale the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("scaling whole model....please wait\n");
		  float scale;
		  printf("please input the scale magnitude:");
		  scanf("%f",&scale);printf("\n");
		Mesh *tempmesh=new Mesh(mesh);
		for(int i=0;i<mesh->nverts;i++){
			Vector3 temp=mesh->verts[i]-mesh->centriod;
			temp=temp*scale;
			tempmesh->verts[i]=mesh->centriod+temp;
		}
		tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);
		return tempmesh;
	}
	Mesh *Attack::DownSample1(Mesh *mesh,float prop){/*downsample mesh2*/
		printf("downsampling.....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		int i,j,cutvertex;
		int orignnum=tempmesh->nverts;
		srand((unsigned)time(NULL));
		for (i=0;i<orignnum*prop;i++){
			//tools.DeleteVertex(tempmesh,rand()%tempmesh->nverts);//cut the vertex
			cutvertex=rand()%tempmesh->nverts;
			if(cutvertex!=tempmesh->nverts-1){
				for (j=cutvertex;j<tempmesh->nverts-1;j++){
					tempmesh->verts[j]=tempmesh->verts[j+1];
				}
			}
			tempmesh->nverts--;
		}
		return tempmesh;
	}
Mesh *Attack::DownSample2(Mesh *mesh){/*downsample mesh2*/
		printf("downsample the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
	float rate;
	printf("please input the downsample rate:");
	scanf("%f",&rate);printf("\n");
	int i,j,k;bool m_bFound=FALSE;
	int iNeighb=0;//neighbourhood index
	for(i=0;i<mesh->nverts;i++)	mesh->VertexExcluded[i]=FALSE;//all vertices can be cut
	int countvert=0;//vertices cut number
	int countface=0;//faces cut number
	int num=(int)(mesh->nverts*rate);
	int *newface=new int[MaxNeighb];
	int nverts=mesh->nverts;
	Face *lastface=NULL;
	for (i=0;i<nverts;i++){
		if((!mesh->VertexExcluded[i-countvert])&&(countvert<num)){//if the vertex can be cut
			printf("%d\t",i);
			for(iNeighb=0;iNeighb<MaxNeighb;iNeighb++) {
				newface[iNeighb]=mesh->neighbourhood[i-countvert][iNeighb];//restore the new face vertices
			}
			tools.DeleteVertex(mesh,i-countvert);//delete a vertex
			//append the new face whose edge number is no less than 3
			mesh->nfaces++;
			//mesh->faces[mesh->nfaces-1]=new Face();//????????
			iNeighb=-1;
			while(newface[++iNeighb]!=-1) {
				if(newface[iNeighb]>i-countvert) newface[iNeighb]--;
				mesh->VertexExcluded[newface[iNeighb]]=TRUE;//vertices connected to the vertex cut should not be cut
				mesh->faces[mesh->nfaces-1].verts[iNeighb]=newface[iNeighb];//(newface[iNeighb]>i-countvert)?(newface[iNeighb]-1):newface[iNeighb];//newface[iNeighb];
			}
			tools.FaceNormal(mesh,mesh->faces[mesh->nfaces-1]);
			mesh->faces[mesh->nfaces-1].nverts=iNeighb;
			countvert++;
			//tools.ReCompNeighb(mesh);//recompute neighbourhoods
			Face *lastface=&mesh->faces[mesh->nfaces-1];
			for(j=0;j<lastface->nverts;j++){
				for(k=0;k<lastface->nverts;k++){
					iNeighb=-1;m_bFound=FALSE;
					while(mesh->neighbourhood[lastface->verts[j]][++iNeighb]!=-1 && !m_bFound){
						if(mesh->neighbourhood[lastface->verts[j]][iNeighb]==lastface->verts[k])
							m_bFound=TRUE;//if mesh->faces[j].verts[kk] is found in vertex i's neighbourhood,do nothing
					}
					if(!m_bFound && j!=k){
						iNeighb=-1;while(mesh->neighbourhood[lastface->verts[j]][++iNeighb]!=-1);//trace the last iNeighb
						mesh->neighbourhood[lastface->verts[j]][iNeighb]=lastface->verts[k];
						iNeighb=-1;while(mesh->neighbourhood[lastface->verts[k]][++iNeighb]!=-1);//trace the last iNeighb
						mesh->neighbourhood[lastface->verts[k]][iNeighb]=lastface->verts[j];
					}
				}
			}
		}
	}
	printf("\n%d vertices downsampled,%f%s\n",countvert,(float)countvert/nverts*100,"%");
	delete []newface;
	return mesh;
}
	Mesh *Attack::DownSample(Mesh *mesh){/*downsample mesh2*/
		printf("downsampling.....please wait\n");
		float rate;
		printf("please input the downsample rate:");
		scanf("%f",&rate);printf("\n");
		Mesh *tempmesh=new Mesh(mesh);
		int num=(int)(tempmesh->nverts*rate);
		int i,j,k1,t1,t2;
		int iNeighb=0;//neighbourhood index
		for(i=0;i<tempmesh->nverts;i++)	tempmesh->VertexExcluded[i]=FALSE;//all vertices can be cut
		int countvert=0;//vertices cut number
		int countface=0;//faces cut number
		for (i=0;i<tempmesh->nverts;i++){
			if((!tempmesh->VertexExcluded[i])&&(countvert<num)){//if the vertex can be cut
				printf("%d\t",i);
				countface=0;//neighbourhood face number of i
				for (t1=0;t1<tempmesh->nfaces;t1++){
					for(t2=0;t2<tempmesh->faces[t1].nverts;t2++){
						if (i==tempmesh->faces[t1].verts[t2])
							countface++;
					}
				}
				static int (*a)[2];//neighbourhood face vertices of i
				a=new int[countface][2];
				tools.NeighbFaces(tempmesh,i,a,&countface);//compute neighbourhood face vertices of i
				iNeighb=0;
				while(tempmesh->neighbourhood[i][iNeighb]!=-1) 
					tempmesh->VertexExcluded[tempmesh->neighbourhood[i][iNeighb++]]=TRUE;//neighbourhoods should not be cut
				for(j=0;j<tempmesh->nfaces;j++){//when the vertex is cut,faces containing it should be cut
					for(k1=0;k1<tempmesh->faces[j].nverts;k1++){
						if(i==tempmesh->faces[j].verts[k1]) {
							k1=tempmesh->faces[j].nverts;
							tools.DeleteFace(tempmesh,j--);
						}
					}
				}
				tools.AppendFaces(tempmesh,i-countvert,a,countface);//append faces
				delete []a;
				tools.DeleteVertex(tempmesh,i-countvert++);//cut the vertex
			}
		}
		tools.DeleteOverlapFaces(tempmesh);
		//tools.DeleteAbnormal(tempmesh);
		/*printf("\n");
		for(i=0;i<tempmesh->nverts;i++){
			  iNeighb=0;
			  printf("%s %d %s","Neighbourhoods of Vertex",i,"is:");
			  while(tempmesh->neighbourhood[i][iNeighb]!=-1) printf(" %d",tempmesh->neighbourhood[i][iNeighb++]);
			  printf("\n");	
		}*/
		printf("\n%d vertices cut\n",countvert);
		return tempmesh;
	}


Mesh *Recover::resort(Mesh *mesh1,Mesh *mesh2){//mesh1 is orignal,mesh2 is disturbed
		printf("resort the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh2;
		printf("resorting....please wait\n");
		Mesh *tempmesh=new Mesh(mesh2);
		int i,j,k1,k2;
		/*reshape mesh2 according to mesh1*/
		Vector3 dCentriod=tempmesh->centriod-mesh1->centriod;
		float SizeMesh1(0),SizeMesh2(0);
		for(i=0;i<tempmesh->nverts;i++){
			SizeMesh1+=vector.Distance(mesh1->verts[i],mesh1->centriod);
			SizeMesh2+=vector.Distance(tempmesh->verts[i],tempmesh->centriod);
			tempmesh->verts[i]=tempmesh->verts[i]-dCentriod;
		}
		tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);
		SizeMesh1/=mesh1->nverts;SizeMesh2/=tempmesh->nverts;
		for(i=0;i<tempmesh->nverts;i++){
			Vector3 temp=tempmesh->verts[i]-tempmesh->centriod;
			temp=temp*(float)SizeMesh1/SizeMesh2;
			tempmesh->verts[i]=tempmesh->centriod+temp;
		}

		/*resort vertices of mesh2*/
	    bool *VertexUsed=new bool[mesh1->nverts];assert(VertexUsed);//if the match vertex is found,never use it
		for(i=0;i<mesh1->nverts;i++){VertexUsed[i]=FALSE;}
		float alpha1=0.005;float alpha2=alpha1/10;bool bFound=FALSE;int changenum=0;
		for(i=0;i<tempmesh->nverts;i++){
			if(bFound&&changenum) {alpha1-=alpha2*changenum;changenum=0;}
			bFound=FALSE;//if the match vertex is not found,change the finding range
			for(j=0;j<mesh1->nverts;j++){
				if (!VertexUsed[j]){
					float distance=vector.Distance(tempmesh->verts[i],mesh1->verts[j]);
					if(distance<alpha1*mesh1->weight[j]){
						bFound=TRUE;VertexUsed[j]=TRUE;
						tempmesh->verts[j]=tempmesh->verts[i];
						//change the index of each face
						for(k1=0;k1<tempmesh->nfaces;k1++){
							for(k2=0;k2<tempmesh->faces[k1].nverts;k2++){
								if (tempmesh->faces[k1].verts[k2]==i)
									tempmesh->faces[k1].verts[k2]=j;
							}
						}
					}					
				}
			}
			if(!bFound) {alpha1+=alpha2;changenum++;i--;}
		}
		delete []VertexUsed;
		return tempmesh;
}
Mesh *Recover::resample(Mesh *mesh1,Mesh *mesh2){
	printf("resampling....please wait\n");
	Mesh *tempmesh=new Mesh(mesh1);
	int i,matchi;float th1,th2,th3,d;
	int ikey,iNeighb; 
    for (i=0;i<mesh1->nverts;i++){mesh1->VertexExcluded[i]=FALSE;}
    for (i=0;i<mesh1->nverts;i++){
	    ikey=(i+key)%mesh1->nverts;
	    if(!mesh1->VertexExcluded[ikey]){
			d=tools.FindND(mesh1,mesh2,ikey,&matchi);//find the nearest vertex
		    th2=mesh1->weight[ikey]*1.5;
			th1=mesh1->weight[ikey]*0.5;
			if(d<th2 && d>th1)
				tempmesh->verts[ikey]=mesh2->verts[matchi];
			iNeighb=-1;
			while(mesh1->neighbourhood[ikey][++iNeighb]!=-1) {//exclude the neighbourhoods
				mesh1->VertexExcluded[mesh1->neighbourhood[ikey][iNeighb]]=TRUE;
			}
		}
	}
    for(i=0;i<mesh1->nverts;i++){
	  if(mesh1->VertexExcluded[i]){
		d=tools.FindND(mesh1,mesh2,i,&matchi);
		th3=0.003*vector.Distance(mesh1->verts[i],mesh1->centriod);
		if(d<th3) tempmesh->verts[i]=mesh2->verts[matchi];
	  }
	}
	return tempmesh; 
}


void Matrix::nrerror(char error_text[]) { 
		 fprintf(stderr,"Numerical Recipes Run-Time Error:\n"); 
		 fprintf(stderr,"%s\n",error_text); 
		 fprintf(stderr,"Now Exiting to system\n"); 
		 exit(1); 
		} 

double *Matrix::vectorm(int nl,int nh) { 
		 double *v; 
  
		 v = (double *) malloc((unsigned) (nh-nl+1) * sizeof(double)); 
		 if (!v) nrerror("allocation failure in vectorm()"); 
		 return v-nl; 
		} 

double **Matrix::matrix(int nrl,int nrh,int ncl,int nch) { 
		 int i; 
		 double **m; 

		 m = (double **) malloc((unsigned) (nrh-nrl+1) * sizeof(double*)); 
		 if (!m) nrerror("row allocation failure in matrix(r,c)"); 
		 m -= nrl; 
		 for (i=nrl;i<=nrh;i++) 
			{ 
			 m[i] = (double *) malloc((unsigned) (nch-ncl+1)*sizeof(double)); 
			 if (!m[i]) nrerror("column allocation failure in matrix(r,c)"); 
			 m[i] -= ncl; 
			} 
		 return m; 
		} 

void Matrix::free_vector(double *v,int nl,int nh){ 
		 free((char *) (v+nl)); 
		} 

void Matrix::free_matrix(double **m,int nrl,int nrh,int ncl,int nch){ 
		 int i; 

		 for (i=nrh;i>=nrl;i--) free((char *) (m[i] + ncl)); 
		 free((char *) (m + nrl)); 
		} 
void Matrix::inverse(double **mat, int dim) { 
	 int i,j,*indx; 
	 double **y,d,*col; 

	 y = matrix(0,dim-1,0,dim-1); 
	 indx = (int *)malloc((unsigned)(dim*sizeof(int))); 
	 col = vectorm(0,dim-1); 
	 ludcmp(mat,dim,indx,&d); 
	 for (j=0;j<dim;j++) 
		{ 
		 for (i=0;i<dim;i++) col[i] = 0.0; 
		 col[j] = 1.0; 
		 lubksb(mat,dim,indx,col); 
		 for (i=0;i<dim;i++) y[i][j] = col[i]; 
		} 
	 for (i=0;i<dim;i++) 
		for (j=0;j<dim;j++) 
		   mat[i][j] = y[i][j]; 
	 free_matrix(y,0,dim-1,0,dim-1); 
	 free_vector(col,0,dim-1); 
	 free(indx); 
	} 

void Matrix::ludcmp(double **a, int n, int *indx, double *d) { 
		 int i,imax,j,k; 
		 double   big,dum,sum,temp; 
		 double   *vv; 

		 vv = (double*)malloc((unsigned)(n*sizeof(double))); 
		 if (!vv)  
		   { 
			fprintf(stderr,"Error Allocating Vector Memory\n"); 
			exit(1); 
		   } 
		 *d = 1.0; 
		 for (i=0;i<n;i++) 
			{ 
			 big = 0.0; 
			 for (j=0;j<n;j++) 
			 { 
			  if ((temp=fabs(a[i][j])) > big) big = temp; 
			 } 
			 if (big == 0.0) 
			   { 
				fprintf(stderr,"Singular Matrix in Routine LUDCMP\n"); 
			 for (j=0;j<n;j++) printf(" %f ",a[i][j]); printf("/n"); 
			 exit(1); 
			   } 
			 vv[i] = 1.0/big; 
			} 
		 for (j=0;j<n;j++) 
			{ 
			 for (i=0;i<j;i++) 
			 { 
			  sum = a[i][j]; 
			  for (k=0;k<i;k++) sum -= a[i][k] * a[k][j]; 
			  a[i][j] = sum; 
			 } 
			 big = 0.0; 
			 for (i=j;i<n;i++) 
				{ 
			  sum = a[i][j]; 
			  for (k=0;k<j;k++) sum -= a[i][k] * a[k][j]; 
			  a[i][j] = sum; 
			  if ((dum=vv[i]*fabs(sum)) >= big) 
				{ 
				 big = dum; 
				 imax = i; 
				} 
			 } 
			 if (j != imax) 
			   { 
			 for (k=0;k<n;k++) 
				{ 
				 dum = a[imax][k]; 
				 a[imax][k] = a[j][k]; 
				 a[j][k] = dum; 
				} 
				*d = -(*d); 
			 vv[imax] = vv[j]; 
			   } 
			 indx[j] = imax; 
			 if (a[j][j] == 0.0) a[j][j] = TINY; 
			 if (j != n-1) 
			   { 
			 dum = 1.0 / a[j][j]; 
			 for (i=j+1;i<n;i++) a[i][j] *= dum; 
			   } 
			} 
		 free(vv); 
		} 

void Matrix::lubksb(double **a, int n, int *indx, double *b) { 
		 int i,ip,j,ii=-1; 
		 double   sum; 

		 for (i=0;i<n;i++) 
			{ 
			 ip = indx[i]; 
			 sum = b[ip]; 
			 b[ip] = b[i]; 
			 if (ii>=0) 
			   for (j=ii;j<i;j++) sum -= a[i][j] * b[j]; 
			 else if (sum) ii = i; 
			 b[i] = sum; 
			} 
		 for (i=n-1;i>=0;i--) 
			{ 
			 sum = b[i]; 
			 for (j=i+1;j<n;j++) sum -= a[i][j] * b[j]; 
			 b[i] = sum / a[i][i]; 
			} 
		} 
void Matrix::print(double **m,int nrl,int nrh,int ncl,int nch){
		for(int i=nrl;i<nrh;i++){
			for(int j=ncl;j<nch;j++){
				printf("%4.4e\t",m[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
Matrix3 Matrix::AddMatrix3(Matrix3 &M1,Matrix3 &M2){
	Matrix3 matrix3;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			matrix3.element[i][j]=M1.element[i][j]+M2.element[i][j];
		}
	}
		return matrix3;
	}
Matrix3 Matrix::Mat3Div(Matrix3 &M,float n){
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			M.element[i][j]/=n;
		}
	}
	return M;
}
float Matrix::Mat3lemda(Matrix3 &M){
	Vector3 v;
	v.x=v.y=v.z=1;
	static int iteration=0;
	float lemda,lemda1;
	lemda=ComMat3lemda(M,&v,iteration);
	lemda1=ComMat3lemda(M,&v,iteration);
	if(lemda<0||lemda1<0){
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	while(fabs(lemda-lemda1)>(0.001*fabs(lemda))){
		lemda=lemda1;
		lemda1=ComMat3lemda(M,&v,iteration);
		if (iteration>10000){
			printf("warning:iteration may not be convergent\n");
			return ((fabs(lemda1)<1e-4)?1e-4:lemda1);
		}
	}
	float result=(fabs(lemda1)<1e-4)?1e-4:fabs(lemda1);
	result=result>1e+8?1e+8:result;
	return result;
}
float Matrix::ComMat3lemda(Matrix3 &M,Vector3 *v,int iteration){
	iteration++;
	*v=tools.Mat3MulV3(M,*v);
	float MaxLemda=vector.Max(*v);
	*v=*v/MaxLemda;
	return MaxLemda;
}
////////////////////////////////////////////////////////////
// GLOBAL VARIABLES
////////////////////////////////////////////////////////////

// Program arguments

// Mesh variables

static Mesh *mesh1 = new Mesh();
static Mesh *mesh2 = new Mesh();
static Mesh *mesh3 = new Mesh();
 bool formatright=FALSE;

// GLUT variables 

static int GLUTwindow = 0;
static int GLUTwindow_height = 1500;
static int GLUTwindow_width = 3600;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTarrows[4] = { 0, 0, 0, 0 };
static int GLUTmodifiers = 0;



// Display variables

static int scaling = 0;
static int translating = 0;
static int rotating = 0;
static float scale = 1.0;
static float center[3] = { 0.0, 0.0, 0.0 };
static float rotation[3] = { 0.0, 0.0, 0.0 };
static float translation[3] = { 0.0, 0.0, -4.0 };
////////////////////////////////////////////////////////////
// READ NOISE DATA
////////////////////////////////////////////////////////////
void ReadNoiseData(float *noise,int length,const char *noisefile){
	printf("reading noise data....\n");
	int i;
	int count=0;
  // Open file
  FILE *fpin;
  if (!(fpin = fopen(noisefile, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", noisefile);
	exit(0);
  }
  
  for(i=0;i<10000,count<length;i++){
	  fscanf(fpin,"%f",&(noise[count++]));
  }
  for(i=0;i<length;i++){
	  if(i>10000)
		  noise[i]=noise[i%10000];
  }
}
void LoadWatermark(unsigned char *Watermark,const char *watermarkfile){
  printf("loading watermark data....\n");
  FILE *fpin;int data;
  if (!(fpin = fopen(watermarkfile, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", watermarkfile);
	exit(0);
  }  
  for(int i=0;i<WatermarkLength;i++){
	  fscanf(fpin,"%d",&data);
	  Watermark[i]=data;
  }
}



////////////////////////////////////////////////////////////
// OFF FILE READING CODE
////////////////////////////////////////////////////////////

Mesh *ReadOffFile(const char *filename)
{
  int i;

  // Open file
  FILE *fpin;
  if (!(fpin = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }
  // Allocate mesh structure
  Mesh *mesh = new Mesh();
  if (!mesh) {
    fprintf(stderr, "Unable to allocate memory for file %s\n", filename);
    fclose(fpin);
    return 0;
  }

  // Read file
  int nverts = 0;
  int nfaces = 0;
  int nedges = 0;
  int line_count = 0;
  char buffer[1024];
  printf("reading vertices and faces of model.....");
  while (fgets(buffer, 1023, fpin)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Check section
    if (nverts == 0) {//in the line 1
      // Read header 
		if (strstr(buffer,"OFF")&&(strlen(buffer)==4)) {
			formatright=TRUE;
			continue;
		}
		
		if (!formatright) {
          fprintf(stderr, "The format of the file %s is %s,not of OFF format\n", filename,buffer);
          fclose(fpin);
          return NULL;
        }
        // Read mesh counts
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) {
          fprintf(stderr, "Syntax error reading header on line %d in file %s\n", line_count, filename);
          fclose(fpin);
          return NULL;
        }

        // Allocate memory for mesh
        mesh->verts = new Vector3 [nverts];
        assert(mesh->verts);
        mesh->faces = new Face [nfaces];
        assert(mesh->faces);
      }

    else if (mesh->nverts < nverts) {
      // Read vertex coordinates
      Vector3& vert = mesh->verts[mesh->nverts++];
      if (sscanf(bufferp, "%f%f%f", &(vert.x), &(vert.y), &(vert.z)) != 3) {
        fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
        fclose(fpin);
        return NULL;
      }
    }
    else if (mesh->nfaces < nfaces) {
      // Get next face
      Face& face = mesh->faces[mesh->nfaces++];

      // Read number of vertices in face 
      bufferp = strtok(bufferp, " \t");
      if (bufferp) face.nverts = atoi(bufferp);
      else {
        fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
        fclose(fpin);
        return NULL;
      }

      // Allocate memory for face vertices
      face.verts = new int[face.nverts];
      assert(face.verts);

      // Read vertex indices for face
      for (i = 0; i < face.nverts; i++) {
        bufferp = strtok(NULL, " \t");
        if (bufferp) face.verts[i] = atoi(bufferp);
        else {
          fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
          fclose(fpin);
          return NULL;
        }
      }

      // Compute normal for face
	  tools.FaceNormal(mesh,face);
      // Normalize normal for face
      face.normal=vector.Normalize(face.normal);
    }



    else {
      // Should never get here
      fprintf(stderr, "Found extra text starting at line %d in file %s\n", line_count, filename);
      break;
    }

  }//get buffer end

  // Check whether read all faces
  if (nfaces != mesh->nfaces) 
    fprintf(stderr, "Expected %d faces, but read only %d faces in file %s\n", nfaces, mesh->nfaces, filename);
  // Close file
  fclose(fpin);


  mesh->centriod=tools.mean(mesh->verts,mesh->nverts);    //compute centriod of the mesh
  mesh->magnitude=0;
  for(i=0;i<mesh->nverts;i++) {//compute magnitude of the mesh
	  float d=vector.Distance(mesh->verts[i],mesh->centriod);
	  if(mesh->magnitude<d) mesh->magnitude=d;
  }
    printf("end\n");
  printf("finding neighbourhoods and computing normal of each vertex.....");
    //find the neighbourhoods and compute normal of each vertex
  mesh->neighbourhood=new int[mesh->nverts][MaxNeighb];assert(mesh->neighbourhood);
  for (i=0;i<mesh->nverts;i++)
	for(int j=0;j<MaxNeighb;j++)
		mesh->neighbourhood[i][j]=-1;
  mesh->VertexNormal=new Vector3[mesh->nverts];assert(mesh->VertexNormal);
  mesh->weight=new float[mesh->nverts];assert(mesh->weight);  
  mesh->VertexExcluded=new bool[mesh->nverts];assert(mesh->VertexExcluded);
  mesh->VNDirection=new bool[mesh->nverts];assert(mesh->VNDirection);
  for(i=0;i<mesh->nverts;i++){mesh->VertexExcluded[i]=mesh->VNDirection[i]=FALSE;mesh->weight[i]=0;}
  for (i=0;i<mesh->nverts;i++){
	  int nFace=0;int iNeighb,tempneighb;bool m_bFound=FALSE;
	  for (int j=0;j<mesh->nfaces;j++){
		  for(int k=0;k<mesh->faces[j].nverts;k++){
			  if (i==mesh->faces[j].verts[k]){
				  k=mesh->faces[j].nverts;//end searching
				  //compute normal of each vertex
				  Face& face = mesh->faces[j];
				  mesh->VertexNormal[i]+=vector.Normalize(face.normal);
				  nFace++;
				  //find the neighbourhoods
				  for(int kk=0;kk<mesh->faces[j].nverts;kk++){/*neighbourhood face of vertex i is found*/
					 m_bFound=FALSE;
					iNeighb=-1;while(mesh->neighbourhood[i][++iNeighb]!=-1 && !m_bFound){
						if(mesh->neighbourhood[i][iNeighb]==mesh->faces[j].verts[kk])
							m_bFound=TRUE;/*if mesh->faces[j].verts[kk] is found in vertex i's neighbourhood,do nothing*/
					}
					if(!m_bFound && i!=mesh->faces[j].verts[kk])
					{
						iNeighb=-1;while(mesh->neighbourhood[i][++iNeighb]!=-1);/*trace the iNeighb*/
						mesh->neighbourhood[i][iNeighb]=mesh->faces[j].verts[kk];/*append vertex i's neighbourhood*/
						tempneighb=-1;while(mesh->neighbourhood[mesh->faces[j].verts[kk]][++tempneighb]!=-1);/*trace the tempneighb*/
						mesh->neighbourhood[mesh->faces[j].verts[kk]][tempneighb]=i;/*append vertex mesh->faces[j].verts[kk]'s neighbourhood*/
					}
				  }
			  }
		  }
	  }
	  mesh->VertexNormal[i]=mesh->VertexNormal[i]/nFace;
	  mesh->VertexNormal[i]=vector.Normalize(mesh->VertexNormal[i]);
//printf("%f %f %f\n",mesh->VertexNormal[i].x,mesh->VertexNormal[i].y,mesh->VertexNormal[i].z);
/* printf("%s %d %s","Neighbourhoods of Vertex",i,"is:");
	  for(int t=0;t<MaxNeighb;t++){if(mesh->neighbourhood[i][t]!=-1) printf(" %d",mesh->neighbourhood[i][t]);}
	  printf("\n");*/ 
  }
    printf("end\n");

  // Return mesh 
  return mesh;
}










////////////////////////////////////////////////////////////
// GLUT USER INTERFACE CODE
////////////////////////////////////////////////////////////

void GLUTRedraw(void)
{
  // Setup viewing transformation
  glLoadIdentity();
  glScalef(scale, scale, scale);
  glTranslatef(translation[0], translation[1], 0.0);

  // Set projection transformation
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, (GLfloat) GLUTwindow_width /(GLfloat) GLUTwindow_height, 0.1, 100.0);

  // Set camera transformation
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(translation[0], translation[1], translation[2]);
  glScalef(scale, scale, scale);
  glRotatef(rotation[0], 1.0, 0.0, 0.0);
  glRotatef(rotation[1], 0.0, 1.0, 0.0);
  glRotatef(rotation[2], 0.0, 0.0, 1.0);
  glTranslatef(-center[0], -center[1], -center[2]);

  // Clear window 
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set lights
  static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  // Set material
  static GLfloat material[] = { 0.9, 0.2, 0.0, 1.0 };
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 

  // Draw faces
  for (int i = 0; i < mesh2->nfaces; i++) {
    Face& face = mesh2->faces[i];
    glBegin(GL_TRIANGLE_FAN);
	//glBegin(GL_LINES);
	//glBegin(GL_POLYGON);
	float normal[3]={face.normal.x,face.normal.y,face.normal.z};
	glNormal3fv(normal);
    for (int j = 0; j < face.nverts; j++) {
      Vector3 &vert = mesh2->verts[face.verts[j]];
      glVertex3f(vert.x, vert.y, vert.z);
    }
    glEnd();
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}



void GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse motion event
  if (rotating) {
    // ChangeOrder model
    rotation[0] += -0.5 * (y - GLUTmouse[1]);
    rotation[2] += 0.5 * (x - GLUTmouse[0]);
  }
  else if (scaling) {
    // Scale window
    scale *= exp(2.0 * (float) (x - GLUTmouse[0]) / (float) GLUTwindow_width);
  }
  else if (translating) {
    // Translate window
    translation[0] += 2.0 * (float) (x - GLUTmouse[0]) / (float) GLUTwindow_width;
    translation[1] += 2.0 * (float) (y - GLUTmouse[1]) / (float) GLUTwindow_height;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse button event
  rotating = (button == GLUT_LEFT_BUTTON);
  scaling = (button == GLUT_MIDDLE_BUTTON);
  translating = (button == GLUT_RIGHT_BUTTON);
  if (rotating || scaling || translating) glutIdleFunc(GLUTRedraw);
  else glutIdleFunc(0);

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Process keyboard button event 
  switch (key) {
  case 'Q':
  case 'q':
    GLUTStop();
    break;

  case 27: // ESCAPE
    GLUTStop();
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = GLUTwindow_height - y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();
}



void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH| GLUT_STENCIL); 
  GLUTwindow = glutCreateWindow("Current 3D Model");

  // Initialize GLUT callback functions 
  glutReshapeFunc(GLUTResize);
  glutDisplayFunc(GLUTRedraw);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
  glutIdleFunc(0);
    
  // Initialize lights 
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glEnable(GL_LIGHT0);
  static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glEnable(GL_LIGHT1);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  // Initialize graphics modes 
  glEnable(GL_DEPTH_TEST);
}



void GLUTMainLoop(void)
{
  // Compute bounding box
  float bbox[2][3] = { { 1.0E30F, 1.0E30F, 1.0E30F }, { -1.0E30F, -1.0E30F, -1.0E30F } };
  for (int i = 0; i < mesh2->nverts; i++) {
    Vector3& vert = mesh2->verts[i];
    if (vert.x < bbox[0][0]) bbox[0][0] = vert.x;
    else if (vert.x > bbox[1][0]) bbox[1][0] = vert.x;
    if (vert.y < bbox[0][1]) bbox[0][1] = vert.y;
    else if (vert.y > bbox[1][1]) bbox[1][1] = vert.y;
    if (vert.z < bbox[0][2]) bbox[0][2] = vert.z;
    else if (vert.z > bbox[1][2]) bbox[1][2] = vert.z;
  }

  // Setup initial viewing scale
  float dx = bbox[1][0] - bbox[0][0];
  float dy = bbox[1][1] - bbox[0][1];
  float dz = bbox[1][2] - bbox[0][2];
  scale = 2.0 / sqrt(dx*dx + dy*dy + dz*dz);

  // Setup initial viewing center
  center[0] = 0.5 * (bbox[1][0] + bbox[0][0]);
  center[1] = 0.5 * (bbox[1][1] + bbox[0][1]);
  center[2] = 0.5 * (bbox[1][2] + bbox[0][2]);

  // Run main loop -- never returns 
  glutMainLoop();
}





////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////

char *ParseArgs(int argc, char **argv,char *filename)
{
  // Innocent until proven guilty
  int print_usage = 0;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-help")) { print_usage = 1; }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!filename) filename = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }
  if (print_usage)     printf("Usage: offviewer <filename>\n");

  // Check filename
  if (!filename) 
	  return NULL;
  else
      return filename;
}






///////////////////////////////////////////////////////////
// MODIFY DATA
///////////////////////////////////////////////////////////
Mesh *EmbedWatermark(Mesh *mesh,unsigned char *Watermark){//embedding method of literature [33]
  printf("please input embedding method(0,1,2),-1 for quit\n");
  scanf("%d",&iMethod);
  switch(iMethod) {
  case 0: 
	alpha=0.0450;
    mesh2=embed.EmbedWatermark0(mesh1,Watermark,alpha);
	break;
  case 1:
	alpha=0.2112;
    mesh2=embed.EmbedWatermark1(mesh1,Watermark,alpha);
	break;
  case 2:
	alpha=0.2112;
    mesh2=embed.EmbedWatermark2(mesh1,Watermark,alpha);
  	break;
  case -1:
	  return mesh;
	break;
  default:
	alpha=0.2112;
    mesh2=embed.EmbedWatermark1(mesh1,Watermark,alpha);
  }
  return mesh2;
}
Mesh *Embed::EmbedWatermark0(Mesh *mesh,unsigned char *Watermark,float alpha){//embedding method of literature [33]
	printf("embedding watermark(method 0).....\n");
	int i,j,k;
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,Watermark,WatermarkLength);
	permute.PSEUDO(Watermark1,Watermark,WatermarkColumn,WatermarkColumn,key);
	Vector3 *v=mesh->verts;
	float *weight=mesh->weight;
	Face *f=mesh->faces;
	float mod=MBIG;
	float S,s,a,b,c,d,dis;
	float thr=0.50;
  //compute the weight of each vertex
	for (i=0;i<mesh->nverts;i++){
		weight[i]=MBIG;
		  for (j=0;j<mesh->nfaces;j++){
			  for(k=0;k<f[j].nverts;k++){
				  if (i==f[j].verts[k]) {//compute each weight[i]
					  k=f[i].nverts;//end searching
					  a=vector.Distance(v[f[i].verts[0]],v[f[i].verts[1]]);
					  b=vector.Distance(v[f[i].verts[1]],v[f[i].verts[2]]);
					  c=vector.Distance(v[f[i].verts[2]],v[f[i].verts[0]]);
					  s=(a+b+c)/2;
					  S=sqrt(s*(s-a)*(s-b)*(s-c));
					  if(i==f[j].verts[0]) d=b;
					  else if(i==f[j].verts[1]) d=c;
					  else d=a;
					  dis=2*S/d;
					  mod=vector.Normalize(v[i]-mesh->centriod)/f[i].normal;
					  mod=fabs(1.0/mod*dis*alpha);
				  }
			  }
			  if(weight[i]>mod)	
				  weight[i]=mod;//find the minimum weight[i]
		  }
		  if(weight[i]>thr){
			weight[i]=thr;
			printf("warning\n");
		  }
		//printf("%s %d %s %f\n","Weight of Vertex",i,"is:",weight[i]);
	}
  // Allocate mesh structure
  Mesh *tempmesh = new Mesh(mesh);
  if (!tempmesh) {
    fprintf(stderr, "Unable to allocate memory to embed watermark\n");
    return 0;
  }
  //definite and initialize variables
  int iBit=0;int ikey=0;int iNeighb=0;Vector3 temp; 
  	v=tempmesh->verts;
	weight=tempmesh->weight;
	f=tempmesh->faces;
  for (i=0;i<tempmesh->nverts;i++){tempmesh->VertexExcluded[i]=FALSE;}
  for (i=0;i<tempmesh->nverts;i++){
	  //embed the watermark with the weight
	  iNeighb=0;
	  ikey=(i+key)%tempmesh->nverts;
	  if(!tempmesh->VertexExcluded[ikey]){
		  temp=vector.Direction(v[ikey],mesh->centriod);
		  temp=temp*(weight[ikey]*(Watermark[iBit++%WatermarkLength]==0?-1:1));
		  v[ikey]+=temp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {
			  tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;
		  }
	  }
  }
  tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);//recompute the centriod of mesh
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",tools.SubEnergy(tempmesh,mesh));//display watermark energy
  delete []Watermark1;
//return the watermarked model
  return tempmesh;
}
Mesh *Embed::EmbedWatermark1(Mesh *mesh,unsigned char *Watermark,float alpha){
	printf("embedding watermark(method 1).....\n");
	int i=0;
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,Watermark,WatermarkLength);
	permute.PSEUDO(Watermark1,Watermark,WatermarkColumn,WatermarkColumn,key);
	Vector3 *v=mesh->verts;
	float *weight=mesh->weight;
	Face *f=mesh->faces;
	double theta=0,z=0,dtemp=0;Vector3 vtemp;
  /*compute the weight of each vertex*/
	int iNeighb=-1;
	for(i=0;i<mesh->nverts;i++){
		iNeighb=-1;theta=0;vtemp.x=vtemp.y=vtemp.z=0;
		while(mesh->neighbourhood[i][++iNeighb]!=-1){
			dtemp=vector.Distance(v[i],v[mesh->neighbourhood[i][iNeighb]]);
			weight[i]+=1.0/dtemp;
			theta+=dtemp;
			vtemp+=v[mesh->neighbourhood[i][iNeighb]];
		}
		if(iNeighb==2)
			weight[i]=1.0/weight[i]*alpha;
		else
			weight[i]=1.0/weight[i]*alpha*iNeighb*sin(PI/2-PI/iNeighb);
		//printf("%f\n",weight[i]);
		theta/=iNeighb;
		if(iNeighb<5)
			theta/=2.5;
		else if (iNeighb<8)
			theta/=2;
		else theta/=1.5;
		vtemp=vtemp*(1.0/iNeighb);
		dtemp=vector.Distance(v[i],vtemp);
		theta=atan(dtemp/theta);
		z=fabs((v[i]-mesh->centriod)/mesh->VertexNormal[i]);
		if((1-z*z)*sin(theta)*sin(theta)>pow((1-z)*cos(theta),2)*(pow((1-cos(2*PI/iNeighb))/sin(2*PI/iNeighb),2)+1))
			mesh->VNDirection[i]=TRUE;
		if(mesh->VNDirection[i])
			printf("%f\t%f\t%f\t%f\t%d\n",theta*RTOD,z,(1-z*z)*sin(theta)*sin(theta),pow((1-z)*cos(theta),2)*(pow((1-cos(2*PI/iNeighb))/sin(2*PI/iNeighb),2)+1),mesh->VNDirection[i]);
	} 
  // Allocate mesh structure
  Mesh *tempmesh = new Mesh(mesh);
  if (!tempmesh) {
    fprintf(stderr, "Unable to allocate memory to embed watermark\n");
    return 0;
  }


  /*embed the watermark with the weight*/
  int iBit=0;int ikey=0;iNeighb=0;   //definite and initialize variables
  v=tempmesh->verts;
  weight=tempmesh->weight;
  f=tempmesh->faces;
  for (i=0;i<tempmesh->nverts;i++){tempmesh->VertexExcluded[i]=FALSE;}
  for (i=0;i<tempmesh->nverts;i++){
	  ikey=(i+key)%tempmesh->nverts;iNeighb=0;
	  if(!tempmesh->VertexExcluded[ikey]){
		  if(!mesh->VNDirection[ikey])
		      vtemp=vector.Direction(v[ikey],tempmesh->centriod);
		  else
			  vtemp=tempmesh->VertexNormal[ikey];
		  vtemp=vtemp*weight[ikey]*(Watermark[iBit++%WatermarkLength]==0?-1:1);
		  v[ikey]+=vtemp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;}
	  }
	//printf("%s %d %s %f\n","Weight of Vertex",i,"is:",weight[i]);
  }
    //compute the centriod of mesh
  tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",tools.SubEnergy(tempmesh,mesh));//display watermark energy
  delete []Watermark1;
//return the watermarked model
  return tempmesh;
}
Mesh *Embed::EmbedWatermark2(Mesh *mesh,unsigned char *Watermark,float alpha){
	printf("embedding watermark(method 2).....\n");
	int i=0;
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,Watermark,WatermarkLength);
	permute.PSEUDO(Watermark1,Watermark,WatermarkColumn,WatermarkColumn,key);
 	Vector3 *v=mesh->verts;
	float *weight=mesh->weight;
	Face *f=mesh->faces;
	//compute the weight of each vertex
      int iNeighb=0;
	  for(i=0;i<mesh->nverts;i++){
		     Matrix3 m=tools.Variance(mesh,i);
			 float lemda=mat.Mat3lemda(m);
			 //printf("%f\n",lemda);
			 mesh->weight[i]=1.0/lemda*2;

		  int tempiNeighb=0;float D0=0;
		  while(mesh->neighbourhood[i][tempiNeighb]!=-1){
				D0+=1.0/vector.Distance(mesh->verts[i],mesh->verts[mesh->neighbourhood[i][tempiNeighb++]]);
		  }
	  D0=1.0/D0*tempiNeighb*sin(PI/2-PI/tempiNeighb)*5;
	  D0=10;
		  if (mesh->weight[i]>D0){
			//printf("weight %d is %f,larger than D0 %f\n",i,mesh->weight[i],D0);
			mesh->weight[i]=D0;
		  }
			//printf("%f\n",mesh->weight[i]);

	  } 
  // Allocate mesh structure
  Mesh *tempmesh = new Mesh(mesh);
  if (!tempmesh) {
    fprintf(stderr, "Unable to allocate memory to embed watermark\n");
    return 0;
  }
  //definite and initialize variables
  int iBit=0;int ikey=0;iNeighb=0;Vector3 temp; 
  v=tempmesh->verts;
  weight=tempmesh->weight;
  f=tempmesh->faces;
  for (i=0;i<tempmesh->nverts;i++){tempmesh->VertexExcluded[i]=FALSE;}
  for (i=0;i<tempmesh->nverts;i++){
	  //embed the watermark with the weight
	  ikey=(i+key)%tempmesh->nverts;iNeighb=0;
	  if(!ikey){
		  ikey=0;
	  }
	  if(!tempmesh->VertexExcluded[ikey]){
		  temp=vector.Direction(v[ikey],mesh->centriod);
		  temp=tempmesh->VertexNormal[ikey];
		  temp=temp*weight[ikey]*(Watermark[iBit++%WatermarkLength]==0?-1:1);
		  v[ikey]+=v[ikey],temp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;}
	  }
	//printf("%s %d %s %f\n","Weight of Vertex",i,"is:",tempmesh->weight[i]);
  }
    //compute the centriod of mesh
  tempmesh->centriod=tools.mean(tempmesh->verts,tempmesh->nverts);
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",tools.SubEnergy(tempmesh,mesh));//display watermark energy
  delete []Watermark1;
//return the watermarked model
  return tempmesh;
}
///////////////////////////////////////////////////////////
// MODIFY DATA
///////////////////////////////////////////////////////////
void SaveModified(const char *filename1){
	//open file
	FILE *fpout;
	if (!(fpout = fopen(filename1, "w"))) {
    fprintf(stderr, "Unable to open file %s\n", filename1);
	fclose(fpout);
	exit(1);
	}
	//write header information
	fputs("OFF\n",fpout);
	char temp[20];
	ltoa(mesh2->nverts,temp,10);strcpy(temp+strlen(temp),"\0");strcat(temp," ");fputs(temp,fpout);
	ltoa(mesh2->nfaces,temp,10);strcpy(temp+strlen(temp),"\0");strcat(temp," ");fputs(temp,fpout);
	ltoa(0,temp,10);strcpy(temp+strlen(temp),"\0");strcat(temp,"\n");fputs(temp,fpout);
	//write vertices information
	for (int i=0;i<mesh2->nverts;i++){
		char buffer[1024];strcpy(buffer,"");
		char temp[20];strcpy(temp,"");
		gcvt(mesh2->verts[i].x,10,temp);strcat(buffer,temp);strcat(buffer," ");
		gcvt(mesh2->verts[i].y,10,temp);strcat(buffer,temp);strcat(buffer," ");
		gcvt(mesh2->verts[i].z,10,temp);strcat(buffer,temp);strcat(buffer,"\n");
		fputs(buffer,fpout);
	}
	//write topology information
	for (i=0;i<mesh2->nfaces;i++){
		char buffer[1024];strcpy(buffer,"");
		char temp[20];strcpy(temp,"");
		ltoa(mesh2->faces[i].nverts,temp,10);strcpy(temp+strlen(temp),"\0");strcat(buffer,temp);strcat(buffer," ");
		int nv=0;
		while (nv!=mesh2->faces[i].nverts) {
			ltoa(mesh2->faces[i].verts[nv++],temp,10);strcpy(temp+strlen(temp),"\0");strcat(buffer,temp);
			if (nv==mesh2->faces[i].nverts) 
				strcat(buffer,"\n");
			else strcat(buffer," ");
		}
		fputs(buffer,fpout);
	}
	//close file
	fclose(fpout);
}
/*permute watermark series*/
void Permute::PSEUDO(unsigned char *origwmdata,unsigned char *pmwmdata,int N1,int N2,int seed)
{
	int number=N1*N2;
	memcpy(pmwmdata,origwmdata,number);
  int bit[18];
  int i,j,scale,temp;
  int tempkey;
  int sequence[4096];
  if (number<4||number>(65536*2)) return;

    //取得m序列的阶数
	scale=0;//scale级反馈移位寄存器--用于产生长为2^scale-1的m序列     
    temp=number;
    while(temp!=1)//求阶数
	{    
	  temp=temp/2;
      scale=scale+1; 
	}
	
   tempkey=seed;//置随机序列的种子

   for (i=0;i<4095;i++) //周期为64*64,产生4095个伪随机值（种子给定）
	{
      for(j=0;j<scale;j++) //取各位0/1值-->bit[i]
	  {
        temp=tempkey>>j;
	    bit[j]=temp&1;
	  }
	  switch(scale)//作模2加法
	  {
	    case 2:
		  temp=bit[0]+bit[1];
		break;
		case 3:
          temp=bit[0]+bit[2];
		break;
		case 4:
		  temp=bit[0]+bit[3];
		break;
		case 5:
		  temp=bit[0]+bit[3];
		break;
		case 6:
		  temp=bit[0]+bit[5];
		break;
		case 7:
		  temp=bit[0]+bit[4];
		break;
		case 8:
		  temp=bit[0]+bit[4]+bit[5]+bit[6];
		break;
		case 9:
		  temp=bit[0]+bit[5];
		break;
		case 10:
		  temp=bit[0]+bit[7];
		break;
		case 11:
		  temp=bit[0]+bit[9];
		break;
		case 12:
		  temp=bit[0]+bit[6]+bit[8]+bit[11];
		break;
		case 13:
		  temp=bit[0]+bit[9]+bit[10]+bit[12];
		break;
		case 14:
          temp=bit[0]+bit[4]+bit[8]+bit[13];
	    break;
		case 15:
		  temp=bit[0]+bit[14];
		break;
		case 16:
		  temp=bit[0]+bit[4]+bit[13]+bit[15];
		break;
		case 17:
		  temp=bit[0]+bit[14];
		break;
		default:
		break;
	  }
	  bit[scale]=temp&1;
      tempkey=(int)(bit[scale]*(pow(2,(scale-1)))+(tempkey>>1));
	  sequence[i]=tempkey;

  }
     sequence[4095]=0;
	for(i=0;i<N1*N2;i++)
	{
		j=sequence[i];
		pmwmdata[j]=origwmdata[i];
	}
}
/*diverse permute watermark series*/
void Permute::DPSEUDO(unsigned char *pmwmdata,unsigned char *rewmdata,int N1,int N2,int seed)
{
  int number=N1*N2;
  int bit[18];
  int i,j,scale,temp;
  int tempkey;
  int sequence[4096];
  if (number<4||number>(65536*2)) return;

    //取得m序列的阶数
	scale=0;//scale级反馈移位寄存器--用于产生长为2^scale-1的m序列     
    temp=number;
    while(temp!=1)//求阶数
	{    
	  temp=temp/2;
      scale=scale+1; 
	}
	
   tempkey=seed;//置随机序列的种子

   for (i=0;i<4095;i++) //周期为64*64,产生4095个伪随机值（种子给定）
	{
      for(j=0;j<scale;j++) //取各位0/1值-->bit[i]
	  {
        temp=tempkey>>j;
	    bit[j]=temp&1;
	  }
	  switch(scale)//作模2加法
	  {
	    case 2:
		  temp=bit[0]+bit[1];
		break;
		case 3:
          temp=bit[0]+bit[2];
		break;
		case 4:
		  temp=bit[0]+bit[3];
		break;
		case 5:
		  temp=bit[0]+bit[3];
		break;
		case 6:
		  temp=bit[0]+bit[5];
		break;
		case 7:
		  temp=bit[0]+bit[4];
		break;
		case 8:
		  temp=bit[0]+bit[4]+bit[5]+bit[6];
		break;
		case 9:
		  temp=bit[0]+bit[5];
		break;
		case 10:
		  temp=bit[0]+bit[7];
		break;
		case 11:
		  temp=bit[0]+bit[9];
		break;
		case 12:
		  temp=bit[0]+bit[6]+bit[8]+bit[11];
		break;
		case 13:
		  temp=bit[0]+bit[9]+bit[10]+bit[12];
		break;
		case 14:
          temp=bit[0]+bit[4]+bit[8]+bit[13];
	    break;
		case 15:
		  temp=bit[0]+bit[14];
		break;
		case 16:
		  temp=bit[0]+bit[4]+bit[13]+bit[15];
		break;
		case 17:
		  temp=bit[0]+bit[14];
		break;
		default:
		break;
	  }
	  bit[scale]=temp&1;
      tempkey=(int)(bit[scale]*(pow(2,(scale-1)))+(tempkey>>1));
	  sequence[i]=tempkey;

  }
     sequence[4095]=0;
	for(i=0;i<N1*N2;i++)
	{
		j=sequence[i];
		rewmdata[i]=pmwmdata[j];
	}
}
Mesh *Rotate(Mesh *mesh2,Mesh *mesh3){//rotate the model
		printf("rotate the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh2;
  float a,b,r,dx,dy,dz;
  printf("please input the rotate parameters(a,b,r is degree and dx,dy,dz is percent):a,b,r,dx,dy,dz\n");
  scanf("%f,%f,%f,%f,%f,%f",&a,&b,&r,&dx,&dy,&dz);printf("\n");
  float mag=mesh2->magnitude;
  mesh2=attack.ChangeView(mesh2,mesh3,DTOR*a,DTOR*b,DTOR*r,mag*dx/100,mag*dy/100,mag*dz/100);
  return mesh2;
}

Mesh *Resample(Mesh *mesh1,Mesh *mesh2){//resample the attacked model
	printf("resample the model? 'y' for yes and 'q' for quit\t");
	char option;scanf("%s",&option);
	if (option=='q' || option=='Q') return mesh2;
    mesh2=recover.resample(mesh1,mesh2);
    return mesh2;
}
////////////////////////////////////////////////////////////
// Attact Watermark
////////////////////////////////////////////////////////////
void AttractWatermark(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
		printf("attract watermarked from the attacked model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return;
	switch(iMethod) {
	  case 0:
		attract.AttractWatermark0(mesh1,mesh2,Watermark);
  		break;
	  case 1:
		attract.AttractWatermark1(mesh1,mesh2,Watermark);
		break;
	  case 2:
		attract.AttractWatermark1(mesh1,mesh2,Watermark);
  		break;
	  default:
		attract.AttractWatermark1(mesh1,mesh2,Watermark);
  }
}
void Attract::AttractWatermark0(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
  printf("attractting watermark <method 0>...\n");
  unsigned char *AttracttedWatermark=new unsigned char[WatermarkLength];
  float *sum=new float[WatermarkLength];assert(sum);
	//make sure that when there is no watermark,NC is approximately zero
	srand((unsigned)time(NULL));
	for (int i=0;i<WatermarkLength;i++)	sum[i]=(float)(rand()-16383.5)/32767*1e-9;
	//attract watermark
	int iBit=0; int iNeighb=0;int ikey=0;Vector3 sub;
	for (i=0;i<mesh2->nverts;i++){mesh2->VertexExcluded[i]=FALSE;}
	for (i=0;i<mesh2->nverts;i++) {
	  ikey=(i+key)%mesh2->nverts;
	  if(!mesh2->VertexExcluded[ikey]){
		sub=mesh2->verts[ikey]-mesh1->verts[ikey];
		sum[iBit%WatermarkLength]+=vector.Magnitude(sub)*(sub/(mesh1->verts[ikey]-mesh1->centriod)>0?1:-1);
		  //exclude the neighbourhoods
		  iNeighb=0;
		  while(mesh2->neighbourhood[ikey][iNeighb]!=-1) {
			mesh2->VertexExcluded[mesh2->neighbourhood[ikey][iNeighb++]]=TRUE;
		  }
		  iBit++;
	  }	
	}
	for (i=0;i<WatermarkLength;i++){
		if (sum[i]>=0)
			AttracttedWatermark[i]=1;
		else AttracttedWatermark[i]=0;
	}
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,AttracttedWatermark,WatermarkLength);
	permute.DPSEUDO(Watermark1,AttracttedWatermark,WatermarkColumn,WatermarkColumn,key);
	delete []Watermark1;
	//display orignal watermark and attractted watermark
	int temp=0;
	for (i=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		if (Watermark[i]<0) Watermark[i]=0;
		printf("  %d",Watermark[i]);
	}
	printf("\n");
	for (i=0,temp=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		printf("  %d",AttracttedWatermark[i]);
	}
	printf("\n");
	//display the NC between the watermarks (method 1)
	/*int right=0;int wrong=0;
	for (i=0;i<WatermarkLength;i++){
		printf("%f\t",sum[i]);
		if (Watermark[i]==AttracttedWatermark[i]) right++;
		else wrong++;
	}
	printf("%s : %f\n","NC",(float)(right-wrong)/WatermarkLength);
	delete []sum;*/

	//display the NC between the watermarks (method 2)
	float mean1,mean2,sum1=0,sum2=0,sum3=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=AttracttedWatermark[i];
		sum2+=Watermark[i];
	}
	mean1=sum1/WatermarkLength;mean2=sum2/WatermarkLength;
	sum1=sum2=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=(AttracttedWatermark[i]-mean1)*(Watermark[i]-mean2);
		sum2+=(AttracttedWatermark[i]-mean1)*(AttracttedWatermark[i]-mean1);
		sum3+=(Watermark[i]-mean2)*(Watermark[i]-mean2);
	}
	float NC=sum1/sqrt(sum2*sum3);
	printf("\n%s : %f\n","NC",NC);
	delete []sum;
	delete []AttracttedWatermark;
}
void Attract::AttractWatermark1(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
	printf("attractting watermark...\n");
  unsigned char *AttracttedWatermark=new unsigned char[WatermarkLength];
  float *sum=new float[WatermarkLength];assert(sum);
	//make sure that when there is no watermark,NC is approximately zero
	srand((unsigned)time(NULL));
	for (int i=0;i<WatermarkLength;i++)	sum[i]=(float)(rand()-16383.5)/32767*1e-9;
	//attract watermark
	int iBit=0; int iNeighb=0;int ikey=0;Vector3 sub;
	for (i=0;i<mesh2->nverts;i++){mesh2->VertexExcluded[i]=FALSE;}
	for (i=0;i<mesh2->nverts;i++) {
	  ikey=(i+key)%mesh2->nverts;
	  if(!mesh2->VertexExcluded[ikey]){
		sub=mesh2->verts[ikey]-mesh1->verts[ikey];
		if(mesh2->VNDirection[ikey])
			sum[iBit%WatermarkLength]+=vector.Magnitude(sub)*(sub/mesh1->VertexNormal[ikey]>0?1:-1);
		else
			sum[iBit%WatermarkLength]+=vector.Magnitude(sub)*(sub/(mesh1->verts[ikey]-mesh1->centriod)>0?1:-1);			
		  //exclude the neighbourhoods
		  iNeighb=0;
		  while(mesh2->neighbourhood[ikey][iNeighb]!=-1) {
			mesh2->VertexExcluded[mesh2->neighbourhood[ikey][iNeighb++]]=TRUE;
		  }
		  iBit++;
	  }	
	}
	for (i=0;i<WatermarkLength;i++){
		if (sum[i]>=0)
			AttracttedWatermark[i]=1;
		else AttracttedWatermark[i]=0;
	}
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,AttracttedWatermark,WatermarkLength);
	permute.DPSEUDO(Watermark1,AttracttedWatermark,WatermarkColumn,WatermarkColumn,key);
	delete []Watermark1;
	//display orignal watermark and attractted watermark
	int temp=0;
	for (i=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		if (Watermark[i]<0) Watermark[i]=0;
		printf("  %d",Watermark[i]);
	}
	printf("\n");
	for (i=0,temp=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		printf("  %d",AttracttedWatermark[i]);
	}
	printf("\n");
	//display the NC between the watermarks (method 1)
	/*int right=0;int wrong=0;
	for (i=0;i<WatermarkLength;i++){
		printf("%f\t",sum[i]);
		if (Watermark[i]==AttracttedWatermark[i]) right++;
		else wrong++;
	}
	printf("%s : %f\n","NC",(float)(right-wrong)/WatermarkLength);
	delete []sum;*/

	//display the NC between the watermarks (method 2)
	float mean1,mean2,sum1=0,sum2=0,sum3=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=AttracttedWatermark[i];
		sum2+=Watermark[i];
	}
	mean1=sum1/WatermarkLength;mean2=sum2/WatermarkLength;
	sum1=sum2=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=(AttracttedWatermark[i]-mean1)*(Watermark[i]-mean2);
		sum2+=(AttracttedWatermark[i]-mean1)*(AttracttedWatermark[i]-mean1);
		sum3+=(Watermark[i]-mean2)*(Watermark[i]-mean2);
	}
	float NC=sum1/sqrt(sum2*sum3);
	printf("\n%s : %f\n","NC",NC);
	delete []sum;
	delete []AttracttedWatermark;
}




  

////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	/*filename is the orignal,filename1 is the watermarked
	  filename2 is the attacked,filename3 is the recovered*/
	 char filename[1024]="model\\";	 char tempfilename[1024];
	 printf("please input the model filename:");
	 scanf("%s",tempfilename);printf("\n");
	 strcat(filename,tempfilename);
	 //char *filename=NULL;filename=ParseArgs(argc,argv,filename);
     char filename1[1024];strcpy(filename1,"");
	 char filename2[1024];strcpy(filename2,"");
	 char filename3[1024];strcpy(filename3,"");
     strcat(filename1,filename);strcat(filename1,"watermarked.off");
	 strcat(filename2,filename);strcat(filename2,"attacked.off");
	 strcat(filename3,filename);strcat(filename3,"recovered.off");
     strcat(filename,".off");
  /* load the watermark*/
  unsigned char *Watermark=new unsigned char[WatermarkLength];
  LoadWatermark(Watermark,"watermark.dat");
  /* Read original model and allocate the detected model*/
//mesh1 = ReadOffFile("model\\m309.off");
//mesh2 = ReadOffFile("model\\m309watermarkedwatermarked.off");
  mesh1 = ReadOffFile(filename);//orignal model
  mesh2=new Mesh(mesh1);//the present model and is to be rendered
  mesh3=new Mesh(mesh1);//temporary mesh,use this mesh can reduce the time of allocating memory
  /*read noise data*/
  float *noise=new float[mesh1->nverts];
  ReadNoiseData(noise,mesh1->nverts,"noise.dat");

  mesh2=EmbedWatermark(mesh1,Watermark);  //Embed Watermark
  SaveModified(filename1);

  /*attack the watermarked model*/
  mesh2=attack.AddNoise(mesh2,noise);//always robust
  //mesh2=attack.Disturb(mesh2);//can be recovered by resort
  //mesh2=attack.ChangeOrder(mesh2,600);//can be recovered by resort
  //mesh2=attack.AddPulseNoise(mesh2,0.005);//always robust
  mesh2=attack.Cut(mesh2);//cut the model  
  //mesh2=attack.Move(mesh2);//move the model
  //mesh2=attack.Scale(mesh2);
  //mesh2=Rotate(mesh2,mesh3);
  mesh2=attack.DownSample2(mesh2);
 SaveModified(filename2);
  /*recover the attacked watermarked model*/
  mesh2=Resample(mesh1,mesh2);
  //mesh2=recover.resort(mesh1,mesh2);
  //mesh2=anneal.registration(mesh1,mesh2);
  SaveModified(filename3);
  AttractWatermark(mesh1,mesh2,Watermark);//Attract the Watermark



  /* Initialize GLUT*/
  GLUTInit(&argc, argv);
  /* Run GLUT interface */
  if(formatright){
	  mesh3=new Mesh(mesh2);
	  mesh2=attack.ChangeView(mesh2,mesh3,DTOR*0,DTOR*0,DTOR*0);
	  GLUTMainLoop();  //draw mesh2 
  }
  delete []noise;
  delete []Watermark;
  /* Return success */
  return 0;
}