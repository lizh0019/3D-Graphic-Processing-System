 #include "anneal.h"
extern Anneal anneal;
extern Matrix mat;
extern Vector3 vector;
extern Mesh model;
extern Attack attack;
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
	 cost=model.SubEnergy(mesh2,tempmesh);
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
	Vector3 dcentroid=mesh2->centroid-mesh1->centroid;
	for(i=0;i<mesh1->nverts;i++){
		mesh2->verts[i]=mesh2->verts[i]-dcentroid;
	}
	mesh2->centroid=mesh1->centroid;
	float num;
	printf("please input the anneal number:");
	 scanf("%f",&num);printf("\n");
	 printf("registrating using anneal....\n");
	Mesh *tempmesh=new Mesh(mesh1);
	float x[6],t[6];
	float a,b,r,dx,dy,dz;
  printf("please input the approximate rotate parameters(a,b,r is degree and dx,dy,dz is percent,relative error less than 20%s):a,b,r,dx,dy,dz\n","%");
  scanf("%f,%f,%f,%f,%f,%f",&a,&b,&r,&dx,&dy,&dz);printf("\n");
  float mag=mesh2->amplitude;
  float mx[6]={a*DTOR,b*DTOR,r*DTOR,0,0,0};//initial condition
    float c[6],mincost=MBIG,cost,dc;
    float up[6]={mx[0]*1.2f,mx[1]*1.2f,mx[2]*1.2f,mx[3]*1.2f,mx[4]*1.2f,mx[5]*1.2f};//up bound
	float low[6]={mx[0]*0.8f,mx[1]*0.8f,mx[2]*0.8f,mx[3]*0.8f,mx[4]*0.8f,mx[5]*0.8f};//low bound
    int m=0,j,n,s,k;
	for(i=0;i<6;i++){
		t[i]=0.005f;c[i]=25;
	}
		for(k=0;k<6;k++){
			mx[k]=low[k]+(int)((up[k]-low[k])*anneal.ran(&anneal.idum)/FAC)*FAC;
			    if (mx[k]>up[k]) mx[k]=up[k];
			    if(mx[k]<low[k])mx[k]=low[k];
		}
		anneal.T=1.0f;mincost=MBIG;
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
			Vector3 outv3(vector.Mat4MulV4(iRotMat,inv4));
			tempmesh->verts[i]=outv3;//实现重定位
		}
		mat.free_matrix(m1,0,4,0,4); //释放矩阵。 
		return tempmesh;
}

