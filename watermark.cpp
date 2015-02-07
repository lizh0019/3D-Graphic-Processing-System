        
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#include "vector.h"
#include "global.h"
#include "mesh.h"
#include "watermark.h"
extern Mesh *mesh1;
extern Mesh *mesh2;
extern Mesh *mesh3;
extern Mesh model;
extern Vector3 vector;
extern Matrix mat;
extern Embed embed;
extern Extract extract;
extern Permute permute;
extern Attack attack;
extern Recover recover;
extern int key;
extern int iMethod;
extern float alpha;
////////////////////////////////////////////////////////////
// CLASSES IMPLEMENTATIONS
////////////////////////////////////////////////////////////
Mesh *Recover::resort(Mesh *mesh1,Mesh *mesh2){//mesh1 is orignal,mesh2 is disturbed
		printf("resort the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh2;
		printf("resorting....");
		Mesh *tempmesh=new Mesh(mesh2);
		int i,j,k1,k2;
		/*reshape mesh2 according to mesh1*/
		Vector3 dcentroid=tempmesh->centroid-mesh1->centroid;
		float SizeMesh1(0),SizeMesh2(0);
		for(i=0;i<tempmesh->nverts;i++){
			SizeMesh1+=vector.Distance(mesh1->verts[i],mesh1->centroid);
			SizeMesh2+=vector.Distance(tempmesh->verts[i],tempmesh->centroid);
			tempmesh->verts[i]=tempmesh->verts[i]-dcentroid;
		}
		tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);
		SizeMesh1/=mesh1->nverts;SizeMesh2/=tempmesh->nverts;
		for(i=0;i<tempmesh->nverts;i++){
			Vector3 temp=tempmesh->verts[i]-tempmesh->centroid;
			temp=temp*(float)SizeMesh1/SizeMesh2;
			tempmesh->verts[i]=tempmesh->centroid+temp;
		}

		/*resort vertices of mesh2*/
	    bool *VertexUsed=new bool[mesh1->nverts];assert(VertexUsed);//if the match vertex is found,never use it
		for(i=0;i<mesh1->nverts;i++){VertexUsed[i]=FALSE;}
		float alpha1=0.005f;float alpha2=alpha1/10;bool bFound=FALSE;int changenum=0;
		float threashold=(float)sqrt((4*PI*mesh1->amplitude*mesh1->amplitude)/mesh1->nverts);
		for(i=0;i<mesh2->nverts;i++){
			if(bFound&&changenum) {
					alpha1-=alpha2*changenum;
					if(alpha1<=0) alpha1=alpha2;
				changenum=0;
			}
			bFound=FALSE;//if the match vertex is not found,change the finding range
			for(j=0;j<mesh1->nverts;j++){
				if (!VertexUsed[j]){
					float distance=vector.Distance(mesh2->verts[i],mesh1->verts[j]);
					if(distance<alpha1*threashold){
						bFound=TRUE;VertexUsed[j]=TRUE;
						tempmesh->verts[j]=mesh2->verts[i];
						//change the index of each face
						for(k1=0;k1<mesh2->nfaces;k1++){
							for(k2=0;k2<mesh2->faces[k1].nverts;k2++){
								if (mesh2->faces[k1].verts[k2]==i)
									tempmesh->faces[k1].verts[k2]=j;
							}
						}
					}	
				}
			}
			if(!bFound) {alpha1+=alpha2;changenum++;i--;}
		}
		delete []VertexUsed;
		printf("end\n");
		return tempmesh;
}
Mesh *Recover::resample(Mesh *mesh1,Mesh *mesh2){
	printf("resample the model? 'y' for yes and 'q' for quit\t");
	char option;scanf("%s",&option);
	if (option=='q' || option=='Q') return mesh2;
	printf("resampling....");
	Mesh *tempmesh=new Mesh(mesh1);
	int i,matchi;float th1,th2,d;
	int ikey,iNeighb; 
    for (i=0;i<mesh1->nverts;i++){mesh1->VertexExcluded[i]=FALSE;}
    for (i=0;i<mesh1->nverts;i++){
	    ikey=(i+key)%mesh1->nverts;
	    if(!mesh1->VertexExcluded[ikey]){
			d=model.FindND(mesh1,mesh2,ikey,&matchi);//find the nearest vertex
		    th2=mesh1->weight[ikey]*1.5f;
			th1=mesh1->weight[ikey]*0.5f;
			if(d<th2 && d>th1)
				tempmesh->verts[ikey]=mesh2->verts[matchi];
			iNeighb=-1;
			while(mesh1->neighbourhood[ikey][++iNeighb]!=-1) {//exclude the neighbourhoods
				mesh1->VertexExcluded[mesh1->neighbourhood[ikey][iNeighb]]=TRUE;
			}
		}
	}
	printf("end\n");
	return tempmesh; 
}




///////////////////////////////////////////////////////////
// MODIFY DATA
///////////////////////////////////////////////////////////
Mesh *EmbedWatermark(Mesh *mesh,unsigned char *Watermark){//embedding method of literature [33]
  printf("please input embedding method(0,1,2),-1 for quit\n");
  scanf("%d",&iMethod);
  switch(iMethod) {
  case 0: 
	alpha=0.04492f;
    mesh2=embed.EmbedWatermark0(mesh2,Watermark,alpha);
	break;
  case 1:
	alpha=0.2112f;
    mesh2=embed.EmbedWatermark1(mesh2,Watermark,alpha);
	break;
  case 2:
	alpha=2.9153e-3f;
    mesh2=embed.EmbedWatermark2(mesh2,Watermark,alpha);
  	break;
  case -1:
	  return mesh;
	break;
  default:
	alpha=0.2112f;
    mesh2=embed.EmbedWatermark1(mesh2,Watermark,alpha);
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
	float mod=MAXDATA;
	float S,s,a,b,c,d,dis;
	float thr=0.50;
  //compute the weight of each vertex
	for (i=0;i<mesh->nverts;i++){
		weight[i]=MAXDATA;
		  for (j=0;j<mesh->nfaces;j++){
			  for(k=0;k<f[j].nverts;k++){
				  if (i==f[j].verts[k]) {//compute each weight[i]
					  k=f[i].nverts;//end searching
					  a=vector.Distance(v[f[i].verts[0]],v[f[i].verts[1]]);
					  b=vector.Distance(v[f[i].verts[1]],v[f[i].verts[2]]);
					  c=vector.Distance(v[f[i].verts[2]],v[f[i].verts[0]]);
					  s=(a+b+c)/2;
					  S=(float)sqrt(s*(s-a)*(s-b)*(s-c));
					  if(i==f[j].verts[0]) d=b;
					  else if(i==f[j].verts[1]) d=c;
					  else d=a;
					  dis=2*S/d;
					  mod=vector.Normalize(v[i]-mesh->centroid)/f[i].normal;
					  mod=(float)fabs(1.0/mod*dis*alpha);
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
 	for(i=0;i<mesh->nverts;i++){mesh1->weight[i]=mesh->weight[i];}
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
		  temp=vector.Direction(v[ikey],mesh->centroid);
		  temp=temp*(weight[ikey]*(Watermark[iBit++%WatermarkLength]==0?-1:1));
		  v[ikey]+=temp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {
			  tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;
		  }
	  }
  }
  tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);//recompute the centroid of mesh
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",model.SubEnergy(tempmesh,mesh));//display watermark energy
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
			weight[i]+=1.0f/(float)dtemp;
			theta+=dtemp;
			vtemp+=v[mesh->neighbourhood[i][iNeighb]];
		}
		if(iNeighb==2)
			weight[i]=1.0f/weight[i]*alpha;
		else
			weight[i]=1.0f/weight[i]*alpha*iNeighb*(float)sin(PI/2-PI/iNeighb);
		//printf("%f\n",weight[i]);
		theta/=iNeighb;
		if(iNeighb<5)
			theta/=2.5;
		else if (iNeighb<8)
			theta/=2;
		else theta/=1.5;
		vtemp=vtemp*(1.0f/iNeighb);
		dtemp=vector.Distance(v[i],vtemp);
		theta=atan(dtemp/theta);
		z=fabs((v[i]-mesh->centroid)/mesh->VertexNormal[i]);
		if((1-z*z)*sin(theta)*sin(theta)>pow((1-z)*cos(theta),2)*(pow((1-cos(2*PI/iNeighb))/sin(2*PI/iNeighb),2)+1))
			mesh->VNDirection[i]=TRUE;
		//if(mesh->VNDirection[i])
			//printf("%f\t%f\t%f\t%f\t%d\n",theta*RTOD,z,(1-z*z)*sin(theta)*sin(theta),pow((1-z)*cos(theta),2)*(pow((1-cos(2*PI/iNeighb))/sin(2*PI/iNeighb),2)+1),mesh->VNDirection[i]);
	} 
	for(i=0;i<mesh->nverts;i++){mesh1->weight[i]=mesh->weight[i];}
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
		      vtemp=vector.Direction(v[ikey],tempmesh->centroid);
		  else
			  vtemp=tempmesh->VertexNormal[ikey];
		  vtemp=vtemp*weight[ikey]*(Watermark[iBit++%WatermarkLength]==0?-1.0f:1.0f);
		  v[ikey]+=vtemp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;}
	  }
	//printf("%s %d %s %f\n","Weight of Vertex",i,"is:",weight[i]);
  }
    //compute the centroid of mesh
  tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",model.SubEnergy(tempmesh,mesh));//display watermark energy
  delete []Watermark1;
//return the watermarked model
  return tempmesh;
}
Mesh *Embed::EmbedWatermark2(Mesh *mesh,unsigned char *Watermark,float alpha){
	printf("embedding watermark(method 2).....\n");
	int i,j,k;
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,Watermark,WatermarkLength);
	permute.PSEUDO(Watermark1,Watermark,WatermarkColumn,WatermarkColumn,key);
	Vector3 *v=mesh->verts;
	float *weight=mesh->weight;
	Face *f=mesh->faces;
	double theta=0,z=0,dtemp=0;Vector3 vtemp;Matrix m;float lemda;
	//compute the weight of each vertex
	  for(i=0;i<mesh->nverts;i++){
		     Matrix3 m3=model.Variance(mesh,i);
			 m=m.Add(m3,TINY);			 
			 mat.inverse(m.element,3);
			 for(j=0;j<3;j++)
				 for(k=0;k<3;k++)
					 m3.element[j][k]=m.element[j][k];
			 int iNeighb=-1;while(mesh->neighbourhood[i][++iNeighb]!=-1);
			 if(iNeighb>2){
				 lemda=mat.Mat3lemda(m3);
				 mesh->weight[i]=alpha*(1-exp(-1.0f/((1e-2f)*lemda)));
			 }
			 //printf("%f\n",lemda);
			 //printf("%f\n",mesh->weight[i]);

	  } 
	  mat.free_matrix(m.element,0,3,0,3);
	for(i=0;i<mesh->nverts;i++){mesh1->weight[i]=mesh->weight[i];}
  // Allocate mesh structure
  Mesh *tempmesh = new Mesh(mesh);
  if (!tempmesh) {
    fprintf(stderr, "Unable to allocate memory to embed watermark\n");
    return 0;
  }


  /*embed the watermark with the weight*/
  int iBit=0;int ikey=0;int iNeighb=0;   //definite and initialize variables
  v=tempmesh->verts;
  weight=tempmesh->weight;
  f=tempmesh->faces;
  for (i=0;i<tempmesh->nverts;i++){tempmesh->VertexExcluded[i]=FALSE;}
  for (i=0;i<tempmesh->nverts;i++){
	  ikey=(i+key)%tempmesh->nverts;iNeighb=0;
	  if(!tempmesh->VertexExcluded[ikey]){
		  if(!mesh->VNDirection[ikey])
		      vtemp=vector.Direction(v[ikey],tempmesh->centroid);
		  else
			  vtemp=tempmesh->VertexNormal[ikey];
		  vtemp=vtemp*weight[ikey]*(float)(Watermark[iBit++%WatermarkLength]==0?-1.0f:1.0f);
		  v[ikey]+=vtemp;
		  //exclude the neighbourhoods
		  while(tempmesh->neighbourhood[ikey][iNeighb]!=-1) {tempmesh->VertexExcluded[tempmesh->neighbourhood[ikey][iNeighb++]]=TRUE;}
	  }
	//printf("%s %d %s %f\n","Weight of Vertex",i,"is:",weight[i]);
  }
    //compute the centroid of mesh
  tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);
  memcpy(Watermark,Watermark1,WatermarkLength);
  printf("watermark energy is:%f\n",model.SubEnergy(tempmesh,mesh));//display watermark energy
  delete []Watermark1;
//return the watermarked model
  return tempmesh;
}


void Extract::ExtractWatermark0(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
  printf("extractting watermark <method 0>...\n");
  unsigned char *extracttedWatermark=new unsigned char[WatermarkLength];
  float *sum=new float[WatermarkLength];assert(sum);
	//make sure that when there is no watermark,NC is approximately zero
	int iBit=0,ikey=0,iNeighb=0,i;Vector3 sub;
	//make sure that when there is no watermark,NC is approximately zero
	for (i=0;i<WatermarkLength;i++){if(i%2) sum[i]=(float)TINY;else sum[i]=-(float)TINY;}
	//extract watermark
	for (i=0;i<mesh2->nverts;i++){mesh2->VertexExcluded[i]=FALSE;}
	for (i=0;i<mesh2->nverts;i++) {
	  ikey=(i+key)%mesh2->nverts;
	  if(!mesh2->VertexExcluded[ikey]){
		sub=mesh2->verts[ikey]-mesh1->verts[ikey];
		sum[iBit%WatermarkLength]+=vector.amplitude(sub)*(sub/(mesh1->verts[ikey]-mesh1->centroid)>0?1:-1);
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
			extracttedWatermark[i]=1;
		else extracttedWatermark[i]=0;
	}
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,extracttedWatermark,WatermarkLength);
	permute.DPSEUDO(Watermark1,extracttedWatermark,WatermarkColumn,WatermarkColumn,key);
	delete []Watermark1;
	//display orignal watermark and extractted watermark
	int temp=0;
	for (i=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		if (Watermark[i]<0) Watermark[i]=0;
		printf("  %d",Watermark[i]);
	}
	printf("\n");
	for (i=0,temp=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		printf("  %d",extracttedWatermark[i]);
	}
	printf("\n");
	//display the NC between the watermarks (method 1)
	/*int right=0;int wrong=0;
	for (i=0;i<WatermarkLength;i++){
		printf("%f\t",sum[i]);
		if (Watermark[i]==extracttedWatermark[i]) right++;
		else wrong++;
	}
	printf("%s : %f\n","NC",(float)(right-wrong)/WatermarkLength);
	delete []sum;*/

	//display the NC between the watermarks (method 2)
	float mean1,mean2,sum1=0,sum2=0,sum3=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=extracttedWatermark[i];
		sum2+=Watermark[i];
	}
	mean1=sum1/WatermarkLength;mean2=sum2/WatermarkLength;
	sum1=sum2=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=(extracttedWatermark[i]-mean1)*(Watermark[i]-mean2);
		sum2+=(extracttedWatermark[i]-mean1)*(extracttedWatermark[i]-mean1);
		sum3+=(Watermark[i]-mean2)*(Watermark[i]-mean2);
	}
	float NC=(float)((sum1+TINY)/sqrt((sum2+TINY)*sum3));
	printf("\n%s : %f\n","NC",NC);
	delete []sum;
	delete []extracttedWatermark;
}
void Extract::ExtractWatermark1(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
	printf("extractting watermark...\n");
  unsigned char *extracttedWatermark=new unsigned char[WatermarkLength];
  double *sum=new double[WatermarkLength];assert(sum);
	int iBit=0,ikey=0,iNeighb=0,i;Vector3 sub;
	//make sure that when there is no watermark,NC is approximately zero
	for (i=0;i<WatermarkLength;i++){if(i%2) sum[i]=TINY;else sum[i]=-TINY;}
	//extract watermark
	for (i=0;i<mesh2->nverts;i++){mesh2->VertexExcluded[i]=FALSE;}
	for (i=0;i<mesh2->nverts;i++) {
	  ikey=(i+key)%mesh2->nverts;
	  if(!mesh2->VertexExcluded[ikey]){
		sub=mesh2->verts[ikey]-mesh1->verts[ikey];
		if(mesh2->VNDirection[ikey])
			sum[iBit%WatermarkLength]+=vector.amplitude(sub)*(sub/mesh1->VertexNormal[ikey]>0?1:-1);
		else
			sum[iBit%WatermarkLength]+=vector.amplitude(sub)*(sub/(mesh1->verts[ikey]-mesh1->centroid)>0?1:-1);			
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
			extracttedWatermark[i]=1;
		else extracttedWatermark[i]=0;
	}
	unsigned char *Watermark1=new unsigned char[WatermarkLength];
	memcpy(Watermark1,extracttedWatermark,WatermarkLength);
	permute.DPSEUDO(Watermark1,extracttedWatermark,WatermarkColumn,WatermarkColumn,key);
	delete []Watermark1;
	//display orignal watermark and extractted watermark
	int temp=0;
	for (i=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		if (Watermark[i]<0) Watermark[i]=0;
		printf("  %d",Watermark[i]);
	}
	printf("\n");
	for (i=0,temp=0;i<WatermarkLength;i++){
		if((temp++ % WatermarkColumn)==0 && temp!=0) printf("\n"); 
		printf("  %d",extracttedWatermark[i]);
	}
	printf("\n");
	//display the NC between the watermarks (method 1)
	/*int right=0;int wrong=0;
	for (i=0;i<WatermarkLength;i++){
		printf("%f\t",sum[i]);
		if (Watermark[i]==extracttedWatermark[i]) right++;
		else wrong++;
	}
	printf("%s : %f\n","NC",(float)(right-wrong)/WatermarkLength);
	delete []sum;*/

	//display the NC between the watermarks (method 2)
	float mean1,mean2,sum1=0,sum2=0,sum3=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=extracttedWatermark[i];
		sum2+=Watermark[i];
	}
	mean1=sum1/WatermarkLength;mean2=sum2/WatermarkLength;
	sum1=sum2=0;
	for(i=0;i<WatermarkLength;i++){
		sum1+=(extracttedWatermark[i]-mean1)*(Watermark[i]-mean2);
		sum2+=(extracttedWatermark[i]-mean1)*(extracttedWatermark[i]-mean1);
		sum3+=(Watermark[i]-mean2)*(Watermark[i]-mean2);
	}
	float NC=(float)((sum1+TINY)/sqrt((sum2+TINY)*sum3));
	printf("\n%s : %f\n","NC",NC);
	delete []sum;
	delete []extracttedWatermark;
}


/*permute watermark series*/
void Permute::PSEUDO(unsigned char *origwmdata,unsigned char *pmwmdata,int N1,int N2,int key)
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
	
   tempkey=key;//置随机序列的种子

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
void Permute::DPSEUDO(unsigned char *pmwmdata,unsigned char *rewmdata,int N1,int N2,int key)
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
	
   tempkey=key;//置随机序列的种子

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

		int *indexcut=new int[VertsCut];//store the indexes cut
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
				model.DeleteVertex(tempmesh,indexcut[i]-displacement++);
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
			Noise=(float)rand();
			if (Noise>=16383.5) Noise=AlphaNoise;
			else Noise =-AlphaNoise;
			temp=mesh->verts[i]-mesh->centroid;
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
		printf("please input the noise amplitude:");
	    scanf("%f",&noisep);printf("\n");
		Mesh *tempmesh=new Mesh(mesh);
		int i;float max=0;Vector3 temp;
	    ReadNoiseData(noise,mesh1->nverts,"noise.dat");//read noise data
		for(i=0;i<mesh->nverts;i++){
			temp=mesh->verts[i]-mesh->centroid;
			float t=vector.amplitude(temp);
			if(max<t)
				max=t;
		}
		for(i=0;i<mesh->nverts;i++){
			noise[i]*=max*noisep/100;
		    //printf("noise[%3d]=%f\n",i,noise[i]);
			temp=mesh->verts[i]-mesh->centroid;
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
			tempmesh->verts[i]=vector.Rotate(mesh->verts[i],a,b,r,dx,dy,dz);
		}
		return tempmesh;
	}
	Mesh *Attack::Disturb(Mesh *mesh){
		printf("disturb the order of the model randomly? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("disturbing....");
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
		printf("end\n");
		return tempmesh;
	}
	Mesh *Attack::Move(Mesh *mesh){
		printf("move the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		Vector3 dcentroid(1,2,3);
		printf("moving whole model....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		for(int i=0;i<mesh->nverts;i++){
			tempmesh->verts[i]=mesh->verts[i]+dcentroid;
		}
		tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);
		return tempmesh;
	}
	Mesh *Attack::Scale(Mesh *mesh){
		printf("scale the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh;
		printf("scaling whole model....please wait\n");
		  float scale;
		  printf("please input the scale amplitude:");
		  scanf("%f",&scale);printf("\n");
		Mesh *tempmesh=new Mesh(mesh);
		for(int i=0;i<mesh->nverts;i++){
			Vector3 temp=mesh->verts[i]-mesh->centroid;
			temp=temp*scale;
			tempmesh->verts[i]=mesh->centroid+temp;
		}
		tempmesh->centroid=vector.mean(tempmesh->verts,tempmesh->nverts);
		return tempmesh;
	}
	Mesh *Attack::DownSample1(Mesh *mesh,float prop){/*downsample mesh2*/
		printf("downsampling.....please wait\n");
		Mesh *tempmesh=new Mesh(mesh);
		int i,j,cutvertex;
		int orignnum=tempmesh->nverts;
		srand((unsigned)time(NULL));
		for (i=0;i<orignnum*prop;i++){
			//model.DeleteVertex(tempmesh,rand()%tempmesh->nverts);//cut the vertex
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
				newface[iNeighb]=mesh->neighbourhood[i-countvert][iNeighb];//store the new face vertices
			}
			model.DeleteVertex(mesh,i-countvert);//delete a vertex
			//append the new face whose edge number is no less than 3
			mesh->nfaces++;
			//mesh->faces[mesh->nfaces-1]=new Face();
			iNeighb=-1;
			while(newface[++iNeighb]!=-1) {
				if(newface[iNeighb]>i-countvert) newface[iNeighb]--;
				mesh->VertexExcluded[newface[iNeighb]]=TRUE;//vertices connected to the vertex cut should not be cut
				mesh->faces[mesh->nfaces-1].verts[iNeighb]=newface[iNeighb];//(newface[iNeighb]>i-countvert)?(newface[iNeighb]-1):newface[iNeighb];//newface[iNeighb];
			}
			model.FaceNormal(mesh,mesh->faces[mesh->nfaces-1]);
			mesh->faces[mesh->nfaces-1].nverts=iNeighb;
			countvert++;
			//model.ReCompNeighb(mesh);//recompute neighbourhoods
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
				model.NeighbFaces(tempmesh,i,a,&countface);//compute neighbourhood face vertices of i
				iNeighb=0;
				while(tempmesh->neighbourhood[i][iNeighb]!=-1) 
					tempmesh->VertexExcluded[tempmesh->neighbourhood[i][iNeighb++]]=TRUE;//neighbourhoods should not be cut
				for(j=0;j<tempmesh->nfaces;j++){//when the vertex is cut,faces containing it should be cut
					for(k1=0;k1<tempmesh->faces[j].nverts;k1++){
						if(i==tempmesh->faces[j].verts[k1]) {
							k1=tempmesh->faces[j].nverts;
							model.DeleteFace(tempmesh,j--);
						}
					}
				}
				model.AppendFaces(tempmesh,i-countvert,a,countface);//append faces
				delete []a;
				model.DeleteVertex(tempmesh,i-countvert++);//cut the vertex
			}
		}
		model.DeleteOverlapFaces(tempmesh);
		//model.DeleteAbnormal(tempmesh);
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