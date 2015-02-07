#include "mesh.h"
#include "globalfunc.h"

extern Mesh *mesh1;
extern Mesh *mesh2;
extern Mesh *mesh3;
extern Mesh model;
extern Vector3 vector;
extern Matrix mat;

Mesh::Mesh(void) : nverts(0), verts(0), nfaces(0), faces(0) ,amplitude(0),\
		VertexNormal(0),neighbourhood(0),VertexExcluded(0),weight(0),VNDirection(0){};
Mesh::Mesh(Mesh *mesh):nverts(mesh->nverts),nfaces(mesh->nfaces),\
			centroid(mesh->centroid),amplitude(mesh->amplitude){
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
  delete []neighbourhood;
  delete []verts;
  delete []VertexExcluded;
  delete []weight;
  delete []VNDirection;
}
	float Mesh::SubEnergy(Mesh *mesh2,Mesh *mesh3){
		float energy=0;
		for (int i=0;i<mesh2->nverts;i++){
			energy+=vector.Distance(mesh2->verts[i],mesh3->verts[i])*vector.Distance(mesh2->verts[i],mesh3->verts[i]);
		}
		return energy;
	}

void Mesh::DeleteVertex(Mesh *mesh,int iVertex){
		int i,j,k,k1,temp;bool endsearch=FALSE;

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
		}
		mesh->nverts--;


}
	void Mesh::ReCompNeighb(Mesh *mesh){//recompute neighbourhoods
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
	void Mesh::DeleteFace(Mesh *mesh,int iFace){
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
	void Mesh::NeighbFaces(Mesh *mesh,int cut,int (*a)[2],int *countface){
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
	void Mesh::AppendFaces(Mesh *mesh,int cut,int (*a)[2],int countface){
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
	void Mesh::DeleteOverlapFaces(Mesh *mesh){
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
							model.DeleteFace(mesh,t2);
						}
					}
				}
	}
	void Mesh::DeleteAbnormal(Mesh *mesh){
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
					model.DeleteFace(mesh,deleteface);
					model.DeleteVertex(mesh,i-countvert1++);
					abnormalend=FALSE;
				}
			}
		}while(!abnormalend);
	}
	bool Mesh::IsOdd(int vertex,int (*a)[2],int countface){
		bool isodd=FALSE;
		for(int i=0;i<countface;i++){
			if((vertex==a[i][0])||(vertex==a[i][1]))
				isodd=!isodd;
		}
		return isodd;
	}
	bool Mesh::bIsNeighb(Mesh *mesh,int a,int b){
		int iNeighb1=-1;
		while(mesh->neighbourhood[a][++iNeighb1]!=-1){
			if((b!=a)&&(mesh->neighbourhood[a][iNeighb1]==b))
				return TRUE;
			}
		return FALSE;
	}
	Matrix3 Mesh::Variance(Mesh *mesh,int vertex){
	  int iNeighb=0;Vector3 average,temp;
	  while(mesh->neighbourhood[vertex][iNeighb]!=-1){
		  temp+=mesh->verts[mesh->neighbourhood[vertex][iNeighb++]];
	  }
	  temp=temp+mesh->verts[vertex];
	  average=temp/(float)(iNeighb+1);
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
	  matrix3=mat.Mat3Div(matrix3,(float)(iNeighb+1));
	  return matrix3;
	}
	float Mesh::FindND(Mesh *mesh1,Mesh *mesh2,int ivertex,int *matchi){
		long int i;float min=MAXDATA;float distance;
		for(i=0;i<mesh2->nverts;i++){
			distance=vector.Distance(mesh1->verts[ivertex],mesh2->verts[i]);
			if(min>distance){
				min=distance;
				*matchi=i;
			}
		}
		return min;
	}
	Vector3 Mesh::FaceNormal(Mesh *mesh,Face &face){
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


