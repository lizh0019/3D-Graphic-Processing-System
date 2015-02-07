#ifndef _TRIANGLE_MESH
#define _TRIANGLE_MESH 4

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <MATH.H>
#include "GraphicsGems.h"
#include "global.h"
#include "matrix.h"
#include "vector.h"
#include "face.h"


#ifndef MaxNeighb
#define MaxNeighb 250
#endif


class Mesh {
public:  
	Mesh(void);
	Mesh(Mesh *mesh);
	~Mesh();
  Vector3 centroid;
  int nverts;
  Vector3 *verts;
  int nfaces;
  Face *faces;
  Vector3 *VertexNormal;
  int (*neighbourhood)[MaxNeighb];
  bool *VertexExcluded;
  float *weight;
  bool *VNDirection;
  float amplitude;
	Vector3 FaceNormal(Mesh *mesh,Face &face);
	float SubEnergy(Mesh *mesh1,Mesh *mesh2);
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


};



#endif
