#include "face.h"



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