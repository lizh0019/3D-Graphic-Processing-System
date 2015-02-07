#ifndef _TRIANGLE_FACE
#define _TRIANGLE_FACE 2

#include "global.h"
#include "vector.h"


class Face {
public:
  Face(void);
  Face(Face &face);
  Face(Face *face);
  ~Face();
  int nverts;
  int *verts;
  Vector3 normal;
};


#endif