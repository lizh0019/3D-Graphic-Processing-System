#ifndef _GLOBAL_FUNC
#define _GLOBAL_FUNC 0

#include "global.h"
#include "mesh.h"

////////////////////////////////////////////////////////////
// functions in main() declarations
////////////////////////////////////////////////////////////
void ExtractWatermark(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark);
Mesh *EmbedWatermark(Mesh *mesh,unsigned char *Watermark);
void GLUTInit(int *argc, char **argv);
void GLUTKeyboard(unsigned char keyboard, int x, int y);
void GLUTMainLoop(void);
void GLUTMotion(int x, int y);
void GLUTMouse(int button, int state, int x, int y);
void GLUTRedraw(void);
void GLUTResize(int w, int h);
void GLUTSpecial(int key, int x, int y);
void GLUTStop(void);
void LoadWatermark(unsigned char *Watermark,const char *watermarkfile);
char *ParseArgs(int argc, char **argv,char *filename);
void ReadNoiseData(float *noise,int length,const char *noisefile);
Mesh *ReadOffFile(const char *filename);
Mesh *Rotate(Mesh *mesh2,Mesh *mesh3);
void SaveModified(const char *filename1);




#endif
