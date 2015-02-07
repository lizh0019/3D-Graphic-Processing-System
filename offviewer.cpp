////////////////////////////////////////////////////////////
// Source file for the OFF file viewer
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////
//Windows include files 
#ifdef _WIN32
#include <windows.h>
#endif
// OpenGL include files 
#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"
// include class headers
#include "global.h"
#include "globalfunc.h"
#include "vector.h"
#include "face.h"
#include "matrix.h"
#include "mesh.h"
#include "anneal.h"
#include "watermark.h"
////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////
Vector3 vector;
Face face;
Matrix mat;
Mesh *mesh1=new Mesh();
Mesh *mesh2=new Mesh();
Mesh *mesh3=new Mesh();
Mesh model;
Anneal anneal;
Embed embed;
Extract extract;
Permute permute;
Attack attack;
Recover recover;
int key=66;
int iMethod;
float alpha;
bool formatright=FALSE;
float *noise;
unsigned char *Watermark;
char filename[1024]="model\\";//the orignal
char filename1[1024];//the watermarked
char filename2[1024];//the attacked
char filename3[1024];//the recovered

////////////////////////////////////////////////////////////
// READ NOISE DATA
////////////////////////////////////////////////////////////
void ReadNoiseData(float *noise,int length,const char *noisefile){
  int i,count=0;
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
////////////////////////////////////////////////////////////
// LOAD WATERMARK SERIES
////////////////////////////////////////////////////////////
void LoadWatermark(unsigned char *Watermark,const char *watermarkfile){
  printf("loading watermark data....");
  FILE *fpin;int data;
  if (!(fpin = fopen(watermarkfile, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", watermarkfile);
	exit(0);
  }  
  for(int i=0;i<WatermarkLength;i++){
	  fscanf(fpin,"%d",&data);
	  Watermark[i]=data;
  }
  printf("end\n");
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
	  model.FaceNormal(mesh,face);
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


  mesh->centroid=vector.mean(mesh->verts,mesh->nverts);    //compute centroid of the mesh
  mesh->amplitude=0;
  for(i=0;i<mesh->nverts;i++) {//compute amplitude of the mesh
	  float d=vector.Distance(mesh->verts[i],mesh->centroid);
	  if(mesh->amplitude<d) mesh->amplitude=d;
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
    //glBegin(GL_TRIANGLE_FAN);
	glBegin(GL_LINES);
	//glBegin(GL_POLYGON);
	//float normal[3]={face.normal.x,face.normal.y,face.normal.z};
	//glNormal3fv(normal);
    for (int j = 0; j < face.nverts; j++) {
      Vector3 &vert = mesh2->verts[face.verts[j]];
      glVertex3f(vert.x, vert.y, vert.z);
    }
    glEnd();
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTStop(void){
    // Destroy window 
    glutDestroyWindow(GLUTwindow);
    // Exit
    exit(0);
}



void GLUTResize(int w, int h){
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



void GLUTKeyboard(unsigned char keyboard, int x, int y)
{
  // Process keyboard button event 
  switch (keyboard) {
  case 27: 
    GLUTStop();// ESCAPE
    break;
  case 'c':case 'C':
	  mesh2=attack.Cut(mesh2);SaveModified(filename2);
  break;
  case 'N':case 'n':
	  mesh2=attack.AddNoise(mesh2,noise);SaveModified(filename2);
  break;
  case 'D':case 'd':
	  mesh2=attack.DownSample2(mesh2);SaveModified(filename2);
  break;
  case 'R':case 'r':
	  mesh2=recover.resample(mesh1,mesh2);SaveModified(filename3);
  break;
  case 'E':case 'e':
	  mesh2=EmbedWatermark(mesh2,Watermark);SaveModified(filename1);
  break;
  case 'x':case 'X':
	  ExtractWatermark(mesh1,mesh2,Watermark);
  break; 
  case 32:
	  mesh2 = new Mesh(mesh1);
  break; 
  case 'S':case 's':
	  mesh2=attack.Scale(mesh2);SaveModified(filename2);
  break; 
  case 'O':case 'o':
	  mesh2=Rotate(mesh2,mesh3);SaveModified(filename2);
  break; 
  case 'G':case 'g':
	  mesh2=anneal.registration(mesh1,mesh2);SaveModified(filename3);
  break; 
  case 'U':case 'u':
      mesh2=attack.Disturb(mesh2);SaveModified(filename2);
  break;
  case 'V':case 'v':
	  mesh2=recover.resort(mesh1,mesh2);SaveModified(filename3);
  break;
  case 'M':case 'm':
	  mesh2=attack.Move(mesh2);SaveModified(filename2);
  break;  
  case 'K':case 'k':
	  printf("please input the key(0~255):");
	  int *tempkey=&key;
	  scanf("%d",tempkey);printf("%d\n",key);
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
  //glutCreateSubWindow(GLUTwindow,10,10,200,200);

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
//glutFullScreen();
  // Run main loop -- never returns 
	glutMainLoop();
}


////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////

char *ParseArgs(int argc, char **argv,char *filename){
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
// SAVE MODIFIED MODEL
///////////////////////////////////////////////////////////
void SaveModified(const char *filename1){
	//open file
	ofstream fpout(filename1);
	//write header information
	fpout<<"OFF"<<endl;
	fpout<<mesh2->nverts<<" "<<mesh2->nfaces<<" "<<0<<endl;
	
	//write vertices information
	for (int i=0;i<mesh2->nverts;i++){
		fpout<<mesh2->verts[i].x<<" "<<mesh2->verts[i].y<<" "<<mesh2->verts[i].z<<endl;
	}
	//write topology information
	for (i=0;i<mesh2->nfaces;i++){
		fpout<<mesh2->faces[i].nverts<<" ";
		int nv=0;
		while (nv!=mesh2->faces[i].nverts) {
			fpout<<mesh2->faces[i].verts[nv++];
			if (nv==mesh2->faces[i].nverts) 
				fpout<<endl;
			else fpout<<" ";
		}
	}
	fpout.close();
}

Mesh *Rotate(Mesh *mesh2,Mesh *mesh3){//rotate the model
		printf("rotate the model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return mesh2;
  float a,b,r,dx,dy,dz;
  printf("please input the rotate parameters(a,b,r is degree and dx,dy,dz is percent):a,b,r,dx,dy,dz\n");
  scanf("%f,%f,%f,%f,%f,%f",&a,&b,&r,&dx,&dy,&dz);printf("\n");
  float mag=mesh2->amplitude;
  mesh3=new Mesh(mesh2);
  mesh2=attack.ChangeView(mesh2,mesh3,DTOR*a,DTOR*b,DTOR*r,mag*dx/100,mag*dy/100,mag*dz/100);
  return mesh2;
}

////////////////////////////////////////////////////////////
//  Watermark RETRIEVING
////////////////////////////////////////////////////////////
void ExtractWatermark(Mesh *mesh1,Mesh *mesh2,unsigned char *Watermark){
		printf("extract watermarked from the attacked model? 'y' for yes and 'q' for quit\t");
		char option;scanf("%s",&option);
		if (option=='q' || option=='Q') return;
		printf("key:%d\n",key);
	switch(iMethod) {
	  case 0:
		extract.ExtractWatermark0(mesh1,mesh2,Watermark);
  		break;
	  case 1:
		extract.ExtractWatermark1(mesh1,mesh2,Watermark);
		break;
	  case 2:
		extract.ExtractWatermark1(mesh1,mesh2,Watermark);
  		break;
	  default:
		extract.ExtractWatermark1(mesh1,mesh2,Watermark);
  }
}


////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

int main(int argc, char **argv){
	 /*input the filename*/
	 printf("please input the model filename:");
	 char tempfilename[1024];scanf("%s",tempfilename);printf("\n");
	 strcat(filename,tempfilename);
	 //char *filename=NULL;filename=ParseArgs(argc,argv,filename);
	 /*generate filenames*/
     strcpy(filename1,"");strcpy(filename2,"");strcpy(filename3,"");
     strcat(filename1,filename);strcat(filename1,"watermarked.off");
	 strcat(filename2,filename);strcat(filename2,"attacked.off");
	 strcat(filename3,filename);strcat(filename3,"recovered.off");
     strcat(filename,".off");
	  /* Read original model and allocate the detected and temporary model*/
	  mesh1 = ReadOffFile(filename);//orignal model
	  mesh2=new Mesh(mesh1);//the present model and is to be rendered
	  mesh3=new Mesh(mesh1);//temporary mesh,use this mesh can reduce the time of allocating memory
	  /*load noise and watermark*/
	  noise=new float[mesh1->nverts];//allocate the noise series
	  Watermark=new unsigned char[WatermarkLength];
	  LoadWatermark(Watermark,"watermark.dat");//load the watermark
	  /*Initialize GLUT*/
	  GLUTInit(&argc, argv);
	  /*Run GLUT interface and draw the current model,never returns*/
	  GLUTMainLoop();
	  /*recharge memory*/
	  delete []noise;
	  delete []Watermark;
	  /*return*/
	  return 0;
}

////////////////////////////////////////////////////////////
// opertations manual
////////////////////////////////////////////////////////////
//E or e  mesh2=EmbedWatermark(mesh1,Watermark);  //Embed Watermark

/*attack the watermarked model*/
//N or n  mesh2=model.AddNoise(mesh2,noise);//always robust
//U or u  mesh2=model.Disturb(mesh2);//can be recovered by resort
//  mesh2=attack.ChangeOrder(mesh2,600);//can be recovered by resort
//  mesh2=attack.AddPulseNoise(mesh2,0.005);//always robust
//C or c  mesh2=model.Cut(mesh2);//cut the model  
//M or m  mesh2=model.Move(mesh2);//move the model
//S or s  mesh2=model.Scale(mesh2);
//O or o  mesh2=Rotate(mesh2,mesh3);
//D or d  mesh2=model.DownSample2(mesh2);

/*recover the attacked watermarked model*/
//R or r  mesh2=recover.resample(mesh1,mesh2);
//V or v  mesh2=recover.resort(mesh1,mesh2);
//G or g  mesh2=anneal.registration(mesh1,mesh2);
//X or x  ExtractWatermark(mesh1,mesh2,Watermark);//extract the Watermark
//Space   mesh1 = ReadOffFile(filename);//reload orignal model
//K or k  input secure key