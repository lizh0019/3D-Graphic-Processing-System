OVERVIEW
========

offviewer is a platform for watermarking 3D models
described in the Geometry Center's OFF file format.  
Information about OFF can be found at the following Web sites: 
http://astronomy.swin.edu.au/~pbourke/geomformats/oogl/#OFF
http://www.dcs.ed.ac.uk/home/mxr/gfx/3d/OOGL.spec
http://www.neuro.sfc.keio.ac.jp/~aly/polygon/format/off.html

offviewer provides code for: 
(1) using GLUT and OpenGL to build an interactive viewing 
program and (2) watermarking OFF models.


The program depends on the GLUT user interface library. 
If your system does not already have GLUT installed, you can download it at:
http://www.opengl.org/developers/documentation/glut/index.html



EXECUTION
=========

Open a C++ compiler and execute it, or double-click the offviewer.exe directly.
////////////////////////////////////////////////////////////
// opertations manual
////////////////////////////////////////////////////////////
At the beginning, please input the model file name,e.g. snake 
then it pops up a viewing window.
Dragging the left/middle/right mouse in the window rotates/scales/translates the view.  Hitting the ESC key
exits the program.
//E or e  mesh2=EmbedWatermark(mesh1,Watermark);  //render the watermark
method 0 is the algorithm proposed in "A robust watermarking scheme for 3D triangular mesh models"
method 1 is the algorithm proposed in my thesis "A Geometry-based Robust Watermarking Scheme for 3D Model"
method 2 is the algorithm we are investigating in.
/*impose attacks*/
Type N or n  mesh2=attack.AddNoise(mesh2,noise);//add Gaussian noise
Type U or u  mesh2=attack.Disturb(mesh2);//scramble the vertices
Type C or c  mesh2=attack.Cut(mesh2);//insection
Type M or m  mesh2=attack.Move(mesh2);//move the model
Type S or s  mesh2=attack.Scale(mesh2);// uniform scaling
Type O or o  mesh2=Rotate(mesh2,mesh3);//rotation
Type D or d  mesh2=attack.DownSample2(mesh2);//simplification

/*recover the attacked model*/
Type R or r  mesh2=recover.resample(mesh1,mesh2);// resample the model against insection or downsampling
Type V or v  mesh2=recover.resort(mesh1,mesh2);//resort the vertices
Type G or g  mesh2=anneal.registration(mesh1,mesh2);//registration against scaling or rotation(including traslation)
Type X or x  ExtractWatermark(mesh1,mesh2,Watermark);//Extract the watermark
Type Space   mesh2=new Mesh(mesh1);//reload orignal model
Type K or k  input key(integer)

Note: Type the keyboard when the current window is the Opengl one titled "Current 3D model".
////////////////////////////////////////////////////////////
// EXAMPLE
////////////////////////////////////////////////////////////
Input commands as follows:
car
e
1
n
y
0.5
c
y
0.2
z
d
y
0.1
r
y
x
k
88
x
Space
e
0
o
y
25,50,80,1,2,3
g
y
10
21,56,89,0.9,2.1,2.8
x



by the above operations
1)the 3D model "car.off" is loaded and rendered.
2)the model is watermarked via algorithm 1
3)the model is added on gaussian noise whose amplitude is 0.5% of the maximum magnitude of the vector that links the model centroid to vertices
4)20% of the vertices are removed from the model along the z axis
5)10% of the vertices are downsampled
6)the model is resampled using original 3D model escrow
7)watermark is extracted using stego-key
8)stego-key is changed to 88
9)watermark is extracted using the new stego-key
10)the model is reloaded
11)watermark the model via algorithm 0
12)rotate the model with 25,50,80 degrees around x,y,z axis respectively and translate the model with the vector (1,2,3)
13)registrate the model ,in which the proximate initial parameters is recommanded to be given, and the relative error is less than 20% to ensure the convergence,
however,other initial parameters may result in correct recovered model also;
The anneal times of the registration is 10.The more the times are ,the higher accuracy is expected .
15)watermark is extracted and the correlation value is caculated

note: the model scene is refreshed whenever you modify the model


DISCLAIMER
==========

This program and source code is provided as-is for non-commercial use.
It is not supported at all. The author is not responsible for any 
bugs/errors/etc. or harm from its use.

