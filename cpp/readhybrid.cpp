#include <iostream>     // for reading the file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
int readhybrid(int iter, float *lx,float *ly,float *lz,float *bx,float *by,float *bz,float *ex,float *ey,float *ez, float *x0hc, float *b0hc, float *v0hc) {
  FILE * pFile;                                 // create file handler
  size_t result;

  //! hybrid code output file

  pFile = fopen ("./data/uniform_grid/Ganymede_uniform_grid_TL28500.dat", "rb");
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  //! in the following lines, each parameter of the hybrid code output needs
  //! to be read with the exact size and variable type; for this, fread is
  //! used, which allows to set the number of items to be read with the
  //! corresponding size

  //! Run name, string of 42 chars  
  char runname[42];
  result = fread(runname,sizeof(char),42,pFile);
  if (result == 0)
    cout << "WARNING: No filename read!" << endl;

  //! Number of meshpoints
  short meshpoints[3];
  result = fread(meshpoints,sizeof(short),3,pFile);
  if (result == 0)
    cout << "WARNING: Zero meshpoints read!" << endl;
  if (iter == 1)
    return meshpoints[0];

  cout << "Mesh x: " << meshpoints[0] << endl
       << "Mesh y: " << meshpoints[1] << endl
       << "Mesy z: " << meshpoints[2] << endl << endl;
  
  //! Time level
  int timelevel[1];                     
  result = fread(timelevel,sizeof(int),1,pFile);
  if (result == 0)
    cout << "WARNING: no TL read!" << endl;

  //! Moon origin (3 floats)
  float origin[3];                      
  result = fread(origin,sizeof(float),3,pFile);
  if (result == 0)
    cout << "WARNING: no origin read!" << endl;
 
  //! Length of the box in normalized units (1 float)
  float length[3];
  result = fread(length,sizeof(float),3,pFile);
  if (result == 0)
    cout << "WARNING: no boxlength read!" << endl;

  //! Radius of moon in normalized units
  float radiushc[1];
  result = fread(radiushc,sizeof(float),1,pFile);
  if (result == 0)
    cout << "WARNING: no radius read!" << endl;
  
  //! Ion inertial length
  float x0[1];
  result = fread(x0,sizeof(float),1,pFile);
  if (result == 0)
    cout << "WARNING: x0 not read!" << endl;
  *x0hc = x0[0];
  
  //! Inertial (Alfven) velocity
  float v0[1];
  result = fread(v0,sizeof(float),1,pFile);
  if (result == 0)
    cout << "WARNING: vA not read!" << endl;
  *v0hc = v0[0];

  //! Magnetic field factor
  float B0[1];
  result = fread(B0,sizeof(float),1,pFile);
  if (result == 0)
    cout << "WARNING: B0 not read!" << endl;
  *b0hc = B0[0];

  //! Density factor
  float n0[1];
  result = fread(n0,sizeof(float),1,pFile);
  if (result == 0)
    cout << "WARNING: rho0 not read!" << endl;

  //! Total number of points  
  const int npoints = meshpoints[0] * meshpoints[1] * meshpoints[2];
  cout << "Total number of points (mesh^3): " << npoints << endl;

  //! Now, memory has to be allocated for arrays containing the grid positions;
  //! assigning arrays of corresponding size like e.g. float grid[npoints]
  //! generates a segmentation fault error since the memory is allocated in the
  //! stack, and the stack is overflown; by using malloc(), the memory is
  //! allocated on the heap

  float * gridx;
  float * gridy;
  float * gridz;

  gridx = (float*) malloc(npoints*sizeof(float));
  gridy = (float*) malloc(npoints*sizeof(float));
  gridz = (float*) malloc(npoints*sizeof(float));

  //! Check if any of the grids were not read
  result = fread(gridx,sizeof(float),npoints,pFile);
  if (result == 0)
    cout << "WARNING: no x-grid information read!" << endl;
  result = fread(gridy,sizeof(float),npoints,pFile);
  if (result == 0)
    cout << "WARNING: no x-grid information read!" << endl;
  result = fread(gridz,sizeof(float),npoints,pFile);
  if (result == 0)
    cout << "WARNING: no x-grid information read!" << endl;

  //! Number of components
  //! (Not used, just to keep reading the file...)
  short components[1];
  result = fread(components,sizeof(short),1,pFile);
  if (result == 0)
    cout << "WARNING: no components(!) read!" << endl;

  //! Type of variable
  //! (Not used)
  int type[1];
  result = fread(type,sizeof(int),1,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  //! Variable name
  //! (Not used)
  char name[50];
  result = fread(name,sizeof(char),50,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  
  //! Similar as above with the grid positions, memory has to be allocated for
  //! the magnetic field data; additionally, the reading flow requires to read
  //! the number of components, the type and the name
  float * bxhc;
  bxhc = (float*) malloc(npoints*sizeof(float));
  result = fread(bxhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  float * byhc;
  byhc = (float*) malloc(npoints*sizeof(float));
  result = fread(byhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  float * bzhc;
  bzhc = (float*) malloc(npoints*sizeof(float));
  result = fread(bzhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  
  result = fread(components,sizeof(short),1,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  
  result = fread(type,sizeof(int),1,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;


  result = fread(name,sizeof(char),50,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  
  //! Again with the electric field
  float * exhc;
  exhc = (float*) malloc(npoints*sizeof(float));
  result = fread(exhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  float * eyhc;
  eyhc = (float*) malloc(npoints*sizeof(float));
  result = fread(eyhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  float * ezhc;
  ezhc = (float*) malloc(npoints*sizeof(float));
  result = fread(ezhc,sizeof(float),npoints,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;
  
  result = fread(components,sizeof(short),1,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;
  
  result = fread(type,sizeof(int),1,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;


  result = fread(name,sizeof(char),50,pFile);
  if (result == 0)
   cout << "Warning: zero lines read!" << endl;

  //! And again with the velocities
  float * vxhc;
  vxhc = (float*) malloc(npoints*sizeof(float));
  result = fread(vxhc,sizeof(float),npoints,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (vxhc)\n");

  float * vyhc;
  vyhc = (float*) malloc(npoints*sizeof(float));
  result = fread(vyhc,sizeof(float),npoints,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (vyhc)\n");

  float * vzhc;
  vzhc = (float*) malloc(npoints*sizeof(float));
  result = fread(vzhc,sizeof(float),npoints,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (vzhc)\n");

  
  result = fread(components,sizeof(short),1,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (rhohc components)\n");

  
  result = fread(type,sizeof(int),1,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (rhohc type)\n");


  result = fread(name,sizeof(char),50,pFile);
  if (result == 0)
    printf("Warning... zero lines read! (rhohc name)\n");

  //! Normalized density  
  float * rhohc;
  rhohc = (float*) malloc(npoints*sizeof(float));
  result = fread(rhohc,sizeof(float),npoints,pFile);
  if (result == 0)
    cout << "WARNING: No rho values read!" << endl;
  
  fclose(pFile);
  
  //! Next, copy the contents of the read arrays to the output ones
  memcpy(lx,gridx,npoints*sizeof(float));
  memcpy(ly,gridy,npoints*sizeof(float));
  memcpy(lz,gridz,npoints*sizeof(float));
  memcpy(bx,bxhc,npoints*sizeof(float));
  memcpy(by,byhc,npoints*sizeof(float));
  memcpy(bz,bzhc,npoints*sizeof(float));
  memcpy(ex,exhc,npoints*sizeof(float));
  memcpy(ey,eyhc,npoints*sizeof(float));
  memcpy(ez,ezhc,npoints*sizeof(float));
  
  //! Output useful constants:

  cout << endl
       << "B-field normalization: " << B0[0] << endl
       << "E-field normalization: " << v0[0] << endl
       << "Ion inertial length:   " << x0[0] << endl;

  //! Free the memory used
  free(gridx);
  free(gridy);
  free(gridz);
  free(bxhc);
  free(byhc);
  free(bzhc);
  free(exhc);
  free(eyhc);
  free(ezhc);
  free(vxhc);
  free(vyhc);
  free(vzhc);
  
  return 0;
  
}
