#include "voronoiMesh2.h"
#include <limits.h>
#include <math.h>

struct voronoiMesh *createMesh ()
{
  struct voronoiMesh *mesh;
  mesh = (struct voronoiMesh*) malloc (sizeof(struct voronoiMesh));
  mesh->npols = 0;
  return mesh;
}

void submkim(float *xcel, float *ycel, float *icel, int ncel, int nx, int ny, float im[nx][ny]) 
{
  struct voronoiMesh *mesh = createMesh();
  double x[2], value;
  int i; 

  int j;
  //  float *pixs = (float *) malloc (nx*ny*sizeof(float));

  /*  printf("TESTING submkim \n");
  printf("%d %d  %d \n",ncel,nx,ny);
      
  printf("test im assign 0 \n");
  im[1][1] = -1; */


  /*  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
	  printf("i %d j %d im %g \n",i,j,im[i][j]);
	}
	} */
   
  for (i = 0; i < ncel; i++)
    {
      x[0] = xcel[i];
      x[1] = ycel[i];
      // printf("CHECK1: %d %g %g %g \n",i,xcel[i],ycel[i],icel[i]);
      addPolygon (mesh, x, icel[i]);
      //      printf("CHECK1 - 2\n");
    }

  // printf("ANCHOR 1\n");

  toImage2 (mesh, nx, ny, im);
  
  //  printf("ANCHOR 2\n");
  
  /*  for (i = 0; i < nx; i++)
    {
      printf("CHECK2: %g \n",im[i][2]);
      } */ 
  
}



void addPolygon (struct voronoiMesh *mesh, double x[2], double value)
{
  int i;
  struct voronoiPolygon *pol_aux;
  for (i = 0; i < mesh->npols; i++)
    {
      if ((mesh->polygons[i].x[0] == x[0]) 
	  && (mesh->polygons[i].x[1] == x[1]))
      {
	mesh->polygons[i].value = value;
	return;
      }
    }
  if (mesh->npols < NV)
  {
    mesh->npols++;
  }
  mesh->polygons[mesh->npols - 1].value = value;
  mesh->polygons[mesh->npols - 1].x[0] = x[0];
  mesh->polygons[mesh->npols - 1].x[1] = x[1];
}

void printMesh (struct voronoiMesh *mesh)
{
  int i;
  for (i = 0; i < mesh->npols; i++)
    {
      printf("(%g, %g) = %g\n", mesh->polygons[i].x[0],
	     mesh->polygons[i].x[1], mesh->polygons[i].value);
    }
}

double getValue(struct voronoiMesh *mesh, int x, int y)
{
  int i;
  double d = INT_MAX, value;

  for (i = 0; i < mesh->npols; i++)
    {
      double d2;
      d2  = pow((x - mesh->polygons[i].x[0]), 2) + pow(y - mesh->polygons[i].x[1], 2);
      if (d2 < d)
	{
	  d = d2;
	  value = mesh->polygons[i].value;
	}
    }
  //  printf ("out of getValue \n");
  return value;
}

/* void toImage(struct voronoiMesh *mesh, int nx, int ny, float im[nx][ny])
{
  int x, y, k = 0;
  
  for (y = 0; y < ny; y++)
    {
      for (x = 0; x < nx; x++)
	{
	  im[k] = (float) getValue(mesh, x, y);
	  k++;
	}
    }
} */

void toImage2(struct voronoiMesh *mesh, int nx, int ny, float im[nx][ny])
{
  int x, y;
  float value;
  
  for (y = 0; y < ny; y++)
    {
      for (x = 0; x < nx; x++)
	{
	  value = (float) getValue(mesh, x, y);
	  //	  printf("x = %d y = %d valor = %g\n",x,y, value);
	  // printf("test im assign 1 \n");
	  // im[1][1] = -1;
	  // printf("test im assign 2 \n");
	  // im[x][y] = -1;
	  // printf("test im assign 3 \n");
	  im[x][y] = (float) getValue(mesh, x, y);
	}
    }
}
