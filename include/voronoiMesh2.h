/****************************************
* Representation of a 2D Voronoi        *
* mesh.	                              	*
*                                       *
* Created by Guillermo Cabrera 14/10/04 *
****************************************/
//#include "image_routines.h"

#define NV 1000

struct voronoiPolygon
{
  double x[2];    // Position of the point (x, y) that defines the polygon.
  double value;   // Value of the polygon.
};

struct voronoiMesh
{
  struct voronoiPolygon polygons[NV];	// Array with the polygons of the mesh.
  int npols;		// Number of polygons.
};

struct voronoiMesh *createMesh ();
void addPolygon (struct voronoiMesh *mesh, double x[2], double value);
void printMesh (struct voronoiMesh *mesh);
double getValue(struct voronoiMesh *mesh, int x, int y);
//void toImage(struct voronoiMesh *mesh, int nx, int ny, float im[nx][ny]);
void toImage2(struct voronoiMesh *mesh, int nx, int ny, float im[nx][ny]);

void submkim(float *xcel, float *ycel, float *icel, int ncel, int nx, int ny, float im[nx][ny]) ;
