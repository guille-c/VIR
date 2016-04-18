#ifndef VIS_RASTERIZADOR
#define VIS_RASTERIZADOR

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mallaVoronoi.h"

#ifdef __cplusplus
extern "C" {
#endif  

#define NULL_PIX -1e300

double **toImage(MallaVoronoi *m, int nx, int ny, int **mask);
void linea1 (double **im, int nx, int ny,
	     PuntoVoronoi *p1, 
	     PuntoVoronoi *p2,
	     PoligonoVoronoi *polD,
	     PoligonoVoronoi *polI, 
	     int **mask);
void linea2 (double **im, int nx, int ny,
	     PuntoVoronoi *p1, 
	     PuntoVoronoi *p2,
	     PoligonoVoronoi *polD,
	     PoligonoVoronoi *polI, 
	     int **mask);
void linea3 (double **im, int nx, int ny,
	     PuntoVoronoi *p1, 
	     PuntoVoronoi *p2,
	     PoligonoVoronoi *polD,
	     PoligonoVoronoi *polI,
	     int **mask);
void linea4 (double **im, int nx, int ny,
	     PuntoVoronoi *p1, 
	     PuntoVoronoi *p2,
	     PoligonoVoronoi *polD,
	     PoligonoVoronoi *polI,
	     int **mask);
int interseccionCuadrado (PuntoVoronoi *pini, PuntoVoronoi *pfin,
			  PuntoVoronoi *qini, PuntoVoronoi *qfin,
			  int tipoLinea);
void raster (double **im, int **mask, int nx, int ny, MallaVoronoi *m);
int bordeV (double **im, int **mask, int i, int ny);
int bordeH (double **im, int **mask, int nx, int j);

//void submkim(float *xcel, float *ycel, float *icel, int ncel, 
//int nx, int ny, float im[nx][ny]) ;
double **toImageSinRaster(MallaVoronoi *malla, int nx, int ny, int **mask);
int comprobarImagen (MallaVoronoi *m, int **im, int nx, int ny);

#ifdef __cplusplus
}
#endif  

#endif
