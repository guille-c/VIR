#ifndef VIS_RECONSTRUCTOR
#define VIS_RECONSTRUCTOR

#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "funcL.h"
#include "frprmn.h"
#include "gauss.h"

#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.
#define MIN_PIX 1e-10

typedef struct{
  funcL *fL;
  int iter;
  float *p;
  char *nombreImagen, **nombresVis;
  int nVis;
  int entropia;
  //double cuantaSize;
  double Acutoff;
}Reconstructor;

Reconstructor *newReconstructor (char * nombreImagen, 
				 char **nombresVis, int nVis, int n,
				 double init_value, int init_gauss,
				 int entropia, double cuantaSize, 
				 double cutoff);
double minimizador (Reconstructor *r, double ftol);
double run (Reconstructor *r, double ftol);

void calcularNuevaPosicion(Reconstructor *r, double *x, double *y,
			   double *valor);
double calcularPromedio (MallaVoronoi *m, double x, double y);
void truncarPorAtenuacion (Reconstructor *r);
void encontrarPosicionAtenuacion (Reconstructor *r, int *i_ori, int *j_ori,
				  int nivel);
void reinicializar (Reconstructor *r, char *nombreArchivo);
void eliminarReconstructor (Reconstructor *r);

Reconstructor *leerArchivoEntrada (char *nombreArchivo);
void imprimirLog (char* s, ...);

#endif
