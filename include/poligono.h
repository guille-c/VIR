#ifndef VIS_POLIGONO
#define VIS_POLIGONO

#include "estructuras_malla.h"
#include "arista.h"
#include "punto.h"

#ifdef __cplusplus
extern "C" {
#endif  

#define TRUE 1
#define FALSE 0

#define BORDE_INF 1
#define BORDE_DER 2
#define BORDE_SUP 3
#define BORDE_IZQ 4

struct poligonoVoronoi {
  double x, y, valor;
  AristaVoronoi *a;
  NodoLista *n;
  int id;
};

PoligonoVoronoi *newPoligonoVoronoi(double x, double y, double valor, int id);
int dentro (PoligonoVoronoi *pol, double x, double y);
double area (PoligonoVoronoi *pol);
double areaCuadrado (PoligonoVoronoi *pol);
//private int bordeCuadrado (double u1, double v1, double u2, double v2,
//			   double *x, double *y);
//int bordeCuadrado (PuntoVoronoi *pini, PuntoVoronoi *pfin,
//		   PuntoVoronoi *qini, PuntoVoronoi *qfin);
int bordeCuadrado (PuntoVoronoi *pini, PuntoVoronoi *pfin, 
		   PuntoVoronoi *qini, PuntoVoronoi *qfin);
int dentroCuadrado(AristaVoronoi *a);
int tope (PuntoVoronoi *p);
double sumarBordes(AristaVoronoi **a, 
		   PuntoVoronoi *p1, PuntoVoronoi *p2, 
		   PuntoVoronoi *q1, PuntoVoronoi *q2,
		   PoligonoVoronoi *pol);

#ifdef __cplusplus
}
#endif  

#endif
