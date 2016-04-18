#ifndef VIS_MALLA_VORONOI
#define VIS_MALLA_VORONOI

#include "estructuras_malla.h"
#include "arista.h"
#include "lista.h"

#ifdef __cplusplus
extern "C" {
#endif  
  
#define TRUE 1
#define FALSE 0
#define PRECISION_DOUBLE 1e-15
  
  struct mallaVoronoi{
    Lista *aristas;
    Lista *poligonos;
    Lista *puntos;
    int nPols;
  };
  
  MallaVoronoi *newMallaVoronoi();
  PuntoVoronoi *circuncentro(PoligonoVoronoi *p1,
			     PoligonoVoronoi *p2,
			     PoligonoVoronoi *p3);
  void eliminarMallaVoronoi(MallaVoronoi *m);
  
  PoligonoVoronoi *encontrarPoligono(MallaVoronoi *m, 
				     double x, double y);
  PoligonoVoronoi *encontrarPoligono2(PoligonoVoronoi *polIni, 
				      double x, double y);
  void insertarSitio (MallaVoronoi *m, double *x, double *y, double valor);
  double H (PoligonoVoronoi *pi, PoligonoVoronoi *pj,
	    PoligonoVoronoi *pk, PoligonoVoronoi *p);
  int buscarPuntos (PuntoVoronoi *q, PuntoVoronoi *qPadre,
		    PoligonoVoronoi *pol, Lista *TPuntos,
		    Lista *aristasDentro, Lista *aristasBorde);
  void crearVertices (MallaVoronoi *m, Lista *aristasBorde, 
		      PoligonoVoronoi *pol, PoligonoVoronoi *pol1);
  AristaVoronoi *siguienteArista(AristaVoronoi *a1, 
				 PuntoVoronoi *q,
				 Lista *aristasBorde);
  void cambiarPoligonosAristasAEliminar (Lista *aristas);
  void cambiarAristaPoligono (PoligonoVoronoi *pol, 
			      Lista *aristas);
  PuntoVoronoi *puntoDentroPoligono (AristaVoronoi *a,
				     PoligonoVoronoi *pol);
  
  void imprimirAristasMalla (MallaVoronoi *m);
  void imprimirAristasMallaArchivo (MallaVoronoi *m, char *nombre,
				    int nx, int ny);
  void imprimirMallaArchivo (MallaVoronoi *m, char *nombre);
  
  int validarMallaPoligonosAristas (MallaVoronoi *m);
  int checkAristas (MallaVoronoi *m, double tol);
  void ajustarPoligonoCuadrado (double *x, double *y);
  //int round (double x);
  
#ifdef __cplusplus
}
#endif

#endif
