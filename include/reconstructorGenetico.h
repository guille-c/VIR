#ifndef VIRPP_RECONSTRUCTORGENETICO
#define VIRPP_RECONSTRUCTORGENETICO

#include "reconstructor.h"
#include "funcion.h"
#include <string.h>
#include <iostream.h>

class ReconstructorGenetico: public Reconstructor {
 protected:
  Funcion *f;
  double individuos, generaciones;
  double *Nmax, *Nmin, **seleccion;
  float **generacion;
  int nSeleccion;

 public:
  ReconstructorGenetico (Funcion *f, Inicializador *ini, double Nmin, double Nmax,
			 double individuos, double generaciones);
  ReconstructorGenetico (Funcion *f, Inicializador *ini, double *Nmin, double *Nmax,
			 double individuos, double generaciones, 
			 int nSelec = 0, char **nombres = NULL);
  virtual double run ();
  Funcion * getFunc() {return f;}

  double getNMax(int i) {return Nmax[i];}
  double getNMin(int i) {return Nmin[i];}
  double getIndividuos() {return individuos;}
  double getParGeneracion (int indiv, int par) {return generacion[indiv][par];}
  int getNSeleccion () {return nSeleccion;}

  void setIndividuo(int indiv, int i, float valor) {generacion[indiv][i] = valor;}
  void setIter (int i) {iter = i;}

  void leerSeleccion (int nSelec, char **nombres);
  void leerGeneracion (char *nombreArchivo);
  void guardarGeneracion (char *nombreArchivo);
};

#endif
