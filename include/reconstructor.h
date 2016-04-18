#ifndef VIRPP_RECONSTRUCTOR
#define VIRPP_RECONSTRUCTOR

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <stdarg.h>
#include <string>
#include "frprmn.h"
#include "gauss.h"
#include "inicializador.h"
using namespace std;

//#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.
#define MIN_PIX 1e-10

string i2s (int numero, int largo);

class Reconstructor {
 protected:
  int iter, n_pars;
  double *p;
  Inicializador *ini;

 public:
  Reconstructor ();
  Reconstructor (int n_pars, Inicializador *ini, double xini = 0.0);
  virtual double run () = 0;
  void reinicializar (char *nombreArchivo, int n_iter);
  void reinicializar (char *nombreArchivo);

  double getPar (int i) {
    if (i >= n_pars) {
      cerr << "ERROR: en getPar " << i << " >= " << n_pars << "\n";
      exit (1);
    }
    return p[i];
  }
  double *getPars () {return p;}
  int getNPars () {return n_pars;}
  virtual void setNPars (int n_pars);
  void setPars (double *pars);
  void setPar (int i, double valor);
  int getIter() {return iter;}
  void addIter() {iter++;};

  virtual void guardarInfo (string info) {}

  virtual ~Reconstructor();

  Reconstructor& operator= (const Reconstructor&);

  void imprimirLog (char *s, ...);
  void imprimirLog (char * nombreArchivo, char *s, ...);
  
};

#endif
