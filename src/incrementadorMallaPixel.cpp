#include "incrementadorMallaPixel.h"

IncrementadorMallaPixel::IncrementadorMallaPixel (int nx, int ny, double valor_ini) {
  this->nx = nx;
  this->ny = ny;
  this->valor_ini = valor_ini;
}

double IncrementadorMallaPixel::incrementar (double **p, Funcion *f) {
  int i, j;
  double x, y, *p_aux, *p_mejor, func, func_mejor;
  MallaVoronoi *malla = newMallaVoronoi();
  
  // Creamos la malla antigua
  for (i = 0; i < f->getNPars() / 3; i++) {
    insertarSitio (malla, &((*p)[3 * i]), &((*p)[3 * i + 1]), (*p)[3 * i + 2]);
  }

  p_aux = new double [f->getNPars() + 3];
  p_mejor = new double [f->getNPars() + 3];
  func_mejor = f->f(*p);

  for (i = 0; i < f->getNPars(); i++) {
    p_aux[i] = (*p)[i];
  }

  f->setNPars (f->getNPars() + 3);

  for (i = 0; i < f->getNPars(); i++) {
    p_mejor[i] = 0;
  }

  for (i = 0; i < nx; i++) {
    x = i / (nx - 1.0);
    for (j = 0; j < ny; j++) {
      y = j / (ny - 1.0);
      p_aux [f->getNPars() - 3] = x;
      p_aux [f->getNPars() - 2] = y;
      if (valor_ini >= 0) {
	p_aux [f->getNPars() - 1] = valor_ini;
      }
      else{
	PoligonoVoronoi *pol = encontrarPoligono (malla, x, y);
	p_aux [f->getNPars() - 1] = pol->valor;
      }
      cout << "f->npars = " << f->getNPars() << "\n";

      func = f->f(p_aux);
      if (func < func_mejor) { // Minimiza
	for (int k = 0; k < f->getNPars(); k++) {
	  p_mejor[k] = p_aux[k];
	  func_mejor = func;
	}
      }
    }
  }
  delete [] p_aux;

  for (i = 0; i < f->getNPars(); i++) {
    if (p_aux[i] != 0) {
      delete [] (*p);
      (*p) = p_aux;
      return func_mejor;
    }
  }
  
  delete [] p_aux;
  f->setNPars (f->getNPars() - 3);
  return func_mejor;
}
