#include "incrementadorMallaAristas.h"

IncrementadorMallaAristas::IncrementadorMallaAristas (double valor_ini) {
  this->valor_ini = valor_ini;
}

double IncrementadorMallaAristas::incrementar (double **p, Funcion *f) {
  int i, iter = 0;
  double *p_aux, func, x, y, valor;
  MallaVoronoi *malla = newMallaVoronoi();
  
  // Creamos la malla antigua
  for (i = 0; i < f->getNPars() / 3; i++) {
    //cout << "insertando {" << (*p)[3 * i] << ", " << (*p)[3 * i + 1] << ", "<< (*p)[3 * i + 2] << "}\n ";
    insertarSitio (malla, &((*p)[3 * i]), &((*p)[3 * i + 1]), (*p)[3 * i + 2]);
  }

  do {
    if (iter > 2000) {
      cerr << "ERROR en IncrementadorMallaAristas::incrementar, no se encuentra nueva posicion\n";
      cerr << "      Los poligonos de la malla tendran la misma intensidad?\n";
      exit (1);
    }
    iter++;
    calcularNuevaPosicion (malla, &x, &y, &valor);
    //cout << "(x, y) = (" << x << ", " << y << ")\n";
  } while (x < 0 || x > 1 || y < 0 || y > 1);

  p_aux = new double [f->getNPars() + 3];
  for (i = 0; i < f->getNPars(); i++) {
    p_aux[i] = (*p)[i];
  }
  p_aux[f->getNPars()] = x;
  p_aux[f->getNPars() + 1] = y;
  if (valor_ini >= 0) {
    p_aux[f->getNPars() + 2] = valor_ini;
  }
  else {
    p_aux[f->getNPars() + 2] = valor;
  }
  delete [] (*p);
  f->setNPars (f->getNPars() + 3);
  (*p) = new double [f->getNPars()];
  for (i = 0; i < f->getNPars(); i++) {
    (*p)[i] = p_aux[i];
  }
  delete [] p_aux;

  eliminarMallaVoronoi(malla);

  //cout << "\nRetorno:\n";
  //for (i = 0; i < f->getNPars() / 3; i++) {
  //cout << "{" << (*p)[3 * i] << ", " << (*p)[3 * i + 1] << ", "<< (*p)[3 * i + 2] << "}\n ";
  //}
  double ret = f->f(*p);
  //for (i = 0; i < f->getNPars() / 3; i++) {
  //cout << "Luego de ret {" << (*p)[3 * i] << ", " << (*p)[3 * i + 1] << ", "<< (*p)[3 * i + 2] << "}\n ";
  //}

  return ret;
}

void IncrementadorMallaAristas::calcularNuevaPosicion (MallaVoronoi *m, double *x,
						       double *y, double *valor) {
  double dITotal = 0, I;
  NodoLista *n;
  AristaVoronoi *a;
  
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    //imprimirArista (a);
    if (a->poliDer != NULL && a->poliIzq != NULL &&
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      dITotal += fabs (a->poliDer->valor - a->poliIzq->valor);
      //cout << "dITotal = " << dITotal << "\n";
    }
  }

  I = dITotal * rand () / RAND_MAX;
  //cout << "dITotal = " << dITotal << ", I = " << I << "\n";
  dITotal = 0;
  for (n = m->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    if (a->poliDer != NULL && a->poliIzq != NULL &&
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      dITotal += fabs (a->poliDer->valor - a->poliIzq->valor);
      
      if (dITotal > I) {
	(*x)     = (a->ptoIni->x + a->ptoFin->x) / 2;
	(*y)     = (a->ptoIni->y + a->ptoFin->y) / 2;
	(*valor) = (a->poliDer->valor + a->poliIzq->valor) / 2;
	return;
      }
    }
  }
  (*x)     = (a->ptoIni->x + a->ptoFin->x) / 2;
  (*y)     = (a->ptoIni->y + a->ptoFin->y) / 2;
  (*valor) = (a->poliDer->valor + a->poliIzq->valor) / 2;
  cout << "calcularNuevaPosicion: Se llego al final.\n";
}
