#include "incrementadorMallaImagen.h"

IncrementadorMallaImagen::IncrementadorMallaImagen (char * nombreImagen, 
						    double init_value) {
  FILE *archivo = fopen ("chi2.dat", "w");
  fclose (archivo);
  im = do_read (nombreImagen);
  if (im == NULL) {
    cerr << "ERROR: " << nombreImagen << " no corresponde a una imagen valida\n";
    exit (1);
  }
  valor_ini = init_value;
}

/* Busca el pixel imax, jmax con mallor diferencia con la malla */

void IncrementadorMallaImagen::encontrarMayorError (struct image *im, MallaVoronoi *m, 
						    int *imax, int *jmax, double *valor){
  int i, j, nx, ny;
  double x, y, errorMax = 0, error, chi2 = 0;
  PoligonoVoronoi *pol;
  FILE *archivo = fopen ("chi2.dat", "a");
  
  nx = im->size[0];
  ny = im->size[1];
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (m, x, y);
      error = fabs (pol->valor - im->pixels[i + j * nx]);
      chi2 += error * error;
      if (error > errorMax) {
	errorMax = error;
	(*imax) = i;
	(*jmax) = j;
	(*valor) = im->pixels[i + j*nx];
      }
    }
  }
  //chi2 = chi2;
  fprintf (archivo, "%d\t%g\n", m->nPols - 3, chi2);
  fclose (archivo);
}

/* Busca el poligono que tenga mayor error en promedio y retorna en
   imax, jmax el punto de este poligono mas distinto a la malla*/

void IncrementadorMallaImagen::encontrarMayorError2 (struct image *im, MallaVoronoi *m, 
						    int *imax, int *jmax, double *valor){
  int i, j, nx, ny, *conts, id, cont = 0;
  double x, y, errorMax = 0, error, chi2 = 0, *errores;
  PoligonoVoronoi *pol;
  FILE *archivo = fopen ("chi2.dat", "a");
  
  errores = new double [m->nPols - 3];
  conts = new int [m->nPols - 3];
  nx = im->size[0];
  ny = im->size[1];
  for (i = 0; i < m->nPols - 3; i++) {
    errores [i] = 0;
    conts [i] = 0;
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (m, x, y);
      id = pol->id - 3;
      conts [id] ++;
      error = fabs (pol->valor - im->pixels[i + j * nx]);
      errores [id] += error * error;
    }
  }

  for (i = 0; i < m->nPols - 3; i++) {
    //errores [i] /= (conts [i] * conts [i]);
    if (errores [i] > errorMax) {
      errorMax = errores [i];
      id = i;
    }
  }
  
  errorMax = 0;
  (*valor) = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (m, x, y);
      error = fabs (pol->valor - im->pixels[i + j * nx]);
      chi2 += error * error;
      if (pol->id - 3 == id) {
	//(*valor) += im->pixels[i + j*nx];
	//cont++;
	if (error > errorMax) {
	  errorMax = error;
	  (*imax) = i;
	  (*jmax) = j;
	  (*valor) = im->pixels[i + j*nx];
	}
      }
    }
  }
  //(*valor) /= cont;
  //chi2 = chi2;
  fprintf (archivo, "%d\t%g\n", m->nPols - 3, chi2);
  fclose (archivo);
  delete [] errores;
  delete [] conts;
}

double IncrementadorMallaImagen::incrementar (double **p, Funcion *f) {
  MallaVoronoi *malla = newMallaVoronoi();
  int i, j, nx = im->size[0], ny = im->size[1], cont = 0;
  double valor, *p_aux, x, y;
  PoligonoVoronoi *pol;

  // Creamos la malla antigua
  for (i = 0; i < f->getNPars() / 3; i++) {
    //cout << "Insertando {" << (*p)[3 * i] << ", " << (*p)[3 * i + 1] << ", " << (*p)[3 * i + 2] << "}\n";
    insertarSitio (malla, &((*p)[3 * i]), &((*p)[3 * i + 1]), (*p)[3 * i + 2]);
  }

  p_aux = new double [f->getNPars() + 3];

  for (i = 0; i < f->getNPars(); i++) {
    p_aux[i] = (*p)[i];
  }

  if (f->getNPars() == 0) {
    encontrarMayorError (im, malla, &i, &j, &valor);
  }
  else {
    encontrarMayorError2 (im, malla, &i, &j, &valor);
  }
  cout << "(i, j) = (" << i << ", " << j << ") = " << valor << "\n";
  p_aux [f->getNPars()]     = i / (nx - 1.0);
  p_aux [f->getNPars() + 1] = j / (ny - 1.0);
  p_aux [f->getNPars() + 2] = valor;

  //ActualizarIntensidades (p_aux, f->getNPars() / 3 + 1); 
  delete [] (*p);
  (*p) = p_aux;

  /*/calculamos el nuevo valor
  x = i / (nx - 1.0);
  y = j / (ny - 1.0);
  i = f->getNPars();
  //insertarSitio (malla, &((*p)[i]), &((*p)[i + 1]), (*p)[i + 2]);
  insertarSitio (malla, &x, &y, (*p)[i + 2]);
  
  pol = encontrarPoligono (malla, x, y);
  int id = pol->id;
  //int id = ((PoligonoVoronoi *)malla->poligonos->primero->info)->id;
  cont = 0;
  valor = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (malla, x, y);
      if (pol->id  == id) {
	valor += im->pixels[i + j*nx];
	cont++;
      }
    }
  }
  valor /= cont;
  cout << "-----------------------------valor = " << valor << "\n";
  (*p) [f->getNPars() + 2] = valor;
  */

  eliminarMallaVoronoi(malla);
  f->setNPars(f->getNPars() + 3);
  //ActualizarIntensidades ((*p), f->getNPars() / 3); 
  return f->f(*p);
  return 0;
}

void IncrementadorMallaImagen::ActualizarIntensidades (double *p, int nPols) {
  MallaVoronoi *malla = newMallaVoronoi();
  int i, j, nx = im->size[0], ny = im->size[1], *conts;
  double valor, *valores, x, y;
  PoligonoVoronoi *pol;

  // Creamos la malla antigua
  for (i = 0; i < nPols; i++) {
    insertarSitio (malla, &(p[3 * i]), &(p[3 * i + 1]), p[3 * i + 2]);
  }

  valores = new double [nPols];
  conts = new int [nPols];
  for (i = 0; i < nPols; i++) {
    valores[i] = 0;
    conts[i] = 0;
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (malla, x, y);
      valores [pol->id - 3]  += im->pixels[i + j*nx];
      conts [pol->id - 3] ++;
      }
    }
  
  eliminarMallaVoronoi(malla);  
  for (i = 0; i < nPols; i++) {
    if (conts[i] != 0) {
      p[3*i + 2] = valores[i] / conts[i];
    }
    else {
      p[3*i + 2] = 0;
    }
  }
  delete [] valores;
  delete [] conts;
}

IncrementadorMallaImagen::~IncrementadorMallaImagen() {
  delete_map (im);
}
