#include "reconstructorGenetico.h"

// ---------Funciones para pikaia-------------------------------
ReconstructorGenetico *r_aux;
float func_max;
int iter_pikaia, iter_aux;
/*extern void pikaia_(float *(*ff)(int *, float []), 
		    int *, float ctrl[12], 
		    float x[], float *, int *);*/

float ff_aux(int *n, float x[]){
  int i, iter, nSeleccion = r_aux->getNSeleccion();
  char fileout [256];
  double *parsD = new double [*n];
  float ret_aux;
  iter = r_aux->getIter();

  if (iter % 100 == 0) {
    cout << "iter: " << r_aux->getIter() << "\n";
  }

  if (iter < nSeleccion){
    if (iter >= r_aux->getIndividuos()) {
      cerr << "ERROR en ReconstructorGenetico::ff_aux: iter = " << iter 
	   << ", individuos = " << r_aux->getIndividuos() 
	   << ", nSeleccion = " << nSeleccion << "\n";
      exit (1);
    }
    for (i = 0; i < (*n); i++) {
      x[i] = r_aux->getParGeneracion (iter, i);
    }
  }
  if (iter == nSeleccion) {
    r_aux->setIter(iter_aux);
  }

  if (r_aux->getIter() % 100 == 0) {
    r_aux->imprimirLog ("reconstructor.log", 
			"        Calculando funcion para la %d iteracion con parametros:\n",
			r_aux->getIter());
    for (i = 0; i < (*n) / 3; i++) {
      r_aux->imprimirLog ("reconstructor.log", "        %g\t%g\t%g\n", 
			  parsD[3 * i], parsD[3 * i + 1], parsD[3 * i + 2]);
    }
  }
  for (i = 0; i < (*n); i++) {
    r_aux->setIndividuo(r_aux->getIter() % (int) r_aux->getIndividuos(), i, x[i]);
    parsD[i] = r_aux->getNMin(i) + (double) x[i] * (r_aux->getNMax(i) - r_aux->getNMin(i));
  }
  
  if (iter % 10000 == 0) {
    r_aux->guardarGeneracion ("generacion.dat");
  }

  r_aux->addIter();

  ret_aux = r_aux->getFunc()->f (parsD);

  r_aux->imprimirLog ("funcion.log", "%d\t%g\n", iter_pikaia, ret_aux);
  if (ret_aux > func_max) {
    r_aux->getFunc()->guardarInfo ("Gen" + i2s (iter_pikaia, 3));
    /*guardarFits (r_out->fL->imagen, r_out->fL->mask, 
      fileout, r_out->fL->nombreFits);*/
    r_aux->imprimirLog ("pikaia.log", "%d\t%g\n", iter_pikaia, ret_aux);
    r_aux->imprimirLog ("reconstructor.log", 
			"        Mejor funcion para la %d iteracion con parametros:\n",
			r_aux->getIter());
    //for (i = 0; i < (*n) / 3; i++) {
    //r_aux->imprimirLog ("reconstructor.log", "        %g\t%g\t%g\n", 
    //		  parsD[3 * i], parsD[3 * i + 1], parsD[3 * i + 2]);
    for (i = 0; i < (*n); i++) {
      r_aux->imprimirLog ("reconstructor.log", "        p[%d] = %g\n", i, parsD[i]);
    }
    r_aux->imprimirLog ("reconstructor.log", 
			"        Valor funcion: %g\n", ret_aux);
    func_max = ret_aux;
    iter_pikaia++;
  }
  
  free (parsD);
  //cout << "ret_aux = " << ret_aux << "\n";
  return ret_aux;
}

#ifdef __cplusplus
extern "C" {
#endif  
extern void pikaia_(float *(*ff_)(int *, float *), 
		    int *, float *, 
		    float *, float *, int *);

extern void rninit_(int *seed);

float ret;
float *ff_(int *n, float x[]){
  float aux;
  ret = ff_aux(n, x);
  //aux = ret;
  printf ("ret = %g\n", ret);
  return &ret;
}
#ifdef __cplusplus
}
#endif

//---------------------------------------------------------

ReconstructorGenetico::ReconstructorGenetico (Funcion *f, Inicializador *ini, double Nmin, double Nmax,
		       double individuos, double generaciones)
  : Reconstructor (f->getNPars(), ini){
  this->f = f;
  this->Nmin = new double [f->getNPars()];
  this->Nmax = new double [f->getNPars()];
  for (int i = 0; i < f->getNPars(); i++) {
    this->Nmin[i] = Nmin;
    this->Nmax[i] = Nmax;
  }
  this->individuos = individuos;
  this->generaciones = generaciones;
  this->generacion = new float * [(int) individuos];
  for (int i = 0; i < individuos; i++) {
    this->generacion[i] = new float[f->getNPars()];
  }
  iter_aux = 0;
  this->nSeleccion = 0;
}

ReconstructorGenetico::ReconstructorGenetico (Funcion *f, Inicializador *ini, double *Nmin, double *Nmax,
		       double individuos, double generaciones, int nSelec, char **nombres)
  : Reconstructor (f->getNPars(), ini){
  this->f = f;
  this->Nmin = new double [f->getNPars()];
  this->Nmax = new double [f->getNPars()];
  for (int i = 0; i < f->getNPars(); i++) {
    this->Nmin[i] = Nmin[i];
    this->Nmax[i] = Nmax[i];
  }
  this->individuos = individuos;
  this->generaciones = generaciones;
  this->generacion = new float * [(int) individuos];
  for (int i = 0; i < individuos; i++) {
    this->generacion[i] = new float[f->getNPars()];
  }
  iter_aux = 0;
  this->nSeleccion = 0;
}

double ReconstructorGenetico::run (){
  int n_pars = f->getNPars();
  float ctrl[12], *parsF = new float [n_pars], func;
  int status;
  int seed = 8;

  ctrl [0] = individuos;   // individuos (100)
  ctrl [1] = generaciones; // generaciones (500)
  ctrl [2] = -1;  // digitos significativos (6)
  ctrl [3] = -1;  // pbb de cruza (0.85)
  ctrl [4] =  5;  // modo de mutacion (2)
  ctrl [5] = -1;  // rango de mutacion inicial (0.005)
  ctrl [6] = -1;  // rango de mutacion minimo (0.0005)
  ctrl [7] = -1;  // rango de mutacion maximo (0.25)
  ctrl [8] = -1;  // diferencial relativo de ajuste (1)
  ctrl [9] =  1;  // plan de reproduccion (3)
  ctrl [10] = 1; // bandera de elitismo (0)
  ctrl [11] = -1; // imprimir output (0)
  func_max = -1e15;
  iter_pikaia = 0;

  for (int i = 0; i < n_pars; i++) {
    //parsF [i] = (float) p[i];
  }
  r_aux = this;
  
  rninit_(&seed);
  pikaia_(ff_, &n_pars, ctrl, parsF, &func, &status);

  for (int i = 0; i < n_pars; i++) {
    p [i] = Nmin [i] + (double) parsF[i] * (Nmax[i] - Nmin[i]);
  }
  f->f(p);
  f->guardarInfo ("GenFin" + i2s (iter_pikaia, 3));
  return (double) func;
}

void ReconstructorGenetico::leerSeleccion (int nSelec, char **nombres) {
}

void ReconstructorGenetico::leerGeneracion (char *nombreArchivo){
  std::ifstream archivo(nombreArchivo);
  char aux [256];

  archivo >> aux >> iter_aux;
  //archivo >> aux >> aux;

  for (int i = 0; i < individuos; i++) {
    for (int j = 0; j < f->getNPars(); j++) {
      archivo >> generacion[i][j];
    }
  }

  nSeleccion = individuos;

  archivo.close();
}

void ReconstructorGenetico::guardarGeneracion (char *nombreArchivo){
  std::ofstream archivo(nombreArchivo);

  archivo << "iter:\t" << iter << "\n";

  for (int i = 0; i < individuos; i++) {
    for (int j = 0; j < f->getNPars(); j++) {
      archivo << generacion[i][j] << "\n";
      //cout << "Guardando generacion [" << i << "][" << j << "] = " << generacion[i][j]  << "\n";
    }
    archivo << "\n";
  }
  archivo.close();
}
 
