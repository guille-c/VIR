#include "funcionBayesVoronoiCBII.h"

FuncionBayesVoronoiCBII::FuncionBayesVoronoiCBII (char * nombreImagen, char **nombresVis, 
						  int nVis, int n,
						  double init_value, int init_gauss, 
						  int entropia, double cuantaSize,
						  double expanded_lx, double expanded_ly,
						  double beamSize) :
  FuncionCBI (nombreImagen, nombresVis, nVis, n) {
  this->n_pars = n;

  pos = new double [this->n_pars * 2];
  for (int i = 0; i < this->n_pars * 2; i++) {
    pos[i] = rand() / (RAND_MAX + 1.0);
    //cout << i << ") " << pos[i] << "\n";
  }
  fL = newFuncL(nombreImagen, nombresVis, nVis, n, entropia, cuantaSize, expanded_lx, expanded_ly, beamSize);
  
  this->malla = fL->malla;
}

FuncionBayesVoronoiCBII::FuncionBayesVoronoiCBII (char* nombre_archivo, Inicializador *ini) {
  std::ifstream archivo (nombre_archivo);
  char *nombreImagen, **nombresVisibilidades, *dato;
  int nPols, entropia, nVis, i, init_gauss;
  double init_value, cuantaSize, cutoff, expanded_lx, expanded_ly;

  
  if (archivo == NULL) {
    cerr <<  "ERROR al intentar leer " << nombre_archivo << "\n";
    exit (1);
  }
  cout << nombre_archivo << " abierto\n";
  
  nombreImagen = new char [512];
  dato = new char [512];
  
  archivo >> dato >> nombreImagen;
  cout << dato << "\t" << nombreImagen << "\n";
  archivo >> dato >> nPols;
  cout << dato << "\t" <<  nPols << "\n";
  archivo >> dato >> entropia;
  cout << dato << "\t" << entropia << "\n";
  archivo >> dato >> init_value;
  cout << dato << "\t" << init_value << "\n";
  archivo >> dato >> init_gauss;
  cout << dato << "\t" << init_gauss << "\n";
  archivo >> dato >> cuantaSize;
  cout << dato << "\t" << cuantaSize << "\n";
  archivo >> dato >> cutoff;
  cout << dato << "\t" << cutoff << "\n";
  archivo >> dato >> expanded_lx;
  cout << dato << "\t" << expanded_lx << "\n";
  archivo >> dato >> expanded_ly;
  cout << dato << "\t" << expanded_ly << "\n";
  archivo >> dato >> nVis;
  cout << dato << "\t" << nVis << "\n";
  
  nombresVisibilidades = new char *[nVis];
  for (i = 0; i < nVis; i++) {
    nombresVisibilidades[i] = new char[512];
    archivo >> nombresVisibilidades[i];
    cout << nombresVisibilidades[i] << "\n";
  }
  
  cout << "Construyendo FuncionBayesVoronoi\n";
  (*this) = FuncionBayesVoronoiCBII (nombreImagen, nombresVisibilidades, 
				     nVis, nPols, init_value, init_gauss, 
				     entropia, cuantaSize,
				     expanded_lx, expanded_ly);
  cout << "Construido\n";
  cout << "npols:" << this->n_pars << " vs " << nPols << " \n";
  if (this->fL->malla == NULL) {
    cout << "2) r->fL->malla == NULL\n";
  }
  delete [] nombreImagen;
  delete [] dato;
  for (i = 0; i < nVis; i++) {
    delete [] nombresVisibilidades[i];
  }
  delete [] nombresVisibilidades;
  cout << "leer: r->nVis = " << this->nVis << "\n";
  cout << "entropia: " << this->fL->entropia << "\n";
  
  // Inicializamos pos
  double *pos_aux = new double [3 * this->n_pars];
  ini->inicializar (pos_aux, this->n_pars * 3, 0.0);
  for (int i = 0; i < this->n_pars; i++) {
    pos [2 * i] = pos_aux [3* i];
    pos [2 * i + 1] = pos_aux [3* i + 1];
  }
  delete [] pos_aux;
}

double FuncionBayesVoronoiCBII::f (double pars[]) {
  double ret, *parsF;
  parsF = new double [this->n_pars * 3];
  for (int i = 0; i < this->n_pars; i++) {
    parsF [3 * i] = pos[2 * i];
    parsF [3 * i + 1] = pos[2 * i + 1];
    parsF [3 * i + 2] = pars[i];
    //cout << "pos[" << 2*i << "] = " << pos[2 * i] 
    // << ", pos[" << 2*i + 1 << "] = " << pos[2 * i + 1] << "\n"; 
    //cout << "Poligono: {" << parsF [3 * i] << ", " << parsF [3 * i + 1] << ", " << parsF [3 * i + 2] << "}\n";
  }
  ret = L (fL, parsF, fL->n_pols, MOCKCBI);
  malla = fL->malla;
  delete [] parsF;
  return ret;
}

void FuncionBayesVoronoiCBII::df (double pars[], double grad[]) {
  double ret, *parsF, *gradF;

  parsF = new double [this->n_pars * 3];
  gradF = new double [this->n_pars * 3];
  for (int i = 0; i < this->n_pars; i++) {
    parsF [3 * i] = pos[2 * i];
    parsF [3 * i + 1] = pos[2 * i + 1];
    parsF [3 * i + 2] = pars[i];
  }
  dL (fL, parsF, gradF, fL->n_pols, 0, MOCKCBI);
  malla = fL->malla;

  for (int i = 0; i < this->n_pars; i++) {
   grad[i] = gradF [3 * i + 2];
  }

  delete [] parsF;
  delete [] gradF;
}

void FuncionBayesVoronoiCBII::guardarInfo (string info){
  // Para este caso usaremos info como el sufijo de los archivos a
  // guargar.
  string nombre_archivo;
  
  nombre_archivo = "!FuncionBayesVoronoiCBI" + info + ".fits";

  char *nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  

  cout << "guardando " << nombreC << "\n";
  do_write_fits (fL->fg_image, nombreC);

  nombre_archivo = "MallaBayesVoronoiCBI" + info + ".dat";
  delete [] nombreC;
  nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;

  if (fL->malla != NULL) {
    imprimirMallaArchivo (fL->malla, nombreC);
  }

  delete [] nombreC;
}

void FuncionBayesVoronoiCBII::setNPars (int n_pars) {
  this->n_pars = n_pars;
  fL->n_pols = n_pars / 3;
  fL->malla->nPols = n_pars / 3;
  this->malla->nPols = n_pars / 3;
}

FuncionBayesVoronoiCBII::~FuncionBayesVoronoiCBII () {
  if (fL != NULL) {
    cout << "eliminando funcL en funcionBayesVoronoi\n";
    eliminarFuncL(fL);
  }
  delete [] pos;
}

FuncionBayesVoronoiCBII& FuncionBayesVoronoiCBII::operator= 
(const FuncionBayesVoronoiCBII& func) {
  int i;
  this->nVis = func.nVis;
  this->n_pars = func.n_pars;
  this->pos = new double [this->n_pars * 2];

  for (i = 0; i < this->n_pars * 2; i++) {
    this->pos[i] = func.pos[i];
  }
  
  this->nombreImagen = new char[256];
  strcpy (this->nombreImagen, func.nombreImagen);

  this->nombresVis = new char*[nVis];
  for (i = 0; i < nVis; i++){
    this->nombresVis[i] = new char[256];
    strcpy(this->nombresVis[i], func.nombresVis[i]);
  }
  
  this->fL = copiarFuncL (func.fL);
  this->malla = fL->malla;
  
  return *this;
}
