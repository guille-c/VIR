#include "funcionBayesVoronoiCBI.h"

FuncionBayesVoronoiCBI::FuncionBayesVoronoiCBI (char * nombreImagen, char **nombresVis, 
						int nVis, int n,
						int entropia, double cuantaSize,
						double expanded_lx, double expanded_ly,
						double beamSize) :
FuncionCBI (nombreImagen, nombresVis, nVis, n){
  std::fstream archivo ("valoresFuncionBVCBI.log");
  archivo.close();

  this->n_pars = n * 3;
  fL = newFuncL(nombreImagen, nombresVis, nVis, n, entropia, -1, 
		expanded_lx, expanded_ly, beamSize);
  //fL = newFuncL(nombreImagen, nombresVis, nVis, n, entropia, 6e8);

  this->ruido = fL->difmapNoise;
  cout << "ruido = " << ruido << " [K]\n";
  cout << "difmapNoise = " << fL->difmapNoise << " [K]\n";
  this->malla = fL->malla;
}

FuncionBayesVoronoiCBI::FuncionBayesVoronoiCBI (char* nombre_archivo) {
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
  (*this) = FuncionBayesVoronoiCBI (nombreImagen, nombresVisibilidades, 
				    nVis, nPols, entropia, cuantaSize,
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
  if (ruido != fL->difmapNoise) {
    cerr << "ERROR en FuncionBayesVoronoiCBI: ruido != difmapNoise\n";
    cerr << "                                 " << ruido << " != " << fL->difmapNoise << "\n";
    exit (1);
  }
  
}

double FuncionBayesVoronoiCBI::f (double pars[]) {
  double ret = L (fL, pars, fL->n_pols, MOCKCBI);
  malla = fL->malla;
  if (ruido != fL->difmapNoise) {
    cerr << "ERROR en f: ruido != difmapNoise\n";
    cerr << "            " << ruido << " != " << fL->difmapNoise << "\n";
    exit (1);
  }
  //cout << "func = " << ret << "\n";
  return ret;
}

void FuncionBayesVoronoiCBI::df (double pars[], double grad[]) {
  //cout << "calculando df/n";
  dL (fL, pars, grad, fL->n_pols, EXACTA, MOCKCBI);
  malla = fL->malla;
}

void FuncionBayesVoronoiCBI::guardarAtenuaciones () {
  struct image *im = do_read (fL->nombreFits);
  char nombre [255];

  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for(int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int k = 1; k <= im->npixels; k++){
	im->pixels[k-1] = fL->atten[arch][k][iff];
      }
      double freq = fL->header_obs[arch]->iffreq[iff] * 1e-9;
      sprintf (nombre, "!atenuacion%d_%g.fits", arch, freq);
      do_write_fits (im, nombre);
    }
  }
  
  delete_map (im);
}

void FuncionBayesVoronoiCBI::guardarInfo (string info){
  // Para este caso usaremos info como el sufijo de los archivos a
  // guardar.
  string nombre_archivo;
  std::ofstream archivo;
  archivo.open ("valoresFuncionBVCBI.log", ios::app);

  archivo << info << "\t" << fL->chi2/2 - fL->S << "\t" << fL->chi2 << "\t"
	  << fL->S << "\n";
  
  archivo.close();

  nombre_archivo = "!FuncionBayesVoronoiCBI" + info + ".fits";

  char *nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  

  //cout << "guardando " << nombreC << "\n";
  do_write_fits (fL->fg_image, nombreC);

  nombre_archivo = "MallaBayesVoronoiCBI" + info + ".dat";
  delete [] nombreC;

  nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;

  if (fL->malla != NULL) {
    imprimirMallaArchivo (fL->malla, nombreC);
  }

  // Visibilidades modelo
  for (int i = 0; i < nVis; i++) {
    nombre_archivo = "!VIS_MOD_" + info + "_" + i2s(i, 2) + ".sub";
    
    nombreC = new char [nombre_archivo.length() + 1];
    nombre_archivo.copy(nombreC, string::npos);
    nombreC [nombre_archivo.length()] = 0;
    
    int status = do_write("write", nombreC, header_obs [i], samples_mod [i], nombresVis[i]);
    if(status != SUCCESS) {
      cerr << "Error al escribir archivo uvf " << i << " en FuncionCBI::guardarResiduos\n";
      //exit(1);
    }
    delete [] nombreC;
  }

  // residuos
  for (int i = 0; i < nVis; i++) {
    nombre_archivo = "!VIS_RES_" + info + "_" + i2s(i, 2) + ".sub";
    
    char *nombreC = new char [nombre_archivo.length() + 1];
    nombre_archivo.copy(nombreC, string::npos);
    nombreC [nombre_archivo.length()] = 0;
    
    struct uvf_sample *samples_res; 
    
    samples_res = (struct uvf_sample*) malloc((header_obs[i]->nsamp) * sizeof(struct uvf_sample));
    residual_vis (header_obs[i], samples_obs[i], samples_mod[i], samples_res);
    
    int status = do_write("write", nombreC, header_obs[i], samples_res, nombresVis[i]);
    if(status != SUCCESS) {
      cerr << "Error al escribir archivo uvf al guardar residuos\n";
      //exit(1);
    }
    free (samples_res);
    delete [] nombreC;
  }
}

void FuncionBayesVoronoiCBI::setNPars (int n_pars) {
  this->n_pars = n_pars;
  fL->n_pols = n_pars / 3;
  //fL->malla->nPols = n_pars / 3;
  //this->malla->nPols = n_pars / 3;
}

FuncionBayesVoronoiCBI::~FuncionBayesVoronoiCBI () {
  if (fL != NULL) {
    cout << "eliminando funcL en funcionBayesVoronoi\n";
    eliminarFuncL(fL);
  }
}

FuncionBayesVoronoiCBI& FuncionBayesVoronoiCBI::operator= 
(const FuncionBayesVoronoiCBI& func) {
  int i;
  this->nVis = func.nVis;
  this->n_pars = func.n_pars;
  
  this->nombreImagen = new char[256];
  strcpy (this->nombreImagen, func.nombreImagen);

  this->nombresVis = new char*[nVis];
  for (i = 0; i < nVis; i++){
    this->nombresVis[i] = new char[256];
    strcpy(this->nombresVis[i], func.nombresVis[i]);
  }
  
  this->fL = copiarFuncL (func.fL);
  this->malla = fL->malla;
  this->ruido = func.ruido;
  
  return *this;
}
