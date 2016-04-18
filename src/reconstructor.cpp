#include "reconstructor.h"

string i2s (int numero, int largo){
  char *numeroC = new char [largo + 1];
  int i;

  for (i = 0; i < largo; i++){
    numeroC [i] = '0' + (int) ((numero % (int) pow (10.0, (double) largo - i))
			       / pow (10.0, largo - i - 1.0));
  }
  numeroC [i] = '\0';

  string s;
  s += numeroC;
  delete [] numeroC;
  return s;
}

Reconstructor::Reconstructor () {
  p = NULL;
}

Reconstructor::Reconstructor (int n_pars, Inicializador *ini, double xini) {
  int i, j, k, nx, ny, imax = 0, jmax = 0;
  FILE* archivoLog = fopen ("reconstructor.log", "w");
  fclose (archivoLog);

  this->iter = 0;
  this->n_pars = n_pars;
  p = new double [n_pars];
  this->ini = ini;
  //this->ini->inicializar (this->p, n_pars);
  if (ini != NULL) this->ini->inicializar (this->p, n_pars, xini);
  //cout << "construido Reconstructor, valor de la funcion: " << f->f(p) << "\n";
}

Reconstructor::~Reconstructor() {
  //delete f;
  //delete ini;
  delete [] p;
}

void Reconstructor::reinicializar (char *nombre_archivo, int n_iter) {
  this->iter = n_iter;
  this->reinicializar (nombre_archivo);
}

void Reconstructor::reinicializar (char *nombre_archivo) {
  std::ifstream archivo (nombre_archivo);
  int n;

  archivo >> n;
  cout << "n = " << n << "\n";
  if (n != n_pars) {
    cerr << "ERROR: No se pudo reinicializar el reconstructor,\n"
	 << "       distinto numero de parametros: " << n << " vs " 
	 << n_pars << "\n";
    exit(1);
  }
  if (p != NULL){
    delete [] p;
  }
  
  p = new double [n];
  for (int i = 0; i < n; i++) {
    archivo >> p[i];
  }
}

void Reconstructor::setNPars (int n_pars) {
  double *p_aux = new double [n_pars];
  for (int i = 0; i < n_pars && i < this->n_pars; i++) {
    p_aux[i] = p[i];
  }
  delete [] p;
  p = p_aux;
  this->n_pars = n_pars;
}

void Reconstructor::setPars (double *pars) {
  delete [] p;
  p = pars;
}

void Reconstructor::setPar (int i, double valor) {
  if (i >= n_pars) {
    cerr << "ERROR en Reconstructor::setPar: i = " << i << " >= n_pars = " << n_pars << ".\n";
    exit (0); 
  }
  p[i] = valor;
}

void Reconstructor::imprimirLog (char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen ("reconstructor.log", "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
}

void Reconstructor::imprimirLog (char *nombreArchivo, char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen (nombreArchivo, "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
}

Reconstructor& Reconstructor::operator= (const Reconstructor& r) {

  this->iter = r.iter;
  this->n_pars = r.n_pars;
 
  this->p = new double[n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->p[i] = r.p[i];
  }

  this->ini = r.ini;
  return *this;
}
