#include "inicializadorArchivo.h"

InicializadorArchivo::InicializadorArchivo (char * nombre_archivo) {
  strcpy (this->nombre_archivo, nombre_archivo);
}

void InicializadorArchivo::inicializar (double pars[], int n, double init) {
  std::ifstream archivo (nombre_archivo);
  
  archivo >> n;
  cout << "n = " << n << "\n";
  
  for (int i = 0; i < n; i++) {
    archivo >> pars [i];
  }
  
  archivo.close();
}
