#include "reconstructorIncremental.h"

ReconstructorIncremental::ReconstructorIncremental (Funcion *f, Inicializador *ini, 
						    Incrementador *inc, int n_iter)
  : Reconstructor (f->getNPars(), ini) {
  this->f = f;
  this->inc = inc;
  this->n_iter = n_iter;
}

double ReconstructorIncremental::run () {
  double valor;
  fstream archivo;
  archivo.open ("ReconstructorIncremental.log", ios::out);
  archivo.close();
  archivo.open ("ReconstructorIncremental.log", ios::app);

  for (int i = 0; i < n_iter; i++) {
    f->guardarInfo (i2s (i, 4));
    valor = inc->incrementar (&p, f);
    archivo << i << "\t" << valor << "\n";
    cout << i << "\t" << valor << "\n";
  }
  this->n_pars = f->getNPars();
}

ReconstructorIncremental::~ReconstructorIncremental() {
  //delete inc;
}
