#include "reconstructorIncrementalDoble.h"

ReconstructorIncrementalDoble::ReconstructorIncrementalDoble (Funcion *f, Inicializador *ini, 
							      Incrementador *inc, int n_iter, 
							      Reconstructor *r)
  : ReconstructorIncremental (f, ini, inc, n_iter) {
  this->r = r;
}

double ReconstructorIncrementalDoble::run () {
  double valor;
  fstream archivo;
  archivo.open ("ReconstructorIncrementalDoble.log", ios::out);
  archivo.close();
  archivo.open ("ReconstructorIncrementalDoble.log", ios::app);

  for (int i = 0; i < n_iter; i++) {
    f->guardarInfo (i2s (i, 4));
    valor = inc->incrementar (&p, f);
    archivo << i << "\t" << valor << "\n";
    n_pars = f->getNPars();

    // Segundo reconstructor
    r->setNPars(f->getNPars());
    for (int i = 0; i < n_pars; i++){
      r->setPar(i, p[i]);
    }
    r->run();
  }
  this->n_pars = f->getNPars();

}

ReconstructorIncrementalDoble::~ReconstructorIncrementalDoble() {
  delete inc;
}
