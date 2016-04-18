#include "reconstructorBayesiano.h"

ReconstructorBayesiano::ReconstructorBayesiano (DistribucionPbb *d, 
						Inicializador *ini, int n_iter)
  : Reconstructor (d->getNPars(), ini) {
  this->n_iter = n_iter;
  this->distribucion = d;
  iter = 0;
}

double ReconstructorBayesiano::run () {
  double f_mejor = -1e300, func;
  double *p_aux = new double [distribucion->getNPars()];
  std::fstream arc ("reconstructor.log");
  arc.close();

  for (int j = 0; j < distribucion->getNPars(); j++) {
    p_aux[j] = p[j];
    //cout << "p_aux[" << j << "] = " << p_aux[j] << "\n";
  }

  for (iter = iter; iter < n_iter; iter++) {
    cout << "---->Iteracion: " << iter << "\n";
    if ((func = distribucion->f(p_aux)) > f_mejor) {
      cout << "Mejor en iteracion " << iter << ", " << p_aux[0] << "\n";
      cout << "valor: " << func << "\n";
      distribucion->guardarInfo ("Mejor" + i2s (iter, 9));
      f_mejor = func;
      imprimirLog ("MEJOR: iteracion %d: %g en\n", iter, f_mejor);
      for (int j = 0; j < distribucion->getNPars(); j++) {
	imprimirLog ("  p[%d] = %g\n", j, p_aux[j]);
	p[j] = p_aux[j];
      }
    }
    else {
      //distribucion->guardarInfo (i2s (iter, 9));
    }
    //imprimimos en reconstructor.log
    //if (iter % 100 == 0) {
      cout << "iter: " << iter << "\n";
      imprimirLog ("Iteracion %d, calculando funcion con parametros:\n", iter);
      for (int j = 0; j < distribucion->getNPars(); j++) {
	imprimirLog ("%g ", j, p_aux[j]);
	//imprimirLog ("  p[%d] = %g\n", j, p_aux[j]);
      }
      imprimirLog ("\n valor funcion: %g, funcion mejor: %g\n", func, f_mejor);
      //}

    distribucion->generar (p_aux);
  }
  delete [] p_aux;
  return 0.0;
}

ReconstructorBayesiano::~ReconstructorBayesiano() {
  //delete distribucion;
  //delete ini;
  //delete [] p;
  //delete f;
}
