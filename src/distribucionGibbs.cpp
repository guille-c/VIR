#include "distribucionGibbs.h"

DistribucionGibbs::DistribucionGibbs (DistribucionCondicional *Pcond)
  : DistribucionPbb (Pcond->getNPars()) {
  this->Pcond = Pcond;
}

double DistribucionGibbs::f (double *p) {
  /*
  double ret = 1;
  
  for (int i = 0; i < n_pars; i++) {
    ret *= Pcond->evaluar (p, p[i], i);
  }
  return ret;
  */
  //cout << "DistribucionGibbs::f\n";
  return Pcond->evaluar (p);
}

string i2s2 (int numero, int largo){
  char *numeroC = new char [largo + 1];
  int i;

  for (i = 0; i < largo; i++){
    numeroC [i] = '0' + (int) ((numero % (int) pow (10.0, (double) largo - i))
			       / pow (10.0, largo - i - 1.0));
  }
  numeroC [i] = 0;

  string s;
  s += numeroC;
  delete [] numeroC;
  return s;
}

void DistribucionGibbs::generar (double *x) {
  for (int i = 0; i < n_pars; i++) {
    cout << "--->Generando " << i << "\n";
    x[i] = Pcond->generar (x, i);
    //Pcond->guardarInfo("Condicional" + i2s2 (i, 4));
  }
}
