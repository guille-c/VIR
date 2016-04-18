#include "distribucionCondicional.h"

void DistribucionCondicional::busquedaCotas (double *p, int i, double *a, double *b) {
  double delta_ini = 1e-10, ret = 0, ret_old,Lmax, Lmin, delta, p_i;
  Lmax = DOUBLE_MAX;
  //Lmin = log (DOUBLE_MIN);
  Lmin = 1e-100;

  // Buscamos b;
  p_i = p[i];
  do {
    ret = this->evaluar (p, p_i, i);
    for (delta = delta_ini; ret > Lmin && ret < Lmax; delta *= 2) {
      p_i += delta;
      ret = this->evaluar (p, p_i, i);
    }
    p_i -= delta / 2;
    delta_ini = delta / 4;
  } while (delta > 1e-9);
  (*b) = p_i;

  // Buscamos a;
  p_i = p[i];
  do {
    ret = this->evaluar (p, p_i, i);
    for (delta = delta_ini; ret > Lmin && ret < Lmax; delta *= 2) {
      p_i -= delta;
      ret = this->evaluar (p, p_i, i);
    }
    p_i += delta / 2;
    delta_ini = delta / 4;
  } while (delta > 1e-9);
  (*a) = p_i;
}
