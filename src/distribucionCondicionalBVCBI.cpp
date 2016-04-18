#include "distribucionCondicionalBVCBI.h"

DistribucionCondicionalBVCBI::DistribucionCondicionalBVCBI(char * nombreImagen, char **nombresVis, 
							   int nVis, int n,
							   int entropia, double cuantaSize,
							   double expanded_lx, double expanded_ly, 
							   int N) {
  f = new FuncionBayesVoronoiCBI(nombreImagen, nombresVis, nVis, n, entropia, cuantaSize,
				 expanded_lx, expanded_ly);
  n_pars = f->getNPars();
  this->N = N;
}

DistribucionCondicionalBVCBI::DistribucionCondicionalBVCBI(char *nombreArchivoInput, int N) {
  f = new FuncionBayesVoronoiCBI(nombreArchivoInput);
  n_pars = f->getNPars();
  this->N = N;
}

double DistribucionCondicionalBVCBI::generar (double *p, int i) {
  if (i % 3 != 2) {
    return generarPosicion (p, i);
  }
  return generarIntensidad (p, i);
}

double DistribucionCondicionalBVCBI::generarPosicion (double *p, int i) {
  double x, a, b, d, fTotal = 0, f_rand;
  double fs [N];
  int cont;

  a = 0;
  b = 1;

  d = (b - a) / N;
  cout << "Cotas: " << a << " -> " << b <<"\n";

  cont = 0;
  for (x = a + d/2; x < b; x += d) {
    fs[cont] = this->evaluar (p, x, i);
    cout << "x pos = " << x << "\n";
    fTotal += fs[cont++];
  }
  
  f_rand = fTotal * ((double) rand() /  RAND_MAX);
  cout << "fTotal = " << fTotal << ", frand = " << f_rand << "\n";
  fTotal = 0;
  cont = 0;
  for (x = a + d/2; fTotal < f_rand && cont < N; x += d) {
    cout  << "x = " << x << "\n";
    fTotal += fs[cont++]; 
  }

  return x - d;
}

double DistribucionCondicionalBVCBI::generarIntensidad (double *p, int i) {
  double x, a, b, fTotal = 0, f_rand;
  int cont, n;
  this->busquedaCotas (p, i, &a, &b);
  
  a = ceil (a);
  b = floor (b);
  n = (int) ceil (b - a);
  double fs [n];

  cout << "Cotas I: " << a << " -> " << b <<"\n";

  cont = 0;
  for (x = a; x <= b; x += 1) {
    cout  << "x1 = " << x << "\n";
    fs [cont] = this->evaluar (p, x, i);
    fTotal += fs[cont++];
  }
  
  f_rand = fTotal * ((double) rand() /  RAND_MAX);
  cout << "fTotal I = " << fTotal << ", frand = " << f_rand << "\n";
  fTotal = 0;
  cont = 0;
  for (x = a; fTotal < f_rand && cont < n; x += 1) {
    cout  << "x2 = " << x << "\n";
    fTotal += fs[cont++];
  }

  return x - 1;
  
}

double DistribucionCondicionalBVCBI::evaluar (double *p) {
  //return exp (-f->f(p));
  return  1 / f->f(p);
}

double DistribucionCondicionalBVCBI::evaluar (double *p, double x, int i) {
  double ret = f->f(p);
  double x_old = p[i];

  //cout << "DistribucionCondicionalBVCBI::evaluar 1\n";
  p[i] = x;
  ret -= f->f(p);
  //cout << "DistribucionCondicionalBVCBI::evaluar 2\n";
  ret = exp (ret);
  p[i] = x_old;
  return ret;
}

void DistribucionCondicionalBVCBI::busquedaCotas (double *p, int i, double *a, double *b) {
  double delta_ini = 1, ret = 0, ret_old,Lmax, Lmin, delta, p_i;
  Lmax = log (DOUBLE_MAX);
  //Lmin = log (DOUBLE_MIN);
  Lmin = log (1e-100);
  cout << "Lmax = " << Lmax << ", Lmin = " << Lmin << "\n";
  // = log (DOUBLE_MAX), LMin = log (1/DOUBLE_MAX);

  // Buscamos b;
  p_i = p[i];
  do {
    //ret = this->evaluar (p, p_i, i);
    double p_old = p[i];
    p[i] = p_i;
    ret = -f->f(p);
    p[i] = p_old;
    ret += f->f(p);
    //cout << "1b-----delta_ini = " << delta_ini << ", p_i = " << p_i <<", ret = " << ret << "\n";

    for (delta = delta_ini; ret > Lmin && ret < Lmax; delta *= 2) {
      p_i += delta;
      //ret = this->evaluar (p, p_i, i);
      p_old = p[i];
      p[i] = p_i;
      ret = -f->f(p);
      if (p[i] != p_i) {
	//(*b) = p[i];
	p_i = p[i];
	cout << "DISTINTO b " << p[i] << " != " << p_i << ", delta = " << delta << "\n";
	delta = 0;
	continue;
      }
      p[i] = p_old;
      ret += f->f(p);
      //cout << "b) delta = " << delta << ", ret = " << ret << ", p_i = " << p_i << "\n";
    }
    p_i -= delta / 2;
    delta_ini = delta / 4;
  } while (delta > 1e-9);
  (*b) = p_i;

  // Buscamos a;
  delta_ini = 1;
  p_i = p[i];
  do {
    //ret = this->evaluar (p, p_i, i);
    double p_old = p[i];
    p[i] = p_i;
    ret = -f->f(p);
    p[i] = p_old;
    ret += f->f(p);
    //cout << "1-----delta = " << delta << ", ret = " << ret << "\n";

    for (delta = delta_ini; ret > Lmin && ret < Lmax; delta *= 2) {
      p_i -= delta;
      //ret = this->evaluar (p, p_i, i);
      p_old = p[i];
      p[i] = p_i;
      ret = -f->f(p);
      if (p[i] != p_i) {
	(*a) = p[i];
	cout << "DISTINTO a " << p[i] << " != " << p_i << ", delta = " << delta << "\n";
	return;
      }
      p[i] = p_old;
      ret += f->f(p);
      //cout << "delta = " << delta << ", ret = " << ret << "\n";
    }
    p_i += delta / 2;
    delta_ini = delta / 4;
  } while (delta > 1e-9);
  (*a) = p_i;
}

DistribucionCondicionalBVCBI::~DistribucionCondicionalBVCBI() {
  delete f;
}
