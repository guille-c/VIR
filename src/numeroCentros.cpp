#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <cstdlib>
#include <fstream.h>

#include "funcionBayesVoronoiCBI.h"
#include "funcionChi2VoronoiImagen.h"
#include "reconstructorGC.h"
#include "reconstructorIncremental.h"
#include "incrementadorMallaImagen.h"
#include "inicializadorVoronoiUniforme.h"

void numeroCentrosGC(int argc, char **argv);
double calculoFuncAjuste(int argc, char **argv);
void numeroCentrosAjusteGC(int argc, char **argv);

int main (int argc, char **argv){
  //double d = calculoFuncAjuste(argc, argv);
  numeroCentrosAjusteGC(argc, argv);
  //numeroCentrosGC(argc, argv);
}

double calculoFuncAjuste(int argc, char **argv) {
  double ret;
  Funcion *func = new FuncionChi2VoronoiImagen (argv[1], 5);
  Reconstructor *r1 = new ReconstructorIncremental (func, new InicializadorVoronoiUniforme (),
						    new IncrementadorMallaImagen (argv[1], -1), 45);
  ret = r1->run();
  delete func;

  double *p = new double [func->getNPars()];
  double p_max = 0;
  for (int i = 0; i < func->getNPars(); i++) {
    p[i] = r1->getPar(i);
    if (i % 3 == 2 && p[i] > p_max) {
      p_max = p[i];
    }
  }
  
  FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[2]);
  double Irms = f->getRuido(), factor = 10 / p_max;
  cout << "Irms = " << Irms << "\n";

  for (int i = 0; i < f->getNPars() / 3; i++) {
    //cout << "Cambiando " << p[3 * i + 2];
    p[3 * i + 2] *= factor;
    //cout << " por " << p[3 * i + 2] << "\n";
  }
  
  ret = f->f(p);
  funcL *fL = f->getFL();

  double n_dat = 2*fL->n_archivos * fL->header_obs[0]->nif * fL->header_obs[0]->nsamp;
  cout << "n_dat = " << n_dat << "\n";
  cout << "L = " << ret << "\n";
  cout << "chi2 = " << f->getFL()->chi2 * 2 << "\n";
  cout << "chi2 reducido = " << f->getFL()->chi2 / (n_dat - f->getNPars()) << "\n";
  cout << "S = " << f->getFL()->S << "\n";
  f->guardarInfo("L");
  delete r1;

  r1 = new ReconstructorGC (f, new InicializadorVoronoiUniforme (), 1e-10);
  for (int i = 0; i < f->getNPars(); i++) {
     r1->setPar(i, p[i]);
  }
  r1->run();
  for (int i = 0; i < f->getNPars(); i++) {
    p[i] = r1->getPar(i);
  }
  ret = f->f(p);
  cout << "n_dat = " << n_dat << "\n";
  cout << "L = " << ret << "\n";
  cout << "chi2 = " << f->getFL()->chi2 * 2 << "\n";
  cout << "chi2 reducido = " << f->getFL()->chi2 / (n_dat - f->getNPars()) << "\n";
  cout << "S = " << f->getFL()->S << "\n";
  delete f;
  delete r1;
  return ret;
}

void numeroCentrosGC(int argc, char **argv) {
  int n;
  std::ofstream archivo ("L_vs_n.dat");
  archivo.close();
  archivo.open ("chi2_vs_n.dat");
  archivo.close();
  archivo.open ("S_vs_n.dat");
  archivo.close();
  for (n = 10; n <= 1e2; n += 10) {
    FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[1], argv + 2, argc - 2, n, 1, -1, 250.0, 250.0);
    double pars [f->getNPars()];
    ReconstructorGC *r = new ReconstructorGC (f, new InicializadorVoronoiUniforme (), 1e-10);
    r->run();

    for (int i = 0; i < f->getNPars(); i++) {
      pars[i] = r->getPar(i);
    }
    double ret = f->f(pars);
    f->guardarInfo ("NC" + i2s(n, 4));
    archivo.open ("L_vs_n.dat", ios::app);
    archivo << n << "\t" << ret << "\n";
    archivo.close();
    archivo.open ("chi2_vs_n.dat", ios::app);
    archivo << n << "\t" << f->getFL()->chi2 << "\n";
    archivo.close();
    archivo.open ("S_vs_n.dat", ios::app);
    archivo << n << "\t" << f->getFL()->S << "\n";
    archivo.close();

    delete f;
    delete r;
  }
}

void numeroCentrosAjusteGC(int argc, char **argv) {
  int n;
  
  std::ofstream archivo ("L_vs_n.dat");
  archivo.close();
  archivo.open ("chi2_vs_n.dat");
  archivo.close();
  archivo.open ("chi2reducido_vs_n.dat");
  archivo.close();
  archivo.open ("S_vs_n.dat");
  archivo.close();
  
  for (n = 1000; n < 1e4; n += 1000) {
    Funcion *func = new FuncionChi2VoronoiImagen (argv[1], 5);
    Reconstructor *r1 = new ReconstructorIncremental (func, new InicializadorVoronoiUniforme (),
						      new IncrementadorMallaImagen (argv[1], -1), n - 5);
    double ret = r1->run();
    delete func;
  
    double *p = new double [func->getNPars()];
    double p_max = 0;
    for (int i = 0; i < func->getNPars(); i++) {
      p[i] = r1->getPar(i);
      if (i % 3 == 2 && p[i] > p_max) {
	p_max = p[i];
      }
    }
    
    //FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[2]);
    FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[1], argv + 2, argc - 2, n, 1, -1, 250.0, 250.0);
    double Irms = f->getRuido(), factor = 10 / p_max;
    cout << "Irms = " << Irms << "\n";
    
    for (int i = 0; i < f->getNPars() / 3; i++) {
      p[3 * i + 2] *= factor;
    }
    
    ret = f->f(p);
    funcL *fL = f->getFL();
    
    double n_dat = 2*fL->n_archivos * fL->header_obs[0]->nif * fL->header_obs[0]->nsamp;
    f->guardarInfo("Ajuste" + i2s(n, 5));
    delete r1;
    
    // Gradiente Conjugado
    double pars [f->getNPars()];
    ReconstructorGC *r = new ReconstructorGC (f, new InicializadorVoronoiUniforme (), 1e-10);
    for (int i = 0; i < f->getNPars(); i++) {
      r->setPar(i, p[i]);
    }
    r->run();

    for (int i = 0; i < f->getNPars(); i++) {
      pars[i] = r->getPar(i);
    }
    ret = f->f(pars);
    f->guardarInfo ("NC" + i2s(n, 4));
    archivo.open ("L_vs_n.dat", ios::app);
    archivo << n << "\t" << ret << "\n";
    archivo.close();
    archivo.open ("chi2_vs_n.dat", ios::app);
    archivo << n << "\t" << f->getFL()->chi2 << "\n";
    archivo.close();
    archivo.open ("chi2reducido_vs_n.dat", ios::app);
    archivo << n << "\t" << f->getFL()->chi2 / (n_dat - f->getNPars()) << "\n";
    archivo.close();
    archivo.open ("S_vs_n.dat", ios::app);
    archivo << n << "\t" << f->getFL()->S << "\n";
    archivo.close();

    delete f;
    delete r;
  }
}

