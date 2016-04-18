#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"

#include "funcionBayesVoronoiCBI.h"
#include "funcionPbbVoronoiCBI.h"
#include "funcionExp.h"
#include "funcionMemCBI.h"
#include "funcionCero.h"
#include "funcionEntropiaNatural.h"
#include "funcionReconstructor.h"
#include "funcionChi2VoronoiImagen.h"
#include "funcionChi2VoronoiImagenN.h"

#include "mcheck.h"

#include "funcionDummy.h"

void comprobacion_gradiente (Funcion *f, double delta);
void comprobacion_gradiente_Voronoi (int argc, char **argv);
void igualar (double *p, double *p_ori, int nPols);

int main (int argc, char **argv) {
  //Funcion *f = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionCero (64 * 64));
  //Funcion *f = new FuncionEntropiaNatural (64 * 64);
  //Funcion *f  = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionEntropiaNatural(64*64));
  //comprobacion_gradiente (f, 1e-3);
  comprobacion_gradiente_Voronoi (argc, argv);
}

void comprobacion_gradiente (Funcion *f, double delta) {
  int n = f->getNPars();
  double *p = new double [n];
  double *g = new double [n];
  double *g_aprox = new double [n];
  double f1, f2, p_old;
  ofstream out;

  for (int i = 0; i < n; i++) {
    p[i] = (double) rand() / RAND_MAX * 10;
  }
  f1 = f->f(p);
  f->df (p, g);
  out.open ("gradiente_exacto.dat");
  for (int i = 0; i < n; i++) {
    out << g[i] << "\n";
  }  
  out.close();

  out.open ("gradiente_aprox.dat");
  for (int i = 0; i < n; i++) {
    if (i % 100 == 0) cout << i << "\n";
    p_old = p[i];
    p[i] += delta;
    f2 = f->f(p);
    p[i] = p_old;
    g_aprox[i] = (f2 - f1) / delta;
    out << g_aprox[i] << "\n";
  }  
  out.close();

  out.open ("gradientes.dat");
  for (int i = 0; i < n; i++) {
    out << g[i] << "\t" << g_aprox[i] << "\n";
  }  
  out.close();

  delete [] p;
  delete [] g;
  delete [] g_aprox;
}

void comprobacion_gradiente_Voronoi (int argc, char **argv) {
  int nPols, i, nIter = 200;
  double ftol = 1e-10;
  char ** vis = new char *[1];
  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (argv[1]);
  double *p, *g, *p_ori;
  ofstream out;

  funcL *fL = func->getFL();
  func->guardarAtenuaciones();
  
  mtrace ();

  nPols = func->getNPars() / 3;
  p = new double [nPols * 3];
  p_ori = new double [nPols * 3];
  g = new double [nPols * 3];
  
  for (i = 0; i < nPols; i++) {
    p[3 * i]     = (double) rand() / RAND_MAX;
    p[3 * i + 1] = (double) rand() / RAND_MAX;
    p[3 * i + 2] = (double) rand() / RAND_MAX * 10.5;
    p_ori[3 * i]     = p[3 * i];
    p_ori[3 * i + 1] = p[3 * i + 1];
    p_ori[3 * i + 2] = p[3 * i + 2];
  }
  
  func->df(p, g);
  func->guardarInfo("");
  //exit(1);

  out.open ("gradientex_mockcbi.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i] << "\n";
  }
  out.close();
  out.open ("gradientey_mockcbi.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i + 1] << "\n";
  }
  
  out.close();
  out.open ("gradienteI_mockcbi.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i + 2] << "\n";
  }
  out.close();

  dL (fL, p, g, fL->n_pols, EXACTA, MOCK_EXACTO);

  out.open ("gradientex_exacto.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i] << "\n";
  }
  out.close();
  out.open ("gradientey_exacto.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i + 1] << "\n";
  }
  
  out.close();
  out.open ("gradienteI_exacto.dat");
  for (i = 0; i < nPols; i++) {
    out << i << "\t" << g[3 * i + 2] << "\n";
  }
  out.close();

  //exit(0);

  for (int indice = 0; indice < nPols; indice++) {
    char nombre[256];
    double delta;

    cout << "Comenzando iteracion " << indice << "\n";
    
    sprintf(nombre, "funcionx_%d.dat", indice);
    out.open (nombre);
    for (delta = -0.1; delta < 0.1; delta += 0.2 / nIter) {
      double f;
      
      p[indice * 3] -= delta;
      f = func->f (p);
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();
    
    sprintf(nombre, "funciony_%d.dat", indice);
    out.open (nombre);
    for (delta = -0.1; delta < 0.1; delta += 0.2 / nIter) {
      double f;
      
      p[indice * 3 + 1] += delta;
      f = func->f (p);
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();  
    
    sprintf(nombre, "funcionI_%d.dat", indice);
    out.open (nombre);
    for (delta = -1; delta < 1; delta += 2.0 / nIter) {
      double f;
      
      p[indice * 3 + 2] += delta;
      f = func->f (p);
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();
  }

  for (int indice = 0; indice < nPols; indice++) {
    char nombre[256];
    double delta;

    cout << "Comenzando iteracion " << indice << "\n";
    
    sprintf(nombre, "funcionx_exacto_%d.dat", indice);
    out.open (nombre);
    for (delta = -0.1; delta < 0.1; delta += 0.2 / nIter) {
      double f;
      
      p[indice * 3] -= delta;
      //f = func->f (p);
      f = L (fL, p, fL->n_pols, MOCK_EXACTO);
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();
    
    sprintf(nombre, "funciony_exacto_%d.dat", indice);
    out.open (nombre);
    for (delta = -0.1; delta < 0.1; delta += 0.2 / nIter) {
      double f;
      
      p[indice * 3 + 1] += delta;
      //f = func->f (p);
      f = L (fL, p, fL->n_pols, MOCK_EXACTO); 
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();  
    
    sprintf(nombre, "funcionI_exacto_%d.dat", indice);
    out.open (nombre);
    for (delta = -1; delta < 1; delta += 2.0 / nIter) {
      double f;
      
      p[indice * 3 + 2] += delta;
      //f = func->f (p);
      f = L (fL, p, fL->n_pols, MOCK_EXACTO);
      igualar (p, p_ori, nPols);
      out << delta << "\t" << f << "\n";
      //printf ("%g\t%.16g\n", delta, f);
    }
    out.close ();
  }

}

void igualar (double *p, double *p_ori, int nPols) {
  for (int i = 0; i < nPols; i++) {
    p[3 * i]     = p_ori[3 * i];
    p[3 * i + 1] = p_ori[3 * i + 1];
    p[3 * i + 2] = p_ori[3 * i + 2];
  }
}
