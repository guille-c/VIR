#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <cstdlib>
#include <fstream.h>

#include "funcionMemCBI.h"
#include "funcionEntropiaNatural.h"
#include "funcionBayesVoronoiCBI.h"
#include "funcionChi2VoronoiImagen.h"
#include "incrementadorMallaImagen.h"
#include "inicializadorVoronoiUniforme.h"

#include "mcheck.h"

double calcularFuncionBVCBI(char *nombreImagen, char **nombresVis, char *nombreArchivo, 
			    int n_archivos, char *nombre_ori);
double calcularFuncionAjuste (char *nombreImagen, char *nombreArchivo);
double calcularFuncionMemCBI (char *nombreImagenOri, char **nombreVis, char *nombreImagen, int nVis);
double calcularChi2Im (char *nombre_ori, char *nombre_mod);

int main (int argc, char **argv){
  mtrace();
  if (argc < 5) {
    cerr << "uso: nombreMalla.dat nombreImagen.fits nombreVis.sub nombreOri.fits\n";
    //exit (1);
  }

  //calcularFuncionMemCBI(argv[1], argv + 2, argv[3], argc - 3);
  //exit(0);

  //for (int i = 0; i < 1e5; i++) {
  //cout << "========================================\ni = " << i << "\n";
  //FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[2], argv + 3, argc - 3, 1000, 1, -1);
  //delete f;
  //}
  calcularFuncionBVCBI(argv[2], argv + 3, argv[1], argc - 4, argv[4]);
  //calcularFuncionAjuste (argv[2], argv[3]);
}

double calcularFuncionMemCBI (char *nombreImagenOri, char **nombreVis, char *nombreImagen, int nVis) {
  struct image *im = do_read (nombreImagen);
  cout << "Im = " << nombreImagen << "\n";
  FuncionMemCBI *f = new FuncionMemCBI (nombreImagenOri, nombreVis, nVis, 
					new FuncionEntropiaNatural(im->npixels), 0.0, 1.0, 1e100, -1);

  //f->setRuido(0.000958957);

  double *p = new double [f->getNPars()];
  for (int i = 0; i < f->getNPars(); i++) {
    p[i] = im->pixels[i] / f->getRuido();
  }
  
  double ret = f->f(p);

  cout << "ndat = " << f->getNDat() << "\n";
  double n_dat = f->getNDat();
  
  cout << "n_dat = " << n_dat << "\n";
  cout << "L = " << ret << "\n";
  cout << "chi2 = " << f->getChi2() << "\n";
  cout << "chi2 reducido = " << f->getChi2() / (n_dat - f->getNPars()) << "\n";
  cout << "chi2 reducido = " << f->getChi2() / (n_dat) << "\n";
  cout << "S = " << f->getS() << "\n";
  cout << "P = " << exp (-ret) << "\n";

  return ret;
}

double calcularFuncionBVCBI(char *nombreImagen, char **nombresVis, char *nombreArchivo, 
			    int n_archivos, char *nombre_ori) {
  std::ifstream archivo (nombreArchivo);
  int n;
  double *p, ret;

  archivo >> n;
  n /= 3;
  cout << "n = " << n << "\n";
  p = new double [3*n];
  FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, n_archivos, 
							  n, 1, -1, 250.0, 250.0);
  for (int i = 0; i < n; i++) {
    archivo >> p[3*i] >> p[3*i + 1] >> p[3*i + 2];
    p[3*i + 2] /= f->getRuido();
    //cout << "p[" << i << "] = " << p[i] << "\n";
  }

  ret = f->f(p);
  f->guardarInfo("CalculoFuncion");
  funcL *fL = f->getFL();

  double n_dat = 2*fL->n_archivos * fL->header_obs[0]->nif * fL->header_obs[0]->nsamp;
  cout << "n_dat = " << n_dat << "\n";
  cout << "L = " << ret << "\n";
  cout << "chi2 = " << f->getFL()->chi2 << "\n";
  cout << "chi2 reducido = " << f->getFL()->chi2 / (n_dat - f->getNPars()) << "\n";
  cout << "chi2 reducido = " << f->getFL()->chi2 / (n_dat) << "\n";
  cout << "S = " << f->getFL()->S << "\n";
  cout << "P = " << exp (-ret) << "\n";
  cout << "Chi2Im = " << calcularChi2Im (nombreImagen, nombre_ori) << "\n";

  delete [] p;
  return ret;
}

double calcularFuncionAjuste (char *nombreImagen, char *nombreArchivo) {
  MallaVoronoi *m = newMallaVoronoi();
  struct image *im =  do_read (nombreImagen);
  std::ifstream archivo (nombreArchivo);
  int n, nx = im->size[0], ny = im->size[1];
  double ret, x, y, valor, **imagen;

  archivo >> n;
  for (int i = 0; i < n; i++) {
    archivo >> x >> y >> valor;
    insertarSitio (m, &x, &y, valor);
  }

  imagen = toImage (m, nx, ny, NULL);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      ret += pow (imagen [i][j] - im->pixels[i + j * nx], 2);
    }
    delete [] (imagen[i]);
  }
  delete [] imagen;

  double n_dat = nx * ny;
  cout << "n_dat = " << n_dat << "\n";
  cout << "chi2 = " << ret << "\n";
  cout << "chi2 reducido = " << ret / (n_dat - 3 * n) << "\n";

  delete_map (im);
  return ret;
}

double calcularChi2Im (char *nombre_ori, char *nombre_mod) {
  struct image *ori = do_read (nombre_ori);
  struct image *mod = do_read (nombre_mod);
  double ret = 0;

  if (ori->npixels != mod->npixels) {
    cerr << "ERROR en calcularChi2Im: != npixels";
    exit (1);
  }
  
  for (int i = 0; i < ori->npixels; i++) {
    ret += (ori->pixels[i] - mod->pixels[i]) * (ori->pixels[i] - mod->pixels[i]);
  }
  return ret;
}
