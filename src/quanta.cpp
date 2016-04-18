#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <cstdlib>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "mockcbiRoutines.h"

#include "funcionBayesVoronoiCBI.h"
#include "funcionMemCBI.h"
#include "funcionEntropiaNatural.h"
#include "funcionVoronoiI.h"
#include "funcionPbbVoronoiCBI.h"
#include "funcionNegativa.h"
#include "funcionDividida.h"
#include "funcionReconstructor.h"

#include "reconstructorGC.h"
#include "reconstructorGCMemCBI.h"

#include "inicializadorArchivo.h"
#include "inicializadorConstante.h"

#include "mcheck.h"

#define POW(x) ((x)*(x))

void comparacionFuncion (FuncionBayesVoronoiCBI *func);
void comparacionReconstructor (FuncionBayesVoronoiCBI *func);
void comparacionEntropia (char *nombreFunc, char *nombreMalla);
void comparacionQuantaVis (char *nombreFunc, char *nombreMalla);
void add_noise (struct uvf_header **header,   /* add noise to visibilities  */
		struct uvf_sample **samples, int n_archivos,
		long *seed);
double entropia (double *Ni, int n);

double cross_correlate_factor (struct uvf_header **header_obs,
			       struct uvf_sample **samples_obs, 
			       struct uvf_sample **samples_simu,
			       int n_archivos);

double quantaOptimo (struct image *fg_image,
		     struct uvf_header **header, 
		     struct uvf_sample **samples_obs,
		     int n_vis,
		     double quanta);
void noise_scale (struct uvf_header **header, struct uvf_sample **samples, 
		  int n_archivs, float ascale);
void copy_vis (struct uvf_header **header,
	       struct uvf_sample **samples_in, 
	       struct uvf_sample **samples_out,
	       int n_archivos);
void guardar (FuncionBayesVoronoiCBI *func, double *pars, 
	      double quanta, double noise_vis);
void reconstruir (FuncionMemCBI *funcMemCBI, int iter);
void funcvsQuanta (FuncionMemCBI *funcMemCBI);

void SvsQuanta (int argc, char **argv);
void DerivadasvsQuanta (int argc, char **argv);
void moverYComparaGradientes (FuncionBayesVoronoiCBI *f, double *p);
double modulo (double *vec, int n);

void imprimirinfo (char *nombreArchivo, double quanta, double valor);

int main (int argc, char **argv) {
  //DerivadasvsQuanta (argc, argv);
  //exit(0);
  SvsQuanta (argc, argv);
  exit (0);

  FuncionMemCBI *f = new FuncionMemCBI (argv[1], argv + 2, argc - 2, 
					new FuncionEntropiaNatural(64*64), 0.0, 0.0, 1e100);
  f->resizeImage (128, 128, 32, 32);

  funcvsQuanta (f);
  /*
  struct uvf_header **header_obs;  // Arreglo con los encabezados.
  struct uvf_sample **samples_mod; // Arreglo con las visibilidades
				   // modeladas.
  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (argv [1]);
  mtrace ();
  int i;

  funcL *fL = func->getFL();
  double quanta = Irms (fL->header_obs, fL->samples_obs, fL->n_archivos);

  cout << "Irms = " << quanta << "\n";

  comparacionQuantaVis (argv[1], argv[2]);

  //comparacionEntropia (argv [1], argv [2]);
  */
}

void SvsQuanta (int argc, char **argv) {
  std::ifstream archivo (argv[argc - 1]);
  int n;
  double *p, *I, q_ori, Itotal = 0, q, L;

  archivo >> n;
  FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[1], argv + 2, argc - 3, n/3, 1, -1, 250.0, 250.0);
  cout << "n = " << n << "\n";
  p = new double [n];
  I = new double [n/3];
  for (int i = 0; i < n / 3 ; i++) {
    archivo >> p[3*i] >> p[3*i + 1] >> I[i];
    p[3*i + 2] = I[i] / f->getRuido();
    Itotal += I[i];
  }
  
  double n_dat = f->getNDat();

  q_ori = f->getRuido();

  std::ofstream archivoOut ("L_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("Chi2_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("S_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("Chi2red1_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("Chi2red2_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("L_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("Chi2_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("S_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("Chi2red1_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("Chi2red2_vs_N.dat");
  archivoOut.close();

  for (double N = 1e-3; N < 2e3; N++) {
    if (((int) N) % 10 == 0) cout << "N = " << N << "\n";
    q = Itotal / N;
    f->setRuido (q);

    for (int i = 0; i < n / 3; i++) {
      p[3*i + 2] = I[i] / q;
    }
    L = f->f(p);
    if (((int) N) % 100 == 0) f->guardarInfo (i2s(N, 4));
    
    imprimirinfo ("L_vs_q.dat", q, L);
    imprimirinfo ("Chi2_vs_q.dat", q, f->getFL()->chi2);
    imprimirinfo ("S_vs_q.dat", q, f->getFL()->S);
    imprimirinfo ("Chi2red1_vs_q.dat", q, f->getFL()->chi2 / (n_dat - f->getNPars()));
    imprimirinfo ("Chi2red2_vs_q.dat", q, f->getFL()->chi2 / (n_dat));

    imprimirinfo ("L_vs_N.dat", N, L);
    imprimirinfo ("Chi2_vs_N.dat", N, f->getFL()->chi2);
    imprimirinfo ("S_vs_N.dat", N, f->getFL()->S);
    imprimirinfo ("Chi2red1_vs_N.dat", N, f->getFL()->chi2 / (n_dat - f->getNPars()));
    imprimirinfo ("Chi2red2_vs_N.dat", N, f->getFL()->chi2 / (n_dat));


    if (N == 1e-3) N = 0;
  }

  delete [] p;
  delete [] I;
  delete f;
}

void DerivadasvsQuanta (int argc, char **argv) {
  std::ifstream archivo (argv[argc - 1]);
  int n;
  double *p, *I, q_ori, Itotal = 0, q, L, *dchi2, *DS, 
    moddS, moddchi2, moddL, moddLdx = 0, moddLdy = 0, moddLdI = 0, moddchi2dI = 0, moddSdNi = 0;

  archivo >> n;
  FuncionBayesVoronoiCBI *f = new FuncionBayesVoronoiCBI (argv[1], argv + 2, argc - 3, n/3, 1, -1, 250.0, 250.0);
  cout << "n = " << n << "\n";
  p = new double [n];
  I = new double [n/3];
  dchi2 = new double [n];
  DS = new double [n];
  for (int i = 0; i < n / 3 ; i++) {
    archivo >> p[3*i] >> p[3*i + 1] >> I[i];
    p[3*i + 2] = I[i] / f->getRuido();
    Itotal += I[i];
  }
  
  moverYComparaGradientes (f, p);
  double n_dat = f->getNDat();

  q_ori = f->getRuido();

  printf ("N = %g, q = %g\n", Itotal/q_ori, q_ori);
  //exit (0);
  std::ofstream archivoOut ("dchi2_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dchi2_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dchi2dI_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dS_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dS_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dL_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dL_vs_N.dat");
  archivoOut.close();

  archivoOut.open ("dLdx_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dLdx_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dLdy_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dLdy_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dLdNi_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dLdNi_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dchi2dNi_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dchi2dNi_vs_N.dat");
  archivoOut.close();
  archivoOut.open ("dSdNi_vs_q.dat");
  archivoOut.close();
  archivoOut.open ("dSdNi_vs_N.dat");
  archivoOut.close();

  archivoOut.open ("I1_vs_I2.dat");
  archivoOut.close();

  funcL *fL = f->getFL();

  int cont = 0;
  for (q = q_ori / 10; q < q_ori * 10; q += q_ori / 10) {
    double N = Itotal / q, N2 = 0;
    if (((int) N) % 10 == 0) cout << "N = " << N << "\n";
    //for (double N = 1e-3; N < 2e2; N++) {
    //if (((int) N) % 10 == 0) cout << "N = " << N << "\n";
    //q = Itotal / N;
    f->setRuido (q);

    for (int i = 0; i < n / 3; i++) {
      p[3*i + 2] = I[i] / q;
      N2 += p[3*i + 2];
    }
    imprimirinfo ("I1_vs_I2.dat", N, N2);

    dS (fL, p, DS, n/3);
    f->df(p, dchi2);
    moddL = modulo (dchi2, n);

    imprimirinfo ("dL_vs_N.dat", N, moddL);
    imprimirinfo ("dL_vs_q.dat", q, moddL);

    for (int i = 0; i < n; i++) {
      if (i%3 == 0) moddLdx += dchi2[i]*dchi2[i];
      if (i%3 == 1) moddLdy += dchi2[i]*dchi2[i];
      if (i%3 == 2) moddLdI += dchi2[i]*dchi2[i];
      dchi2[i] -= DS[i];
      if (i%3 == 2) {moddchi2dI += dchi2[i]*dchi2[i]; moddSdNi += DS[i]*DS[i];}
    }
    moddLdx = sqrt(moddLdx);
    moddLdy = sqrt(moddLdy);
    moddLdI = sqrt(moddLdI);
    moddchi2dI = sqrt(moddchi2dI);
    moddSdNi = sqrt(moddSdNi);
    imprimirinfo ("dLdx_vs_N.dat", N, moddLdx);
    imprimirinfo ("dLdx_vs_q.dat", q, moddLdx);
    imprimirinfo ("dLdy_vs_N.dat", N, moddLdy);
    imprimirinfo ("dLdy_vs_q.dat", q, moddLdy);
    imprimirinfo ("dLdNi_vs_N.dat", N, moddLdI);
    imprimirinfo ("dLdNi_vs_q.dat", q, moddLdI);
    imprimirinfo ("dchi2dNi_vs_N.dat", N, moddchi2dI);
    imprimirinfo ("dchi2dNi_vs_q.dat", q, moddchi2dI);
    imprimirinfo ("dSdNi_vs_N.dat", N, moddSdNi);
    imprimirinfo ("dSdNi_vs_q.dat", q, moddSdNi);

    moddS = modulo (DS, n);
    moddchi2 = modulo (dchi2, n);

    if (cont % 10 == 0) f->guardarInfo (i2s(cont, 4)); cont++;
    //if (((int) N) % 100 == 0) f->guardarInfo (i2s(N, 4));
    
    imprimirinfo ("dchi2_vs_N.dat", N, moddchi2);
    imprimirinfo ("dchi2_vs_q.dat", q, moddchi2);
    imprimirinfo ("dchi2dI_vs_q.dat", q, moddchi2 * q);
    imprimirinfo ("dS_vs_N.dat", N, moddS);
    imprimirinfo ("dS_vs_q.dat", q, moddS);

    if (N == 1e-3) N = 0;
  }

  delete [] DS;
  delete [] dchi2;
  delete [] p;
  delete [] I;
  delete f;
}

void moverYComparaGradientes (FuncionBayesVoronoiCBI *f, double *p){
  std::ofstream archivoOut ("comparacionGradientes.dat");
  int n = f->getFL()->n_pols * 3, nx = f->getFL()->fg_image->size[0];
  double dL1, dL2, dL3;
  double *p2 = new double [n], *p3 = new double [n];
  double *g1 = new double [n], *g2 = new double [n], *g3 = new double [n];

  for (int i = 0; i < n / 3; i++) {
    p2[3 * i]     = p[3 * i] +  1.0 / (nx - 1);
    p2[3 * i + 1] = p[3 * i + 1] +  1.0 / (nx - 1);
    p2[3 * i + 2] = p[3 * i + 2] + 1;
  }
  for (int i = 0; i < n; i++) {
    p3[i] = p[i];
  }
  for (int i = 0; i < n; i++) {
    p3[i] = p2[i];
    g3[i] = (f->f(p3) - f->f(p)) / (p3[i] - p[i]);
    //cout << i << " " << f->f(p3) - f->f(p) << " vs " << p3[i] - p[i] << "\n";
    p3[i] = p[i];
  }
  
  f->df(p, g1);
  f->df(p2, g2);
  
  dL1 = modulo (g1, n);
  dL2 = modulo (g2, n);
  dL3 = modulo (g3, n);
  
  archivoOut << "dL1 = " << dL1 << ", dL2 = " << dL2 << ", dL3 = " << dL3 << "\n\n";
  for (int i = 0; i < n / 3; i++) {
    archivoOut << g1[i] << "\t" << g2[i] << "\t" << g3[i] << "\n";
  }
  archivoOut.close();
}

double modulo (double *vec, int n){
  double suma = 0;
  for (int i = 0; i < n; i++) {
    suma += vec[i] * vec[i];
  }
  return sqrt (suma);
}

void funcvsQuanta (FuncionMemCBI *f) {
  std::ofstream archivo ("L_vs_q_media.dat");
  archivo.close();
  archivo.open ("Chi2_vs_q_media.dat");
  archivo.close();
  archivo.open ("S_vs_q_media.dat");
  archivo.close();
  archivo.open ("Chi2red1_vs_q_media.dat");
  archivo.close();
  archivo.open ("Chi2red2_vs_q_media.dat");
  archivo.close();

  archivo.open ("L_vs_q_escalado.dat");
  archivo.close();
  archivo.open ("Chi2_vs_q_escalado.dat");
  archivo.close();
  archivo.open ("S_vs_q_escalado.dat");
  archivo.close();
  archivo.open ("Chi2red1_vs_q_escalado.dat");
  archivo.close();
  archivo.open ("Chi2red2_vs_q_escalado.dat");
  archivo.close();

  archivo.open ("L_vs_q_rerec.dat");
  archivo.close();
  archivo.open ("Chi2_vs_q_rerec.dat");
  archivo.close();
  archivo.open ("S_vs_q_rerec.dat");
  archivo.close();
  archivo.open ("Chi2red1_vs_q_rerec.dat");
  archivo.close();
  archivo.open ("Chi2red2_vs_q_rerec.dat");
  archivo.close();

  double ruido_ori = f->getRuido();
  //for (double factor = 1e-3; factor <= 1e3; factor *= sqrt (10.0)) {
  for (int i = -4; i <= 6; i++) {
    double factor = pow (sqrt(10.0), i);
    cout << "factor = " << factor << "\n";
    f->setRuido (ruido_ori * factor);
    reconstruir (f, 4 + i);
  }
}

void reconstruir (FuncionMemCBI *funcMemCBI, int iter) {
    ReconstructorGCMemCBI *r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, 0.0008628 / funcMemCBI->getRuido()); /* Para usar Nmax */

    cout << "npars rec = " << funcMemCBI->getNPars() << "\n";
    r1->run();
    double n_dat = funcMemCBI->getNDat();

    imprimirinfo ("L_vs_q_media.dat", funcMemCBI->getRuido(), funcMemCBI->f(r1->getPars()));
    imprimirinfo ("Chi2_vs_q_media.dat", funcMemCBI->getRuido(), funcMemCBI->getChi2());
    imprimirinfo ("S_vs_q_media.dat", funcMemCBI->getRuido(), funcMemCBI->getS());
    imprimirinfo ("Chi2red1_vs_q_media.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / (n_dat - funcMemCBI->getNPars()));
    imprimirinfo ("Chi2red2_vs_q_media.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / n_dat);

    funcMemCBI->guardarInfo ("Final1_" + i2s (iter, 3));
    double pMin = 1e300;
    double *pr1 = r1->getPars();
    for (int i = 0; i < r1->getNPars(); i++) {
      if (pr1[i] < pMin) {
	pMin = pr1[i];
      }
    }
    
    for (int i = 0; i < r1->getNPars(); i++) {
      pr1[i] = pr1[i] - pMin;
    }

    funcMemCBI->scalePars(r1->getPars());
    funcMemCBI->guardarInfo ("Escalado" + i2s (iter, 3));

    imprimirinfo ("L_vs_q_escalado.dat", funcMemCBI->getRuido(), funcMemCBI->f(r1->getPars()));
    imprimirinfo ("Chi2_vs_q_escalado.dat", funcMemCBI->getRuido(), funcMemCBI->getChi2());
    imprimirinfo ("S_vs_q_escalado.dat", funcMemCBI->getRuido(), funcMemCBI->getS());
    imprimirinfo ("Chi2red1_vs_q_escalado.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / (n_dat - funcMemCBI->getNPars()));
    imprimirinfo ("Chi2red2_vs_q_escalado.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / n_dat);

    r1->run ();

    imprimirinfo ("L_vs_q_rerec.dat", funcMemCBI->getRuido(), funcMemCBI->f(r1->getPars()));
    imprimirinfo ("Chi2_vs_q_rerec.dat", funcMemCBI->getRuido(), funcMemCBI->getChi2());
    imprimirinfo ("S_vs_q_rerec.dat", funcMemCBI->getRuido(), funcMemCBI->getS());
    imprimirinfo ("Chi2red1_vs_q_rerec.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / (n_dat - funcMemCBI->getNPars()));
    imprimirinfo ("Chi2red2_vs_q_rerec.dat", funcMemCBI->getRuido(), 
		  funcMemCBI->getChi2() / n_dat);
    
    funcMemCBI->guardarInfo ("ReRec" + i2s (iter, 3));
    delete r1;
}

void imprimirinfo (char *nombreArchivo, double quanta, double valor) {
  std::ofstream archivo;

  archivo.open (nombreArchivo, ios::app);
  archivo << quanta << "\t" << valor << "\n";

  archivo.close();
}

void comparacionQuantaVis (char *nombreFunc, char *nombreMalla) {
  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (nombreFunc);
  funcL *fL = func->getFL();
  double quanta = 0, ret, quanta_ori;
  struct uvf_sample **samples_ori; // Arreglo con las visibilidades
  Inicializador *ini = new InicializadorArchivo(nombreMalla);
  ReconstructorGC *r;
  std::ofstream archivo ("L_vs_quanta.dat");
  archivo.close();
  archivo.open ("L_vs_noise_vis.dat");
  archivo.close();
  archivo.open ("L_vs_N.dat");
  archivo.close();
  
  archivo.open ("chi2_vs_quanta.dat");
  archivo.close();
  archivo.open ("chi2_vs_noise_vis.dat");
  archivo.close();
  archivo.open ("chi2_vs_N.dat");
  archivo.close();

  archivo.open ("S_vs_quanta.dat");
  archivo.close();
  archivo.open ("S_vs_noise_vis.dat");
  archivo.close();
  archivo.open ("S_vs_N.dat");
  archivo.close();

  archivo.open ("dchi2_vs_quanta.dat");
  archivo.close();
  archivo.open ("dchi2_vs_noise_vis.dat");
  archivo.close();
  archivo.open ("dchi2_vs_N.dat");
  archivo.close();

  archivo.open ("dS_vs_quanta.dat");
  archivo.close();
  archivo.open ("dS_vs_noise_vis.dat");
  archivo.close();
  archivo.open ("dS_vs_N.dat");
  archivo.close();

  archivo.open ("quanta_vs_noise_vis.dat");
  archivo.close();

  double *pars = new double [func->getNPars()];
  struct image *im_ori = do_read (fL->nombreFits);

  samples_ori = (struct uvf_sample **) malloc(fL->n_archivos * 
					      sizeof(struct uvf_sample *));
  quanta_ori = func->getRuido();

  for(int i = 0; i < fL->n_archivos; i++) {
    // abrimos los archivos de visibilidades
    int status = do_read_uvdata(fL->infile[i], &(fL->header_obs[i]), 
				&(samples_ori[i]));
    if (status != SUCCESS) {
      printf("Error con el archivo uvf\n");
      exit(1);
    }
  }

  do_write_fits (im_ori, "!fLfgImage.fits");

  for (int arch = 0; arch < fL->n_archivos; arch++) {
    mockcbi_sub (fL->cmb_image, im_ori,
		 fL->header_obs[arch], fL->samples_obs[arch], samples_ori[arch],
		 fL->beam);
  }

  long seed = 0;
  int cont = 0;
  for (double noiseScale = 1e-2; noiseScale < 1; noiseScale += 1e-2) {
    cont++;
    r = new ReconstructorGC (func, ini, 1e-10);
    copy_vis (fL->header_obs, samples_ori, fL->samples_mod, fL->n_archivos);
    noise_scale (fL->header_obs, fL->samples_mod, fL->n_archivos, noiseScale);
    add_noise (fL->header_obs, fL->samples_mod, fL->n_archivos, &seed);
    quanta = Irms (fL->header_obs, fL->samples_mod, fL->n_archivos);
    func->setRuido (quanta);
    for (int i = 0; i < func->getNPars() / 3; i++) {
      r->setPar (3 * i + 2, r->getPar (3 * i + 2) / quanta_ori);
    }
    ret = r->run ();
    for (int i = 0; i < func->getNPars(); i++) {
      pars[i] = r->getPar(i);
    }
    guardar (func, pars, quanta, noiseScale);
    func->guardarInfo ("Quanta" + i2s (cont, 5));
    delete r;
  }
}

void guardar (FuncionBayesVoronoiCBI *func, double *pars, 
	      double quanta, double noise_vis) {
  double f = func->f (pars), S = func->getFL()->S, chi2 = func->getFL()->chi2, N = 0;
  funcL *fL = func->getFL();

  for (int i = 0; i < func->getNPars(); i++) {
    N += pars[3 * i + 2] / func->getRuido();
  }
  
  // --------------- L ------------------

  std::ofstream archivo;
  archivo.open ("L_vs_quanta.dat", ios::app);
  archivo << quanta << "\t" << f << "\n";    
  archivo.close();    
  
  archivo.open ("L_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << f << "\n";    
  archivo.close();    
  
  archivo.open ("L_vs_N.dat", ios::app);
  archivo << N << "\t" << f << "\n";    
  archivo.close();    
    
  // --------------- chi2 ------------------

  archivo.open ("chi2_vs_quanta.dat", ios::app);
  archivo << quanta << "\t" << chi2 << "\n";
  archivo.close();
  
  archivo.open ("chi2_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << chi2 << "\n";
  archivo.close();
  
  archivo.open ("chi2_vs_N.dat", ios::app);
  archivo << N << "\t" << chi2 << "\n";
  archivo.close();

  // --------------- S ------------------

  archivo.open ("S_vs_quanta.dat", ios::app);
  archivo << quanta << "\t" << S << "\n";
  archivo.close();
  
  archivo.open ("S_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << S << "\n";
  archivo.close();
  
  archivo.open ("S_vs_N.dat", ios::app);
  archivo << N << "\t" << S << "\n";
  archivo.close();
    
  // --------------- dchi2 ------------------

  double *grad = new double [func->getNPars()];
  double norma;
  for (int i = 0; i < func->getNPars(); i++) {
    grad[i] = 0;
  }
  fL->entropia = 0;
  dL (fL, pars, grad, fL->n_pols, 0, MOCKCBI);
  for (int i = 0; i < fL->n_pols; i++) {
    norma += grad[3 * i + 2] * grad[3 * i + 2];
  }
  //norma /= fL->n_pols;
  norma = sqrt (norma);
  fL->entropia = 1;
  
  archivo.open ("dchi2_vs_quanta.dat", ios::app);
  archivo << quanta << "\t" << norma << "\n";
  archivo.close();
  archivo.open ("dchi2_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << norma << "\n";
  archivo.close();
  archivo.open ("dchi2_vs_N.dat", ios::app);
  archivo << N << "\t" << norma << "\n";
  archivo.close();

  // --------------- dS ------------------

  for (int i = 0; i < func->getNPars(); i++) {
    grad[i] = 0;
  }
  dS (fL, pars, grad, fL->n_pols);
  for (int i = 0; i < fL->n_pols; i++) {
    norma += grad[3 * i + 2] * grad[3 * i + 2];
  }
  //norma /= fL->n_pols;
  norma = sqrt (norma);

  archivo.open ("dS_vs_quanta.dat", ios::app);
  archivo << quanta << "\t" << norma << "\n";
  archivo.close();
  archivo.open ("dS_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << norma << "\n";
  archivo.close();
  archivo.open ("dS_vs_N.dat", ios::app);
  archivo << N << "\t" << norma << "\n";
  archivo.close();

  archivo.open ("quanta_vs_noise_vis.dat", ios::app);
  archivo << noise_vis << "\t" << quanta << "\n";
  archivo.close();

}

void comparacionEntropia (char *nombreFunc, char *nombreMalla) {
  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (nombreFunc);
  Inicializador *ini = new InicializadorArchivo(nombreMalla);
  int n = func->getNPars() / 3;
  double *Ni = new double [n];
  double *p = new double [func->getNPars()];
  ini->inicializar (p, n);
  std::ofstream archivo ("S_vs_quanta.dat");
  archivo.close();
  archivo.open ("S_vs_N.dat");
  archivo.close();
  double quanta, S, quanta_old, N;

  for (int i = 0; i < n; i++) {
    Ni[i] = p[3 * i + 2] / func->getRuido();
  }

  quanta_old = func->getRuido();
  for (quanta = 1e-6; quanta < 1e3; quanta *= 2) {
    N = 0;
    for (int i = 0; i < n; i++) {
      Ni[i] = Ni[i] * quanta_old / quanta;
      N += Ni[i];
    }
    S = entropia (Ni, n);
    archivo.open ("S_vs_quanta.dat", ios::app);
    archivo << quanta << "\t" << -S << "\n";    
    archivo.close();    
    archivo.open ("S_vs_N.dat", ios::app);
    archivo << N << "\t" << -S << "\n";    
    archivo.close();    
    quanta_old = quanta;
  }
}

double entropia (double *Ni, int n) {
  double N = 0, S = 0;

  for (int i = 0; i < n; i++) {
    N += Ni [i];
    S -= lgamma (Ni[i] + 1);
  }
  S += lgamma (N + 1) - N * log ((double)n);
  return S;
}

void comparacionFuncion (FuncionBayesVoronoiCBI *func) {
  funcL *fL;
  double ftol = 1e-10, quanta, ret;
  double *p = new double [func->getNPars()];
  std::ofstream archivo ("f_vs_quanta.dat");
  archivo.close();
  archivo.open ("chi_vs_quanta.dat");
  archivo.close();
  archivo.open ("S_vs_quanta.dat");
  archivo.close();
  int i = 0;

  func->f(p);
  fL = func->getFL();
  quanta =  quantaOptimo (fL->fg_image, fL->header_obs, fL->samples_obs, 
			  fL->n_archivos, func->getRuido());
  cout << "quanta = " << quanta << ", quanta func = " << func->getRuido() << "\n";

  for (quanta = 1e-6; quanta < 1e-3; quanta += 1e-5) {
  //for (quanta = 1e-6; quanta < 1e3; quanta *= 2 ) {
    func->setRuido (quanta);
    ret = func->f(p);
    //ret = r->run();
    func->guardarInfo ("Quanta" + i2s (i, 5));

    archivo.open ("f_vs_quanta.dat", ios::app);
    archivo << quanta << "\t" << ret << "\n";    
    cout << "quanta = " << quanta << ", func = " << ret << "\n";
    archivo.close();

    archivo.open ("chi_vs_quanta.dat", ios::app);
    archivo << quanta << "\t" << fL->chi2 << "\n";    
    archivo.close();

    archivo.open ("S_vs_quanta.dat", ios::app);
    archivo << quanta << "\t" << fL->S << "\n";    
    archivo.close();
    i++;
  }
  delete [] p;
}

void comparacionReconstructor (FuncionBayesVoronoiCBI *func, char *nombreMalla) {
  funcL *fL;
  int i;
  double ret, quanta;
  //Inicializador *ini = new InicializadorArchivo(argv[2]);
  Inicializador *ini = new InicializadorArchivo(nombreMalla);
  ReconstructorGC *r = new ReconstructorGC (func, ini, 1e-10);
  std::ofstream archivo ("f_vs_quanta.dat");
  archivo.close();
  double *p = new double [func->getNPars()];

  fL = func->getFL();
  for (i = 0; i < func->getNPars() / 3; i++) {
    //cout << r->getPar (i) << " / "  << func->getRuido() << " = " << r->getPar (i) / func->getRuido() << "\n";
    r->setPar (3 * i + 2, r->getPar (3 * i + 2) / func->getRuido());
    p[3 * i] = r->getPar(3 * i);
    p[3 * i + 1] = r->getPar(3 * i + 1);
    p[3 * i + 2] = r->getPar(3 * i + 2);
  }


  i = 0;
  r->run();
  for (int i = 0; i < 1e4; i++) {
    fL = func->getFL();
    quanta =  quantaOptimo (fL->fg_image, fL->header_obs, fL->samples_obs, 
			    fL->n_archivos, func->getRuido());
    func->setRuido (quanta);
    ret = r->run();
    
    func->guardarInfo ("Quanta" + i2s (i, 5));

    archivo.open ("f_vs_quanta.dat", ios::app);
    archivo << quanta << "\t" << ret << "\n";    
    cout << "quanta = " << quanta << ", func = " << ret << "\n";
    archivo.close();
  }
  delete r;
  delete [] p;
}

double quantaOptimo (struct image *fg_image,
		     struct uvf_header **header, 
		     struct uvf_sample **samples_obs,
		     int n_vis,
		     double quanta) {
  double c1 = 0, c2 = 0;
  struct image *cmb_image, *fg_image_aux;
  struct uvf_sample *samples_mod;
  struct pbeam beam;
  
  init_beam (&beam);
  beam.type = CBI;

  cmb_image = new_map();
  fg_image_aux = new_map();
  copy_empty_map (cmb_image, fg_image);
  copy_empty_map (fg_image_aux, fg_image);
  for (int i = 0; i < fg_image->npixels; i++){
    fg_image_aux->pixels[i] = fg_image->pixels[i] / quanta;
  }
  
  for (int arch = 0; arch < n_vis; arch++) {
    samples_mod = (struct uvf_sample*) malloc((header[arch]->nsamp) * 
					      sizeof(struct uvf_sample));
    if (samples_mod == NULL) {
      fprintf (stderr, "ERROR en quantaOptimo, fL->samples_mod[%d] = NULL.\n", arch);
      exit (1);
    }
    
    mockcbi_sub (cmb_image, fg_image_aux,
		 header[arch], samples_obs[arch], 
		 samples_mod, beam);
    for (int iff = 0; iff < header[arch]->nif; iff++) {
      for (int samp = 0; samp < header[arch]->nsamp; samp++) {
	c1 += ((samples_mod[samp].rdata[iff * 3] * 
		samples_obs[arch][samp].rdata[iff * 3] + 
		samples_mod[samp].rdata[iff * 3 + 1] * 
		samples_obs[arch][samp].rdata[iff * 3 + 1]) *
	       samples_obs[arch][samp].rdata[iff * 3 + 2]);
	       
	c2 += ((POW(samples_mod[samp].rdata[iff * 3]) +
		POW(samples_mod[samp].rdata[iff * 3 + 1])) *
	       samples_obs[arch][samp].rdata[iff * 3 + 2]);
//	ret += ((SQR(fL->samples_mod[arch][samp].rdata[iff * 3] 
//		     - fL->samples_obs[arch][samp].rdata[iff * 3])
//		 + SQR(fL->samples_mod[arch][samp].rdata[iff * 3 + 1] 
//		       - fL->samples_obs[arch][samp].rdata[iff * 3 + 1]))
//		* fL->samples_mod[arch][samp].rdata[iff * 3 + 2]);
      }
    }
    free (samples_mod);
  }


  delete_map (cmb_image);
  delete_map (fg_image_aux);

  //cout << "c1 = " << c1 << ", c2 = " << c2 << "\n";

  return c1 / c2;
}

/*--------------------------------------------------------------------
 * Calculo del factor de correlacion entre visibilidades de samples_31
 * y samples_simu.
 *--------------------------------------------------------------------*/

double cross_correlate_factor (struct uvf_header **header_31,
			       struct uvf_sample **samples_31, 
			       struct uvf_sample **samples_simu,
			       int n_archivos) {
  int arch, samp, chan;
  double c1 = 0, c2 = 0;

  for (arch = 0; arch < n_archivos; arch++) {
    for (chan = 0; chan < header_31[arch]->nchan; chan++) {
      for (samp = 0; samp < header_31[arch]->nsamp; samp++) {
	c1 += ((samples_simu[arch][samp].rdata[chan * 3] * 
		samples_31[arch][samp].rdata[chan * 3] + 
		samples_simu[arch][samp].rdata[chan * 3 + 1] * 
		samples_31[arch][samp].rdata[chan * 3 + 1]) *
	       samples_31[arch][samp].rdata[chan * 3 + 2]);
	       
	c2 += ((POW(samples_simu[arch][samp].rdata[chan * 3]) +
		POW(samples_simu[arch][samp].rdata[chan * 3 + 1])) *
	       samples_31[arch][samp].rdata[chan * 3 + 2]);
      }
    }
  }
  cout << "c1 = " << c1 << ", c2 = " << c2 << "\n";
  return c1 / c2;
}

/*--------------------------------------------------------------------
 * scale input  visibilities noise
 *--------------------------------------------------------------------*/

void noise_scale (struct uvf_header **header, struct uvf_sample **samples, 
		 int n_archivos, float ascale) {
  int arch, i, k, nsamp, nchan;
  printf ("scaling vis by %g \n",ascale);
  
  for (arch = 0; arch < n_archivos; arch++) {
    nsamp = header[arch]->nsamp;
    nchan = header[arch]->nchan;

    for (i = 0; i < nsamp; i++) {      /* Loop through samples */
      for (k = 0; k < nchan; k++) {    /* Loop through channels */
	samples[arch][i].rdata[3*k + 2] /= (ascale * ascale); 
      }
    }
  }
}

/*--------------------------------------------------------------------
 * Copia Visibilidades de in a out.
 *--------------------------------------------------------------------*/

void copy_vis (struct uvf_header **header,
	       struct uvf_sample **samples_in, 
	       struct uvf_sample **samples_out,
	       int n_archivos) {
  int arch, samp, chan;

  for (arch = 0; arch < n_archivos; arch++) {
    for (chan = 0; chan < header[arch]->nchan; chan++) {
      for (samp = 0; samp < header[arch]->nsamp; samp++) {
	samples_out[arch][samp].rdata[chan * 3] = 
	  samples_in[arch][samp].rdata[chan * 3];

	samples_out[arch][samp].rdata[chan * 3 + 1] = 
	  samples_in[arch][samp].rdata[chan * 3 + 1];

	samples_out[arch][samp].rdata[chan * 3 + 2] = 
	  samples_in[arch][samp].rdata[chan * 3 + 2];
      }
    }
  }
}

/*--------------------------------------------------------------------
 * add noise to visibilities.
 *--------------------------------------------------------------------*/

void add_noise(struct uvf_header **header, struct uvf_sample **samples,
		      int n_archivos, long *seed) {
  int arch, i, k, nsamp, nchan;
  float sum =0; 
  //float target_Chi2 = 31291. ;
  double sigma;
  double visnoise;
  float scale_factor;
  double sum_VR = 0, sum_VI = 0, sum_noise = 0, noise;

  for (arch = 0; arch < n_archivos; arch++) {
    nsamp = header[arch]->nsamp;
    nchan = header[arch]->nchan;
    for (i = 0; i < nsamp; i++) {      /* Loop through samples */
      for (k = 0; k < nchan; k++) {    /* Loop through channels */
	if (samples[arch][i].rdata[3*k+2]  != 0) {
	  visnoise = sqrt(1/samples[arch][i].rdata[3*k+2]);
	  noise = gasdev(seed)*visnoise;
	  
	  sum_noise += fabs(noise);
	  sum_VR += fabs(samples[arch][i].rdata[3*k]);
	  sum_VI += fabs(samples[arch][i].rdata[3*k+1]);
	  
	  samples[arch][i].rdata[3*k] += gasdev(seed)*visnoise;   // DEVTEST
	  samples[arch][i].rdata[3*k+1] += gasdev(seed)*visnoise;   // DEVTEST
	}    /* weights */
      }
    }
  }
  printf ("suma = %g\t%g,\t noise = %g\n", 
	  sum_VR, sum_VI, sum_noise);
}

