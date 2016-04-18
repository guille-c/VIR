#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "image_routines.h"
#include "uvsubs.h"
#include "mockcbiRoutines.h"

#include "funcionCBI.h"
#include "funcionMemCBI.h"
#include "funcionCero.h"
#include "funcionEntropiaNatural.h"
#include "funcionChi2VoronoiImagen.h"
#include "funcionBayesVoronoiCBI.h"
#include "reconstructorGCMemCBI.h"
#include "reconstructorIncremental.h"

#include "incrementadorMallaImagen.h"

#include "inicializadorConstante.h"
#include "inicializadorVoronoiUniforme.h"

#include "mcheck.h"

double *optimizar (char *nombreImagen, char **nombresVis, int n_archivos, int nIni, int nFin, int dn,
		   int iter, double qFactor, long *seed);
static void add_noise(struct uvf_header *header,
		      struct uvf_sample *samples, long *seed);
void escalarImagenRuido (struct image *im, struct uvf_header **header, 
			 struct uvf_sample **samples, int n_archivos, double n_ruido);
void leerVisibilidades (char **nombresVis, struct uvf_header ***header, struct uvf_sample ***samples_obs, 
			struct uvf_sample ***samples_mod, int n_archivos);

int main (int argc, char **argv) {
  char *nombreImagen = new char [255], *dato = new char [512], **nombresVis, **nombresVisAux;
  double nsr, qIni, qFin, dq, q;
  int nVis, nIni, nFin, iterIni, nIter, dn;
  mtrace ();

  std::ifstream archivoIn (argv[1]);
  archivoIn >> dato >> nombreImagen;
  cout << dato << "\t" << nombreImagen << "\n";
  archivoIn >> dato >> nsr;
  cout << dato << "\t" << nsr << "\n";

// MTRACE
//   struct image *im1 = do_read (nombreImagen);
//   resize_map (im1, 1024, 1024, im1->crpix[0], im1->crpix[1]); 
//   delete_map(im1);
//   delete [] nombreImagen;
//   delete [] dato;
//   exit(0);

  archivoIn >> dato >> nIni;
  cout << dato << "\t" << nIni << "\n";
  archivoIn >> dato >> nFin;
  cout << dato << "\t" << nFin << "\n";
  archivoIn >> dato >> dn;
  cout << dato << "\t" << dn << "\n";
  archivoIn >> dato >> iterIni;
  cout << dato << "\t" << iterIni << "\n";
  archivoIn >> dato >> nIter;
  cout << dato << "\t" << nIter << "\n";

  archivoIn >> dato >> qIni;
  cout << dato << "\t" << qIni << "\n";
  archivoIn >> dato >> qFin;
  cout << dato << "\t" << qFin << "\n";
  archivoIn >> dato >> dq;
  cout << dato << "\t" << dq << "\n";

  archivoIn >> dato >> nVis;
  cout << dato << "\t" << nVis << "\n";
  
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    archivoIn >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  archivoIn.close();

  struct image *im = do_read (nombreImagen);
  struct uvf_header **header;
  struct uvf_sample **samples_obs, **samples_mod;
  leerVisibilidades (nombresVis, &header, &samples_obs, &samples_mod, nVis);

  //escalarImagenRuido (im, header, samples_obs, nVis, nsr);
  nombreImagen = "!imagenAux.fits";
  do_write_fits (im, nombreImagen);
  nombreImagen = "imagenAux.fits";
  nombresVisAux = new char *[nVis];
  for (int i = 0;i < nVis; i++) {
    nombresVisAux[i] = new char[512];
    sprintf (nombresVisAux[i], "!visAux%d.sub", i);
    do_write ("write uvf", nombresVisAux[i], header[i], samples_obs[i], nombresVis[i]);
    sprintf (nombresVisAux[i], "visAux%d.sub", i);
  }

  delete_map (im);
  for (int i = 0; i < nVis; i++){
    delete [] nombresVis[i]; 
    free (samples_obs[i]);
    free (header[i]);
    free (samples_mod[i]);
  }
  delete [] nombresVis;

  free (header);
  free (samples_obs);
  free (samples_mod);


  char *nombreArchivo = new char[255];
  std::ofstream archivoOut;
  double *L, *Lprom = new double [nFin - nIni];
  //long seed = -labs((long)time(0));
  long seed = 0;

  for (q = qIni; q <= qFin; q *= dq) {
    for (int j = nIni; j < nFin; j++) {
      Lprom[j - nIni] = 0;
    }
    int i = 0;

    // Leemos archivos en caso de que se comienza de una iteracion != 0
    if (iterIni > 0) {
      for (i = 0; i < iterIni; i++) {
	int jLeido;
	L = new double [nFin - nIni];
	sprintf (nombreArchivo, "LvsnPols%d_q%g.dat", i, q);
	cout << "Leyendo " << nombreArchivo << "\n";
	archivoIn.open(nombreArchivo);
	for (int j = nIni; j < nFin; j += dn) {
	  archivoIn >> jLeido >> L[j - nIni]; 
	  if (jLeido != j) {
	    cerr << "jLeido != j: " << jLeido << " != " << j << "\n";
	    exit (1);
	  }
	  Lprom[j - nIni] += L [j - nIni];
	}
	archivoIn.close();
	delete [] L;
      
	// Imprimimos el promedio parcial.
	sprintf (nombreArchivo, "LvsnPols_parcial_Prom_q%g_i%d.dat", q, i);
	archivoOut.open (nombreArchivo);
	//archivoOut.open ("LvsnPolsProm.dat");
	for (int j = nIni; j < nFin; j += dn) {
	  archivoOut << j << "\t" <<  Lprom[j - nIni] / (i + 1) << "\n";
	}
	archivoOut.close ();
      
      }
      iterIni = 0;
    }
    
    for (; i < nIter; i++) {
      L = optimizar (nombreImagen, nombresVisAux, nVis, nIni, nFin, dn, i, q, &seed);
      sprintf (nombreArchivo, "LvsnPols%d_q%g.dat", i, q);
      cout << "Imprimiendo en " << nombreArchivo << "\n";
      archivoOut.open(nombreArchivo);
      for (int j = nIni; j < nFin; j += dn) {
	Lprom[j - nIni] += L [j - nIni];
	archivoOut << j << "\t" << L [j - nIni] << "\n";
      }
      archivoOut.close();
      delete [] L;

      // Imprimimos el promedio parcial.
      sprintf (nombreArchivo, "LvsnPols_parcial_Prom_q%g_i%d.dat", q, i);
      archivoOut.open (nombreArchivo);
      //archivoOut.open ("LvsnPolsProm.dat");
      for (int j = nIni; j < nFin; j += dn) {
	archivoOut << j << "\t" <<  Lprom[j - nIni] / (i + 1) << "\n";
      }
      archivoOut.close ();
    }
    
    sprintf (nombreArchivo, "LvsnPolsProm_q%g.dat", q);
    archivoOut.open (nombreArchivo);
    //archivoOut.open ("LvsnPolsProm.dat");
    for (int j = nIni; j < nFin; j += dn) {
      Lprom[j - nIni] /= nIter;
      archivoOut << j << "\t" <<  Lprom[j - nIni] << "\n";
    }
    archivoOut.close ();
  }
  for (int i = 0; i < nVis; i++){
    delete [] nombresVisAux[i]; 
  }
  delete [] nombresVisAux;
  delete [] Lprom;
  delete [] nombreArchivo;
}

double *optimizar (char *nombreImagen, char **nombresVis, int n_archivos, int nIni, int nFin, 
		   int dn, int iter, double qFactor, long *seed) {
  if (nFin <= nIni) {
    cerr << "ERROR en optimizar: nFin = " << nFin << " <= nIni = " << nIni << "\n";
    exit (0);
  }

  double *ret = new double [nFin - nIni];

  struct image *im = do_read (nombreImagen);
  struct image *cmb_im = do_read (nombreImagen);
  struct uvf_header **header;
  struct uvf_sample **samples_obs, **samples_mod;
  int nx, ny, n1x, n1y;
  double lx = 32, ly = 32;
  n1x = round (fabs (lx / (im->cdelt[0] * 60)));
  n1y = round (fabs(ly / (im->cdelt[1] * 60)));
  nx = round (pow (2, floor (log ((double) n1x) / log (2.0)) + 1));
  ny = round (pow (2, floor (log ((double) n1y) / log (2.0)) + 1));

  leerVisibilidades (nombresVis, &header, &samples_obs, &samples_mod, n_archivos);

  for (int i = 0; i < cmb_im->npixels; i++) {
    cmb_im->pixels[i] = 0;
  }

  resize_map (im, nx, ny, im->crpix[0], im->crpix[1]); 
  resize_map (cmb_im, nx, ny, cmb_im->crpix[0], cmb_im->crpix[1]);

  struct pbeam beam;
  init_beam (&beam);
  beam.type = CBI;

  char **nombresVisAux = new char *[n_archivos];
  for (int i = 0; i < n_archivos; i++) {
    mockcbi_sub (cmb_im, im, header[i], samples_obs[i], samples_mod[i], beam);
    add_noise (header[i], samples_mod[i], seed);

    nombresVisAux[i] = new char[512];
    sprintf (nombresVisAux[i], "!visRuidoIter%d_%d.sub", iter, i);
    do_write ("write uvf", nombresVisAux[i], header[i], samples_mod[i], nombresVis[i]);
    sprintf (nombresVisAux[i], "visRuidoIter%d_%d.sub", iter, i);
  }
  do_write_fits (im, "!mock.fits");

  for (int i = 0; i < n_archivos; i++) {
    free(samples_obs[i]);
    free(header[i]);
    free(samples_mod[i]);
  }
  free (header);
  free (samples_obs);
  free (samples_mod);
  delete_map (im);
  delete_map (cmb_im);

  //// TEST DIRECTO
  //funcL * fL1 = newFuncL (nombreImagen, nombresVis, n_archivos, 
  //			  nIni, 1, -1, 128.0, 128.0);
  //fL1->fg_image = do_read (nombreImagen);
  //
  //resize_map (fL1->fg_image,  nx, ny, fL1->fg_image->crpix[0],  fL1->fg_image->crpix[1]); 
  //resize_map (fL1->cmb_image, nx, ny, fL1->cmb_image->crpix[0], fL1->cmb_image->crpix[1]);
  //do_write_fits (fL1->fg_image, "!fg_image.fits");
  //
  //fL1 = newFuncL ("fg_image.fits", nombresVis, n_archivos, nIni, 1, -1, 128.0, 128.0);
  //fL1->fg_image = do_read ("fg_image.fits");
  //
  ////fL1->atten = atenuacion(fL1->fg_image, fL1->header_obs, fL1->n_archivos, fL1->beam);
  ////calculardxdy (fL1);
  //do_write_fits (fL1->fg_image, "!fg_image.fits");
  //mockCBI(fL1);
  //add_noise (fL1->header_obs[0], fL1->samples_mod[0], seed);
  //
  //do_write ("write uvf", "!visRuidoExacto.sub", fL1->header_obs[0], fL1->samples_mod[0], 
  //	    nombresVis[0]);
  //cout << nombreImagen << "\n";
  //exit (0);


  // MEM
  //FuncionEntropiaNatural *funcEN = new FuncionEntropiaNatural(64*64);
  Funcion *funcEN = new FuncionCero(64*64);
  FuncionMemCBI *funcMemCBI = new FuncionMemCBI (nombreImagen, nombresVisAux, n_archivos, funcEN, 1e-5, 0.0, 1e100);
  //FuncionMemCBI *funcMemCBI = new FuncionMemCBI (nombreImagen, nombresVisAux, n_archivos, funcEN, 0.0, 0.0, 1e100);
  //FuncionMemCBI *funcMemCBI = new FuncionMemCBI (nombreImagen, nombresVisAux, n_archivos, new FuncionEntropiaNatural(128*128), 0.0, 0.0);
  funcMemCBI->setRuido (funcMemCBI->getRuido() * qFactor);
  //ReconstructorGCMemCBI rMem (funcMemCBI, new InicializadorConstante (), 1e-10, 1, 0.0);
  Inicializador *ini = new InicializadorConstante ();
  ReconstructorGCMemCBI rMem (funcMemCBI, ini, 1e-6, 2, 0.0, 6);

  rMem.run();

  char nombreArchivo [255];
  sprintf (nombreArchivo, "q%g_iter", qFactor);
  //funcMemCBI->guardarInfo (nombreArchivo + i2s (iter, 3));

  delete ini;
  delete funcMemCBI;

  // Ajustar hasta nIni
  //Funcion *func = new FuncionChi2VoronoiImagen ("MEM_CBIFinal.fits", 0);
  Funcion *func = new FuncionChi2VoronoiImagen ("MEM_CBIGC005.fits", 0);
  //IncrementadorMallaImagen *inc = new IncrementadorMallaImagen ("MEM_CBIFinal.fits", -1);
  IncrementadorMallaImagen *inc = new IncrementadorMallaImagen ("MEM_CBIGC005.fits", -1);
  InicializadorVoronoiUniforme *iniVU = new InicializadorVoronoiUniforme ();
  ReconstructorIncremental  rInc (func, iniVU, inc , nIni);
  rInc.run();
  
  // Ajustar + GC hasta nFin
  FuncionCBI *funcCBI = new FuncionBayesVoronoiCBI (nombreImagen, nombresVisAux, 
						    n_archivos, nIni, 1, -1, 128.0, 128.0);
  funcCBI->setRuido (funcCBI->getRuido() * qFactor);
  ReconstructorGC recGC (funcCBI, iniVU, 1e-10);
  
//  for (int i = 0; i < n_archivos; i++) {
//    delete [] nombresVisAux[i];
//  }
//  delete [] nombresVisAux;
//  delete funcEN;
//  delete func;
//  delete inc;
//  delete iniVU;
//  delete funcCBI;
//  cout << "FIN\n";
//  return ret;

  double *p = new double [func->getNPars()];
  for (int i = 0; i < func->getNPars(); i++) {
    p[i] = rInc.getPar(i);
  }
  
  //double valor = inc->incrementar (&p, func);
  double valor = func->f(p);
  for (int i = nIni; i < nFin; i += dn) {
    //double valor = inc->incrementar (&p, func);
    if (3*i != func->getNPars()) {
      cerr << "ERROR: 3*i != func->getNPars()\n";
      cerr << "       " << 3*i << " != " << func->getNPars() << "\n";
      exit (1);
    }

    recGC.setNPars (func->getNPars());
    
    //cout << "---------------->func->getNPars() = " << func->getNPars() << "\n";
    
    if (recGC.getNPars() != func->getNPars()) {
      cerr << "ERROR: recGC.getNPars() != func->getNPars()\n";
      cerr << "       " << recGC.getNPars() << " != " << func->getNPars() << "\n";
      exit (1);
    }

    //archivo << i << "\t" << valor << "\n";
    for (int j = 0; j < recGC.getNPars() / 3; j++) {
      //cout << "cambiando " << j << "): " << recGC.getPar(j) << " por " << p[j] << "\n";
      recGC.setPar(3 * j, p[3 * j]);
      recGC.setPar(3 * j + 1, p[3 * j + 1]);
      recGC.setPar(3 * j + 2, p[3 * j + 2] / funcCBI->getRuido());
    }
    ret [i - nIni] = recGC.run();
    cout << ">>>>>>>>>>>funcion = " << funcCBI->f(recGC.getPars()) << " vs " 
	 << ret [i - nIni] << "\n";
    sprintf (nombreArchivo, "Busqueda_q%g_iter", qFactor);
    //funcMemCBI->guardarInfo (nombreArchivo + i2s (iter, 3));
    funcCBI->guardarInfo (nombreArchivo + i2s (iter, 3) + "_" + i2s (i, 4));

    for (int j = 0; j < dn; j++) {
      func->guardarInfo (i2s (i + j, 4));
      valor = inc->incrementar (&p, func);
    }
  }
  cout << "OK\n";
  //exit (0);
  delete func;
  delete funcCBI;
  delete inc;
  delete iniVU;
  delete funcEN;
  delete [] p;

  for (int i = 0; i < n_archivos; i++) {
    delete [] nombresVisAux[i];
  }
  delete [] nombresVisAux;

  return ret;
}

/*--------------------------------------------------------------------
 * add noise to visibilities.
 *--------------------------------------------------------------------*/

static void add_noise(struct uvf_header *header,
		      struct uvf_sample *samples, long *seed) {
  int i, k, nsamp, nchan;
  float sum =0; 
  //float target_Chi2 = 31291. ;
  double sigma;
  double visnoise;
  float scale_factor;
  nsamp = header->nsamp;
  nchan = header->nchan;
  double sum_VR = 0, sum_VI = 0, sum_noise = 0, noise;

  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    for (k = 0; k < nchan; k++) {    /* Loop through channels */
      //if (!ch_list[k])
      //continue;
      if (samples[i].rdata[3*k+2]  != 0) {
	visnoise = sqrt(1/samples[i].rdata[3*k+2]);
	noise = gasdev(seed)*visnoise;
	//printf ("V = %g\t%g,\t noise = %g\n", 
	//samples[i].rdata[3*k], samples[i].rdata[3*k+1],
	//noise);
	sum_noise += fabs(noise);
	sum_VR += fabs(samples[i].rdata[3*k]);
	sum_VI += fabs(samples[i].rdata[3*k+1]);

	//sigma = visnoise * scale_factor;
	//printf("%g %g %g \n",gasdev(seed),visnoise,sigma);
	//printf("%g --> ",samples[i].rdata[3*k]);
	samples[i].rdata[3*k] += gasdev(seed)*visnoise;   // DEVTEST
	//printf("%g \n ",samples[i].rdata[3*k]);
	samples[i].rdata[3*k+1] += gasdev(seed)*visnoise;   // DEVTEST
	}    /* weights */
    }
  }
  printf ("suma = %g\t%g,\t noise = %g\n", 
	  sum_VR, sum_VI, sum_noise);
}

void escalarImagenRuido (struct image *im, struct uvf_header **header, 
			 struct uvf_sample **samples, int n_archivos, double n_ruido) {
  double s = Irms (header, samples, n_archivos);
  double Imax = 0, factor, xf, g, ref_freq = 30.0e9, bmaj, bmin;
  int i;

  bmaj = im->bmaj * RPDEG; // en radianes.
  bmin = im->bmin * RPDEG; // en radianes.
  xf = (PLANCK_H * ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
  s *= (4 * log(2.0) / (PI * bmaj * bmin));
  s /= (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(ref_freq);

  for (i = 0; i < im->npixels; i++) {
    if (im->pixels[i] > Imax) {
      Imax = im->pixels[i];
    }
  }
  factor = n_ruido * s / Imax;
  for (i = 0; i < im->npixels; i++) {
    //printf ("cambiando %g ", im->pixels[i]);
    im->pixels[i] *= factor;
    //printf ("por %g\n", im->pixels[i]);
  }
}

void leerVisibilidades (char **nombresVis, struct uvf_header ***header, struct uvf_sample ***samples_obs, 
			struct uvf_sample ***samples_mod, int n_archivos) {

  (*header) = (struct uvf_header**) malloc(n_archivos * sizeof(struct uvf_header *));
  if ((*header) == NULL) {
    cerr << "ERROR en leerVisibilidades, header = NULL.\n";
    exit (1);
  }

  (*samples_obs) = (struct uvf_sample **) malloc(n_archivos * sizeof(struct uvf_sample *));
  if ((*samples_obs) == NULL) {
    cerr << "ERROR en leerVisibilidades, samples_obs = NULL.\n";
    exit (1);
  }

  (*samples_mod) = (struct uvf_sample **) malloc(n_archivos * sizeof(struct uvf_sample *));
  if ((*samples_mod) == NULL) {
    cerr << "ERROR en leerVisibilidades, samples_mod = NULL.\n";
    exit (1);
  }

  for(int i = 0; i < n_archivos; i++) {
    // abrimos los archivos de visibilidades
    int status = do_read_uvdata (nombresVis[i], &(*header[i]), &(*samples_obs[i]));
    if(status != SUCCESS) {
      cerr << "Error con el archivo uvf\n";
      exit(1);
    }
  }

  for (int i = 0; i < n_archivos; i++) {
    (*samples_mod)[i] = (struct uvf_sample*) malloc(((*header)[i]->nsamp) * sizeof(struct uvf_sample));
    if ((*samples_mod)[i] == NULL) {
      cerr << "ERROR en newFuncL, fL->samples_mod[" << i << "] = NULL.\n";
      exit (1);
    }
  }
  // ya tenemos abiertos los archivos de visibilidades
}
