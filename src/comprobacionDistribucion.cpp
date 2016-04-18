#include <iostream.h>
#include <fstream>
#include "funcL.h"
#include "distribucionPbb.h"
#include "distribucionGibbs.h"
#include "distribucionCondicional.h"
#include "distribucionCondicionalBVCBI.h"
#include "distribucionCondicionalMemCBI.h"
#include "distribucionCondGaussiana.h"

#include "funcionEntropiaNatural.h"

#include "inicializadorArchivo.h"

void escalarIntensidades (double *p, int n, double factor);
void generarMuestras (DistribucionPbb *d, char *nombreArchivo, int N);

void comprobarDistribucionMemCBI (int argc, char ** argv);
void comprobarVisibilidades (int argc, char ** argv);

int main (int argc, char **argv) {
  //comprobarVisibilidades (argc, argv);
  //exit (0);
  comprobarDistribucionMemCBI (argc, argv);
  exit (0);

  //DistribucionCondicionalBVCBI *d = new DistribucionCondicionalBVCBI (argv[1], 10);
  double *medias = new double[2];
  double *sigmas = new double[2];
  medias[0] = 1;
  medias[1] = 2;
  sigmas[0] = 0.5;
  sigmas[1] = 0.25;
  DistribucionCondicional *d = new DistribucionCondGaussiana (2, medias, sigmas);
  DistribucionPbb *dist = new DistribucionGibbs (d);
  double *p = new double [d->getNPars()];
  InicializadorArchivo *ini = new InicializadorArchivo (argv[2]);
  ini->inicializar (p, d->getNPars());

  generarMuestras(dist, "muestras.dat", 1e7);
  exit(0);
  /*
  double quanta = d->getFunc()->getRuido();
  cout << "quanta = " << quanta << "\n";

  escalarIntensidades (p, d->getNPars(), 8.204e-3 / 2.2561);

  for (int i = 0; i < d->getNPars()/3; i++) {
    p[3 * i + 2] /= quanta;
    //p[3 * i + 2] /= 1e-4;
    cout << "p[" << 3 * i + 2 << "] = " << p[3 * i + 2] * quanta << " -> " << p[3 * i + 2] << "\n";
  }
  */
  double a, b;
  int i = 5;

  d->busquedaCotas (p, i, &a, &b);
  cout << "p[2] = " << p[i] << ", a = " << a << ", b = " << b << "\n";
  cout << "P(a) = " << d->evaluar(p, a, i) << ", P(b) = " << d->evaluar(p, b, i) << "\n";
  //cout << p[2] << " -> " << a << "\t" << d->evaluar (p, a, 2) << "\n";
  //cout << p[2] - a << "\n";
  //exit (0);

  
  std::ofstream archivo ("PbbCondicional.dat");
  archivo.close();
  double ret;
  double p_i;
  archivo.precision (8);
  for (p_i = a; p_i < b; p_i += (b - a) * 1e-3) {
    //for (double delta = -1e-5; delta < 1e-5; delta += 1e-7) {
    //p_i = p[2] * (1+delta);
    
    ret = d->evaluar (p, p_i, i);
    archivo.open ("PbbCondicional.dat", ios::app);
    archivo << p_i << "\t" << ret << "\n";
    archivo.close();
    cout << p_i << "\t" << ret << "\n";
  }
  
  archivo.open ("muestraCond.dat");
  archivo.close();
  int N = (int) 1e6;
  archivo.precision (8);
  archivo << N << "\n";
  for (int cont = 0; cont < N; cont++) {
    p_i = d->generar (p, i);
    archivo.open ("muestrasCond.dat", ios::app);
    archivo << p_i <<  "\n";
    archivo.close();
  }
}

void generarMuestras (DistribucionPbb *d, char *nombreArchivo, int N) {
  double *p = new double [d->getNPars()];
  std::ofstream archivo (nombreArchivo);
  archivo.close();
  archivo.precision (8);
  archivo << N << "\n";
  for (int cont = 0; cont < N; cont++) {
    if (cont % 1000 == 0) cout << "cont = " << cont << "\n";
    d->generar (p);
    archivo.open (nombreArchivo, ios::app);
    for (int i = 0; i < d->getNPars(); i++) {
      archivo << p[i] << "\t";
    }
    archivo <<"\n";
    archivo.close();
  }
}


void escalarIntensidades (double *p, int n, double factor) {

  for (int i = 0; i < n / 3; i++) {
    cout << "Cambiando " << p [3 * i + 2];
    p [3 * i + 2] *= factor;
    cout << " por " << p [3 * i + 2] << "\n";
  }
}

void comprobarDistribucionMemCBI (int argc, char ** argv) {
  InicializadorArchivo *ini = new InicializadorArchivo (argv [3]);
  DistribucionCondicionalMemCBI *dcMem1 = new DistribucionCondicionalMemCBI (argv[1], argv + 2, 1/*argc - 2*/, 
									    new FuncionEntropiaNatural(64*64), 0, 51,
									    ini);
  DistribucionCondicionalMemCBI *dcMem2 = new DistribucionCondicionalMemCBI (argv[1], argv + 2, 1/*argc - 2*/, 
									    new FuncionEntropiaNatural(64*64), 0, 51,
									    ini);
  dcMem1->guardarInfo("dcMem1");
  double *p = new double [64*64];
  ini->inicializar (p, 64*64);
  for (int i = 0; i < 64*64; i++) {
    p[i] /= dcMem1->getQuanta();
  }

  std::ofstream archivo ("ComprobacionDistribucion.dat");
  archivo.close();
  archivo.open ("distribCond1.dat");
  archivo.close();
  archivo.open ("distribCond1_x.dat");
  archivo.close();
  archivo.open ("distribCond2.dat");
  archivo.close();
  archivo.open ("distribCond2_x.dat");
  archivo.close();
  double dc1, dc2;
  for (int i = 0; i < 64*64; i++) {
    dc1 = dcMem1->evaluar (p, p[i] * 1.1, i);
    dc2 = dcMem2->evaluar_old (p, p[i] * 1.1, i);

    archivo.open ("ComprobacionDistribucion.dat", ios::app);
    archivo << dc1 << "\t" << dc2 << "\n";
    archivo.close();

    archivo.open ("distribCond1.dat", ios::app);
    archivo << p[i] << "\t" << dc1 << "\n";
    archivo.close();

    archivo.open ("distribCond2.dat", ios::app);
    archivo << p[i] << "\t" << dc2<< "\n";
    archivo.close();

    archivo.open ("distribCond1_x.dat", ios::app);
    archivo << i << "\t" << dc1 << "\n";
    archivo.close();

    archivo.open ("distribCond2_x.dat", ios::app);
    archivo << i << "\t" << dc2 << "\n";
    archivo.close();
  }
}

void comprobarVisibilidades (int argc, char ** argv) {
  funcL *fL = newFuncL(argv[1], argv + 2, argc - 2, 50, 1, -1, 250.0, 250.0, -1);
  int nx =  fL->fg_image->size[0], ny =  fL->fg_image->size[1];

  for (int i = 0; i < fL->fg_image->npixels; i++) {
    fL->fg_image->pixels[i] = 0;
  }

  int ix = 32, iy = 32;
  fL->fg_image->pixels [ix + ny * iy] = fL->difmapNoise * 100;
  
  mockcbi_sub (fL->cmb_image, fL->fg_image, 
	       fL->header_obs[0], fL->samples_obs[0], 
	       fL->samples_mod[0], fL->beam);

  do_write("write", "!Vis_mockcbi.sub", fL->header_obs[0], fL->samples_mod[0], argv[2]);

  double raimage, decimage, dx, dy, kx, seno, coseno, x0, y0, x, y;
  double obsra, obsdec, lobs, mobs, ra_offset = 0.0, raizxy, k, VR, VI;
  std::ofstream archivoMock ("Vis_mockcbi.dat");
  std::ofstream archivoMockR ("Vis_mockcbiR.dat");
  std::ofstream archivoMockI ("Vis_mockcbiI.dat");
  std::ofstream archivoIter ("Vis_iter.dat");
  std::ofstream archivoIterR ("Vis_iterR.dat");
  std::ofstream archivoIterI ("Vis_iterI.dat");
  
  
  dx = fL->fg_image->cdelt[0] * RPDEG;   /* radians */
  dy = fL->fg_image->cdelt[1] * RPDEG;   /* radians */
  //double pixel = fabs(fg_image->cdelt[0] * RPDEG * fg_image->cdelt[1] * RPDEG);    /* pixel solid angle in ster */
  
  /* coordinates of the image reference pixel in radians */
  raimage = fL->fg_image->crval[0] * RPDEG;
  decimage = fL->fg_image->crval[1] * RPDEG;
  
  for (int arch = 0; arch < fL->n_archivos; arch++) {                                   /* filenames loop */
    obsra = ra_offset + fL->header_obs[arch]->obsra*RPDEG;
    obsdec = fL->header_obs[arch]->obsdec*RPDEG;  
    /* direction cosines of the phase center in the image coordinate system */
    direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
    /* Find the phase center in the pixel coordinate system of the image (where pixel numbers are zero-based) */
    x0 = (fL->fg_image->crpix[0] - 1.0) + lobs / dx;
    y0 = (fL->fg_image->crpix[1] - 1.0) + mobs / dy;
    
    x = (ix - x0) * dx;          /* radians */
    y = (iy - y0) * dy;          /* radians */
    
    raizxy = sqrt (1 - x*x - y*y);
    //cout << "raizxy = " << raizxy << "\n";
    
    for (int chan = 0; chan < fL->header_obs[arch]->nchan; chan++) {        /* chan: loop through channels */
      
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {      /* i: loop through samples */
	
	/* u*freq : radians , x,y : radians */
	if (fL->samples_mod[arch][samp].rdata[3 * chan + 2] != 0) {
	  k = sqrt (pow (fL->samples_mod[arch][samp].u, 2) + pow (fL->samples_mod[arch][samp].u, 2)) * fL->header_obs[arch]->iffreq[chan];
	  kx = (fL->samples_mod[arch][samp].u * x + fL->samples_mod[arch][samp].v * y) * fL->header_obs[arch]->iffreq[chan];
	  seno = sin(2 * PI * kx);
	  coseno = cos(2 * PI * kx);
	  
	  //cout << " Cambiando " << samples_mod[arch][samp].rdata[3 * chan] << " por ";
	  VR = fL->samples_mod[arch][samp].rdata[3 * chan];
	  VI = fL->samples_mod[arch][samp].rdata[3 * chan + 1];
	  
	  archivoMock << k << "\t" << sqrt (VR * VR + VI * VI) << "\n";
	  archivoMockR << k << "\t" << VR << "\n";
	  archivoMockI << k << "\t" << VI << "\n";
	  
	  fL->samples_mod[arch][samp].rdata[3 * chan]     = (fabs (dx * dy) * coseno * fL->atten[arch][ix + nx * iy + 1][chan] * 
							   100 * fL->difmapNoise * fL->fg_scale / raizxy);
	  fL->samples_mod[arch][samp].rdata[3 * chan + 1] = (fabs (dx * dy) * seno * fL->atten[arch][ix + nx * iy + 1][chan] * 
							   100 * fL->difmapNoise * fL->fg_scale / raizxy);
	  
	  VR = fL->samples_mod[arch][samp].rdata[3 * chan];
	  VI = fL->samples_mod[arch][samp].rdata[3 * chan + 1];
	  
	  archivoIter << k << "\t" << sqrt (VR * VR + VI * VI) << "\n";
	  archivoIterR << k << "\t" << VR << "\n";
	  archivoIterI << k << "\t" << VI << "\n";
	}
      }
    }
  }

  do_write("write", "!Vis_iter.sub", fL->header_obs[0], fL->samples_mod[0], argv[2]);
  
  
}
