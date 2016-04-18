#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "random.h"

#define POW(x) ((x)*(x))

double cross_correlate_factor (struct uvf_header **header_31,
			       struct uvf_sample **samples_31, 
			       struct uvf_sample **samples_simu,
			       int n_archivos);
double multi_factor (struct uvf_header **header_31, struct uvf_sample **samples_31, 
		     struct uvf_sample **samples_simu, int n_archivos);
//static void scale (struct uvf_header *header,             /* scale visibilities  */
//	   struct uvf_sample *samples, float ascale);
static void add_noise (struct uvf_header *header,             /* add noise to visibilities  */
		      struct uvf_sample *samples,
		      long *seed);
void escalarImagenRuido (struct image *im, struct uvf_header **header, 
			 struct uvf_sample **samples, int n_archivos, double n_ruido);

void insertarGaussiana (double mx, double my, double sigma, struct image *im, double factor);
void insertarRectangulo (int i1, int j1, int i2, int j2, struct image *im, double valor);

void imagen2ascii (char *nombreImagen, char *nombreAscii);

int main(int argc, char **argv) {
  struct image *im = do_read (argv[1]);
  struct uvf_header *header [1];
  struct uvf_sample *samples [1];
  funcL *fL;
  long seed;

  //resize_map(im , 128, 128, 32, 32);
  //do_write_fits (im, argv[2]);
  //exit (0);

  do_read_uvdata(argv[2], &(header[0]), &(samples[0]));

  struct pbeam beam;               // Beam del CBI
  init_beam (&(beam));
  beam.type = CBI;

  double ***aten = atenuacion(im, header, 1, beam);
  dirtyMap (im, header, samples, 1, aten);

  do_write_fits (im, "!dirtyMap.fits");

  fL = newFuncL (argv[1], argv + 2, argc - 4, 150, 0, 0, 500.0, 500.0, -1);

  mockcbi_sub (fL->cmb_image, fL->fg_image, 
	       fL->header_obs[0], fL->samples_obs[0], fL->samples_mod[0], 
	       fL->beam);

  seed =  -labs((long)time(0));
  add_noise (fL->header_obs[0],
	     fL->samples_mod[0], &seed);
  do_write ("write uvf", argv[4], fL->header_obs[0],
	    fL->samples_mod[0], argv[2]);

  //imagen2ascii (argv[1], argv[2]);

  /*
  int n = 150, cx, cy, arch, chan, samp, i;
  float ftol = 1e-20;
  funcL *fL;
  struct image *punt = do_read (argv[1]);
  double ascale;
  FILE *archivo = fopen ("correlacion.dat", "w");
  FILE *output  = fopen ("crear_imagen_out.dat", "w");
  struct image *imagen_in = do_read (argv[1]);
  //struct image *imagen_out = do_read (argv[2]);
  long seed;
  
  printf ("%s, %s, %s, %s\n\n", argv[1], argv[2], argv[3], argv[4]);
  fprintf (output, "crear_imagen %s, %s, %s, %s\n\n", 
	   argv[1], argv[2], argv[3], argv[4]);

  fL = newFuncL (argv[1], argv + 2, argc - 4, n, 0, 0);
  //do_write_fits (fL->fg_image, argv[3]);
  
  copy_empty_map (punt, fL->fg_image);
  delete_map(fL->fg_image);
  //for (i = 0; i < MAXDIM; i++) {
    //punt->cdelt[i] *= 1.5;
    //fL->cmb_image->cdelt[i] *= 1.5;
  //}

  fL->fg_image = punt;

  //punt->pixels[64 + 64 * fL->fg_image->size[0]] = 0.1;
  
  insertarGaussiana (30, 30, 10, punt, 50);
  //insertarGaussiana (25, 27,  5, punt, 1.5);
  //insertarGaussiana (40, 40,  3, punt, 1);
  insertarGaussiana (45, 45,  4, punt, 10);
  insertarGaussiana (32, 15,  4, punt, 10);

  insertarRectangulo (32, 27, 50, 37, punt, 0.3);

  escalarImagenRuido (punt, fL->header_obs, fL->samples_obs, fL->n_archivos, 100);

  punt->datamin = 1e10;
  punt->datamax = -1e10;
  for (i = 0; i < punt->npixels; i++) {
    if (punt->pixels[i] < punt->datamin) {
      punt->datamin = punt->pixels[i];
    }
    if (punt->pixels[i] > punt->datamax) {
      punt->datamax = punt->pixels[i];
    }
  }

  do_write_fits (fL->cmb_image, "!cmb.fits");
  do_write_fits (punt, argv[3]);
  
  mockcbi_sub (fL->cmb_image, fL->fg_image, 
	       fL->header_obs[0], fL->samples_obs[0], fL->samples_mod[0], 
	       fL->beam);

  seed =  -labs((long)time(0));
  add_noise (fL->header_obs[0],
	     fL->samples_mod[0], &seed);
  do_write ("write uvf", argv[4], fL->header_obs[0],
	    fL->samples_mod[0], argv[2]);
  */
  
  return 0;
}

void imagen2ascii (char *nombreImagen, char *nombreAscii) {
  struct image *im = do_read (nombreImagen);
  
  FILE *archivo = fopen (nombreAscii, "w");
  double Imax = 0;
  int i, j;

  for (i = 0; i < im->npixels; i++) {
    if (im->pixels[i] > Imax) Imax = im->pixels[i];
  }
  
  fprintf (archivo, "P2\n# Created by Guille\n");
  fprintf (archivo, "%d %d\n", im->size[0], im->size[1]);
  fprintf (archivo, "255\n");
  for (j = im->size[1]; j > 0; j--) {
    for (i = 0; i < im->size[0]; i++) {
      fprintf (archivo, "%d ", (int) (im->pixels[i + im->size[0] * j] * 255 / Imax));
    }
    fprintf (archivo, "\n");
  }

  close (archivo);
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
  return c1 / c2;
}

/*--------------------------------------------------------------------
 * Calculo del factor de multiplicacion entre visibilidades de
 * samples_31 y samples_simu.
 *--------------------------------------------------------------------*/

double multi_factor (struct uvf_header **header_31,
		     struct uvf_sample **samples_31, 
		     struct uvf_sample **samples_simu, int n_archivos) {
  int arch, samp, chan;
  double c1 = 0, c2 = 0;

  for (arch = 0; arch < n_archivos; arch++) {
    for (chan = 0; chan < header_31[arch]->nchan; chan++) {
      for (samp = 0; samp < header_31[arch]->nsamp; samp++) {
 	c1 += (POW(samples_31[arch][samp].rdata[chan * 3]) +
	       POW(samples_31[arch][samp].rdata[chan * 3 + 1]));

 	c2 += (POW(samples_simu[arch][samp].rdata[chan * 3]) +
	       POW(samples_simu[arch][samp].rdata[chan * 3 + 1]));
      }
    }
  }
  return sqrt(c1 / c2);
}


/*--------------------------------------------------------------------
 * scale input  visibilities
 *--------------------------------------------------------------------*/
/*
static void scale (struct uvf_header *header,             // scale visibilities
	   struct uvf_sample *samples, float ascale) {
  int i, k, nsamp, nchan;

  nsamp = header->nsamp;
  nchan = header->nchan;


  //  printf ("scaling vis by %g \n",*scale);
  printf ("scaling vis by %g \n",ascale);
  
  for (i = 0; i < nsamp; i++) {      // Loop through samples 
  for (k = 0; k < nchan; k++) {      // Loop through channels
      //if (!ch_list[k])
      //continue;
      if (samples[i].rdata[3*k+2]  != 0) {
	samples[i].rdata[3*k]     *=  ascale;   
 	samples[i].rdata[3*k + 1] *=  ascale; 
      }    // weights 
    }
  }
}
*/
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
  s *= (4 * log(2) / (PI * bmaj * bmin));
  s /= (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(ref_freq);

  for (i = 0; i < im->npixels; i++) {
    if (im->pixels[i] > Imax) {
      Imax = im->pixels[i];
    }
  }
  factor = n_ruido * s / Imax;
  for (i = 0; i < im->npixels; i++) {
    printf ("cambiando %g ", im->pixels[i]);
    im->pixels[i] *= factor;
    printf ("por %g\n", im->pixels[i]);
  }
}

/*--------------------------------------------------------------------
 * Inserta un Gaussiana en im centrada en (mx, my) (en pixeles).
 *--------------------------------------------------------------------*/

void insertarGaussiana (double mx, double my, double sigma, struct image *im, double factor){
  int i, j, nx = im->size[0], ny = im->size[1];
  
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      im->pixels[i + j * ny] += factor * (exp (-((i - mx) * (i - mx) + (j - my) * (j - my))
					       / (2 * sigma * sigma))
					  / (sigma* sqrt(2 * PI)));
    }
  }
}

/*--------------------------------------------------------------------
 * Inserta un Rectangulo en im con vertices en (i1, j1) y (i2, j2) (en
 * pixeles)
 *--------------------------------------------------------------------*/

void insertarRectangulo (int i1, int j1, int i2, int j2, struct image *im, double valor) {
  int i, j, ny = im->size[1];
  
  for (i = i1; i <= i2; i++) {
    for (j = j1; j <= j2; j++) {
      im->pixels[i + j * ny] += valor;
    } 
  } 
}
