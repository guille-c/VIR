/*********************
 * Comprobacion dVdI *
 *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "funcionBayesVoronoiCBI.h"
#include "mcheck.h"

double *abrirArchivo (char *nombreArchivo, int *n);
void probarMalla(funcL *fL, double pars[]);
void compararVisibilidades(funcL *fL);
void compararVisibilidades(struct uvf_sample *sin, struct uvf_sample *sout, int n);
void imprimirVisibilidades (funcL *fL, char *nombre_archivo, int tipo);
void imprimirVisibilidadesu (funcL *fL, char *nombre_archivo, int tipo);
void crearImagenDelta(struct image* im, int x0, int y0, double dx, double dy, 
		      int nx, int ny, double F0);
void visFuentePuntual (funcL *fL, int x0, int y0, double F0);
void apodizar (funcL *fL);

enum impresionVis {RE, IM, MOD};

int main (int argc, char **argv) {
  struct uvf_header *header_1;
  struct uvf_header *header_2;
  struct uvf_sample *samples_1;
  struct uvf_sample *samples_2;

//  do_read_uvdata(argv[1], &(header_1), &(samples_1));
//  do_read_uvdata(argv[2], &(header_2), &(samples_2));

//  compararVisibilidades (samples_1, samples_2, header_1->nsamp);
//  exit(0);

  funcL *fL;
  FILE *outR, *outI;
  int x0 = 16, y0 = 16, nx = 512, ny = 512;
  double F0 = 1; // [Jy]
  double dx = -1, dy = 1;

  //x0 = -x0 * 0.75 / dx;
  //y0 = y0 * 0.75 / dy;

  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (argv[1]);
  
  if (argc < 2) {
    printf ("ERROR usage: comparar_visibilidades archivo.in\n");
    exit (1);
  }
  
  fL = func->getFL();

  //cout << "A: " << fL->atten[0][x0 + y0 * fL->fg_image->size[0] + 1][0] << "\n";
  //cout << "A: " << fL->atten[0][256 *128][0] << "\n";
  // func->guardarAtenuaciones ();

  fL->fg_image = do_read (fL->nombreFits);
  //crearImagenDelta (fL->fg_image, x0, y0, dx, dy, nx, ny, F0 / fL->fg_scale);
  func->guardarInfo ("Comparacion");
  //exit (0);
  //apodizar (fL);
  //compararVisibilidades (fL);
  visFuentePuntual (fL, x0, y0, F0);
  do_write("write", "!vis_exacto.sub", fL->header_obs[0],fL->samples_mod[0], fL->infile[0]);
  imprimirVisibilidades (fL, "vis_exacto.dat", MOD);
  imprimirVisibilidades (fL, "vis_exactoR.dat", RE);
  imprimirVisibilidades (fL, "vis_exactoI.dat", IM);
  return 0;
}

void compararVisibilidades(struct uvf_sample *s1, struct uvf_sample *s2, int n) {
  FILE *archivo = fopen ("comparacion.dat", "w");
  FILE *archivoR = fopen ("comparacionR.dat", "w");
  FILE *archivoI = fopen ("comparacionI.dat", "w");
  double mod1, mod2;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < MAX_DAT / 3; j++) {
      mod1 = sqrt (s1[i].rdata[3 * j] * s1[i].rdata[3 * j] + 
		   s1[i].rdata[3 * j + 1] * s1[i].rdata[3 * j + 1]);
      mod2 = sqrt (s2[i].rdata[3 * j] * s2[i].rdata[3 * j] + 
		   s2[i].rdata[3 * j + 1] * s2[i].rdata[3 * j + 1]);
      fprintf (archivo, "%g\t%g\n", mod1, mod2);
      fprintf (archivoR, "%g\t%g\n", s1[i].rdata[3 * j], s2[i].rdata[3 * j]);
      fprintf (archivoI, "%g\t%g\n", s1[i].rdata[3 * j + 1], s2[i].rdata[3 * j + 1]);
    }
  }
}

void compararVisibilidades (funcL *fL) {
  double *visR =  new double [fL->n_archivos * fL->header_obs[0]->nif * 
			     fL->header_obs[0]->nsamp * sizeof(double)];
  double *visI =  new double [fL->n_archivos * fL->header_obs[0]->nif * 
			     fL->header_obs[0]->nsamp * sizeof(double)];
  double u, mockcbiMod, directMod, dx = fL->dx[0];
  int i = 0, arch = 0, nx = fL->fg_image->size[0];

  FILE *archivo = fopen ("comparacionVis.dat", "w");
  FILE *archivoR = fopen ("comparacionVisR.dat", "w");
  FILE *archivoI = fopen ("comparacionVisI.dat", "w");
  FILE *archivoError  = fopen ("erroruVis.dat", "w");
  FILE *archivoErrorR = fopen ("erroruVisR.dat", "w");
  FILE *archivoErrorI = fopen ("erroruVisI.dat", "w");
  
  expanded_mockcbi(fL, arch);
  do_write("write", "!vis_expanded_mockcbi.sub", fL->header_obs[0],fL->samples_mod[0], fL->infile[0]);
  imprimirVisibilidades (fL, "vis_expanded_mockcbi.dat", MOD);
  imprimirVisibilidades (fL, "vis_expanded_mockcbiR.dat", RE);
  imprimirVisibilidades (fL, "vis_expanded_mockcbiI.dat", IM);
  imprimirVisibilidadesu (fL, "vis_expanded_mockcbi_u.dat", MOD);
  imprimirVisibilidadesu (fL, "vis_expanded_mockcbiR_u.dat", RE);
  imprimirVisibilidadesu (fL, "vis_expanded_mockcbiI_u.dat", IM);

  mockcbi_sub (fL->cmb_image, fL->fg_image, 
	       fL->header_obs[arch], fL->samples_obs[arch], 
	       fL->samples_mod[arch], fL->beam);

  do_write("write", "!vis_mockcbi.sub", fL->header_obs[0],fL->samples_mod[0], fL->infile[0]);
  imprimirVisibilidades (fL, "vis_mockcbi.dat", MOD);
  imprimirVisibilidades (fL, "vis_mockcbiR.dat", RE);
  imprimirVisibilidades (fL, "vis_mockcbiI.dat", IM);
  imprimirVisibilidadesu (fL, "vis_mockcbi_u.dat", MOD);
  imprimirVisibilidadesu (fL, "vis_mockcbiR_u.dat", RE);
  imprimirVisibilidadesu (fL, "vis_mockcbiI_u.dat", IM);
  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 2] != 0) {
	  visR[i] = fL->samples_mod[arch][samp].rdata[3 * iff];
	  visI[i] = fL->samples_mod[arch][samp].rdata[3 * iff + 1];
	}
	i++;
      }
    }
  }
  cout << "i = " << i << ", " << fL->n_archivos * fL->header_obs[0]->nif * 
    fL->header_obs[0]->nsamp << "\n";

//  mockCBI (fL);
//  do_write("write", "!vis_directo.sub", fL->header_obs[0],fL->samples_mod[0], fL->infile[0]);
//  imprimirVisibilidades (fL, "vis_directo.dat", MOD);
//  imprimirVisibilidades (fL, "vis_directoR.dat", RE);
//  imprimirVisibilidades (fL, "vis_directoI.dat", IM);
//  imprimirVisibilidadesu (fL, "vis_directo_u.dat", MOD);
//  imprimirVisibilidadesu (fL, "vis_directoR_u.dat", RE);
//  imprimirVisibilidadesu (fL, "vis_directoI_u.dat", IM);
//  i = 0;
//  for (int arch = 0; arch < fL->n_archivos; arch++) {
//    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
//      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
//	if (fL->samples_mod[arch][samp].rdata[3 * iff + 2] != 0) {
//	  mockcbiMod = sqrt (visR[i] * visR[i] + visI[i] * visI[i]);
//	  directMod = sqrt (fL->samples_mod[arch][samp].rdata[3 * iff] * 
//			    fL->samples_mod[arch][samp].rdata[3 * iff] + 
//			    fL->samples_mod[arch][samp].rdata[3 * iff + 1] *
//			    fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
//	  fprintf (archivo, "%g\t%g\n", mockcbiMod, directMod);
//	  fprintf (archivoR, "%g\t%g\n", visR[i], fL->samples_mod[arch][samp].rdata[3 * iff]);
//	  fprintf (archivoI, "%g\t%g\n", visI[i], fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
//
//	  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff] * (nx*dx);
//
//	  fprintf (archivoError, "%g\t%g\n", u, fabs (mockcbiMod - directMod / directMod));
//	  fprintf (archivoErrorR, "%g\t%g\n", u, fabs (visR[i] - fL->samples_mod[arch][samp].rdata[3 * iff]
//						       / fL->samples_mod[arch][samp].rdata[3 * iff]));
//	  fprintf (archivoErrorI, "%g\t%g\n", u, fabs (visI[i] - fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
//						       / fL->samples_mod[arch][samp].rdata[3 * iff + 1]));
//	}
//	i++;
//      }
//    }
//  }
  fclose (archivo);
}

void imprimirVisibilidades (funcL *fL, char *nombre_archivo, int tipo) {
  double vis, radio;
  FILE *archivo = fopen (nombre_archivo, "w");

  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 2] != 0) {
	  radio = sqrt(fL->samples_mod[arch][samp].u * fL->samples_mod[arch][samp].u +
		       fL->samples_mod[arch][samp].v * fL->samples_mod[arch][samp].v)
	    * fL->header_obs[arch]->iffreq[iff];
	  if (tipo == MOD) {
	    vis = sqrt (fL->samples_mod[arch][samp].rdata[3 * iff] * 
			fL->samples_mod[arch][samp].rdata[3 * iff] + 
			fL->samples_mod[arch][samp].rdata[3 * iff + 1] *
			fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
	  }
	  else if (tipo == RE) {
	    vis = fL->samples_mod[arch][samp].rdata[3 * iff];
	  }
	  else if (tipo == IM) {
	    vis = fL->samples_mod[arch][samp].rdata[3 * iff + 1];
	  }
	  fprintf (archivo, "%g\t%g\n", radio, vis);
	}
	else {
	  //printf ("w = 0\n");
	}
      }
    }
  }
  
  fclose (archivo);
}

void imprimirVisibilidadesu (funcL *fL, char *nombre_archivo, int tipo) {
  double vis;
  FILE *archivo = fopen (nombre_archivo, "w");

  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 2] != 0) {
	  if (tipo == MOD) {
	    vis = sqrt (fL->samples_mod[arch][samp].rdata[3 * iff] * 
			fL->samples_mod[arch][samp].rdata[3 * iff] + 
			fL->samples_mod[arch][samp].rdata[3 * iff + 1] *
			fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
	  }
	  else if (tipo == RE) {
	    vis = fL->samples_mod[arch][samp].rdata[3 * iff];
	  }
	  else if (tipo == IM) {
	    vis = fL->samples_mod[arch][samp].rdata[3 * iff + 1];
	  }
	  fprintf (archivo, "%g\t%g\n", fL->samples_mod[arch][samp].u, vis);
	}
	else {
	  //printf ("w = 0\n");
	}
      }
    }
  }
  
  fclose (archivo);
}

void crearImagenDelta(struct image* im, int x0, int y0, double dx, double dy, 
		      int nx, int ny, double F0) {
  int xc = im->crpix[0];
  int yc = im->crpix[0];

  resize_map (im, nx, ny, im->crpix[0], im->crpix[1]);
  //x0 = im->crpix[0] - x0 * (dx * RPARCM) / (im->cdelt[0] * RPDEG);
  //y0 = im->crpix[1] - y0 * im->cdelt[1] * RPDEG / (dy * RPARCM);
  x0 = im->crpix[0] - x0;
  y0 = im->crpix[1] - y0;

  for (int i = 0; i < im->npixels; i++) {
    im->pixels[i] = 0;
  }
  im->cdelt[0] = dx * RPARCM / RPDEG;
  im->cdelt[1] = dy * RPARCM / RPDEG;
  //  int nx = im->size[0], ny = im->size[1];
  //  cout << "imagen cambiada a (" << nx << ", " << ny << ")\n";

  im->pixels [x0 + y0 * nx] = F0 / fabs(im->cdelt[0] * RPDEG * im->cdelt[1] * RPDEG);
}

void apodizar (funcL *fL){
  struct image *im = fL->fg_image;
  int nx = im->size[0], ny = im->size[1], arch = 0, iff = 0;

  for (int i = 0; i < im->npixels; i++) {
    im->pixels[i] *= fL->atten[arch][i + 1][iff] ;
  }
}

void visFuentePuntual (funcL *fL, int x0, int y0, double F0) {
  int i, j, arch, iff, samp;
  double Re, Im, x, y, kx, u, v, dx, dy;
  struct image *im = fL->fg_image;

  dx = im->cdelt[0] * RPDEG;
  dy = im->cdelt[1] * RPDEG;

  cout << "dx = " << im->cdelt[0] << " grad\n";

  x = (x0 - fL->x0[0]) * dx;
  y = (y0 - fL->y0[0]) * dy;

  x0 = im->crpix[0] - x0;
  y0 = im->crpix[1] - y0;  

  struct image *im_aux = do_read (fL->nombreFits);
  for (i = 0; i < im_aux->npixels; i++) {
    im_aux->pixels[i] = 0;
  }
  im_aux->pixels[x0 + y0 * im_aux->size[0]] = F0;
  do_write_fits (im_aux, "!FtePuntual.fits");

  for (arch = 0; arch < fL->n_archivos; arch++) {
    copiar_uvf_samples(fL->samples_obs[arch], fL->samples_mod[0], fL->header_obs[arch]->nsamp);
    do_clear("clear", fL->header_obs[arch], fL->samples_mod[arch]);   /* Visibilidades = 0 */
    for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	//fL->samples_mod[arch][samp].u = fL->samples_obs[arch][samp].u;
	//fL->samples_mod[arch][samp].v = fL->samples_obs[arch][samp].v;
	//fL->samples_mod[arch][samp].w = fL->samples_obs[arch][samp].w;
	u = fL->samples_obs[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
	v = fL->samples_obs[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
	kx = (u * x + v * y);

	cout << "u = " << u << ", v = " << v << "\n";

	Re = F0 * fL->atten[arch][x0 + y0 * fL->fg_image->size[0] + 1][iff] 
	  * cos (-2 * PI * kx);
	Im = F0 * fL->atten[arch][x0 + y0 * fL->fg_image->size[0] + 1][iff] 
	  * sin (-2 * PI * kx);
	
	//cout << "arch = " << arch << ", iff = " << iff << ", samp = " << samp << "\n";
	//cout << "x0 = " << x0 << ", y0 = " << y0 << "\n";
	//cout << "Re: " << Re << ", IM: " << Im << ", A: " << fL->atten[arch][x0 + y0 * fL->fg_image->size[0] + 1][iff] << "\n\n";
	//	if (fL->samples_mod[arch][samp].rdata[3 * iff + 2] != 0) {
 	  fL->samples_mod[arch][samp].rdata[3 * iff] = Re; 
 	  fL->samples_mod[arch][samp].rdata[3 * iff + 1] = Im;
	  //fL->samples_mod[arch][samp].rdata[3 * iff + 2] = fL->samples_obs[arch][samp].rdata[3 * iff + 2];
	  //}
      }
    }
  }
}
