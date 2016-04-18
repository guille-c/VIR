/*****************************
 * Comprobacion integraldVdx *
 *****************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "mcheck.h"

double gradienteNumerico (funcL *fL, double *pars, int i, double delta);
double *abrirArchivo (char *nombreArchivo, int *n);
void probarMalla(funcL *fL, double pars[]);
void integraldVdxHorizontal (funcL *fL, int arch, int iff, int samp, 
			     double y, double *integralR, double *integralI);
int integral (funcL *fL, int arch, int iff, int samp, 
	      PoligonoVoronoi *pol, AristaVoronoi *a,
	      double *integralR, double *integralI);

int main (int argc, char **argv) {
  int n = 2, i, indice = 4;
  funcL *fL;
  FILE *outxR, *outxI, *outyR, *outyI;
  double *pars, *grad, delta;
  double dVdx_R = 0, dVdx_I = 0, dVdy_R = 0, dVdy_I = 0, 
    datamin = 1e300, datamax = -1e300;
  int arch = 0, samp = 1, iff = 1;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  int ix, iy, j, k;
  AristaVoronoi *a;
  double integralR, integralI;
  mtrace();

  outxR = fopen ("integralesR.dat", "w");
  outxI = fopen ("integralesI.dat", "w");

  if (argc < 2) {
    printf ("ERROR usage: comprobacion_gradiente archivo.fits archivo.sub\n");
    exit (1);
  }
  
  //pars = abrirArchivo (argv[argc - 1], &n);
  //pars = (double *) malloc (3 * n * sizeof(double));
  pars = (double *) calloc (3 * n, sizeof(double));

  pars[0] = 0.0;
  pars[1] = 0.0;
  pars[2] = 1;
  pars[3] = 1.0;
  pars[4] = 1.0;
  pars[5] = 2;

  
  //grad = (double *) malloc (3 * n * sizeof(double));
  grad = (double *) calloc (3 * n, sizeof(double));

  fL = newFuncL (argv[1], argv + 2, argc - 2, n);
  L (fL, pars, n);
  
  for (i = 0; i < 1; i++) {
    k = 1;
    for (iy = 0; iy < fL->fg_image->size[1]; iy++) {
      for (ix = 0; ix < fL->fg_image->size[0]; ix++) {
	for(j = 0; j < fL->header_obs[i]->nif; j++) {
	  fL->atten[i][k][j] = 1;
	}
	k++;
      }
    }
  }
  
  //pol = (PoligonoVoronoi *) fL->malla->poligonos->primero->info;

  for (pars[0] = 0.0, pars[3] = 1.0; pars[0] <= 1.0 && pars[3] >= 0.0; 
       pars[0] += 0.1 , pars[3] -= 0.1) {
    double alfa = atan2(pars[4] - pars[1], pars[3] - pars[0]) * 180 / PI;

    L (fL, pars, n);
    for (i = 0, nodo = fL->malla->poligonos->primero; 
	 i < n; i++, nodo = nodo->sgte) {
      
      pol = (PoligonoVoronoi *) nodo->info;
      printf ("\nPoligono %d, (%g, %g)\n", i, pol->x, pol->y);
      a = pol->a;
      do {
	if (integral (fL, arch, iff, samp, pol, a, 
		      &integralR, &integralI)) {
	  //double alfa = a
	  printf ("\na: ");
	  imprimirArista (a);
	  printf ("integralR = %g, integralI = %g vs\n", 
		  integralR, integralI);
	  fprintf (outxR, "%g\t%g\t", alfa, integralR);
	  fprintf (outxI, "%g\t%g\t", alfa, integralI);

	  integraldVdx (fL, arch, iff, samp, pol, a, &integralR, &integralI);
	  fprintf (outxR, "%g\t", integralR);
	  fprintf (outxI, "%g\t", integralI);
	  printf ("integralR = %g, integralI = %g vs\n",
		  integralR, integralI);
	  integraldVdxAprox (fL, arch, iff, samp, pol, a, &integralR, &integralI);
	  fprintf (outxR, "%g\n", integralR);
	  fprintf (outxI, "%g\n", integralI);
	  printf ("integralR = %g, integralI = %g\n",
		  integralR, integralI);
	}
	if (a->poliDer == pol) {
	  a = a->cwPred;
	}
	else if (a->poliIzq == pol) {
	  a = a->cwSucc;
	}
      } while (a != pol->a);
    }
  }

  //x1 = (p1->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  printf ("x1= %g, y1 = %g, x2 = %g, y2 = %g\n",
	  - fL->x0[arch] * fL->dx[arch], - fL->y0[arch] * fL->dx[arch],
	  (fL->fg_image->size[0] - fL->x0[arch]) * fL->dx[arch], 
	  (fL->fg_image->size[1] - fL->y0[arch]) * fL->dx[arch]);
  printf ("sinalfa = %g\n", sin(PI / 4));
  eliminarFuncL (fL);
  free (pars);
  free (grad);
  
  fclose (outxR);
  fclose (outxI);
  
  return 0;
}

void integraldVdxHorizontal (funcL *fL, int arch, int iff, int samp, 
			     double y, double *integralR, double *integralI) {
  int i, nx, ny, j;
  double x1, y1, x2, u, v, I1R, I1I, factor;
  
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  (*integralR) = 0;
  (*integralI) = 0;
  
  j = y * (ny - 1.0);
  y1 = (y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  for (i = 0; i < nx; i++) {
    x1 = (i - fL->x0[arch]) * fL->dx[arch];
    factor = sin (PI * u * fL->dx[arch]) / (PI * u);
    I1R = factor * cos (2 * PI * (x1 * u + y1 * v));
    I1I = factor * sin (2 * PI * (x1 * u + y1 * v));
    (*integralR) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1R;
    (*integralI) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1I;
  }
}

int integral (funcL *fL, int arch, int iff, int samp, 
	      PoligonoVoronoi *pol, AristaVoronoi *a,
	      double *integralR, double *integralI) {
  PuntoVoronoi *p1, *p2;
  int nx, ny, tipoLinea;
  double u, v, m, x1, x2, y1, y2, x, y, dx, sin_alfa, cos_alfa, a1;
  PuntoVoronoi *paux;

  p1 = newPuntoVoronoi (0, 0);
  p2 = newPuntoVoronoi (0, 0);
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  (*integralR) = 0;
  (*integralI) = 0;
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  //printf ("1) u, v = (%g, %g)\n", 
  //fL->samples_mod[arch][samp].u, fL->samples_mod[arch][samp].v);
  
  m = (a->ptoIni->y - a->ptoFin->y) * ny / 
    ((a->ptoIni->x - a->ptoFin->x) * nx);
  if (m > 1) {
    tipoLinea = 3;
  }
  else if (m > 0) {
    tipoLinea = 1;
  }
  else if (m > -1) {
    tipoLinea = 2;
  }
  else {
    tipoLinea = 4;
  }

  if (a->ptoIni->x < a->ptoFin->x) {
    if (!interseccionCuadrado (a->ptoIni, a->ptoFin, p1, p2, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return FALSE;
    }
  }
  else {
    if (!interseccionCuadrado (a->ptoFin, a->ptoIni, p2, p1, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return FALSE;
    }
  }

  if (a->poliDer == pol) {
    paux = p1;
    p1 = p2;
    p2 = paux;
    cos_alfa = a->cosDer;
    sin_alfa = a->sinDer;
  }
  else if (a->poliIzq == pol) {
    cos_alfa = a->cosIzq;
    sin_alfa = a->sinIzq;
  }
  else {
    fprintf (stderr, "ERROR en integraldVdx: arista no coincide con poligono.\n");
    exit (1);
  }

  //printf ("p1 = (%g, %g), p2 = (%g, %g)\n", p1->x, p1->y, p2->x, p2->y);
  //printf ("cosalfa = %g, sinalfa = %g\n", cos_alfa, sin_alfa);

  // Damos vuelta la imagen
    
  paux = p1;
  p1 = p2;
  p2 = paux;
  cos_alfa = -cos_alfa;
  //p1->x = 1 - p1->x;
  //p2->x = 1 - p2->x;

  //x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  //y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  x1 = (p1->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  y1 = (p1->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  x2 = (p2->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  y2 = (p2->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  x = (x1 + x2)/ 2;
  y = (y1 + y2)/ 2;
  dx = x2 - x1;
  a1 = PI * (v * cos_alfa - u * sin_alfa);

  (*integralR) = cos (2 * PI * (u * x + v * y)) * sin (a1 * dx / sin_alfa) / a1;
  (*integralI) = sin (2 * PI * (u * x + v * y)) * sin (a1 * dx / sin_alfa) / a1;

  free(p1);  
  free(p2);
  return TRUE;
}
