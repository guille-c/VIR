#ifndef VIS_FUNC_L
#define VIS_FUNC_L

#include "mallaVoronoi.h"
#include "rasterizador.h"
#include "mockcbiRoutines.h"
#include "image_routines.h" /* For reading FITS images */
#include "uvsubs.h"         /* UVFITS routines */
#include "newstring.h"      /* defines newstring() */
#include "nrutil.h"

#ifdef __cplusplus
extern "C" {
#endif  

enum integral {EXACTA, APROX_DX, APROX_VORONOI};
enum mock {MOCKCBI, MOCK_EXACTO};

#define MIN_PIX 1e-15
#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.
//#define BORDE_INF 1  ver poligono.h
//#define BORDE_DER 2
//#define BORDE_SUP 3
//#define BORDE_IZQ 4

  typedef struct{
    MallaVoronoi *malla;
    //int nx, ny;
    double **imagen, chi2, S;
    struct image *cmb_image;         // Imagen en blanco
    struct image *fg_image;          // Parametros libres.
    struct uvf_header **header_obs;  // Arreglo con los encabezados.
    struct uvf_sample **samples_obs; // Arreglo con las visibilidades
    // observadas.
    struct uvf_sample **samples_mod; // Arreglo con las visibilidades 
    // modeladas.
    int n_pols;                      // Numero de poligonos.
    int n_archivos;                  // Numero de archivos uvf.
    struct pbeam beam;               // Beam del CBI
    double ***atten;
    // Arreglo de 3 dimensiones con las atenuaciones. La 1era es
    // los archivos de visibilidades, la 2da las coordenadas de la
    // imagen y la 3ra los canales.
    double ***sin;
    // Arreglos de 3 dimensiones con el calculo de sin(PI u dx)sin(PI v dy).
    // La 1era es los archivos de visibilidades, la 2da las 
    // lineas de base y la 3ra los canales.
    char **infile;   // Arreglo con los nombres de los archivos de vis.
    int **mask;
    double ref_freq, fg_scale;
    int *x0, *y0;    // Arreglos con x0 e y0 para cada archivo.
    double *dx, *dy; // Arreglos con dx y dy para cada archivo.
    double **fourierI;
    double difmapNoise;
    double expanded_lx, expanded_ly;
    
    // Para debugging
    int nPixIntegrados_I;
    int **PixIntegrados_I;
    int **PixIntegrados_x;
    char *nombreFits;
    int entropia;
  }funcL;
  
  funcL * newFuncL (char * nombreImagen, char **nombresVis, int nVis, 
		    int n, int entropia, double cuantaSize, 
		    double expanded_lx, double expanded_ly, double beamSize);
  double L (funcL *L, double *pars, int n, int tipoMock);
  double S (funcL *L, double *pars, int n);
  void dL (funcL *L, double *pars, double *grad, int n, int aprox, int tipoMock);
  void dS (funcL *fL, double *pars, double *grad, int n);
  
  void calcularSin (funcL *L);
  void calculardxdy (funcL *L);
  
  void dLdI (funcL *fL, double *grad, int n, int aprox);
  void dLdx (funcL *fL, double *grad, int n, int aprox);
  
  void dLdI_Voronoi (funcL *fL, double *grad, int n);
  void dVdx (funcL *fL, int arch, int chan, int samp, PoligonoVoronoi *pol, 
	     double *dVdx_R, double *dVdx_I, double *dVdy_R, double *dVdy_I,
	     int aprox);
  void integraldVdx (funcL *fL, int arch, int chan, int samp, 
		     PoligonoVoronoi *pol,  AristaVoronoi *a,
		     double s0, double c1, double c2, double M, double b,
		     double *integralR, double *integralI);
  void integraldVdxAprox (funcL *fL, int arch, int chan, int samp, 
			  PoligonoVoronoi *pol,  AristaVoronoi *a,
			  double s0, double c1, double c2, double M, double b,
			  double *integralR, double *integralI);
  int encontrarSgtePtoPixel (funcL *fL, PuntoVoronoi *p1, PuntoVoronoi *pFin,
			     PuntoVoronoi *pSgte, int lineaAnterior);
  
  void calcularFourierI (funcL *fL);
  void calcularFourierIAprox (funcL *fL);
  /* void dVdI (funcL *fL, int arch, int chan, int samp,  */
  /* 	   PoligonoVoronoi *pol, double *dVdI_R, double *dVdI_I); */
  /* //void integraldVdI (funcL *fL, int arch, int chan, int samp,  */
  /* //	   PoligonoVoronoi *pol, double *dVdI_R, double *dVdI_I, */
  /* //	   int i, int j); */
  /* void integraldVdIH (funcL *fL, int arch, int chan, int samp,  */
  /* 		    PoligonoVoronoi *pol, double *dVdI_R, double *dVdI_I, */
  /* 		    int i_ini, int j_ini, int *i1, int *i2); */
  void eliminarFuncL (funcL *L);
  funcL *copiarFuncL (funcL *fL_ori);
  void normalizarVisibilidades (struct uvf_header **header,
				struct uvf_sample **samples, int narch);
  void guardarFits (double **im, int **mask, 
		    char *nombre, char *nombreFits);

  void expanded_mockcbi(funcL *fL, int arch);
  void mockCBI (funcL *fL);
  void visMod (funcL *fL, int arch, int iff, int samp, double *visR, double *visI);
  //funcL *copiar (funcL *fL);
  
#ifdef __cplusplus
}
#endif

#endif
