#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "random.h"

#define POW(x) ((x)*(x))

int main(int argc, char **argv) {
  int i, samp, iff;
  struct image *fg = do_read (argv[1]);
  struct image *cmb = do_read (argv[1]);
  //struct image *imagen_out = do_read (argv[2]);
  struct uvf_header *header_obs;  // Arreglo con los encabezados.
  struct uvf_sample *samples_obs; // Arreglo con las visibilidades
                                  // observadas.
  struct uvf_sample *samples_mod; // Arreglo con las visibilidades 
  struct pbeam beam;               // Beam del CBI
  FILE *archivo = fopen (argv[5], "w");
  double radio, u, v;

  init_beam (&beam);
  beam.type = CBI;

  do_read_uvdata(argv[2], &header_obs, &samples_obs);
  samples_mod = (struct uvf_sample*) malloc((header_obs->nsamp) * 
					    sizeof(struct uvf_sample));


  for (i = 0; i < cmb->npixels; i++) {
    cmb->pixels[i] = 0.0;
  }

  
  do_write_fits (cmb, "!cmb.fits");
  do_write_fits (fg, argv[3]);
  

  mockcbi_sub (cmb, fg, header_obs, samples_obs, samples_mod,
	       beam);

  do_write ("write uvf", argv[4], header_obs,
	    samples_mod, argv[2]);
  for (iff = 0; iff < header_obs->nif; iff++) {
    for (samp = 0; samp < header_obs->nsamp; samp++) {
      u = samples_obs[samp].u * header_obs->iffreq[iff];
      v = samples_obs[samp].v * header_obs->iffreq[iff];
      
      radio = sqrt (u * u + v * v);
      fprintf (archivo, "%g\t%g\n", radio,
	       sqrt (pow (samples_obs[samp].rdata[iff * 3], 2) + 
	       pow (samples_obs[samp].rdata[iff * 3 + 1], 2)));
      
    }
  }

  
  fclose (archivo);
  return 0;
}
