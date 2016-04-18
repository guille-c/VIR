#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <cstdlib>
#include <fstream.h>

#include "mockcbiRoutines.h"
#include "mcheck.h"

int main (int argc, char **argv){
  mtrace();
  
  struct pbeam beam;               // Beam del CBI
  init_beam (&(beam));
  beam.type = CBI;

  if (argc < 4) {
    cerr << "uso: nombreImagen.fits nombreVis.sub nombreRes.sub\n";
    exit (1);
  }

  struct image *im = do_read (argv[1]);
  struct image *cmb = do_read (argv[1]);

  resize_map (im, 256, 256, im->crpix[0], im->crpix[1]); 
  resize_map (cmb, 256, 256, cmb->crpix[0], cmb->crpix[1]);

  for (int i = 0; i < cmb->npixels; i++) {
    cmb->pixels[i] = 0;
  }

  struct uvf_header *header_obs;  
  struct uvf_sample *samples_obs; 
  struct uvf_sample *samples_mod; 
  struct uvf_sample *samples_res; 

  int status = do_read_uvdata(argv[2], &header_obs, &samples_obs);
  if(status != SUCCESS) {
    cerr << "Error con el archivo uvf\n";
    exit(1);
  }

  samples_mod = (struct uvf_sample*) malloc((header_obs->nsamp) * sizeof(struct uvf_sample));
  mockcbi_sub (cmb, im, header_obs, samples_obs, samples_mod, beam);

  samples_res = (struct uvf_sample*) malloc((header_obs->nsamp) * sizeof(struct uvf_sample));
  residual_vis (header_obs, samples_obs, samples_mod, samples_res);

  status = do_write("write", argv[3], header_obs, samples_res, argv[2]);
  if(status != SUCCESS) {
    cerr << "Error al escribir archivo uvf al guardar residuos\n";
    exit(1);
  }


}

