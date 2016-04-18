#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcL.h"

int comparar_strings (char *s1, char *s2) {
  if (s1 == NULL) {
    s1 = "NULL";
  }
  if (s2 == NULL) {
    s2 = "NULL";
  }
  return strcmp (s1, s2);
}

void comparar_header (struct image *i1, struct image *i2)
{
  int i;

  if (i1->type != i2->type) {
    printf ("type:\t%d\t%d\n", i1->type, i2->type);
  }
  if (i1->bitpix != i2->bitpix) {
    printf ("bitpix:\t%d\t%d\n", i1->bitpix, i2->bitpix);
  }
  if (i1->naxis != i2->naxis) {
    printf ("naxis:\t%d\t%d\n", i1->naxis, i2->naxis);
  }
  for (i = 0; i < i1->naxis; i++) {
    if (i1->size[i] != i2->size[i]) {
      printf ("size[%d]:\t%d\t%d\n", i, i1->size[i], i2->size[i]);
    }
    if (i1->crval[i] != i2->crval[i]) {
      printf ("crval[%d]:\t%.99g\t%.99g\n", i, i1->crval[i], i2->crval[i]);
    }
    if (i1->crpix[i] != i2->crpix[i]) {
      printf ("crpix[%d]:\t%.99g\t%.99g\n", i, i1->crpix[i], i2->crpix[i]);
    }
    if (i1->cdelt[i] != i2->cdelt[i]) {
      printf ("cdelt[%d]:\t%.99g\t%.99g\n", i, i1->cdelt[i], i2->cdelt[i]);
    }
    if (i1->crota[i] != i2->crota[i]) {
      printf ("crota[%d]:\t%.99g\t%.99g\n", i, i1->crota[i], i2->crota[i]);
    }
    if (comparar_strings (i1->ctype[i], i2->ctype[i])) {
      printf ("ctype[%d]:\t%s\t%s\n", i, i1->ctype[i], i2->ctype[i]);
    }
  }
  if (i1->npixels != i2->npixels) {
    printf ("npixels:\t%d\t%d\n", i1->npixels, i2->npixels);
  }
  if (i1->niter != i2->niter) {
    printf ("niter:\t%d\t%d\n", i1->niter, i2->niter);
  }
  if (i1->bmaj != i2->bmaj) {
    printf ("bmaj:\t%.99g\t%.99g\n", i1->bmaj, i2->bmaj);
  }
  if (i1->bmin != i2->bmin) {
    printf ("bmin:\t%.99g\t%.99g\n", i1->bmin, i2->bmin);
  }
  if (i1->bpa != i2->bpa) {
    printf ("bpa:\t%.99g\t%.99g\n", i1->bpa, i2->bpa);
  }
  if (i1->equinox  != i2->equinox) {
    printf ("equinox:\t%.99g\t%.99g\n", i1->equinox, i2->equinox);
  }
  if (i1->freq != i2->freq) {
    printf ("freq:\t%.99g\t%.99g\n", i1->freq, i2->freq);
  }
  if (i1->noise != i2->noise) {
    printf ("noise:\t%.99g\t%.99g\n", i1->noise, i2->noise);
  }
  if (i1->datamax != i2->datamax) {
    printf ("datamax:\t%.99g\t%.99g\n", i1->datamax, i2->datamax);
  }
  if (i1->datamin != i2->datamin) {
    printf ("datamin:\t%.99g\t%.99g\n", i1->datamin, i2->datamin);
  }
  if (i1->projection != i2->projection) {
        printf ("projection:\t%d\t%d\n", i1->projection, i2->projection);
  }
  if (i1->obsra != i2->obsra) {
      printf ("obsra:\t%.99g\t%.99g\n", i1->obsra, i2->obsra);
  }
  if (i1->obsdec != i2->obsdec) {
    printf ("obsdec:\t%.99g\t%.99g\n", i1->obsdec, i2->obsdec);
  }
  if (comparar_strings (i1->filename, i2->filename)) {
    printf ("filename:\t%s\t%s\n", i1->filename, i2->filename);
  }
  if (comparar_strings (i1->object, i2->object)) {
    printf ("object:\t%s\t%s\n", i1->object, i2->object);
  }
  if (comparar_strings (i1->telescope, i2->telescope)) {
    printf ("telescope:\t%s\t%s\n", i1->telescope, i2->telescope);
  }
  if (comparar_strings (i1->instrument, i2->instrument)) {
    printf ("instrument:\t%s\t%s\n", i1->instrument, i2->instrument);
  }
  if (comparar_strings (i1->observer, i2->observer)) {
    printf ("observer:\t%s\t%s\n", i1->observer, i2->observer);
  }
  if (comparar_strings (i1->date_obs, i2->date_obs)) {
    printf ("date_obs:\t%s\t%s\n", i1->date_obs, i2->date_obs);
  }
  if (comparar_strings (i1->timesys, i2->timesys)) {
    printf ("timesys:\t%s\t%s\n", i1->timesys, i2->timesys);
  }
  if (comparar_strings (i1->date_map, i2->date_map)) {
    printf ("date_map:\t%s\t%s\n", i1->date_map, i2->date_map);
  }
  if (comparar_strings (i1->bunit, i2->bunit)) {
    printf ("bunit:\t%s\t%s\n", i1->bunit, i2->bunit);
  }
  if (comparar_strings (i1->map_type, i2->map_type)) {
    printf ("map_type:\t%s\t%s\n", i1->map_type, i2->map_type);
  }
  if (comparar_strings (i1->origin, i2->origin)) {
    printf ("origin:\t%s\t%s\n", i1->origin, i2->origin);
  }
  if (comparar_strings (i1->stokes, i2->stokes)) {
    printf ("stokes:\t%s\t%s\n", i1->stokes, i2->stokes);
  }
}

int main(int argc, char **argv) {
  struct image *im1 = do_read (argv[1]);
  struct image *im2 = do_read (argv[2]);

  comparar_header (im1, im2);
  return 0;
}
