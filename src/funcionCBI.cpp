#include "funcionCBI.h"

FuncionCBI::FuncionCBI (char * nombreImagen, char **nombresVis, 
			int nVis, int n_pars, double beamSize) : Funcion (n_pars){
  int i;
  double xf, g, bmaj, bmin; //, pixel;
  this->nVis = nVis; 

  this->nombreImagen = new char[256];
  strcpy (this->nombreImagen, nombreImagen);

  this->nombresVis = new char*[nVis];
  for (i = 0; i < nVis; i++){
    this->nombresVis[i] = new char[256];
    strcpy(this->nombresVis[i], nombresVis[i]);
  }
  
  ref_freq = 30.0e9; // 31.5e9;

  // Imagenes
  cmb_image = do_read (nombreImagen);
  if (!cmb_image) {
    fprintf (stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit (1);
  }
  fg_image = do_read (nombreImagen);
  for (i = 0; i < cmb_image->npixels; i++) {
    cmb_image->pixels[i] = 0;
    fg_image->pixels[i] = 0;
  }

  init_beam (&(beam));
  
  if (beamSize == -2) {
    beam.type = CBI2;  
    printf ("Usando beam CBI2.\n");
  }
  else if (beamSize <= 0) {
    beam.type = CBI;
    printf ("Usando beam CBI.\n");
  }
  else {
    beam.type = GAUSS;
    beam.fwhm = beamSize * RPARCM;
    printf ("Usando beam GAUSS de %g arcmin.\n", beamSize);
  }

  // Visibilidades
  cout << "nVis = " << nVis << "\n";
  header_obs = (struct uvf_header**) malloc(nVis * sizeof(struct uvf_header *));
  //header_obs[0] = (struct uvf_header*) malloc(sizeof(struct uvf_header));

  if (header_obs == NULL) {
    cerr << "ERROR en FuncionCBI, header_obs = NULL.\n";
    exit (1);
  }

  samples_obs = (struct uvf_sample **) malloc(nVis * sizeof(struct uvf_sample *));
  if (samples_obs == NULL) {
    cerr << "ERROR en newFuncL, samples_obs = NULL.\n";
    exit (1);
  }
  samples_mod = (struct uvf_sample **) malloc(nVis * sizeof(struct uvf_sample *));
  if (samples_mod == NULL) {
    cerr << "ERROR en newFuncL, samples_mod = NULL.\n";
    exit (1);
  }

  for(i = 0; i < nVis; i++) {
    // abrimos los archivos de visibilidades
    int status = do_read_uvdata(this->nombresVis[i], &(header_obs[i]), &(samples_obs[i]));
    cout << "header[" << i << "]->object[0] = " << header_obs[i]->object[0] << "\n";
    if(status != SUCCESS) {
      cerr << "Error con el archivo uvf\n";
      exit(1);
    }
  }

  for (i = 0; i < nVis; i++) {
    samples_mod[i] = (struct uvf_sample*) malloc((header_obs[i]->nsamp) * sizeof(struct uvf_sample));
    if (samples_mod[i] == NULL) {
      cerr << "ERROR en newFuncL, samples_mod[%d] = NULL.\n";
      exit (1);
    }
  }
  // ya tenemos abiertos los archivos de visibilidades

  xf = (PLANCK_H * ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
  fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(ref_freq);   /* K to Jy sr^-1 */
  double pixel = fabs(fg_image->cdelt[0] * RPDEG * fg_image->cdelt[1] * RPDEG);    /* pixel solid angle in ster */
  fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * pixel * SQR(ref_freq);    /* K to Jy pixel^-1 */

  cout << "fg_scale = " << fg_scale << "\n";
  bmaj = fg_image->bmaj * RPDEG; // en radianes.
  bmin = fg_image->bmin * RPDEG; // en radianes.
  ruido = Irms (header_obs, samples_obs, nVis);
  ruido *= (4 * log(2.0) / (PI * bmaj * bmin)); // [Jy / sr]
  ruido *= pixel;
  ruido /= fg_scale; // en Kelvin.
  cout << "ruido = " << ruido << " [K]\n";

  double Nbeam = (PI * fg_image->bmaj * fg_image->bmin / (4 * log(2.0))/ 
		  fabs (fg_image->cdelt[0] * fg_image->cdelt[1]));
  printf ("Nbeam = %g\n", Nbeam);
  ruido *= sqrt (Nbeam);

  //ruido /= 2.0; //un medio del ruido real.
  //ruido /= 100;

  atten = atenuacion(fg_image, header_obs, nVis, beam);

  noise_image = NULL;
  calcularNoiseImage();

}

void FuncionCBI::resizeImage (int nx, int ny, int cx, int cy) {
  for (int i = 0; i < nVis; i++) {
    for (int j = 0; j < fg_image->npixels; j++) {
      free (atten[i][j + 1]);
    }
    free (atten[i]);
  }
  free (atten);

  resize_map(fg_image , nx, ny, cx, cy);
  resize_map(cmb_image , nx, ny, cx, cy);
  setNPars(nx * ny);

  atten = atenuacion(fg_image, header_obs, nVis, beam);  
}

/*--------------------------------------------------------------------
 * Calculo del factor de correlacion entre visibilidades de samples_31
 * y samples_simu.
 *--------------------------------------------------------------------*/

double FuncionCBI::cross_correlate_factor () {
  int arch, samp, chan;
  double c1 = 0, c2 = 0;

  for (arch = 0; arch < nVis; arch++) {
    for (chan = 0; chan < header_obs[arch]->nchan; chan++) {
      for (samp = 0; samp < header_obs[arch]->nsamp; samp++) {
	c1 += ((samples_mod[arch][samp].rdata[chan * 3] * 
		samples_obs[arch][samp].rdata[chan * 3] + 
		samples_mod[arch][samp].rdata[chan * 3 + 1] * 
		samples_obs[arch][samp].rdata[chan * 3 + 1]) *
	       samples_obs[arch][samp].rdata[chan * 3 + 2]);
	       
	c2 += ((samples_mod[arch][samp].rdata[chan * 3] * 
		samples_mod[arch][samp].rdata[chan * 3] +
		samples_mod[arch][samp].rdata[chan * 3 + 1] * 
		samples_mod[arch][samp].rdata[chan * 3 + 1]) *
	       samples_obs[arch][samp].rdata[chan * 3 + 2]);
      }
    }
  }
  return c1 / c2;
}

void FuncionCBI::guardarResiduos (string info) {
  string nombre_archivo;
  struct uvf_sample *samples_res;
  struct uvf_header *header_res;

  for (int i = 0; i < nVis; i++) {
    int status = do_read_uvdata(nombresVis[i], &(header_res), &(samples_res));
    if(status != SUCCESS) {
      cerr << "Error con el archivo uvf " << i << " en FuncionCBI::guardarResiduos\n";
      exit(1);
    }

    residual_vis (header_obs[i], samples_obs[i], samples_mod[i], samples_res);

    nombre_archivo = "!RES_" + info + "_" + i2s(i, 2) + ".sub";
    
    char *nombreC = new char [nombre_archivo.length() + 1];
    nombre_archivo.copy(nombreC, string::npos);
    nombreC [nombre_archivo.length()] = 0;
    
    status = do_write("write", nombreC, header_res, samples_res, nombresVis[i]);
    if(status != SUCCESS) {
      cerr << "Error al escribir archivo uvf " << i << " en FuncionCBI::guardarResiduos\n";
      exit(1);
    }
    
    free (samples_res);
    free (header_res);
    delete [] nombreC;
  }
}

void FuncionCBI::calcularNoiseImage(){
  int nx = fg_image->size[0], ny = fg_image->size[1];
  double atten_image;

  if (noise_image != NULL) {
    delete_map (noise_image);
  }

  noise_image = new_map();
  copy_changed_map(noise_image, fg_image, fg_image->pixels);

  for (int i = 0; i < noise_image->npixels; i++) {
    noise_image->pixels[i] = 0;
  }
  cout << "ruido : " << ruido << "\n";

  for (int i = 0; i < nVis; i++) {
    for (int k = 0; k < noise_image->npixels; k++) {
      atten_image = 0;
      for (int j = 0; j < header_obs[i]->nchan; j++)	{
	  atten_image += atten[i][k + 1][j];
      }
      atten_image /= header_obs[i]->nchan;
      noise_image->pixels[k] += SQR(atten_image / ruido);
      //cout << "atten_image =  " << atten_image << "\n";
      //noise_image->pixel[i] += SQR(atten_image->pixels[k - 1] / difmap_noise[i]);
      //cout << "noise_image->pixels[i] " << noise_image->pixels[k] << "\n";
    }
  }

  for (int k = 0; k < noise_image->npixels; k++) {
    noise_image->pixels[k] = sqrt (1 / noise_image->pixels[k]);
  }

  do_write_fits (noise_image, "!noise_image.fits");
  //exit(0);
}

FuncionCBI::~FuncionCBI () {
  int i, j;
  
  delete [] nombreImagen;
  for (i = 0; i < nVis; i++){
    delete [] nombresVis[i]; 
    free (samples_obs[i]);
    free (header_obs[i]);
    free (samples_mod[i]);
  }
  delete [] nombresVis;

  free (header_obs);
  free (samples_obs);
  free (samples_mod);

  for (i = 0; i < nVis; i++) {
    for (j = 0; j < fg_image->npixels; j++) {
      free (atten[i][j + 1]);
    }
    free (atten[i]);
  }
  free (atten);

  delete_map (cmb_image);
  delete_map (fg_image);
  delete_map (noise_image);

  cout << "ELIMINANDO FUNCIONCBI!!!!!\n";
}

FuncionCBI& FuncionCBI::operator = (const FuncionCBI& func) {
  int i;
  this->nVis = func.nVis;
  
  this->nombreImagen = new char[256];
  strcpy (this->nombreImagen, func.nombreImagen);

  this->nombresVis = new char*[nVis];
  for (i = 0; i < nVis; i++){
    this->nombresVis[i] = new char[256];
    strcpy(this->nombresVis[i], func.nombresVis[i]);
  }

  return *this;
}
