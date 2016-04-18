#include "distribucionCondicionalCBI.h"

DistribucionCondicionalCBI::DistribucionCondicionalCBI (char *nombreImagen, char **nombresVis, int nVis) {
  double xf, g, bmaj, bmin; //, pixel;
  this->nVis = nVis; 
  int i;

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
  beam.type = CBI;

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

  bmaj = fg_image->bmaj * RPDEG; // en radianes.
  bmin = fg_image->bmin * RPDEG; // en radianes.
  ruido = Irms (header_obs, samples_obs, nVis);
  ruido *= (4 * log(2.0) / (PI * bmaj * bmin)); // [Jy / sr]
  ruido *= pixel;
  ruido /= fg_scale; // en Kelvin.
  //ruido /= 2.0; //un medio del ruido real.

  cout << "ruido = " << ruido << " [K].\n";

  atten = atenuacion(fg_image, header_obs, nVis, beam);
}

void DistribucionCondicionalCBI::guardarInfo (string info) {
      // Para este caso usaremos info como el sufijo de los archivos a
  // guardar.
  string nombre_archivo;
  
  nombre_archivo = "!MEM_CBI" + info + ".fits";

  char *nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  
  cout << "guardando " << nombreC << "\n";
  do_write_fits (fg_image, nombreC);

  delete [] nombreC;

  nombre_archivo = "MEM_CBI" + info + ".dat";

  nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  
  cout << "guardando " << nombreC << "\n";

  ofstream archivo (nombreC);
  
  archivo << n_pars << "\n";
  for (int i = 0; i < n_pars; i++) {
    archivo << fg_image->pixels[i] << "\n";
  }

  delete [] nombreC;
}

DistribucionCondicionalCBI::~DistribucionCondicionalCBI() {
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
}
