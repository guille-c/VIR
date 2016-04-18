#include "funcionChi2VoronoiImagen.h"

FuncionChi2VoronoiImagen::FuncionChi2VoronoiImagen (char * nombreImagen, int n_pols) {
  im = do_read (nombreImagen);
  if (im == NULL) {
    cerr << "ERROR: " << nombreImagen << " no corresponde a una imagen valida\n";
    exit (1);
  }
  n_pars = 3 * n_pols;
  malla = newMallaVoronoi ();
  this->nombreImagen = new char [255];
  strcpy (this->nombreImagen, nombreImagen);
  cout << "nombreImagen " << this->nombreImagen << "\n";
}

double FuncionChi2VoronoiImagen::f (double pars[]) {
  double **imagen, ret = 0;
  int nx = im->size[0], ny = im->size[1];

  // Verificamos que los poligonos esten dentro del cuadrado.
  for (int i = 0; i < n_pars / 3; i++) {
    if (pars[3 * i] < -1e-5 || pars[3 * i] > 1 ||
	pars[3 * i + 1] < -1e-5 ||  pars[3 * i + 1] > 1) {
      cerr << "ERROR en FuncionChi2VoronoiImagen::f\n";
      cerr << "   pol [" << i << "] = {" << pars[3 * i] << ", " << pars[3 * i + 1] << ", "
	   << pars[3 * i + 2] << "}\n";
      cerr << "   n_pars = " << n_pars << "\n";
      //return 1e300;
    }
  }
  
  if (malla != NULL) {
    eliminarMallaVoronoi (malla);
  }

  malla = newMallaVoronoi ();
  for (int i = 0; i < n_pars / 3; i++) {
    //cout << "insertando {" << pars[3 * i] << ", " << pars[3 * i + 1] << ", "
    // << pars[3 * i + 2] << "}\n";
    insertarSitio (malla, &(pars[3 * i]), &(pars[3 * i + 1]), pars[3 * i + 2]);
  }

  imagen = toImage (malla, nx, ny, NULL);
  
  /*
  struct image *fg_image = do_read ("../../gaussianas3.fits");
  for (int i = 0; i < fg_image->size[0]; i++) {
    for (int j = 0; j < fg_image->size[1]; j++) {
      fg_image->pixels[i + j * fg_image->size[0]] = imagen[i][j];
    }
  }
  do_write_fits (fg_image, "!auxF.fits");
  do_write_fits (im, "!auxIm.fits");
  */

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      ret += pow (imagen [i][j] - im->pixels[i + j * nx], 2);
      //cout << "restando " << imagen [i][j] << " - " << im->pixels[i + j * nx] << "\n";
    }
    free (imagen[i]);
    //delete [] (imagen[i]);
  }
  free (imagen);
  //delete [] imagen;
  //cout << "ret f = " << ret << "\n";
  //cout << "ret 1/f = " << 1/ret << "\n";
  return ret;
}

void FuncionChi2VoronoiImagen::df (double pars[], double grad[]) {
  // No se :P
}

void FuncionChi2VoronoiImagen::guardarInfo (string info) {
  // Para este caso usaremos info como el sufijo de los archivos a
  // guargar.
  string nombre_archivo;
  struct image *im_aux = do_read (this->nombreImagen);
  int nx = im_aux->size[0], ny = im_aux->size[1];

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      double x = i / (nx - 1.0);
      double y = j / (ny - 1.0);
      im_aux->pixels[i + j * nx] = (encontrarPoligono (malla, x, y)->valor);
      //cout << im_aux->pixels[i + j * nx] << " ";
    }
  }
  cout << "\n";
  
  nombre_archivo = "!FuncionChi2VorIm" + info + ".fits";

  char *nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  

  cout << "guardando " << nombreC << "\n";
  do_write_fits (im_aux, nombreC);

  nombre_archivo = "MallaChi2VorIm" + info + ".dat";
  delete [] nombreC;
  nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;

  if (malla != NULL) {
    imprimirMallaArchivo (malla, nombreC);
  }

  delete [] nombreC;
  delete_map (im_aux);
}

FuncionChi2VoronoiImagen::~FuncionChi2VoronoiImagen () {
  delete_map (im);
  if (malla != NULL) {
    eliminarMallaVoronoi (malla);
    malla = NULL;
  }
  delete [] nombreImagen;
}

