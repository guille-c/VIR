#include "reconstructorVIR.h"

ReconstructorVIR::ReconstructorVIR (char *nombreImagen, char **nombresVis, int nVis, 
				    int nPolsIni, int nPolsFin, int dn, 
				    double expanded_lx, double expanded_ly, int MEMIter,
				    double beamSize){
  this->nVis = nVis;
  this->nPolsIni = nPolsIni;
  this->nPolsFin = nPolsFin;
  this->dn = dn;
  this->MEMIter = MEMIter;
  this->nombreImagen = new char[256];
  this->expanded_lx = expanded_lx;
  this->expanded_ly = expanded_ly;
  this->beamSize = beamSize;
  strcpy (this->nombreImagen, nombreImagen);
  this->nombresVis = new char* [nVis];

  for (int i = 0; i < nVis; i++) {
    this->nombresVis[i] = new char[255];
    strcpy (this->nombresVis[i], nombresVis[i]);
  }
  //p = new double [nPolsIni];
}

void ReconstructorVIR::realizarMEM () {

  // MEM
  struct image *im = do_read (nombreImagen);
  //FuncionEntropiaNatural *funcEN = new FuncionEntropiaNatural(im->npixels);
  Funcion *funcEN = new FuncionCero(im->npixels);
  delete_map(im);
  FuncionMemCBI *funcMemCBI = new FuncionMemCBI (nombreImagen, nombresVis, nVis, funcEN, 1e-5, 0.0, 1e100, beamSize);
  //funcMemCBI->setRuido (funcMemCBI->getRuido() * qFactor);
  Inicializador *ini = new InicializadorConstante ();
  //ReconstructorGCMemCBI rMem (funcMemCBI, ini, 1e-6, 2, 0.0, MEMIter + 1);
  ReconstructorGCMemCBI rMem (funcMemCBI, ini, 1e-6, 2, 0.0, MEMIter);

  rMem.run();

  delete ini;
  delete funcMemCBI;
  delete funcEN;
}

double ReconstructorVIR::run () {
  //GUILLE 31/052007
  FILE *archivo = fopen ("frprmin.log", "w");
  fclose (archivo);

  std::ofstream archivoOut;
  archivoOut.open ("reconstructorVIR.out");
  archivoOut.close ();

  realizarMEM ();

  // Ajustar hasta nIni
  char nombreArchivo [255];
  string sIter = i2s (MEMIter, 3);
  sprintf (nombreArchivo, "MEM_CBIGC%s.fits", sIter.c_str());
  Funcion *func = new FuncionChi2VoronoiImagen (nombreArchivo, 0);
  IncrementadorMallaImagen *inc = new IncrementadorMallaImagen (nombreArchivo, -1);
//  Funcion *func = new FuncionChi2VoronoiImagen ("MEM_CBIGC002.fits", 0);
//  IncrementadorMallaImagen *inc = new IncrementadorMallaImagen ("MEM_CBIGC002.fits", -1);
  InicializadorVoronoiUniforme *iniVU = new InicializadorVoronoiUniforme ();
  ReconstructorIncremental  rInc (func, iniVU, inc, nPolsIni);
  rInc.run();
  
  // GC 
  FuncionCBI *funcCBI = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, 
						    nVis, nPolsIni, 1, -1,
						    expanded_lx, expanded_ly, beamSize);
  //FuncionCBI *funcCBI = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, 
  //					    nVis, nPolsIni, 0, -1);
  ReconstructorGC recGC (funcCBI, iniVU, 1e-10);

  double ret;
  p = new double [func->getNPars()];

  for (int i = 0; i < func->getNPars(); i++) {
    p[i] = rInc.getPar(i);
  }

  for (int nPols = nPolsIni; nPols <= nPolsFin; nPols += dn) {

    for (int j = 0; j < recGC.getNPars() / 3; j++) {
      //cout << "cambiando " << j << "): " << recGC.getPar(j) << " por " << p[j] << "\n";
      recGC.setPar(3 * j, p[3 * j]);
      recGC.setPar(3 * j + 1, p[3 * j + 1]);
      recGC.setPar(3 * j + 2, p[3 * j + 2] / funcCBI->getRuido());
    }
    
    ret = recGC.run();
    
    archivoOut.open("reconstructorVIR.out", ios::app);
    archivoOut << nPols << "\t" << ret << "\n";
    archivoOut.close();

    //sprintf (nombreArchivo, "VIR_NPols_%g_iter", qFactor);
    funcCBI->guardarInfo ("VIR_NPols_" + i2s (nPols, 4));
    funcCBI->guardarResiduos ("VIR_NPols_" + i2s (nPols, 4));

    for (int j = 0; j < dn; j++) {
      func->guardarInfo (i2s (nPols + j, 4));
      inc->incrementar (&(p), func);
    }
    
    recGC.setNPars (func->getNPars());
  }

  archivoOut.close();
  delete func;
  delete funcCBI;
  delete inc;
  delete iniVU;

  for (int i = 0; i < func->getNPars() / 3; i++) {
    p[3 * i + 2]  /= funcCBI->getRuido();
    printf ("p[%d] = %g\n", 3 * i + 2, p[3 * i + 2] );
  }

  return ret;
}

ReconstructorVIR::~ReconstructorVIR () {
  delete [] this->nombreImagen;

  for (int i = 0; i < nVis; i++) {
    delete [] this->nombresVis[i];
  }
  delete [] this->nombresVis;
}
