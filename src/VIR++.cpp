#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"

#include "funcionMemCBI.h"
#include "funcionBayesVoronoiCBI.h"
#include "funcionBayesVoronoiCBII.h"
#include "funcionVoronoiI.h"
#include "funcionVoronoiX.h"
#include "funcionPbbVoronoiCBI.h"
#include "funcionNegativa.h"
#include "funcionDividida.h"
#include "funcionExp.h"
#include "funcionEjemploH.h"
#include "funcionReconstructor.h"
#include "funcionChi2VoronoiImagen.h"
#include "funcionChi2VoronoiImagenN.h"

#include "funcionCero.h"
#include "funcionEntropiaNatural.h"
#include "funcionEntropiaIlogI.h"

#include "frprmn.h"
#include "gauss.h"

#include "reconstructor.h"
#include "reconstructorGC.h"
#include "reconstructorGCMemCBI.h"
#include "reconstructorBayesiano.h"
#include "reconstructorGenetico.h"
#include "reconstructorIncremental.h"
#include "reconstructorIndividual.h"
#include "reconstructorVIR.h"
#include "reconstructorVIRIterativo.h"

#include "incrementadorMallaPixel.h"
#include "incrementadorMallaAristas.h"
#include "incrementadorMallaImagen.h"

#include "inicializadorConstante.h"
#include "inicializadorArchivo.h"
#include "inicializadorVoronoiUniforme.h"
#include "inicializadorVoronoiPixelizado.h"
#include "inicializadorVoronoiPixelizadoAlternado.h"

#include "distribucionNormal.h"
#include "distribucionMetropolis.h"
#include "distribucionHamiltoniana.h"
#include "distribucionGibbs.h"
#include "distribucionCondicional.h"
#include "distribucionCondicionalBVCBI.h"
#include "distribucionCondicionalVoronoiI.h"
#include "distribucionCondicionalMemCBI.h"
#include "mcheck.h"

#include "funcionDummy.h"

enum reconst {GC, GENETICO, GENETICOGC, BAYES, INCREMENTAL, INCREMENTALPIX, VIR,
	      INCREMENTALARISTAS, AJUSTE, GIBBS, INDIVIDUAL, MEM, GIBBSMEM, NORMGCMEM,
              VIRITER};

enum funcs {MEMCBI, ILOGI, SNATURAL};

Reconstructor *crearMEM (std::ifstream *archivoRec);
Reconstructor *crearIncremental (std::ifstream *archivoRec);
Reconstructor *crearVIR (std::ifstream *archivoRec);
Reconstructor *crearVIRIterativo (std::ifstream *archivoRec);
Reconstructor *crearGenetico (std::ifstream *archivoRec);
ReconstructorGC *crearGC (std::ifstream *archivoRec);
Reconstructor *crearGeneticoGC (std::ifstream *archivoRec);
Reconstructor *crearBayes (std::ifstream *archivoRec);
Reconstructor *crearGibbs (std::ifstream *archivoRec);

Funcion *func = NULL;
FuncionBayesVoronoiCBI *funcBVCBI = NULL;
FuncionMemCBI *funcMemCBI = NULL;
Funcion *fEnt = NULL;
FuncionReconstructor *funcR = NULL;
Inicializador *ini = NULL;
Incrementador *inc = NULL;
ReconstructorGC *recGC = NULL;
DistribucionMetropolis *distMet = NULL;
DistribucionNormal *distNorm = NULL;
DistribucionCondicionalBVCBI *distCond;
DistribucionGibbs *distGibbs;
funcs tipoFunc;
funcs tipoEnt;


int main (int argc, char **argv) {
  reconst tipoRec;
  //char dato[255];
  std::ifstream archivoRec;

  Reconstructor *r;
  archivoRec.open (argv[2]);
  if (strcmp (argv[1], "VIR") == 0) {
    tipoRec = VIR;
    cout << "VIR\n";
    
    r = crearVIR (&archivoRec);
  }
  else if (strcmp (argv[1], "VIRITER") == 0) {
    tipoRec = VIRITER;
    cout << "VIR\n";
    
    r = crearVIRIterativo (&archivoRec);
  }
  else if (strcmp (argv[1], "MEM") == 0) {
    tipoRec = MEM;
    cout << "MEM\n";
    
    r = crearMEM (&archivoRec);
  }
  else if (strcmp (argv[1], "INCREMENTAL") == 0) {
    tipoRec = INCREMENTAL;
    cout << "INCREMENTAL\n";
    
    r = crearIncremental (&archivoRec);
  }
  else if (strcmp (argv[1], "GEN") == 0) {
    tipoRec = GENETICO;
    cout << "GENETICO\n";
    
    r = crearGenetico (&archivoRec);
  }
  else if (strcmp (argv[1], "GC") == 0) {
    tipoRec = GC;
    cout << "GC\n";
    
    r = crearGC (&archivoRec);
  }
  else if (strcmp (argv[1], "GENGC") == 0) {
    tipoRec = GENETICOGC;
    cout << "GENETICOGC\n";
    
    r = crearGeneticoGC (&archivoRec);
  }
  else if (strcmp (argv[1], "BAYES") == 0) {
    tipoRec = BAYES;
    cout << "BAYES\n";
    
    r = crearBayes (&archivoRec);
  }
  else if (strcmp (argv[1], "GIBBS") == 0) {
    tipoRec = GIBBS;
    cout << "GIBBS\n";
    
    r = crearGibbs (&archivoRec);
  }
  else {
    cerr << "ERROR: tipo de reconstructor inexistente\n";
    cerr << "       use: VIR, MEM, INCREMENTAL, GEN, GC, GENGC, BAYES o GIBBS\n";
    exit (1);
  }

  r->run();

  if (strcmp (argv[1], "VIR") == 0) {
    double *dchi2, *DS, *p;
    ReconstructorVIR *rVIR = (ReconstructorVIR *) r;
    int n = rVIR->getNPolsIni() * 3 ;
    cout << "11) n = " << n << " vs " << r->getNPars() << "\n";
    FuncionBayesVoronoiCBI *f = 
      new FuncionBayesVoronoiCBI (rVIR->getNombreImagen(), rVIR->getNombresVis(), 
				  1, n/3, 1, -1, rVIR->getExpandedLx(), rVIR->getExpandedLy());
    funcL *fL = f->getFL();
    dchi2 = new double [n];
    DS = new double [n];
    p = rVIR->getPars();

    dS (fL, p, DS, n/3);
    f->df(p, dchi2);

    std::ofstream archivoL ("dL.dat");
    std::ofstream archivoChi2 ("dchi2.dat");
    std::ofstream archivoS ("dS.dat");
    for (int i = 0; i < n; i++) {
      cout << "derivadas [" << i << "] " << dchi2[i] << " " << DS[i] << "; p = " << p[i] << "\n";
      archivoL << dchi2[i] << "\n";
      archivoS << DS[i] << "\n";
      archivoChi2 << dchi2[i] - DS[i] << "\n";
    }
    archivoL.close();
    archivoS.close();
    archivoChi2.close();
  }  

  delete r;
  if (func) delete func;
  if (funcBVCBI) delete funcBVCBI;
  if (funcMemCBI) delete funcMemCBI;
  if (fEnt) delete fEnt;
  if (funcR) delete funcR;
  if (ini) delete ini;
  if (inc) delete inc;
  if (recGC) delete recGC;
  if (distMet) delete distMet;
  if (distNorm) delete distNorm;
  if (distCond) delete distCond;
  if (distGibbs) delete distGibbs;
}

Reconstructor *crearVIR (std::ifstream *archivoRec) {
  char nombreImagen[255], **nombresVis, dato [255];
  int nVis, nPolsIni, nPolsFin, dn, MEMIter;
  double expanded_lx, expanded_ly, beamSize;

  (*archivoRec) >> dato >> nombreImagen;
  cout << "nombreImagen: " << nombreImagen << "\n";
  (*archivoRec) >> dato >> nPolsIni;
  cout << "nPolsIni: " << nPolsIni << "\n";
  (*archivoRec) >> dato >> nPolsFin;
  cout << "nPolsFin: " << nPolsFin << "\n";
  (*archivoRec) >> dato >> dn;
  cout << "dn: " << dn << "\n";
  (*archivoRec) >> dato >> MEMIter;
  cout << "MEMIter: " << MEMIter << "\n";
  (*archivoRec) >> dato >> expanded_lx;
  cout << "expandedlx: " << expanded_lx << "\n";
  (*archivoRec) >> dato >> expanded_ly;
  cout << "expandedly: " << expanded_ly << "\n";
  (*archivoRec) >> dato >> beamSize;
  cout << "beamSize: " << beamSize << " arcmin\n";
  (*archivoRec) >> dato >> nVis;

  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis [i] = new char [255];
    (*archivoRec) >> nombresVis [i];
  }

  return new ReconstructorVIR (nombreImagen, nombresVis, nVis, 
			       nPolsIni, nPolsFin, dn,
			       expanded_lx, expanded_ly, MEMIter, beamSize);
}

Reconstructor *crearVIRIterativo (std::ifstream *archivoRec) {
  char nombreImagen[255], **nombresVis, dato [255];
  int nVis, nPolsIni, nPolsFin, dn, MEMIter;
  double expanded_lx, expanded_ly, beamSize;

  (*archivoRec) >> dato >> nombreImagen;
  cout << "nombreImagen: " << nombreImagen << "\n";
  (*archivoRec) >> dato >> nPolsIni;
  cout << "nPolsIni: " << nPolsIni << "\n";
  (*archivoRec) >> dato >> nPolsFin;
  cout << "nPolsFin: " << nPolsFin << "\n";
  (*archivoRec) >> dato >> dn;
  cout << "dn: " << dn << "\n";
  (*archivoRec) >> dato >> MEMIter;
  cout << "MEMIter: " << MEMIter << "\n";
  (*archivoRec) >> dato >> expanded_lx;
  cout << "expandedlx: " << expanded_lx << "\n";
  (*archivoRec) >> dato >> expanded_ly;
  cout << "expandedly: " << expanded_ly << "\n";
  (*archivoRec) >> dato >> beamSize;
  cout << "beamSize: " << beamSize << " arcmin\n";
  (*archivoRec) >> dato >> nVis;

  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis [i] = new char [255];
    (*archivoRec) >> nombresVis [i];
  }

  return new ReconstructorVIRIterativo (nombreImagen, nombresVis, nVis, 
			       nPolsIni, nPolsFin, dn,
			       expanded_lx, expanded_ly, MEMIter, beamSize);
}

Reconstructor *crearIncremental (std::ifstream *archivoRec) {
  char nombreFits [255], nombreArc [255], nombreFunc [255], dato [255];
  int nPolsFin;
  Reconstructor *ret;
  
  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPolsFin;
  (*archivoRec) >> dato >> nombreFunc;

  cout << "nombreFits: " << nombreFits << "\n";
  cout << "nPolsFin: " << nPolsFin << "\n";
  cout << "NombreFunc " << nombreFunc << "\n";

  if (strcmp(nombreFunc, "AJUSTE") == 0) {
    char nombreFits [255];
    int nPols;

    (*archivoRec) >> dato >> nombreFits;
    (*archivoRec) >> dato >> nPols;
    
    func = new FuncionChi2VoronoiImagen (nombreFits, nPols);
    inc = new IncrementadorMallaImagen (nombreFits, -1);
    ini = new InicializadorVoronoiUniforme ();
    
    ret = new ReconstructorIncremental (func, ini, inc, nPolsFin);
  }
  else if (strcmp(nombreFunc, "GCARISTAS") == 0) {
    char nombreFits [255];
    int nPols;

    (*archivoRec) >> dato >> nombreArc;
    cout << "nombreArc: " << nombreArc << "\n";
    //(*archivoRec) >> dato >> nPols;
    
    std::ifstream archivoRecGC;
    archivoRecGC.open (nombreArc);
    recGC = crearGC (&archivoRecGC);
    archivoRecGC.close();
    //exit (0);

    func = new FuncionReconstructor (recGC);
    inc = new IncrementadorMallaAristas (-1);
    ini = new InicializadorVoronoiUniforme ();
    
    ret = new ReconstructorIncremental (func, ini, inc, nPolsFin);
    func->f(ret->getPars());
    
  }
  else {
    cerr << "ERROR: tipo de funcion inexistente: " << nombreFunc << "\n";
    cerr << "       use: AJUSTE o GCARISTAS\n";
    exit (1);
  }
  return ret;
}

Reconstructor *crearGenetico (std::ifstream *archivoRec) {
  char nombreFits [255], nombreFunc [255], dato [255], **nombresVis, nombreGeneracion[255];
  int nPols, nVis;
  double *nMax, *nMin, NiMax, individuos, generaciones, expanded_lx, expanded_ly;

  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPols;
  (*archivoRec) >> dato >> NiMax;  
  (*archivoRec) >> dato >> individuos;
  (*archivoRec) >> dato >> generaciones;
  (*archivoRec) >> dato >> expanded_lx;
  (*archivoRec) >> dato >> expanded_ly;
  (*archivoRec) >> dato >> nombreFunc;
  (*archivoRec) >> dato >> nVis;
  cout << "nVis: " << nVis << "\n";
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    (*archivoRec) >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  (*archivoRec) >> dato >> nombreGeneracion;
  
  funcBVCBI = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
					  expanded_lx, expanded_ly);
  if (strcmp(nombreFunc, "DIV") == 0) {
        func = new FuncionDividida (funcBVCBI);
  }
  else if (strcmp(nombreFunc, "NEG") == 0) {
        func = new FuncionNegativa (funcBVCBI);
  }
  else {
    cerr << "ERROR: tipo de funcion inexistente: " << nombreFunc << "\n";
    cerr << "       use: DIV o NEG\n";
    exit (1);
  }

  nMax = new double [func->getNPars()];
  nMin = new double [func->getNPars()];

  for (int i = 0; i < func->getNPars() / 3; i++) {
    nMax [3 * i]     = 1.0;
    nMax [3 * i + 1] = 1.0;
    nMax [3 * i + 2] = NiMax;
    nMin [3 * i]     = 0;
    nMin [3 * i + 1] = 0;
    nMin [3 * i + 2] = 0.0;
  }

  ini = new InicializadorVoronoiUniforme ();
    
  ReconstructorGenetico *r = new ReconstructorGenetico (func, ini, nMin, nMax, 
							individuos, generaciones);
  if (strcmp(nombreGeneracion, "NO") != 0) {
    r->leerGeneracion(nombreGeneracion);
  }
  return r;
}

ReconstructorGC *crearGC (std::ifstream *archivoRec) {
  char nombreFits [256], nombreFunc [256], dato [256], **nombresVis, archivoIn [256], 
    tipoFuncion [256];
  int nPols, nVis;
  double ftol, expanded_lx, expanded_ly;
  ReconstructorGC *ret;

  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPols;
  (*archivoRec) >> dato >> ftol;
  (*archivoRec) >> dato >> expanded_lx;
  (*archivoRec) >> dato >> expanded_ly;
  (*archivoRec) >> dato >> nVis;
  cout << "nVis: " << nVis << "\n";
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    (*archivoRec) >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  (*archivoRec) >> dato >> archivoIn;
  (*archivoRec) >> dato >> tipoFuncion;

  cout << "nombreFits: " << nombreFits << "\n";
  cout << "nPols: " << nPols << "\n";
  cout << "ftol: " << ftol << "\n";
  cout << "exlx: " << expanded_lx << "\n";
  cout << "exly: " << expanded_ly << "\n";
  cout << "NVis: " << nVis << "\n";
  cout << "archivoIn: " << archivoIn << "\n";
  cout << "TipoFuncion: " << tipoFuncion << "\n";

  if (strcmp(archivoIn, "NO") == 0) {
     ini = new InicializadorVoronoiUniforme ();
  }
  else {
    ini = new InicializadorArchivo (archivoIn);
  }
 
  if (strcmp(tipoFuncion, "FUNCIONBAYESVORONOI") == 0) {
    func = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
				       expanded_lx, expanded_ly);
  }
  else if (strcmp(tipoFuncion, "FUNCIONBAYESVORONOII") == 0) {
    funcBVCBI = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
					    expanded_lx, expanded_ly);
    func = new FuncionVoronoiI (funcBVCBI, ini);

    delete ini;
    (*archivoRec) >> dato >> archivoIn;
    cout << "archivoIn: " << archivoIn << "\n";
    if (strcmp(archivoIn, "NO") == 0) {
      ini = new InicializadorVoronoiUniforme ();
    }
    else {
      ini = new InicializadorArchivo (archivoIn);
    }
  }
  else {
    cerr << "ERROR: tipo de funcion inexistente\n";
    cerr << "       use: FUNCIONBAYESVORONOI o FUNCIONBAYESVORONOII\n";
    exit (1);

  }

  ret = new ReconstructorGC (func, ini, ftol);
  /*
  FuncionBayesVoronoiCBI *f_aux = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 
							      1, -1, expanded_lx, expanded_ly);
  for (int i = 0; i < ret->getNPars()/3; i++) {
    ret->setPar(i * 3 + 2, ret->getPar(i * 3 + 2) / f_aux->getRuido() / f_aux->getRuido());
  }
  f_aux->f(ret->getPars());
  f_aux->guardarInfo ("_in");
  delete f_aux;
  exit(0);
  */
  return ret;
}

Reconstructor *crearGeneticoGC (std::ifstream *archivoRec) {
  char nombreFits [255], nombreFunc [255], dato [255], **nombresVis, nombreGeneracion[255];
  int nPols, nVis;
  double *nMax, *nMin, NiMax, individuos, generaciones, expanded_lx, expanded_ly, ftol;

  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPols;
  (*archivoRec) >> dato >> ftol;
  (*archivoRec) >> dato >> NiMax;  
  (*archivoRec) >> dato >> individuos;
  (*archivoRec) >> dato >> generaciones;
  (*archivoRec) >> dato >> expanded_lx;
  (*archivoRec) >> dato >> expanded_ly;
  (*archivoRec) >> dato >> nombreFunc;
  (*archivoRec) >> dato >> nVis;
  cout << "nVis: " << nVis << "\n";
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    (*archivoRec) >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  (*archivoRec) >> dato >> nombreGeneracion;
  
  funcBVCBI = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
					  expanded_lx, expanded_ly);
  recGC = new ReconstructorGC(funcBVCBI, new InicializadorVoronoiUniforme (), ftol);
  funcR = new FuncionReconstructor (recGC);

  if (strcmp(nombreFunc, "DIV") == 0) {
        func = new FuncionDividida (funcR);
  }
  else if (strcmp(nombreFunc, "NEG") == 0) {
        func = new FuncionNegativa (funcR);
  }
  else {
    cerr << "ERROR: tipo de funcion inexistente: " << nombreFunc << "\n";
    cerr << "       use: DIV o NEG\n";
    exit (1);
  }

  nMax = new double [func->getNPars()];
  nMin = new double [func->getNPars()];

  for (int i = 0; i < func->getNPars() / 3; i++) {
    nMax [3 * i]     = 1.0;
    nMax [3 * i + 1] = 1.0;
    nMax [3 * i + 2] = NiMax;
    nMin [3 * i]     = 0;
    nMin [3 * i + 1] = 0;
    nMin [3 * i + 2] = 0.0;
  }

  ini = new InicializadorVoronoiUniforme ();
    
  ReconstructorGenetico *r = new ReconstructorGenetico (func, ini, nMin, nMax, individuos, generaciones);
  if (strcmp(nombreGeneracion, "NO") != 0) {
    r->leerGeneracion(nombreGeneracion);
  }
  return r;
}

Reconstructor *crearBayes (std::ifstream *archivoRec){
  char nombreFits [256], nombreFunc [256], dato [256], **nombresVis, archivoIn [256];
  int nPols, nVis, n_iter;
  double ftol, expanded_lx, expanded_ly, sigma;
  Reconstructor *ret;

  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPols;
  (*archivoRec) >> dato >> ftol;
  (*archivoRec) >> dato >> expanded_lx;
  (*archivoRec) >> dato >> expanded_ly;
  (*archivoRec) >> dato >> nVis;
  cout << "nVis: " << nVis << "\n";
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    (*archivoRec) >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  (*archivoRec) >> dato >> sigma;
  (*archivoRec) >> dato >> n_iter;
  cout << "n_iter = " << n_iter << "\n";
  (*archivoRec) >> dato >> archivoIn;
  
  func = new FuncionPbbVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
				   expanded_lx, expanded_ly);

  if (strcmp(archivoIn, "NO") == 0) {
     ini = new InicializadorVoronoiUniforme ();
  }
  else {
    ini = new InicializadorArchivo (archivoIn);
  }
 
  //ini = new InicializadorVoronoiUniforme ();
  distNorm = new DistribucionNormal (func->getNPars(), sigma);
  distMet = new DistribucionMetropolis (func, distNorm);
  ret =  new ReconstructorBayesiano (distMet, ini, n_iter);
  
  /*
  FuncionBayesVoronoiCBI *f_aux = new FuncionBayesVoronoiCBI (nombreFits, nombresVis, nVis, nPols, 
							  1, -1, expanded_lx, expanded_ly);
  for (int i = 0; i < ret->getNPars()/3; i++) {
    ret->setPar(i * 3 + 2, ret->getPar(i * 3 + 2) / f_aux->getRuido() / f_aux->getRuido());
  }
  f_aux->f(ret->getPars());
  f_aux->guardarInfo ("_in");
  delete f_aux;
  exit(0);
  */
  return ret;
}

Reconstructor *crearGibbs (std::ifstream *archivoRec) {
  char nombreFits [256], nombreFunc [256], dato [256], **nombresVis, archivoIn [256];
  int nPols, nVis, n_iter, N;
  double ftol, expanded_lx, expanded_ly;
  Reconstructor *ret;

  (*archivoRec) >> dato >> nombreFits;
  (*archivoRec) >> dato >> nPols;
  (*archivoRec) >> dato >> ftol;
  (*archivoRec) >> dato >> expanded_lx;
  (*archivoRec) >> dato >> expanded_ly;
  (*archivoRec) >> dato >> nVis;
  cout << "nVis: " << nVis << "\n";
  nombresVis = new char *[nVis];
  for (int i = 0; i < nVis; i++) {
    nombresVis[i] = new char[512];
    (*archivoRec) >> nombresVis[i];
    cout << nombresVis[i] << "\n";
  }
  (*archivoRec) >> dato >> n_iter;
  (*archivoRec) >> dato >> N;
  cout << "n_iter = " << n_iter << "\n";
  (*archivoRec) >> dato >> archivoIn;
  
  distCond = new DistribucionCondicionalBVCBI (nombreFits, nombresVis, nVis, nPols, 1, -1,
					       expanded_lx, expanded_ly, N);
  distGibbs = new DistribucionGibbs (distCond);

  if (strcmp(archivoIn, "NO") == 0) {
    ini = new InicializadorVoronoiUniforme ();
  }
  else {
    ini = new InicializadorArchivo (archivoIn);
  }
  
  ret =  new ReconstructorBayesiano (distGibbs, ini, n_iter);

  // for (int i = 0; i < ret->getNPars()/3; i++) {
  //   ret->setPar(i * 3 + 2, ret->getPar(i * 3 + 2) / distCond->getQuanta());
  // }

  return ret;
}

Reconstructor *crearMEM (std::ifstream *archivoRec){
  double ftol, xini;
  int Nmax;
  Reconstructor *r;
  char nombreFunc [255], dato [255];

  (*archivoRec) >> dato >> ftol;
  (*archivoRec) >> dato >> Nmax;
  (*archivoRec) >> dato >> xini;
  (*archivoRec) >> dato >> nombreFunc;


  if (strcmp(nombreFunc, "MEMCBI") == 0) {
    tipoFunc = MEMCBI;
    cout << "MEMCBI\n";
    
    char *nombreFits, tipoEntropia[255], **nombresVis;
    double pMin, lambda, noise_cut, beamSize;
    int nVis;

    nombreFits = new char [255];

    (*archivoRec) >> dato >> nombreFits;
    struct image *im = do_read (nombreFits);
    int npixels = im->npixels;
    delete_map (im);
    (*archivoRec) >> dato >> pMin;
    (*archivoRec) >> dato >> lambda;
    (*archivoRec) >> dato >> noise_cut;
    (*archivoRec) >> dato >> beamSize;
    (*archivoRec) >> dato >> nVis;

    nombresVis = new char *[nVis];
    for (int i = 0; i < nVis; i++) {
      nombresVis[i] = new char[512];
      (*archivoRec) >> nombresVis[i];
      cout << nombresVis[i] << "\n";
    }
    (*archivoRec) >> dato >> tipoEntropia;

    if (strcmp (tipoEntropia, "ILOGI") == 0) {
      tipoEnt = ILOGI;
      double lambdaE;
      (*archivoRec) >> dato >> lambdaE;

      fEnt = new FuncionEntropiaIlogI(npixels, pMin, lambdaE * 9.58957e-5);
      funcMemCBI = new FuncionMemCBI (nombreFits, nombresVis, nVis, fEnt, pMin, 
				      lambda, noise_cut, beamSize);
      delete fEnt;
      fEnt = new FuncionEntropiaIlogI(npixels, pMin, lambdaE * funcMemCBI->getRuido());
      delete funcMemCBI;
    }
    else if (strcmp (tipoEntropia, "NATURAL") == 0) {
      tipoEnt = SNATURAL;

      fEnt = new FuncionEntropiaNatural(npixels);
    }
    else {
      cerr << "ERROR: tipo de entropia inexistente en " << nombreFunc << " : " 
	   << tipoEntropia << "\n";
      cerr << "       use: ILOGI\n";
      exit (1);
    }

    funcMemCBI = new FuncionMemCBI (nombreFits, nombresVis, nVis, fEnt, pMin, lambda, noise_cut, beamSize);
    
    //funcMemCBI = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionEntropiaIlogI(64*64, pMin, 200 * 9.58957e-5), pMin, 0.0);
  }
  else {
    cerr << "ERROR: tipo de funcion inexistente para reconstructor MEM: " << nombreFunc << "\n";
    cerr << "       use: MEMCBI\n";
    exit (1);
  }

  ini = new InicializadorConstante();
  //r1 = new ReconstructorGCMemCBI (funcMemCBI, ini, 1e-10, 1, pMin); /* Para usar Nmax */
  r = new ReconstructorGCMemCBI (funcMemCBI, ini, ftol, Nmax, xini); /* Para usar Nmax */

  return r;
}

int main_old (int argc, char **argv) {

  FuncionMemCBI *funcMemCBI;
  double ftol = 1e-10;
  Reconstructor *r1;
  char ** vis = new char *[1];
  reconst tipoRec = MEM;
  //reconst rec = NORMGCMEM;
  //reconst rec = GIBBSMEM;
  //reconst rec = GIBBS;
  //reconst rec = GC;
  //reconst rec = GENETICO;
  //reconst rec = GENETICOGC;
  //reconst rec = BAYES;
  //reconst rec = INCREMENTALPIX;
  //reconst rec = INCREMENTALARISTAS;
  //reconst rec = AJUSTE;
  //reconst rec = INDIVIDUAL;
  Funcion *func;
  FuncionCBI *funcCBI;
  Funcion *funcR;
  Inicializador *ini;
  //FuncionMemCBI *funcMemCBI;
  double *nMax, *nMin, pmin, pmax;
  bool igual;

  DistribucionCondicionalBVCBI *dc;
  DistribucionGibbs *dg;
  DistribucionCondicionalMemCBI *dcMem;
  DistribucionGibbs *dgMem;
  double *pX;
  double *pI;
  double pMin;
  mtrace ();

  switch (tipoRec) {
  case GC:
    funcCBI = new FuncionBayesVoronoiCBI (argv[1]);
    //r1 = new ReconstructorGC (func, new InicializadorVoronoiPixelizado (), 1e-10);
    //r1 = new ReconstructorGC (func, new InicializadorVoronoiPixelizadoAlternado (), 1e-10);
    r1 = new ReconstructorGC (funcCBI, new InicializadorVoronoiUniforme (), 1e-10);
    if (argc == 3) {
      r1->reinicializar (argv [2]);
      for (int i = 0; i < r1->getNPars() / 3; i++) {
	r1->setPar(3 * i + 2, r1->getPar(3 * i + 2) / funcCBI->getRuido());
      }
    }

    break; 
  case GIBBS:
    dc = new DistribucionCondicionalBVCBI (argv[1], 50);
    //DistribucionCondicionalVoronoiI *dc = new DistribucionCondicionalVoronoiI (new DistribucionCondicionalBVCBI (argv[1], 50), 
    //new InicializadorVoronoiPixelizado());
    dg = new DistribucionGibbs (dc);
    r1 = new ReconstructorBayesiano (dg, new InicializadorVoronoiUniforme (), 1000000);
    if (argc == 3) {
      r1->reinicializar (argv [2]);
      //r1->reinicializar (argv [2], 130);
      for (int i = 0; i < r1->getNPars() / 3; i++) {
	r1->setPar(3 * i + 2, r1->getPar(3 * i + 2) / dc->getQuanta());
      }
    }
    
    break;
  case GIBBSMEM:
    dcMem = new DistribucionCondicionalMemCBI (argv[1], argv + 2, 1/*argc - 2*/, 
									      new FuncionEntropiaNatural(64*64), 1.0, 11,
									      new InicializadorArchivo (argv [3]));
    dgMem = new DistribucionGibbs (dcMem);
    r1 = new ReconstructorBayesiano (dgMem, new InicializadorConstante (), 1000000);
    if (argc == 4) {
      cout << "Reinicilizando\n";
      r1->reinicializar (argv [3]);
      //r1->reinicializar (argv [2], 130);
      for (int i = 0; i < r1->getNPars(); i++) {
	r1->setPar(i, r1->getPar(i) / dcMem->getQuanta());
	//cout << "p[" << i << "] = " << r1->getPar(i) << "\n";
      }
    }
    
    break;
  case BAYES:
    func = new FuncionPbbVoronoiCBI (argv[1]);
    //func = new FuncionNegativa(new FuncionBayesVoronoiCBII (argv[1]));
    ini = new InicializadorVoronoiUniforme ();
    //func = new FuncionDividida(new FuncionBayesVoronoiCBII (argv[1], ini));
    //func = new FuncionDividida (new FuncionVoronoiI (new FuncionChi2VoronoiImagen (argv[1], 50), ini));
    delete ini;
    //func = new FuncionEjemploH (3);
    //func = new FuncionBayesVoronoiCBI (argv[1]);
    //r1 = new ReconstructorBayesiano (new DistribucionHamiltoniana 
    //			     (func,
    //			      new DistribucionNormal (func->getNPars(), 
    //						      1),
    //			      10, 0.01),
    //			     new InicializadorVoronoiUniforme (), 1000000);
    
    r1 = new ReconstructorBayesiano (new DistribucionMetropolis 
				     (func,
				      new DistribucionNormal (func->getNPars(), 
							      0.03)),
				     new InicializadorVoronoiUniforme (), 1000000);
    
    //new InicializadorVoronoiUniforme (), 1e15);
    break;
  case GENETICO:
    /*
    ini = new InicializadorVoronoiPixelizadoAlternado ();
    func = new FuncionDividida(new FuncionVoronoiI (new FuncionBayesVoronoiCBI (argv[1])
						    , ini));
    delete ini;
    */
    //func = new FuncionPbbVoronoiCBI (argv[1]);
    //func = new FuncionNegativa (new FuncionVoronoiI (new FuncionChi2VoronoiImagen (argv[1], 50), 
    //					     new InicializadorVoronoiUniforme ()));
    //func = new FuncionDividida (new FuncionChi2VoronoiImagen (argv[1], 300));
    func = new FuncionDividida (new FuncionBayesVoronoiCBI (argv[1]));
    //func = new FuncionExp (3);
    
    nMax = new double [func->getNPars()];
    nMin = new double [func->getNPars()];
    if (func->getNPars() % 3 != 0) {
      cerr << "ERROR en GENETICO de VIR++\n";
      exit (0);
    }
    for (int i = 0; i < func->getNPars() / 3; i++) {
      nMax [3 * i]     = 1.0;
      nMax [3 * i + 1] = 1.0;
      nMax [3 * i + 2] = 11;
      nMin [3 * i]     = 0;
      nMin [3 * i + 1] = 0;
      nMin [3 * i + 2] = 0.0;
    }
    r1 = new ReconstructorGenetico (func, new InicializadorVoronoiUniforme (),
				    nMin, nMax, 100, 1000000);
    break;
  case GENETICOGC:
    func = new FuncionBayesVoronoiCBI (argv[1]);
    funcR = new FuncionDividida (new FuncionReconstructor (new ReconstructorGC(func, new InicializadorVoronoiUniforme (), 1e-5)));
    nMax = new double [func->getNPars()];
    nMin = new double [func->getNPars()];
    if (func->getNPars() % 3 != 0) {
      cerr << "ERROR en GENETICO de VIR++\n";
      exit (0);
    }
    for (int i = 0; i < func->getNPars() / 3; i++) {
      nMax [3 * i]     = 1.0;
      nMax [3 * i + 1] = 1.0;
      nMax [3 * i + 2] = 90;
      nMin [3 * i]     = 0;
      nMin [3 * i + 1] = 0;
      nMin [3 * i + 2] = 0.0;
    }
    r1 = new ReconstructorGenetico (funcR, new InicializadorVoronoiUniforme (),
				    nMin, nMax, 100, 1000000);
    break;
  case INCREMENTALPIX:
    func = new FuncionBayesVoronoiCBI (argv[1]);
    funcR = new FuncionReconstructor (new ReconstructorGC (func, new InicializadorVoronoiUniforme (), 1e-5));
    r1 = new ReconstructorIncremental (funcR, new InicializadorVoronoiUniforme (),
				       new IncrementadorMallaPixel (10, 10, -1), 50);
    break; 
  case INCREMENTALARISTAS:
    func = new FuncionBayesVoronoiCBI (argv[1]);
    funcR = new FuncionReconstructor (new ReconstructorGC (func, new InicializadorVoronoiUniforme (), 1e-5));
    r1 = new ReconstructorIncremental (funcR, new InicializadorVoronoiUniforme (),
				       new IncrementadorMallaAristas (-1), 1000);
    break; 
  case AJUSTE:
     func = new FuncionChi2VoronoiImagen (argv[1], 0);
    r1 = new ReconstructorIncremental (func, new InicializadorVoronoiUniforme (),
				       new IncrementadorMallaImagen (argv[1], -1), 1000);
    break; 
  case INDIVIDUAL:
    ini = new InicializadorVoronoiUniforme ();
    //ini = new InicializadorVoronoiPixelizado();
    //func = new FuncionDividida(new FuncionVoronoiI (new FuncionBayesVoronoiCBI (argv[1]), ini));

    /*/moviendo centros e Intensidades
    func = new FuncionDividida(new FuncionBayesVoronoiCBI (argv[1]));
    double *max = new double [func->getNPars()];
    double *min = new double [func->getNPars()];
    double *N = new double [func->getNPars()];
    for (int i = 0; i < func->getNPars() / 3; i++) {
      min[3*i]     = 0;
      min[3*i + 1] = 0;
      min[3*i + 2] = 0;
      max[3*i]     = 1;
      max[3*i + 1] = 1;
      max[3*i + 2] = 11;
      N[3*i]     = 50;
      N[3*i + 1] = 50;
      N[3*i + 2] = 22;
    }
    r1 = new ReconstructorIndividual (func, min, max, N, ini);
    */

    //r1 = new ReconstructorIndividual (func, 0, 11, 55, new InicializadorConstante());
    func = new FuncionDividida(new FuncionVoronoiI (new FuncionBayesVoronoiCBI (argv[1]), ini));
    pX = new double [func->getNPars() * 2];
    pI = new double [func->getNPars()];
    for (int i = 0; i < func->getNPars() * 2; i++) {
      pX [i] = rand()/(RAND_MAX + 1.0);
    }
    igual = false;
    delete func;
    for (int i = 0; !igual; i++) {
      func = new FuncionDividida(new FuncionVoronoiI (new FuncionBayesVoronoiCBI (argv[1]), pX));
      r1 = new ReconstructorIndividual (func, 0, 11, 22, new InicializadorConstante());
      for (int j = 0; j < func->getNPars(); j++) {
	cout << "r1->getPar(" << j << ") = " << r1->getPar(j) << "\n";
      }
      r1->run();

      igual = true;
      for (int j = 0; j < func->getNPars() && igual; j++){
	if (pI[j] != r1->getPar(j)) {
	  igual = false;
	}
      }
      if (igual) {
	exit(0);
      }
      for (int j = 0; j < func->getNPars(); j++) {
	pI[j] = r1->getPar(j);
      }
      func->f(pI);
      func->guardarInfo("FinalI" + i2s(i, 3));
      delete r1;
      delete func;
      continue;

      func = new FuncionDividida(new FuncionVoronoiX (new FuncionBayesVoronoiCBI (argv[1]), pI));
      r1 = new ReconstructorIndividual (func, 0, 1, 50, new InicializadorConstante());
      for (int j = 0; j < func->getNPars(); j++) {
	r1->setPar(j, pX[j]);
      }
      r1->run();

      igual = true;
      for (int j = 0; j < func->getNPars() && igual; j++){
	if (pX[j] != r1->getPar(j)) {
	  igual = false;
	}
      }
      if (igual) {
	exit(0);
      }
      for (int j = 0; j < func->getNPars(); j++){
	pX[j] = r1->getPar(j);
      }
      
      func->f(pX);
      func->guardarInfo("FinalX" + i2s(i, 3));
      delete r1;
      delete func;
    }
    break;
  case MEM:
    pMin = 1.71543e-8/9.58957e-5;
    //pMin = 9.65927e-9/9.58957e-5;
    cout << "pmin = " << pMin << "\n";

    //funcMemCBI = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionEntropiaNatural(64*64), 0.0, 0.0);
    funcMemCBI = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionEntropiaIlogI(64*64, pMin, 200 * 9.58957e-5), pMin, 0.0, 0.0, -1);

    //funcMemCBI->resizeImage (128, 128, 32, 32);
    //funcMemCBI = new FuncionMemCBI (argv[1], argv + 2, argc - 2, new FuncionCero(64*64));
    /*
    double *p = new double [funcCBI->getNPars()];
    double *g = new double [funcCBI->getNPars()];
    for (int j = 0; j < funcCBI->getNPars(); j++) {
      p[j] = 0;
    }
    for (int j = 0; j < 1e2; j++) {
      if (j % 10 == 0) cout << j << "\n";
      funcCBI->f(p);
    }
    for (int j = 0; j < 1e2; j++) {
      if (j % 10 == 0) cout << j << "\n";
      funcCBI->df(p, g);
      }*/
    //r1 = new ReconstructorGC (func, new InicializadorVoronoiPixelizado (), 1e-10);
    //r1 = new ReconstructorGC (func, new InicializadorVoronoiPixelizadoAlternado (), 1e-10);
    
    //r1 = new ReconstructorGC (funcMemCBI, new InicializadorConstante (), 1e-10);
    //r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, Nmax, xini); /* Para usar Nmax */
    //r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, 1); /* Para usar Nmax */
    //r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, 0.0008628 / funcMemCBI->getRuido()); /* Para usar Nmax */
    r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, pMin); /* Para usar Nmax */
    //r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, 52.143); /* Para usar Nmax */

    //std::ostream arcOut ("");
    //for (double xini = 0; xini < 70 ; xini += 3) {
    //r1 = new ReconstructorGCMemCBI (funcMemCBI, new InicializadorConstante (), 1e-10, 1, xini); /* Para usar Nmax */
    //arcOut << xini << "\t" << r1->run() << "\n";
    //}
    /*
    r1->run();
    funcMemCBI->guardarInfo ("Final1");
    pMin = 1e300;
    double *pr1 = r1->getPars();
    for (int i = 0; i < r1->getNPars(); i++) {
      if (pr1[i] < pMin) {
	pMin = pr1[i];
      }
    }
    
    for (int i = 0; i < r1->getNPars(); i++) {
      pr1[i] = pr1[i] - pMin;
    }

    funcMemCBI->scalePars(r1->getPars());
    funcMemCBI->guardarInfo ("Escalado");
    //exit(0);
    */
    break;
  case NORMGCMEM:
    funcMemCBI = new FuncionMemCBI (argv[1], argv + 2, 1, new FuncionEntropiaNatural(64*64), 0.0, 1.0, 0.0, -1);

    r1 = new ReconstructorGC (funcMemCBI, new InicializadorArchivo (argv[3]), 1e-10);

    //Normalizamos las intensidades para que exista Imax
    pmax = 0;
    pmin = 1e300;
    for (int i = 0; i < r1->getNPars() / 3; i++) {
      if (r1->getPar(3 * i + 2) < pmin) {
	pmin = r1->getPar(3 * i + 2);
      }
      if (r1->getPar(3 * i + 2) > pmax) {
	pmax = r1->getPar(3 * i + 2);
      }
    }

    double p2;//, *p = new double [r1->getNPars() / 3];
    for (int i = 0; i < r1->getNPars(); i++) {
      p2 = (r1->getPar(i) - pmin) * (10 - pmin) / (pmax - pmin) + pmin;
      r1->setPar(i, p2);
    }
    
    /*
    for (int j = 0; j < 1e2; j++) {
      cout << "---------------------->" << j << "\n";
      delete r1;
      Funcion * fExp = new FuncionExp (1);
      r1 = new ReconstructorGC (fExp, new InicializadorConstante (), 1e-10);
      r1->run();
      delete fExp;
      }*/
    /*
    if (argc == 3) {
      r1->reinicializar (argv [2]);
      for (int i = 0; i < r1->getNPars() / 3; i++) {
	r1->setPar(3 * i + 2, r1->getPar(3 * i + 2) / funcCBI->getRuido());
      }
    }
    */
    break; 
  default:
    exit(0);
  }

  cout << "Reconstructor creado \n";
  r1->run ();

  //  cout << "Valor Maximo en " << r1->getPar(0) << "\n";
  delete r1;
}

