#include "reconstructorVIRIterativo.h"

ReconstructorVIRIterativo::ReconstructorVIRIterativo (char *nombreImagen, char **nombresVis, 
						      int nVis, int nPolsIni, int nPolsFin, 
						      int dn, double expanded_lx, 
						      double expanded_ly, int MEMIter,
						      double beamSize, int nIterMax) :
ReconstructorVIR (nombreImagen, nombresVis, nVis, nPolsIni, nPolsFin, dn, expanded_lx, 
		  expanded_ly, MEMIter, beamSize){
  this->nIterMax = nIterMax;
}

double ReconstructorVIRIterativo::realizarReconstruccion (FuncionVoronoi *funcCBI, double *p) {
  double ret1, ret2;
  bool salir = false;
  int iter = 0;
  FuncionVoronoiI * funcI;
  FuncionVoronoiX * funcX;
  ReconstructorGC *recGC;
  Inicializador *ini;
  Inicializador *iniI, *iniX;
  double *pI = new double [funcCBI->getNPars() / 3], 
    *pX = new double [2 * funcCBI->getNPars() / 3];

  for (int i = 0; i < funcCBI->getNPars() / 3; i++) {
    pX [2 * i]     = p [3 * i];
    pX [2 * i + 1] = p [3 * i + 1];
    pI [i]         = p [3 * i + 2];
  }

  while (!salir) {
    ini =  new InicializadorParametros (p, funcCBI->getNPars());

    funcX = new FuncionVoronoiX (funcCBI, ini);
    iniX = new InicializadorParametros (pX, 2 * funcCBI->getNPars() / 3);
    recGC = new ReconstructorGC (funcX, iniX, 1e-10);
    recGC->run();
    for (int i = 0; i < funcCBI->getNPars() / 3; i++) {
      p[3 * i]     = (pX [2 * i]     = recGC->getPar(2 * i));
      p[3 * i + 1] = (pX [2 * i + 1] = recGC->getPar(2 * i + 1));
    }
    funcX->f(pX);
    funcX->guardarInfo("VIRIncX_" + i2s(funcCBI->getNPars() / 3, 3) + "_" + i2s(iter, 3));
    delete funcX;
    delete recGC;
    delete ini;

    ini =  new InicializadorParametros (p, funcCBI->getNPars());
    funcI = new FuncionVoronoiI (funcCBI, ini);
    iniI = new InicializadorParametros (pI, funcCBI->getNPars() / 3);
    recGC = new ReconstructorGC (funcI, iniI, 1e-10);
    salir = (ret1 == (ret2 = recGC->run()) || nIterMax == iter++);
    ret1 = ret2;
    for (int i = 0; i < funcCBI->getNPars() / 3; i++) {
      p[3 * i + 2] = (pI [i] = recGC->getPar(i));
    }
    funcI->f(pI);
    funcI->guardarInfo("VIRIncI_" + i2s(funcCBI->getNPars() / 3, 3) + "_" + i2s(iter, 3));
    delete funcI;
    delete recGC;
    delete ini;
    iter ++;
  }
  
  delete [] pI;
  delete [] pX;
  return ret1;
}

double ReconstructorVIRIterativo::run () {
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
  FuncionBayesVoronoiCBI *funcCBI = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, 
						    nVis, nPolsIni, 1, -1,
						    expanded_lx, expanded_ly, beamSize);
  //FuncionCBI *funcCBI = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, 
  //					    nVis, nPolsIni, 0, -1);
  //ReconstructorGC recGC (funcCBI, iniVU, 1e-10);

  double ret;
  //double *pChi2 = new double [func->getNPars()];

  double *pChi2 = new double [func->getNPars()];
  for (int i = 0; i < func->getNPars()/3; i++) {
      pChi2[3 * i]     = rInc.getPar(3 * i);
      pChi2[3 * i + 1] = rInc.getPar(3 * i + 1);
      pChi2[3 * i + 2] = rInc.getPar(3 * i + 2);
    }

  for (int nPols = nPolsIni; nPols <= nPolsFin; nPols += dn) {
    delete [] p;
    p = new double [func->getNPars()];
    for (int i = 0; i < func->getNPars()/3; i++) {
      p[3 * i]     = pChi2 [3 * i];
      p[3 * i + 1] = pChi2 [3 * i + 1];
      p[3 * i + 2] = pChi2 [3 * i + 2] / funcCBI->getRuido();
    }
    /*
    for (int j = 0; j < recGC.getNPars() / 3; j++) {
      //cout << "cambiando " << j << "): " << recGC.getPar(j) << " por " << p[j] << "\n";
      recGC.setPar(3 * j, p[3 * j]);
      recGC.setPar(3 * j + 1, p[3 * j + 1]);
      recGC.setPar(3 * j + 2, p[3 * j + 2] / funcCBI->getRuido());
      }*/
    
    ret = realizarReconstruccion (funcCBI, p);
    
    archivoOut.open("reconstructorVIR.out", ios::app);
    archivoOut << nPols << "\t" << ret << "\n";
    archivoOut.close();

    //sprintf (nombreArchivo, "VIR_NPols_%g_iter", qFactor);
    //funcCBI->guardarInfo ("VIR_NPols_" + i2s (nPols, 4));
    //funcCBI->guardarResiduos ("VIR_NPols_" + i2s (nPols, 4));

    /*for (int i = 0; i < func->getNPars()/3; i++) {
      p[3 * i + 2] *= funcCBI->getRuido();
      }*/
    for (int j = 0; j < dn; j++) {
      func->guardarInfo (i2s (nPols + j, 4));
      inc->incrementar (&(pChi2), func);
    }
    /*for (int i = 0; i < func->getNPars()/3; i++) {
      p[3 * i + 2] /= funcCBI->getRuido();
      }*/

    funcCBI->setNPars(func->getNPars());
    //recGC.setNPars (func->getNPars());
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

  delete [] pChi2;
  return ret;
}

