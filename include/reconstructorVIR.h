#ifndef VIRPP_RECONSTRUCTORVIR
#define VIRPP_RECONSTRUCTORVIR

#include "funcion.h"			 
#include "funcionCBI.h"			 
#include "funcionEntropiaNatural.h"	 
#include "funcionChi2VoronoiImagen.h"	 
#include "funcionBayesVoronoiCBI.h"	 
#include "funcionMemCBI.h"
#include "funcionCero.h"
	  					 
#include "reconstructor.h"
#include "reconstructorGC.h"
#include "reconstructorGCMemCBI.h"
#include "reconstructorIncremental.h"
	  					 
#include "inicializador.h"		 
#include "inicializadorConstante.h"	 
#include "inicializadorVoronoiUniforme.h"

#include "incrementadorMallaImagen.h"

class ReconstructorVIR : public Reconstructor {
 protected:
  char *nombreImagen, **nombresVis;
  int nVis, MEMIter, nPolsIni, nPolsFin, dn;
  double expanded_lx, expanded_ly, beamSize;

  void realizarMEM ();

 public:
  ReconstructorVIR (char *nombreImagen, char **nombresVis, int nVis, 
		    int nPolsIni, int nPolsFin, int dn, 
		    double expanded_lx, double expanded_ly, int MEMIter = 1, double beamSize = -1);
  char *getNombreImagen() {return nombreImagen;}
  char **getNombresVis() {return nombresVis;}
  double getExpandedLx () {return expanded_lx;}
  double getExpandedLy () {return expanded_ly;}
  int getNPolsIni() {return nPolsIni;}
  double run ();
  
  ~ReconstructorVIR();
};

#endif
