# This file is designed to be invoked from a recursive call to make.
#
# The makefile that includes it is expected to define the following
# macros:
#
#  BIN       - The full name of the bin directory.
#  INC       - The full name of the include directory.
#  LIB       - The full name of the lib directory.
#  CC        - The C compiler to be used.
#  CFLAGS    - Flags for the C compiler.
#  RANLIB    - The archiver randomizer command (use : if not needed).
#  PGDIR     - The directory that contains pgplot.h and libpgplot.
#  FORT      - The Fortran compiler.
#  FFLAGS    - Fortran compiler options
#  XLIB      - linking options for X11
#  THLIB     - linking options for threads
#  CFLIB     - linking options for CFITSIO
#-----------------------------------------------------------------------

SHELL=/bin/csh
CPGPLOT=-L$(PGDIR) -lcpgplot -lpgplot $(XLIB)
PGPLOT=-L$(PGDIR) -lpgpl$(XLIB)

#-----------------------------------------------------------------------
# List the default targets.
#-----------------------------------------------------------------------

all: \
     $(BIN)\
     $(BIN)/VIR++\
     $(BIN)/crear_imagen\
     $(BIN)/crear_intensidades\
     $(BIN)/comprobacion_gradiente\
     $(BIN)/comprobacion_gradiente_dVdI\
     $(BIN)/comparar_visibilidades\
     $(BIN)/do_mockcbi\
     $(BIN)/quanta\
     $(BIN)/numeroCentros\
     $(BIN)/comprobacionDistribucion\
     $(BIN)/calcular_funcion\
     #$(BIN)/busquedaNPols\
     $(BIN)/residuos\

#-----------------------------------------------------------------------
# Make the bin directory if necessary.
#-----------------------------------------------------------------------

$(BIN):
	mkdir $(BIN)

#-----------------------------------------------------------------------
# Source file dependencies.
#-----------------------------------------------------------------------

#command_interp.o: command_interp.c $(INC)/libtecla.h \
# $(INC)/history.h $(INC)/command_interp.h

#help.o: help.c \
# $(INC)/help.h

gauss.o: gauss.c\
 $(INC)/gauss.h

complex.o: complex.c\
 $(INC)/complex.h

fourn.o: fourn.c\
 $(INC)/fourn.h

frprmn.o: frprmn.c\
 $(INC)/frprmn.h

brent.o: brent.c\
 $(INC)/brent.h

linmin.o: linmin.c\
 $(INC)/linmin.h

f1dim.o: f1dim.c\
 $(INC)/f1dim.h

mnbrak.o: mnbrak.c\
 $(INC)/mnbrak.h

nodoLista.o: nodoLista.c\
 $(INC)/nodoLista.h

lista.o: lista.c\
 $(INC)/lista.h\
 $(INC)/nodoLista.h

punto.o: punto.c\
 $(INC)/punto.h

arista.o: arista.c\
 $(INC)/arista.h\
 $(INC)/punto.h

poligono.o: poligono.c\
 $(INC)/poligono.h\
 $(INC)/arista.h

mallaVoronoi.o: mallaVoronoi.c\
 $(INC)/mallaVoronoi.h\
 $(INC)/lista.h\
 $(INC)/arista.h

rasterizador.o: rasterizador.c\
 $(INC)/rasterizador.h\
 $(INC)/mallaVoronoi.h

image_routines.o: image_routines.c\
 $(INC)/fitsio.h\
 $(INC)/longnam.h\
 $(INC)/image_routines.h\
 $(INC)/mathconst.h\
 $(INC)/newstring.h\
 $(INC)/slalib.h

funcionCBI.o: funcionCBI.cpp\
 $(INC)/funcionCBI.h\
 $(INC)/funcion.h

funcionBayesVoronoiCBI.o: funcionBayesVoronoiCBI.cpp\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcionVoronoi.h\
# $(INC)/funcL.h\
# $(INC)/mallaVoronoi.h\
 $(INC)/funcion.h

funcionMemCBI.o: funcionMemCBI.cpp\
 $(INC)/funcionMemCBI.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcion.h

funcionEntropiaNatural.o: funcionEntropiaNatural.cpp\
 $(INC)/funcionEntropiaNatural.h\
 $(INC)/funcion.h

funcionEntropiaIlogI.o: funcionEntropiaIlogI.cpp\
 $(INC)/funcionEntropiaIlogI.h\
 $(INC)/funcion.h

funcionBayesVoronoiCBII.o: funcionBayesVoronoiCBII.cpp\
 $(INC)/funcionBayesVoronoiCBII.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcionVoronoi.h\
 $(INC)/funcion.h

funcionPbbVoronoiCBI.o: funcionPbbVoronoiCBI.cpp\
 $(INC)/funcionPbbVoronoiCBI.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcionVoronoi.h\
 $(INC)/funcion.h

funcionReconstructor.o: funcionReconstructor.cpp\
 $(INC)/funcionReconstructor.h\
 $(INC)/funcion.h\
 $(INC)/reconstructor.h

funcionChi2VoronoiImagen.o: funcionChi2VoronoiImagen.cpp\
 $(INC)/funcionChi2VoronoiImagen.h\
 $(INC)/funcionVoronoi.h\
 $(INC)/funcion.h
# $(INC)/image_routines.h

#funcionChi2VoronoiImagenN.o: $(INC)/funcionChi2VoronoiImagenN.h\
# $(INC)/funcionChi2VoronoiImagen.h\
# $(INC)/funcionVoronoi.h\
# $(INC)/funcion.h

funcionExp.o: funcionExp.cpp\
 $(INC)/funcionExp.h\
 $(INC)/funcion.h

funcionEjemploH.o: funcionEjemploH.cpp\
 $(INC)/funcionEjemploH.h\
 $(INC)/funcion.h

funcionNegativa.o: funcionNegativa.cpp\
 $(INC)/funcionNegativa.h\
 $(INC)/funcion.h

funcionDividida.o: funcionDividida.cpp\
 $(INC)/funcionDividida.h\
 $(INC)/funcion.h

funcionVoronoiI.o: funcionVoronoiI.cpp\
 $(INC)/funcionVoronoiI.h\
 $(INC)/funcion.h\
 $(INC)/funcionVoronoi.h\
 $(INC)/inicializador.h

funcionVoronoiX.o: funcionVoronoiX.cpp\
 $(INC)/funcionVoronoiX.h\
 $(INC)/funcion.h\
 $(INC)/funcionVoronoi.h\
 $(INC)/inicializador.h

inicializadorVoronoiUniforme.o: inicializadorVoronoiUniforme.cpp\
 $(INC)/inicializadorVoronoiUniforme.h\
 $(INC)/inicializador.h

inicializadorVoronoiPixelizado.o: inicializadorVoronoiPixelizado.cpp\
 $(INC)/inicializadorVoronoiPixelizado.h\
 $(INC)/inicializador.h

inicializadorVoronoiPixelizadoAlternado.o: inicializadorVoronoiPixelizadoAlternado.cpp\
 $(INC)/inicializadorVoronoiPixelizadoAlternado.h\
 $(INC)/inicializador.h

inicializadorArchivo.o: inicializadorArchivo.cpp\
 $(INC)/inicializadorArchivo.h\
 $(INC)/inicializador.h

inicializadorParametros.o: inicializadorParametros.cpp\
 $(INC)/inicializadorParametros.h\
 $(INC)/inicializador.h

reconstructor.o: reconstructor.cpp\
 $(INC)/reconstructor.h\
 $(INC)/funcion.h\
 $(INC)/frprmn.h\
 $(INC)/inicializador.h
# $(INC)/gauss.h

reconstructorGC.o: reconstructorGC.cpp\
 $(INC)/reconstructorGC.h\
 $(INC)/frprmn.h\
 $(INC)/reconstructor.h

reconstructorGCMemCBI.o: reconstructorGCMemCBI.cpp\
 $(INC)/reconstructorGCMemCBI.h\
 $(INC)/frprmn.h\
 $(INC)/reconstructor.h\
 $(INC)/funcionMemCBI.h

reconstructorBayesiano.o: reconstructorBayesiano.cpp\
 $(INC)/reconstructorBayesiano.h\
 $(INC)/reconstructor.h

reconstructorGenetico.o: reconstructorGenetico.cpp\
 $(INC)/reconstructorGenetico.h\
 $(INC)/reconstructor.h

reconstructorIncremental.o: reconstructorIncremental.cpp\
 $(INC)/reconstructorIncremental.h\
 $(INC)/reconstructor.h\
 $(INC)/incrementador.h

reconstructorIncrementalDoble.o: reconstructorIncrementalDoble.cpp\
 $(INC)/reconstructorIncrementalDoble.h\
 $(INC)/reconstructor.h\
 $(INC)/incrementador.h

reconstructorIndividual.o: reconstructorIndividual.cpp\
 $(INC)/reconstructorIndividual.h\
 $(INC)/reconstructor.h\
 $(INC)/incrementador.h\
 $(INC)/funcion.h

reconstructorVIR.o: reconstructorVIR.cpp\
 $(INC)/reconstructorVIR.h\
 $(INC)/funcion.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcionEntropiaNatural.h\
 $(INC)/funcionChi2VoronoiImagen.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionMemCBI.h\
 $(INC)/funcionCero.h\
 $(INC)/reconstructor.h\
 $(INC)/reconstructorGC.h\
 $(INC)/reconstructorGCMemCBI.h\
 $(INC)/reconstructorIncremental.h\
 $(INC)/inicializador.h\
 $(INC)/inicializadorConstante.h\
 $(INC)/inicializadorVoronoiUniforme.h\
 $(INC)/incrementadorMallaImagen.h

reconstructorVIRIterativo.o: reconstructorVIRIterativo.cpp\
 $(INC)/reconstructorVIRIterativo.h\
 $(INC)/reconstructorVIR.h\
 $(INC)/funcion.h\
 $(INC)/funcionCBI.h\
 $(INC)/funcionEntropiaNatural.h\
 $(INC)/funcionChi2VoronoiImagen.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionMemCBI.h\
 $(INC)/funcionCero.h\
 $(INC)/funcionVoronoiI.h\
 $(INC)/funcionVoronoiX.h\
 $(INC)/reconstructor.h\
 $(INC)/reconstructorGC.h\
 $(INC)/reconstructorGCMemCBI.h\
 $(INC)/reconstructorIncremental.h\
 $(INC)/inicializador.h\
 $(INC)/inicializadorParametros.h\
 $(INC)/inicializadorConstante.h\
 $(INC)/inicializadorVoronoiUniforme.h\
 $(INC)/incrementadorMallaImagen.h

incrementadorMallaPixel.o: incrementadorMallaPixel.cpp\
 $(INC)/incrementadorMallaPixel.h\
 $(INC)/incrementador.h

incrementadorMallaImagen.o: incrementadorMallaImagen.cpp\
 $(INC)/incrementadorMallaImagen.h\
 $(INC)/mallaVoronoi.h\
# $(INC)/image_routines.h\
 $(INC)/incrementador.h

incrementadorMallaAristas.o: incrementadorMallaAristas.cpp\
 $(INC)/incrementadorMallaAristas.h\
 $(INC)/incrementador.h\
 $(INC)/mallaVoronoi.h\
 $(INC)/arista.h

distribucionNormal.o: distribucionNormal.cpp\
 $(INC)/distribucionNormal.h\
 $(INC)/distribucionPbb.h

distribucionMetropolis.o: distribucionMetropolis.cpp\
 $(INC)/distribucionMetropolis.h\
 $(INC)/distribucionGenerador.h\
 $(INC)/distribucionPbb.h

distribucionHamiltoniana.o: distribucionHamiltoniana.cpp\
 $(INC)/distribucionHamiltoniana.h\
 $(INC)/distribucionPbb.h

VIR++.o: VIR++.cpp \
 $(INC)/mallaVoronoi.h\
# $(INC)/image_routines.h\
 $(INC)/funcion.h\
 $(INC)/fourn.h \
 $(INC)/frprmn.h \
 $(INC)/brent.h \
 $(INC)/linmin.h \
 $(INC)/f1dim.h \
 $(INC)/mnbrak.h \
 $(INC)/nrutil.h \
 $(INC)/gauss.h \
 $(INC)/reconstructor.h\
 $(INC)/reconstructorGC.h\
 $(INC)/reconstructorGCMemCBI.h\
 $(INC)/reconstructorGenetico.h\
 $(INC)/reconstructorBayesiano.h\
 $(INC)/reconstructorIncremental.h\
 $(INC)/inicializador.h\
 $(INC)/inicializadorConstante.h\
 $(INC)/distribucionNormal.h\
 $(INC)/distribucionMetropolis.h\
 $(INC)/funcionExp.h\
 $(INC)/funcionEjemploH.h\
 $(INC)/incrementadorMallaPixel.h\
 $(INC)/incrementadorMallaAristas.h\
 $(INC)/incrementadorMallaImagen.h\
 $(INC)/funcionChi2VoronoiImagen.h\
 $(INC)/funcionChi2VoronoiImagenN.h\
 $(INC)/funcionEntropiaNatural.h\
 $(INC)/funcionEntropiaIlogI.h

quanta.o: quanta.cpp \
 $(INC)/mallaVoronoi.h\
# $(INC)/image_routines.h\
 $(INC)/funcion.h\
 $(INC)/fourn.h \
 $(INC)/nrutil.h \
 $(INC)/inicializador.h\
 $(INC)/inicializadorArchivo.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/mockcbiRoutines.h\
 $(INC)/reconstructorGC.h

numeroCentros.o: numeroCentros.cpp\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/reconstructorGC.h\
 $(INC)/inicializadorVoronoiUniforme.h

calcular_funcion.o: calcular_funcion.cpp\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/mallaVoronoi.h

comprobacionDistribucion.o: comprobacionDistribucion.cpp\
 $(INC)/distribucionCondicional.h\
 $(INC)/distribucionCondicionalBVCBI.h\
 $(INC)/distribucionCondGaussiana.h\
 $(INC)/distribucionCondicionalMemCBI.h\
 $(INC)/funcionEntropiaNatural.h

nrutil.o: nrutil.c\
 $(INC)/nrutil.h

newstring.o: newstring.c\
 $(INC)/newstring.h

parse_sexagesimal.o: parse_sexagesimal.c\
 $(INC)/parse_sexagesimal.h

primary_beam.o: primary_beam.c \
 $(INC)/cbi_beam.h \
 $(INC)/mathconst.h \
 $(INC)/primary_beam.h

ra_dec_string.o: ra_dec_string.c\
 $(INC)/ra_dec_string.h\
 $(INC)/slalib.h

random.o: random.c \
 $(INC)/random.h

uvsubs.o: uvsubs.c\
 $(INC)/fitsio.h\
 $(INC)/longnam.h\
 $(INC)/slalib.h\
 $(INC)/uvsubs.h\
 $(INC)/mathconst.h

mockcbiRoutines.o: mockcbiRoutines.c\
 $(INC)/mockcbiRoutines.h\
 $(INC)/complex.h\
 $(INC)/fitsio.h\
 $(INC)/fourn.h\
 $(INC)/image_routines.h\
 $(INC)/mathconst.h\
 $(INC)/parse_sexagesimal.h\
 $(INC)/primary_beam.h\
 $(INC)/ra_dec_string.h\
 $(INC)/random.h\
 $(INC)/slalib.h\
 $(INC)/uvsubs.h

funcL.o: funcL.c\
 $(INC)/funcL.h\
 $(INC)/rasterizador.h\
 $(INC)/mallaVoronoi.h\
 $(INC)/mockcbiRoutines.h\
 $(INC)/image_routines.h\
 $(INC)/uvsubs.h\
 $(INC)/newstring.h

crear_imagen.o: crear_imagen.c \
 $(INC)/mallaVoronoi.h\
 $(INC)/funcL.h\
 $(INC)/fourn.h \
 $(INC)/nrutil.h

crear_intensidades.o: crear_intensidades.cpp \
 $(INC)/mallaVoronoi.h\
 $(INC)/funcL.h\
 $(INC)/fourn.h \
 $(INC)/nrutil.h\
 $(INC)/image_routines.h\
 $(INC)/rasterizador.h

do_mockcbi.o: do_mockcbi.c \
 $(INC)/mallaVoronoi.h\
 $(INC)/funcL.h\
 $(INC)/fourn.h \
 $(INC)/nrutil.h

comprobacion_gradiente.o: comprobacion_gradiente.cpp\
 $(INC)/funcion.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionCero.h\
 $(INC)/funcionMemCBI.h

comprobacion_gradiente_dVdI.o: comprobacion_gradiente_dVdI.cpp\
 $(INC)/funcion.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionCero.h\
 $(INC)/funcionMemCBI.h\
 $(INC)/reconstructor.h

comparar_visibilidades.o: comparar_visibilidades.cpp\
 $(INC)/funcion.h\
 $(INC)/funcionBayesVoronoiCBI.h\
 $(INC)/funcionCero.h\
 $(INC)/funcionMemCBI.h\
 $(INC)/reconstructor.h

# busquedaNPols.o: busquedaNPols.cpp\
#  $(INC)/image_routines.h\
#  $(INC)/uvsubs.h\
#  $(INC)/mockcbiRoutines.h\
#  $(INC)/reconstructorGCMemCBI.h\
#  $(INC)/funcionCero.h\
#  $(INC)/funcionEntropiaNatural.h

residuos.o: residuos.cpp\
 $(INC)/mockcbiRoutines.h

#--------para Gibbs Sampling--------------------------------------------
distribucionCondicional.o: distribucionCondicional.cpp\
 $(INC)/distribucionCondicional.h

distribucionCondicionalBVCBI.o: distribucionCondicionalBVCBI.cpp\
 $(INC)/distribucionCondicionalBVCBI.h\
 $(INC)/distribucionCondicional.h

distribucionCondicionalVoronoiI.o: distribucionCondicionalVoronoiI.cpp\
 $(INC)/distribucionCondicionalVoronoiI.h\
 $(INC)/distribucionCondicional.h

distribucionCondGaussiana.o: distribucionCondGaussiana.cpp\
 $(INC)/distribucionCondGaussiana.h\
 $(INC)/distribucionCondicional.h

distribucionCondicionalCBI.o: distribucionCondicionalCBI.cpp\
 $(INC)/distribucionCondicionalCBI.h\
 $(INC)/distribucionCondicional.h\
 $(INC)/mockcbiRoutines.h

distribucionCondicionalMemCBI.o: distribucionCondicionalMemCBI.cpp\
 $(INC)/distribucionCondicionalMemCBI.h\
 $(INC)/distribucionCondicionalCBI.h\
 $(INC)/distribucionCondicional.h

distribucionGibbs.o: distribucionGibbs.cpp\
 $(INC)/distribucionGibbs.h\
 $(INC)/distribucionPbb.h


#-----------------------------------------------------------------------
# Program generation.
#-----------------------------------------------------------------------

#PRUEBA_OBJ = prueba.o parse_sexagesimal.o uvsubs.o random.o primary_beam.o cbi_beam.o image_routines.o ra_dec_string.o fourn.o
VIR_OBJ = mallaVoronoi.o arista.o punto.o nodoLista.o lista.o poligono.o funcL.o\
 parse_sexagesimal.o uvsubs.o random.o primary_beam.o cbi_beam.o image_routines.o\
 ra_dec_string.o fourn.o newstring.o mockcbiRoutines.o frprmn.o brent.o linmin.o f1dim.o\
 mnbrak.o nrutil.o complex.o gauss.o rasterizador.o

VIRPP_OBJ = VIR++.o reconstructor.o reconstructorGC.o funcionCBI.o funcionBayesVoronoiCBI.o	\
 inicializadorVoronoiUniforme.o distribucionNormal.o distribucionMetropolis.o			\
 reconstructorBayesiano.o funcionExp.o funcionPbbVoronoiCBI.o reconstructorGenetico.o		\
 pikaia.o reconstructorIncremental.o reconstructorIndividual.o reconstructorGCMemCBI.o		\
 incrementadorMallaPixel.o funcionReconstructor.o funcionMemCBI.o funcionEntropiaNatural.o	\
 incrementadorMallaAristas.o funcionChi2VoronoiImagen.o incrementadorMallaImagen.o		\
 inicializadorVoronoiPixelizado.o inicializadorVoronoiPixelizadoAlternado.o			\
 funcionBayesVoronoiCBII.o funcionNegativa.o funcionDividida.o inicializadorArchivo.o		\
 funcionVoronoiI.o funcionVoronoiX.o distribucionHamiltoniana.o	reconstructorIncrementalDoble.o	\
 funcionEjemploH.o distribucionCondicionalBVCBI.o distribucionCondicionalCBI.o			\
 distribucionGibbs.o distribucionCondicional.o distribucionCondicionalMemCBI.o			\
 distribucionCondicionalVoronoiI.o funcionEntropiaIlogI.o					\
 reconstructorVIR.o reconstructorVIRIterativo.o inicializadorParametros.o
# funcionChi2VoronoiImagenN.o

$(BIN)/VIR++: $(VIR_OBJ) $(VIRPP_OBJ) 
	$(CXX) -o $@ $(VIR_OBJ) $(VIRPP_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses -lg2c


CREAR_IMAGEN_OBJ = crear_imagen.o mallaVoronoi.o arista.o punto.o	\
 nodoLista.o lista.o poligono.o funcL.o parse_sexagesimal.o uvsubs.o	\
 random.o primary_beam.o cbi_beam.o image_routines.o ra_dec_string.o	\
 fourn.o newstring.o mockcbiRoutines.o frprmn.o brent.o linmin.o	\
 f1dim.o mnbrak.o nrutil.o complex.o gauss.o rasterizador.o

$(BIN)/crear_imagen: $(CREAR_IMAGEN_OBJ)
	$(CC) -o $@ $(CREAR_IMAGEN_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

COMPROBACION_GRADIENTE_OBJ = comprobacion_gradiente.o funcionBayesVoronoiCBI.o funcionCBI.o\
 funcionMemCBI.o funcionEntropiaNatural.o

$(BIN)/comprobacion_gradiente: $(VIR_OBJ) $(COMPROBACION_GRADIENTE_OBJ)
	$(CXX) -o $@ $(VIR_OBJ) $(COMPROBACION_GRADIENTE_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

COMPROBACION_GRADIENTE_DVDI_OBJ = comprobacion_gradiente_dVdI.o	\
 funcionBayesVoronoiCBI.o funcionCBI.o				\
 funcionMemCBI.o funcionEntropiaNatural.o			

$(BIN)/comprobacion_gradiente_dVdI: $(VIR_OBJ) $(COMPROBACION_GRADIENTE_DVDI_OBJ)
	$(CXX) -o $@ $(VIR_OBJ) $(COMPROBACION_GRADIENTE_DVDI_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

COMPARAR_VISIBILIDADES_OBJ = comparar_visibilidades.o	\
 funcionBayesVoronoiCBI.o funcionCBI.o				\
 funcionMemCBI.o funcionEntropiaNatural.o			

$(BIN)/comparar_visibilidades: $(VIR_OBJ) $(COMPARAR_VISIBILIDADES_OBJ)
	$(CXX) -o $@ $(VIR_OBJ) $(COMPARAR_VISIBILIDADES_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

DO_MOCKCBI_OBJ = do_mockcbi.o mallaVoronoi.o arista.o punto.o		\
 nodoLista.o lista.o poligono.o funcL.o parse_sexagesimal.o uvsubs.o	\
 random.o primary_beam.o cbi_beam.o image_routines.o ra_dec_string.o	\
 fourn.o newstring.o mockcbiRoutines.o frprmn.o brent.o linmin.o	\
 f1dim.o mnbrak.o nrutil.o complex.o gauss.o rasterizador.o

$(BIN)/do_mockcbi: $(DO_MOCKCBI_OBJ)
	$(CC) -o $@ $(DO_MOCKCBI_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

QUANTA_OBJ = quanta.o funcionCBI.o funcionBayesVoronoiCBI.o		\
 funcionPbbVoronoiCBI.o inicializadorArchivo.o funcionNegativa.o	\
 funcionDividida.o funcionVoronoiI.o reconstructorGC.o reconstructor.o	\
 funcionMemCBI.o reconstructorGCMemCBI.o funcionEntropiaNatural.o

$(BIN)/quanta: $(VIR_OBJ) $(QUANTA_OBJ) 
	$(CXX) -o $@ $(VIR_OBJ) $(QUANTA_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses -lg2c


COMPROBACION_DISTRIBUCION_OBJ = comprobacionDistribucion.o  distribucionCondicionalBVCBI.o	\
 funcionBayesVoronoiCBI.o funcionCBI.o inicializadorArchivo.o					\
 distribucionGibbs.o distribucionCondicional.o	distribucionCondicionalMemCBI.o			\
 distribucionCondGaussiana.o distribucionNormal.o funcionEntropiaNatural.o 			\
 distribucionCondicionalCBI.o funcionMemCBI.o

$(BIN)/comprobacionDistribucion: $(COMPROBACION_DISTRIBUCION_OBJ) $(VIR_OBJ)
	$(CXX) -o $@ $(COMPROBACION_DISTRIBUCION_OBJ) $(VIR_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses -lg2c

NUMEROCENTROS_OBJ = numeroCentros.o funcionCBI.o funcionBayesVoronoiCBI.o	\
 funcionPbbVoronoiCBI.o inicializadorArchivo.o funcionNegativa.o		\
 funcionDividida.o funcionVoronoiI.o reconstructorGC.o reconstructor.o		\
 inicializadorVoronoiUniforme.o funcionChi2VoronoiImagen.o			\
 incrementadorMallaImagen.o reconstructorIncremental.o

$(BIN)/numeroCentros: $(VIR_OBJ) $(NUMEROCENTROS_OBJ) 
	$(CXX) -o $@ $(VIR_OBJ) $(NUMEROCENTROS_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses -lg2c

CALCULARFUNCION_OBJ = calcular_funcion.o funcionCBI.o			\
 funcionBayesVoronoiCBI.o funcionMemCBI.o funcionEntropiaNatural.o

$(BIN)/calcular_funcion: $(VIR_OBJ) $(CALCULARFUNCION_OBJ) 
	$(CXX) -o $@ $(VIR_OBJ) $(CALCULARFUNCION_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses -lg2c

CREAR_INTENSIDADES_OBJ = crear_intensidades.o mallaVoronoi.o arista.o	\
 punto.o								\
 nodoLista.o lista.o poligono.o funcL.o parse_sexagesimal.o uvsubs.o	\
 random.o primary_beam.o cbi_beam.o image_routines.o ra_dec_string.o	\
 fourn.o newstring.o mockcbiRoutines.o frprmn.o brent.o linmin.o	\
 f1dim.o mnbrak.o nrutil.o complex.o gauss.o rasterizador.o

$(BIN)/crear_intensidades: $(CREAR_INTENSIDADES_OBJ)
	$(CXX) -o $@ $(CREAR_INTENSIDADES_OBJ)  -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

# BUSQUEDANPOLS_OBJ = busquedaNPols.o reconstructor.o			\
#  reconstructorGCMemCBI.o reconstructorGC.o reconstructorIncremental.o	\
#  funcionCBI.o funcionMemCBI.o funcionBayesVoronoiCBI.o			\
#  funcionChi2VoronoiImagen.o funcionEntropiaNatural.o			\
#  incrementadorMallaImagen.o inicializadorVoronoiUniforme.o
# 
# $(BIN)/busquedaNPols: $(BUSQUEDANPOLS_OBJ) $(VIR_OBJ)
# 	$(CXX) -o $@ $(BUSQUEDANPOLS_OBJ) $(VIR_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses

RESIDUOS_OBJ = residuos.o mockcbiRoutines.o parse_sexagesimal.o uvsubs.o	\
  random.o primary_beam.o cbi_beam.o image_routines.o ra_dec_string.o		\
 fourn.o newstring.o complex.o

$(BIN)/residuos: $(RESIDUOS_OBJ) 
	$(CXX) -o $@ $(RESIDUOS_OBJ) -lm -L$(LIB) -lsla $(CFLIB) -ltecla -lncurses
