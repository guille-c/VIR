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
#  CFLIB     - linker options for CFITSIO
#-----------------------------------------------------------------------

SHELL=/bin/csh
CPGPLOT=-L$(PGDIR) -lcpgplot -lpgplot
PGPLOT=-L$(PGDIR) -lpgplot

targets: $(BIN)/cbisky $(BIN)/cbifake $(BIN)/polfake $(BIN)/image 

SKYOBJ = cbisky.o four2.o len1.o pangle.o usernm.o nodenm.o dattim.o\
      putout.o error.o udate.o sysver.o keyin.o kctor.o kpack.o keyout.o\
      keyhlp.o kword.o bessj0.o keyput.o kexpre.o keyfnd.o krdlin.o\
      kctod.o kctoi.o kskipb.o keycmp.o keymap.o keychk.o symsub.o\
      prognm.o tsttty.o  gasdev.o ran1.o

FAKOBJ = cbifake.o cbibeam.o four2.o keyin.o kctor.o kpack.o keyout.o keyhlp.o kword.o\
      keyput.o kexpre.o keyfnd.o krdlin.o\
      kctod.o kctoi.o kskipb.o keycmp.o keymap.o keychk.o symsub.o\
      prognm.o tsttty.o usernm.o putout.o len1.o pangle.o\
      error.o timec.o gasdev.o ran1.o

POLFAKOBJ = polfake.o cbibeam.o four2.o keyin.o kctor.o kpack.o keyout.o keyhlp.o kword.o\
      keyput.o kexpre.o keyfnd.o krdlin.o\
      kctod.o kctoi.o kskipb.o keycmp.o keymap.o keychk.o symsub.o\
      prognm.o tsttty.o usernm.o putout.o len1.o pangle.o\
      error.o timec.o gasdev.o ran1.o

IMGOBJ = image.o kctor.o usernm.o sysver.o len1.o nodenm.o pangle.o\
      dattim.o kskipb.o error.o udate.o putout.o kctod.o kctoi.o

MAPOBJ = mapplot.o dattim.o error.o frdt.o getin.o keyin.o kctor.o\
      kpack.o keyout.o keyhlp.o kword.o keyput.o kexpre.o keyfnd.o\
      krdlin.o kctod.o kctoi.o kskipb.o keycmp.o keymap.o keychk.o\
      len1.o mkhist.o nodenm.o prognm.o putout.o sspace.o swap.o\
      symsub.o tsttty.o udate.o upcase.o usernm.o

$(BIN)/cbisky: $(SKYOBJ) $(LIB)/libcfitsio.a
	$(FORT) -o $@  $(SKYOBJ) -L$(LIB) $(CFLIB)

$(BIN)/cbifake: $(FAKOBJ) $(LIB)/libsla_f.a $(LIB)/libcfitsio.a
	$(FORT) -o $@  $(FAKOBJ) -L$(LIB) -lsla_f $(CFLIB)

$(BIN)/polfake: $(POLFAKOBJ) $(LIB)/libsla_f.a $(LIB)/libcfitsio.a
	$(FORT) -o $@  $(POLFAKOBJ) -L$(LIB) -lsla_f $(CFLIB)

$(BIN)/image: $(IMGOBJ) $(LIB)/libcfitsio.a
	$(FORT) -o $@  $(IMGOBJ) -L$(LIB) $(CFLIB)

$(BIN)/mapplot: $(MAPOBJ) $(LIB)/libcfitsio.a
	$(FORT) -o $@ $(MAPOBJ) -L$(LIB) $(CFLIB) $(PGPLOT)

.SUFFIXES:
.SUFFIXES: .o .c .f

.c.o:
	$(CC) $(CFLAGS) -c $<
.f.o:
	$(FORT) $(FFLAGS) -c $<


cbifake.o: cbifake.f mapcom.inc cbifake.inc
polfake.o: polfake.f mapcom.inc polfake.inc

clean: 
	rm -f cbisky cbifake polfake image $(SKYOBJ) $(FAKOBJ) $(IMGOBJ)
