# The parent directory.

CODE = ..

# The binary directory.

BIN = $(CODE)/bin

# The library directory.

LIB = $(CODE)/lib

# The include directory.

INC = $(CODE)/include

# Optimisation flag (use 'make OPT=-g' if you want to enable debugging).

OPT = -O

# The C compiler and its run-time library.

CC = gcc -g
CLIB = `gcc -print-libgcc-file-name`
CXX = g++ -g

# C compiler flags.

CFLAGS = $(OPT) -Wall -Wmissing-prototypes -Wmissing-declarations -I$(INC)
CXXFLAGS = -I$(INC) -Wno-deprecated
# -D_POSIX_C_SOURCE=199506L -D__EXTENSIONS__

# Fortran compiler

FORT = g77
#FORT = f77

# Fortran compiler flags

FFLAGS = $(OPT) -Wunused -Wimplicit -fno-backslash -fugly-init
#FFLAGS = -u -xl

# The library randomizer (not required under Solaris).

RANLIB = @:

# The PGPLOT library and include-file directory.

PGDIR = /usr/local/pgplot

# The loader flags needed to link with X11

XLIB = -L/usr/X11R6/lib64 -lX11

# The loader flags needed to link with threads

THLIB = -lpthread -lrt -lnsl

# Extra libraries need with CFITSIO

CFLIB = -lcfitsio

# The loader flags needed to link with Tcl/Tk.

TKLD = -L/usr/local/lib64 -ltk8.0 -ltcl8.0 -lX11 -ldl

# The flags for including Tcl/Tk files.

TK_INC = -I/usr/local/include

# Flags needed by ld when creating a shared library.

DLFLAGS = -G -i

#-----------------------------------------------------------------------
# Default to making all targets.
#-----------------------------------------------------------------------

all: make_program

#-----------------------------------------------------------------------
# Recursive make targets.
#-----------------------------------------------------------------------

make_program:
	@echo ' ';echo 'Building programs in Slalib_f'; echo ' '
	@cd Slalib_f; $(MAKE) "LIB=$(LIB)" "RANLIB=$(RANLIB)"\
 "FORT=$(FORT)" "FFLAGS=$(FFLAGS)"
	@echo ' ';echo 'Building programs in src'; echo ' '
	@cd src; $(MAKE) "BIN=$(BIN)" "INC=$(INC)" "LIB=$(LIB)" \
 "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$(CFLAGS)" "CXXFLAGS=$(CXXFLAGS)"\
 "RANLIB=$(RANLIB)" "PGDIR=$(PGDIR)"\
 "FORT=$(FORT)" "FFLAGS=$(FFLAGS)" "XLIB=$(XLIB)" "THLIB=$(THLIB)" "CFLIB=$(CFLIB)"
	@echo ' ';echo 'Building programs in Simprog'; echo ' '
	@cd Simprog; $(MAKE) "BIN=$(BIN)" "INC=$(INC)" "LIB=$(LIB)" \
 "CC=$(CC)" "CFLAGS=$(CFLAGS)" "RANLIB=$(RANLIB)" "PGDIR=$(PGDIR)"\
 "FORT=$(FORT)" "FFLAGS=$(FFLAGS)" "CFLIB=$(CFLIB)"
