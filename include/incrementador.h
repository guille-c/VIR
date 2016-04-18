#ifndef INCREMENTADOR
#define INCREMENTADOR

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <stdarg.h>
#include <string>
#include "funcion.h"

class Incrementador {
 public:
  virtual double incrementar (double **p, Funcion *f) = 0;
};

#endif
