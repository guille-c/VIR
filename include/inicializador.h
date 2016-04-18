#ifndef VIRPP_INICIALIZADOR
#define VIRPP_INICIALIZADOR

class Inicializador {
 public:
  virtual void inicializar (double pars[], int n, double init = 0.0) = 0;
};

#endif
