#ifndef VIRPP_FUNCION
#define VIRPP_FUNCION

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <string>
using namespace std;

class Funcion {
 protected:
  int n_pars;
  string i2s (int numero, int largo){
    char *numeroC = new char [largo + 1];
    int i;
    
    for (i = 0; i < largo; i++){
      numeroC [i] = '0' + (int) ((numero % (int) pow (10.0, (double) largo - i))
				 / pow (10.0, largo - i - 1.0));
    }
    numeroC [i] = '\0';
    
    string s;
    s += numeroC;
    delete [] numeroC;
    return s;
  }
 public:
  Funcion () {n_pars = 0;}
  Funcion (int n) {this->n_pars = n;}
  virtual double f(double pars[]) = 0;
  virtual void df(double pars[], double grad[]) = 0;
  int getNPars() {return n_pars;}
  virtual void setNPars (int n_pars) {this->n_pars = n_pars;}
  virtual void guardarInfo (string info){}

  virtual ~Funcion() {}
};



#endif
