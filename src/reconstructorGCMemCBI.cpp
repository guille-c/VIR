#include "reconstructorGCMemCBI.h"

//==========WRAPPER=========

ReconstructorGCMemCBI *r1_wrap_mem;

float f_wrap_mem(float pars[]){
  int i, n = r1_wrap_mem->getFunc()->getNPars();
  double *parsD = new double [n];
  float ret;

  for (i = 0; i < n; i++) {
    parsD[i] = (double) pars[i + 1];
    //cout << "parsD[" << i << "] = " << parsD[i] << "\n";
  }
  ret = (float) r1_wrap_mem->getFunc()->f (parsD);
  //ret = (float) rand() / RAND_MAX;
  delete [] parsD;
  return ret;
}

void df_wrap_mem (float pars[], float grad[]) {
  int i, n = r1_wrap_mem->getFunc()->getNPars();
  double *parsD = new double [n], *gradD = new double [n];

  r1_wrap_mem->getFunc()->guardarInfo ("GC" + i2s (r1_wrap_mem->getIter(), 3));

  for (i = 0; i < n; i++) {
    parsD[i] = (double) pars[i + 1];
    //cout << "parsD[" << i << "] = " << parsD[i] << "\n";
  }

  std::ofstream archivo; 
  archivo.open ("GC_L.dat", ios::app);
  archivo << r1_wrap_mem->getIter() << "\t" << r1_wrap_mem->getFunc()->f(parsD) << "\t" 
	  << r1_wrap_mem->getFunc()->getChi2() << "\t" << r1_wrap_mem->getFunc()->getS() 
	  << "\t" << r1_wrap_mem->getFunc()->getLambda () << "\n";
  archivo.close();

  r1_wrap_mem->getFunc()->df (parsD, gradD);

  for (i = 0; i < n; i++) {
    grad[i + 1] = (float) gradD[i];
    //grad[i + 1] = (float) rand() / RAND_MAX;
  }
  
  int iter = r1_wrap_mem->getIter();
  if (iter >= r1_wrap_mem->getNIterMax()) {
    for (i = 0; i < n; i++) {
      grad[i + 1] = 0.0;
    }
  }
  int Nmax = r1_wrap_mem->getNmax();
  if (Nmax > 0 && iter <= Nmax) {
    double lambda = 0;//iter / Nmax;
    //double lambda = sqrt(iter / Nmax);
    //double lambda = iter * iter / (Nmax * Nmax);
    r1_wrap_mem->getFunc()->setLambda (lambda);
    cout << "--------> Cambiando lambda por " << lambda << "\n";
  }
  else {
    r1_wrap_mem->getFunc()->setLambda (1.0);
  }

  delete [] gradD;
  delete [] parsD;
}

//==========================

ReconstructorGCMemCBI::ReconstructorGCMemCBI (FuncionMemCBI *f, Inicializador *ini, 
					      double ftol, int Nmax, double xini, int nIterMax)
  : Reconstructor (f->getNPars(), ini, xini) {
  this->ftol = ftol; this->f = f;
  this->Nmax = Nmax;
  this->nIterMax = nIterMax;
  //r1_wrap_mem = new ReconstructorGCMemCBI (f, )
  r1_wrap_mem = this;
  //this->f->setLambda (0.0);
}


/*********************************************************
   Minimiza utilizando Polak-Ribiere de Numerical recipes.
   Retorna el valor final de la funcion a minimizar.
*********************************************************/

double ReconstructorGCMemCBI::minimizador () {
  //float f_aux (float pars[]) { return f(pars);}
  
  int i, n = f->getNPars();
  float fret, *pF = new float [n];
  long t;
  fstream archivo;
  char fileout[30];
  std::ofstream archivo1 ("GC_L.dat");
  archivo1.close();

  //func = &ReconstructorGCMemCBI::f;
  if (f->getNPars() < 10) {
    sprintf(fileout, "GC_0%d.dat", f->getNPars());
  }
  else {
    sprintf(fileout, "GC_%d.dat", f->getNPars());
  }
  archivo.open (fileout, ios::out);
 
  //r1_wrap_mem = this;

  for (i = 0; i < n; i++) {
    pF[i] = (float) p[i];
  }

  t = time(0);
  frprmn(pF - 1, f->getNPars(), ftol, &(iter), &fret, &f_wrap_mem, df_wrap_mem);
  t = time(0) - t;

  for (i = 0; i < n; i++) {
    p[i] = (double) pF[i];
  }
  delete [] pF;

  printf("Terminado P-R en %d segundos y %d iteraciones\n", 
	 (int) t, iter);
  //free (p_aux);
  /*
  if (f->getNPars() < 10) {
    sprintf(fileout, "!MEM_NR_00%d.fits", f->getNPars());
  }
  else if (fL->n_pols < 100) {
    sprintf(fileout, "!MEM_NR_0%d.fits", f->getNPars());
  }
  else {
    sprintf(fileout, "!MEM_NR_%d.fits", f->getNPars());
  }
  do_write_fits(fL->fg_image, fileout);
  */
  archivo << "Terminado P-R en " << t << " segundos y " 
	  << iter << " iteraciones\n";
  archivo << "Chi2 = " << fret << "\n";
  archivo.close();
  /*
  if (fL->n_pols < 10) {
    sprintf(fileout, "MEM_NR_00%d.dat", fL->n_pols);
  }
  else if (fL->n_pols < 100) {
    sprintf(fileout, "MEM_NR_0%d.dat", fL->n_pols);
  }
  else {
    sprintf(fileout, "MEM_NR_%d.dat", fL->n_pols);
  }
  imprimirMallaArchivo (fL->malla, fileout);
  //free (p);
  */
  return fret;
}

double ReconstructorGCMemCBI::run () {
  //return 0.0;
  int i, j, nVis;
  double func = 1e299, x, y, valor, 
    funcMin = 1e300, init_value = p[2];
  fstream archivo;

  cout << "rec npars = " << f->getNPars() << "\n";

  archivo.open("L_n.dat", ios::out);
  //n = this->f->getNPars() / 3;
  
  //for (i = 0; i < n; i++) {
  //cout << "  p[" << i << "] = (" << p[3 * i] << ", " << p[3 * i + 1] 
  //<< ", " << p[3 * i + 2] << ")\n";
    //imprimirLog (" p[%d] = (%g, %g, %g)\n", i, p[3 * i], 
    //p[3 * i + 1], p[3 * i + 2]);
  //}
  //imprimirLog ("    n = %d\n", n);

  func =  minimizador ();
  //imprimirLog ("      func = %g\n", func);
  string s = "Final";
  cout << r1_wrap_mem << ", " << this << "\n";
  double ret2 = r1_wrap_mem->getFunc()->f(p);
  r1_wrap_mem->getFunc()->guardarInfo (s);

  cout << "func = " << func << ", x = " << x << ", y= " << y 
       << ", valor = " << valor << "\n";
    
  cout << "FIN, func = " << func << "\n";
  archivo << "FIN, func = " << func << "\n";
  archivo << "FIN, func2 = " << ret2 << "\n";
  archivo.close ();
  return func;
}
