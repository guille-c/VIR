#include "distribucionCondicionalMemCBI.h"

DistribucionCondicionalMemCBI::DistribucionCondicionalMemCBI (char *nombreImagen, char **nombresVis,
							      int nVis, Funcion *entropia, double lambda, int N,
							      Inicializador *ini)
  : DistribucionCondicionalCBI (nombreImagen, nombresVis, nVis) {

  this->entropia = entropia;
  this->n_pars = fg_image->npixels;
  double *p = new double[n_pars];
  if (entropia->getNPars() != n_pars) {
    cerr << "ERROR en FuncionMemCBI::FuncionMemCBI, entropia->getNPars() != n_pars\n";
    cerr << "         " << entropia->getNPars() << " != " << n_pars << "\n";
    exit (1);
  }
  chi2 = 0;
  S = 0;
  minpix = fabs (ruido / 10);
  this->lambda = lambda;
  this->N = N;
  if (ini != 0) {
    ini->inicializar (p, n_pars);
    for (int i = 0; i < n_pars; i++) {
      fg_image->pixels[i] = p[i];
      p[i] /= ruido;
      //cout << fg_image->pixels[i] << "\n";
    }
  }
  for (int arch = 0; arch < nVis; arch++) {
    mockcbi_sub (cmb_image, fg_image, header_obs[arch], samples_obs[arch], samples_mod[arch], beam);
  }
  
  S = entropia->f(p);
  chi2 = calcularChi2();
  L_old = chi2 / 2 - lambda * S;

  delete [] p;  
}

double DistribucionCondicionalMemCBI::generar (double *p, int i) {
  
  double x, a, b, fTotal = 0, f_rand;
  int cont, n;
  L_old = chi2 / 2 - lambda * S;
 
  a = 0;
  b = N;
  n = (int) ceil (b - a);
  double fs [n];

  cont = 0;
  for (x = a; x <= b; x += 1) {
    fs [cont] = this->evaluar (p, x, i);
    fTotal += fs[cont++];
    //cout  << "x1 = " << x << ", fTotal = " << fTotal <<  "\n";
  }
  cout  << "x1 = " << x << "\n";
  
  f_rand = fTotal * ((double) rand() /  RAND_MAX);
  cout << "fTotal I = " << fTotal << ", frand = " << f_rand << "\n";
  fTotal = 0;
  cont = 0;
  for (x = a; fTotal < f_rand && cont < n; x += 1) {
    fTotal += fs[cont++];
    //cout  << "x2 = " << x << ", fTotal = " << fTotal << "\n";
  }
  cout  << "x2 = " << x - 1<< "\n";
  if (x - 1 < minpix / ruido) {
    x = minpix / ruido + 1;
  }
  evaluar (p, x - 1, i);

  fg_image->pixels[i] = (x - 1) * ruido;
  
  return x - 1;
}

double DistribucionCondicionalMemCBI::evaluar (double *p) {
  //return exp (-f->f(p));
  //return  1 / f->f(p);

  // FuncionMemCBI f2 (nombreImagen, nombresVis, nVis, entropia, lambda);
  //double ret2 = f2.f(p);

  this->evaluar (p, p[0], 0);

  /*
  if (fabs(ret2 - chi2/2 + S) > 1e-3) {
    cerr << "ERROR en DistribucionCondicionalMemCBI::evaluar\n";
    cerr << "ret1 = " << chi2 / 2 - S << ", ret2 = " << ret2 << "\n";
    cerr << "1) chi2 = " << chi2 << ", S = " << S << "\n";
    cerr << "2) chi2 = " << f2.getChi2() << ", S = " << f2.getS() << "\n";
    exit (1);
  }
  */
  return 1/ (chi2 / 2 - S);
}

double DistribucionCondicionalMemCBI::evaluar (double *p, double p_i, int i) {
  double ret = 0;
  //double ret = f->f(p);
  bool igual = true;
  double pi_old = p[i];

  //double chi2_old = chi2;
  //double S_old = S;
  //cout << "OLD: chi2 = " << chi2 << ", S = " << S << ", ret = " << L_old - chi2 / 2 + lambda * S << ", lambda = " << lambda << "\n";
  //for (int j = 0; j < n_pars && igual; j++) {
  //if (p[j] < minpix / ruido) {
  //p[j] = minpix / ruido;
  //}
  //}

  if (p_i < minpix / ruido) {
    p_i = minpix / ruido;
  }
  for (int j = 0; j < n_pars && igual; j++) {
    if (fabs (p[j] - fg_image->pixels[j]/ruido) > 1e-3 && j != i) {
      igual = false;
      cout << "p[" << j << "] = " << p[j] << " != " << fg_image->pixels[j] / ruido << "\n";
      cout << "diferencia = " << fabs (p[j] - fg_image->pixels[j]/ruido) << "\n";
    }
  }

  p[i] = p_i;
  S = entropia->f(p);
  p[i] = pi_old;

  if (igual) {
    double raimage, decimage, dx, dy, kx, seno, coseno, x0, y0, x, y;
    double obsra, obsdec, lobs, mobs, ra_offset = 0.0, raizxy, dVR, dVI;
    int ix, iy;
    
    ix = i % fg_image->size[0];
    iy = i / fg_image->size[0];
    if (ix + iy * fg_image->size[0] != i) {
      cerr << "ERROR, ix, iy\n" << "ix = " << ix << ", iy = " << iy << ", i = " << i << "\n";
      exit (1);
    }
    //cout << "ix = " << ix << ", iy = " << iy << ", i = " << i << "\n";
      
    dx = fg_image->cdelt[0] * RPDEG;   /* radians */
    dy = fg_image->cdelt[1] * RPDEG;   /* radians */
    //double pixel = fabs(fg_image->cdelt[0] * RPDEG * fg_image->cdelt[1] * RPDEG);    /* pixel solid angle in ster */

    /* coordinates of the image reference pixel in radians */
    raimage = fg_image->crval[0] * RPDEG;
    decimage = fg_image->crval[1] * RPDEG;

    chi2 = 0;
    for (int arch = 0; arch < nVis; arch++) {                                   /* filenames loop */
      obsra = ra_offset + header_obs[arch]->obsra*RPDEG;
      obsdec = header_obs[arch]->obsdec*RPDEG;  
      /* direction cosines of the phase center in the image coordinate system */
      direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
      /* Find the phase center in the pixel coordinate system of the image (where pixel numbers are zero-based) */
      x0 = (fg_image->crpix[0] - 1.0) + lobs / dx;
      y0 = (fg_image->crpix[1] - 1.0) + mobs / dy;

      x = (ix - x0) * dx;          /* radians */
      y = (iy - y0) * dy;          /* radians */

      raizxy = sqrt (1 - x*x - y*y);
      //cout << "raizxy = " << raizxy << "\n";

      for (int chan = 0; chan < header_obs[arch]->nchan; chan++) {        /* chan: loop through channels */
	
	for (int samp = 0; samp < header_obs[arch]->nsamp; samp++) {      /* i: loop through samples */
	  
	  if (samples_obs[arch][samp].rdata[chan * 3 + 2] != 0) {
	    /* u*freq : radians , x,y : radians */
	    kx = (samples_mod[arch][samp].u * x + samples_mod[arch][samp].v * y) * header_obs[arch]->iffreq[chan];
	    seno = sin(2 * PI * kx);
	    coseno = cos(2 * PI * kx);
	    
	    //cout << " Cambiando " << samples_mod[arch][samp].rdata[3 * chan] << " por ";
	    
	    dVR = (coseno * atten[arch][i + 1][chan] * 
		   (p_i * ruido - fg_image->pixels[i]) * fg_scale / raizxy);
	    dVI = (seno * atten[arch][i + 1][chan] * 
		   (p_i * ruido - fg_image->pixels[i]) * fg_scale / raizxy);
	    
	    samples_mod[arch][samp].rdata[3 * chan]     += dVR;
	    samples_mod[arch][samp].rdata[3 * chan + 1] += dVI;
	    
	    chi2 += samples_obs[arch][samp].rdata[chan * 3 + 2] * (SQR(samples_mod[arch][samp].rdata[chan * 3] - samples_obs[arch][samp].rdata[chan * 3])
								   + SQR(samples_mod[arch][samp].rdata[chan * 3 + 1] - samples_obs[arch][samp].rdata[chan * 3 + 1]));
	    samples_mod[arch][samp].rdata[3 * chan]     -= dVR;
	    samples_mod[arch][samp].rdata[3 * chan + 1] -= dVI;
	    
	    //cout <<  samples_mod[arch][samp].rdata[3 * chan] << "\n";
	    //cout << "dx = " << dx << ", dy = " << dy << ", coseno = " << coseno << ", atten = "
	    //     << atten[arch][i + 1][chan] << ", dI = " << (p_i * ruido - fg_image->pixels[i]) 
	    //     << ", fgscale = " << fg_scale << "\n";
	    //cout << "Total = " << (dx * dy * seno * atten[arch][i + 1][chan] * 
	    //			 (p_i * ruido - fg_image->pixels[i]) * fg_scale / raizxy) << "\n";
	  }
	}
      }
    }
  }
  else {
    cout << "DISTINTO\n";
    for (int j = 0; j < n_pars; j++) {
      if (p[j] < minpix / ruido) {
	p[j] = minpix / ruido;
      }
      fg_image->pixels[j] = p[j] * ruido;
    }
    
    fg_image->pixels[i] = p_i * ruido;
    /* build the modeled visibilities */
    for (int arch = 0; arch < nVis; arch++) {
      mockcbi_sub (cmb_image, fg_image, header_obs[arch], samples_obs[arch], samples_mod[arch], beam);
    }
    chi2 = calcularChi2();
  }
  //chi2 = calcularChi2();
  //fg_image->pixels[i] = x * ruido;

  ret = exp (L_old - chi2 / 2 + lambda * S);
  //cout << "L_old = " << L_old << ", chi2 = " << chi2 << ", S = " << S << ", ret = " << L_old - chi2 / 2 + lambda * S << ", lambda = " << lambda << "\n";
  /*
  cout << "L_old = " << L_old << ", chi2 = " << chi2 << ", S = " << S << ", ret = " << L_old - chi2 / 2 + lambda * S << ", lambda = " << lambda << "\n";

  FuncionMemCBI f2 (nombreImagen, nombresVis, nVis, entropia, lambda);
  p[i] = pi_old;
  double ret2 = f2.f(p);
  cout << "ret2 1 = " << ret2 << ", chi2 1 = " << f2.getChi2() << ", S 1 = " << f2.getS() << "\n";
  cout << "delta chi2 old = " << chi2_old - f2.getChi2() << ", delta S_old = " << S_old - f2.getS() << "\n";
  //cout << "DistribucionCondicionalMemCBI::evaluar 1\n";
  p[i] = x;
  ret2 -= f2.f(p);
  cout << "ret2 2 = " << ret2 << ", chi2 1 = " << f2.getChi2() << ", S 1 = " << f2.getS() << "\n";
  cout << "delta chi2 new = " << chi2 - f2.getChi2() << ", delta S_new = " << S - f2.getS() << "\n";
  //cout << "DistribucionCondicionalMemCBI::evaluar 2\n";
  ret2 = exp (ret2);
  //p[i] = x_old;
  
  cout << "ret = " << ret << ", ret2 = " << ret2 << "\n\n";
  */
  //p[i] = x;
  return ret;
}

double DistribucionCondicionalMemCBI::calcularChi2 () {
  double ret = 0;

  for (int arch = 0; arch < nVis; arch++) {
    for (int samp = 0; samp < header_obs[arch]->nsamp; samp++) {
      for (int chan = 0; chan < header_obs[arch]->nchan; chan++) {
	if (samples_mod[arch][samp].u != samples_obs[arch][samp].u) {
	  cout << "SCANDALE in DistribucionCondicionalMemCBI::calcularChi2:  samples_mod[arch][samp].u != samples_obs[arch][samp].u\n";
	  cout << samples_mod[arch][samp].u << " != " << samples_obs[arch][samp].u << "\n";
	  exit (1);
	}
	
	ret += samples_obs[arch][samp].rdata[chan * 3 + 2] * (SQR(samples_mod[arch][samp].rdata[chan * 3] - samples_obs[arch][samp].rdata[chan * 3])
							       + SQR(samples_mod[arch][samp].rdata[chan * 3 + 1] - samples_obs[arch][samp].rdata[chan * 3 + 1]));
      }
    }
  }
  return ret;
}

DistribucionCondicionalMemCBI::~DistribucionCondicionalMemCBI() {
  //delete f;
}

double DistribucionCondicionalMemCBI::evaluar_old (double *p, double x, int i) {
  FuncionMemCBI f (nombreImagen, nombresVis, nVis, entropia, 0.0, lambda, 1e100, -1);
  double ret = f.f(p);
  double x_old = p[i];

  //cout << "DistribucionCondicionalMemCBI::evaluar 1\n";
  p[i] = x;
  ret -= f.f(p);
  //cout << "DistribucionCondicionalMemCBI::evaluar 2\n";
  ret = exp (ret);
  p[i] = x_old;
  return ret;
}
