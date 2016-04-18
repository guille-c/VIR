#include "funcionMemCBI.h"

FuncionMemCBI::FuncionMemCBI (char *nombreImagen, char **nombresVis,
			      int nVis, Funcion *entropia, double lambda) 
  : FuncionCBI (nombreImagen, nombresVis, nVis){

  this->entropia = entropia;
  this->n_pars = fg_image->npixels;
  if (entropia->getNPars() != n_pars) {
    cerr << "ERROR en FuncionMemCBI::FuncionMemCBI, entropia->getNPars() != n_pars\n";
    cerr << "         " << entropia->getNPars() << " != " << n_pars << "\n";
    exit (1);
  }
  chi2 = 0;
  S = 0;
  minpix = fabs (ruido / 10);
  this->lambda = lambda;
}

double FuncionMemCBI::f (double pars[]) {
  for (int i = 0; i < n_pars; i++) {
    if (pars[i] < minpix / ruido) {
      pars[i] = minpix / ruido;
    }
    fg_image->pixels[i] = pars[i] * ruido;
  }

  /* build the modeled visibilities */
  for (int arch = 0; arch < nVis; arch++) {
    mockcbi_sub (cmb_image, fg_image, header_obs[arch], samples_obs[arch], samples_mod[arch], beam);
  }

  chi2 = 0;
  for (int arch = 0; arch < nVis; arch++)
    {
      for(int samp = 0; samp < header_obs[arch]->nsamp; samp++)
	{
	  for(int chan = 0; chan < header_obs[arch]->nif; chan++)
	    {
	      //if (!ch_list[chan] || chi2 > 1e15) 
	      //continue;
	      
	      if (samples_mod[arch][samp].u != samples_obs[arch][samp].u) {
		cout << "SCANDALE in FuncionMemCBI::f:  samples_mod[arch][samp].u != samples_obs[arch][samp].u\n";
		exit (1);
	      }
	      
	      chi2 += samples_obs[arch][samp].rdata[chan * 3 + 2] * (SQR(samples_mod[arch][samp].rdata[chan * 3] - samples_obs[arch][samp].rdata[chan * 3]) 
							    + SQR(samples_mod[arch][samp].rdata[chan * 3 + 1] - samples_obs[arch][samp].rdata[chan * 3 + 1]));	    
	      //if (samples_obs[arch][samp].rdata[chan * 3 + 2] > 0)
	      //bestchi2++;
	    }
	}
    }
  
  S = entropia->f(pars);
  //cout << "chi2 = " << chi2 << "\n";
  return chi2/2 - lambda * S;
}

void FuncionMemCBI::df (double pars[], double grad[]) {
  double raimage, decimage, obsra, obsdec, lobs, mobs, dx, dy, x0, y0, x, y, kx, ra_offset = 0.0;
  double Re, Im, dchi2, seno, coseno;
  int k, ix, iy;
  double *dS = new double [n_pars];

  cout << "Calculando df\n";
  for (int i = 0; i < n_pars; i++) {
    if (pars[i] < minpix / ruido) {
      pars[i] = minpix / ruido;
    }
    fg_image->pixels[i] = pars[i] * ruido;
  }

  /* build the modeled visibilities */
  for (int arch = 0; arch < nVis; arch++) {
    mockcbi_sub (cmb_image, fg_image, header_obs[arch], samples_obs[arch], samples_mod[arch], beam);
  }

  entropia->df (pars, dS);
  for (int i = 0; i < n_pars; i++) {
    grad [i] = -dS[i] * lambda;
  }

  /* coordinates of the image reference pixel in radians */
  raimage = fg_image->crval[0] * RPDEG;
  decimage = fg_image->crval[1] * RPDEG;
  
  /* dchi2 calculations */
  for (int arch = 0; arch < nVis; arch++)                                   // filenames loop 
    {
      // phase center "absolute" coordinates in radians 
      obsra = ra_offset + header_obs[arch]->obsra*RPDEG;
      obsdec = header_obs[arch]->obsdec*RPDEG;  
      // direction cosines of the phase center in the image coordinate system
      direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
      /* infinitesimal steps */
      dx = fg_image->cdelt[0] * RPDEG;   // radians
      dy = fg_image->cdelt[1] * RPDEG;   // radians
      // Find the phase center in the pixel coordinate system of the image (where pixel numbers are zero-based)
      x0 = (fg_image->crpix[0] - 1.0) + lobs / dx;
      y0 = (fg_image->crpix[1] - 1.0) + mobs / dy;
      
      // show pixel size and conversion factor 
      //printf ("pixel size: %g [arcmin^2] -> %g [ste], conversion factor: %g [Jy pix^-1 K^-1]\n", pixel / RPARCM / RPARCM, pixel, fg_scale);
      
      // d/dI 0.5chi2  

      k = 1;                                                     // k recorre los pixeles de la imagen
      for (iy = 0; iy < fg_image->size[1]; iy++) {
	y = (iy - y0) * dy;                                    // radians 
	for (ix = 0; ix < fg_image->size[0]; ix++) {
	  //if (noise_image->pixels[k - 1] > noise_cut)
	  //{
	  //k++;
	  //continue;
	  //}
	  x = (ix - x0) * dx;                                // radians
	  for(int chan = 0; chan < header_obs[arch]->nif; chan++) {        // chan: loop through channels
	    //if (!ch_list[chan])
	    //continue;
	    
	    dchi2 = 0;
	    for(int samp = 0; samp < header_obs[arch]->nsamp; samp++) {      // i: loop through samples
	      
	      // u*freq : radians , x,y : radians
	      kx = (samples_mod[arch][samp].u * x + samples_mod[arch][samp].v * y) * header_obs[arch]->iffreq[chan];
	      seno = sin(2 * PI * kx);
	      coseno = cos(2 * PI * kx);
	      
	      // partial contribution to chi2
	      Re = (samples_mod[arch][samp].rdata[3 * chan] - samples_obs[arch][samp].rdata[3 * chan]) * coseno;
	      Im = (samples_mod[arch][samp].rdata[3 * chan + 1] - samples_obs[arch][samp].rdata[3 * chan + 1]) * seno ;
	      dchi2 += samples_mod[arch][samp].rdata[3 * chan + 2] * (Re + Im);
	    }
	    dchi2 *= fg_scale * atten[arch][k][chan];
	    grad[k - 1] += dchi2 * ruido;
	  }
	  //cout << "grad[" << k - 1<< "] = " << dchi2 << "\n";
	  k++;
	}
      }
    }
  
  cout << "Calculado df\n";
  delete [] dS;

  //for (int i = 0; i < n_pars; i++) {
  //cout << "grad[" << i << "] = " << grad[i] << "\n";
  //}

}

void FuncionMemCBI::guardarInfo (string info) {
    // Para este caso usaremos info como el sufijo de los archivos a
  // guardar.
  string nombre_archivo;
  
  nombre_archivo = "!MEM_CBI" + info + ".fits";

  char *nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  
  cout << "guardando " << nombreC << "\n";
  do_write_fits (fg_image, nombreC);

  delete [] nombreC;

  nombre_archivo = "MEM_CBI" + info + ".dat";

  nombreC = new char [nombre_archivo.length() + 1];
  nombre_archivo.copy(nombreC, string::npos);
  nombreC [nombre_archivo.length()] = 0;
  
  cout << "guardando " << nombreC << "\n";

  ofstream archivo (nombreC);
  
  archivo << n_pars << "\n";
  for (int i = 0; i < n_pars; i++) {
    archivo << fg_image->pixels[i] << "\n";
  }

  delete [] nombreC;
  
}

//FuncionMemCBI& FuncionMemCBI::operator= (const FuncionMemCBI&) {
//}

void FuncionMemCBI::setLambda (double lambda) {
  this->lambda = lambda;
}
