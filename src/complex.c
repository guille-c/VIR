/* Adapted from Numerical Recipes */

#include <math.h>

#include "complex.h"

complex Cadd(complex a, complex b)
{
  complex c;
  c.r = a.r+b.r;
  c.i = a.i+b.i;
  return c;
}

complex Csub(complex a, complex b)
{
  complex c;
  c.r = a.r-b.r;
  c.i = a.i-b.i;
  return c;
}


complex Cmul(complex a, complex b)
{
  complex c;
  c.r = a.r*b.r-a.i*b.i;
  c.i = a.i*b.r+a.r*b.i;
  return c;
}

complex Complex(double re, double im)
{
  complex c;
  c.r = re;
  c.i = im;
  return c;
}

complex Conjg(complex z)
{
  complex c;
  c.r = z.r;
  c.i  =  -z.i;
  return c;
}

complex Cdiv(complex a, complex b)
{
  complex c;
  double r,den;
  if (fabs(b.r) >= fabs(b.i)) {
    r = b.i/b.r;
    den = b.r+r*b.i;
    c.r = (a.r+r*a.i)/den;
    c.i = (a.i-r*a.r)/den;
  } else {
    r = b.r/b.i;
    den = b.i+r*b.r;
    c.r = (a.r*r+a.i)/den;
    c.i = (a.i*r-a.r)/den;
  }
  return c;
}

double Cabs(complex z)
{
  double x,y,ans,temp;
  x = fabs(z.r);
  y = fabs(z.i);
  if (x == 0.0)
    ans = y;
  else if (y == 0.0)
    ans = x;
  else if (x > y) {
    temp=y/x;
    ans = x*sqrt(1.0+temp*temp);
  } else {
    temp = x/y;
    ans = y*sqrt(1.0+temp*temp);
  }
  return ans;
}

double Carg(complex z)
{
  if (z.r == 0.0 && z.i == 0.0)
    return 0.0;
  else
    return atan2(z.i, z.r);
}

complex Csqrt(complex z)
{
  complex c;
  double x,y,w,r;
  if ((z.r == 0.0) && (z.i == 0.0)) {
    c.r = 0.0;
    c.i = 0.0;
    return c;
  } else {
    x = fabs(z.r);
    y = fabs(z.i);
    if (x >= y) {
      r = y/x;
      w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    } else {
      r = x/y;
      w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
    }
    if (z.r >= 0.0) {
      c.r = w;
      c.i = z.i/(2.0*w);
    } else {
      c.i = (z.i >=  0) ? w : -w;
      c.r = z.i/(2.0*c.i);
    }
    return c;
  }
}

complex RCmul(double x, complex a)
{
  complex c;
  c.r = x*a.r;
  c.i = x*a.i;
  return c;
}

