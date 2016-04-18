typedef struct COMPLEX {double r,i;} complex;

complex Cadd(complex a, complex b);
complex Csub(complex a, complex b);
complex Cmul(complex a, complex b);
complex Complex(double re, double im);
complex Conjg(complex z);
complex Cdiv(complex a, complex b);
double Cabs(complex z);
double Carg(complex z);
complex Csqrt(complex z);
complex RCmul(double x, complex a);
