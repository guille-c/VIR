#include <math.h>
#define NRANSI
#include "nrutil.h"
#include "mnbrak.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	printf("mnbrak fa, a = %g\n", (*ax));
	*fa=(*func)(*ax);
	printf("mnbrak fb, b = %g\n", (*bx));
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	printf("mnbrak fc, c = %g\n", (*cx));
	*fc=(*func)(*cx);
	printf("mnbrak OK, a, b, c\n");
	//while (*fb >= *fc) { // Guille 17/05/06
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
		  printf("mnbrak fu 1, u = %g\n", u);
		  fu=(*func)(u);
		  if (fu < *fc) {
		    *ax=(*bx);
		    *bx=u;
		    *fa=(*fb);
		    *fb=fu;
		    return;
		  } else if (fu > *fb) {
		    *cx=u;
		    *fc=fu;
		    return;
		  }
		  u=(*cx)+GOLD*(*cx-*bx);
		  printf("mnbrak fu 2, u = %g\n", u);
		  fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
		  printf("mnbrak fu 3, u = %g\n", u);
		  fu=(*func)(u);
		  if (fu < *fc) {
		    SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
		    SHFT(*fb,*fc,fu,(*func)(u))
		  }
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
		  u=ulim;
		  printf("mnbrak fu 4, u = %g\n", u);
		  fu=(*func)(u);
		} else {
		  u=(*cx)+GOLD*(*cx-*bx);
		  printf("mnbrak fu 5, u = %g\n", u);
		  fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	      }
	printf("mnbrak OK\n");
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
