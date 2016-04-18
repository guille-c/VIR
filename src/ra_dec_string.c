#include <stdio.h>
#include "ra_dec_string.h"
#include "slalib.h"

/* --------------------------------------------------------------------
 * Convert RA (rad) to string; return pointer to (static) string
 * -------------------------------------------------------------------- */

char *ra_string(double ra)
{
  static char s[80];
  char sign;
  int dmsf[4];
  slaDr2tf(2, ra, &sign, dmsf);
  sprintf(s, "%s%2.2d:%2.2d:%2.2d.%2.2d", sign == '-' ? "-":"",
	  dmsf[0], dmsf[1], dmsf[2], dmsf[3]);
  return s;
}

/* --------------------------------------------------------------------
 * Convert Dec (rad) to string; return pointer to (static) string
 * -------------------------------------------------------------------- */

char *dec_string(double dec)
{
  static char s[80];
  char sign;
  int dmsf[4];
  slaDr2af(1, dec, &sign, dmsf);
  sprintf(s, "%c%2.2d:%2.2d:%2.2d.%1.1d", sign, dmsf[0], dmsf[1],
	   dmsf[2], dmsf[3]);
  return s;
}

