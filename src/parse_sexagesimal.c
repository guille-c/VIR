/*--------------------------------------------------------------------
 * Parse a sexagesimal string [+-]int:int:real
 *--------------------------------------------------------------------*/

#include <stdlib.h>
#include "parse_sexagesimal.h"

int parse_sexagesimal(char *s, double *result)
{
  char *end;
  double sign=1.0, t1, t2, t3;

  *result = 0.0;
  if (*s == '+') {
    s++;
  } else if (*s == '-') {
    s++; sign = -1.0;
  }
  t1 = (double) strtol(s, &end, 10);
  if (*end == '\0') {
    *result = sign*t1;
    return EXIT_SUCCESS;
  } else if (*end != ':') {
    return EXIT_FAILURE;
  }
  end++;
  t2 = (double) strtol(end, &end, 10);
  if (*end == '\0') {
    *result = sign*(t1 + t2/60.0);
    return EXIT_SUCCESS;
  } else if (*end != ':') {
    return EXIT_FAILURE;
  }
  end++;
  t3 = strtod(end, &end);
  if (*end != '\0') {
    return EXIT_FAILURE;
  }
  *result = sign*(t1 + t2/60.0 + t3/3600.0);
  return EXIT_SUCCESS;
}
