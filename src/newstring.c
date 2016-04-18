#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "newstring.h"

/*--------------------------------------------------------------------
 * Return a pointer to a dynamically allocated new string that is a
 * duplicate of the string pointed to by "s", or NULL on failure
 *--------------------------------------------------------------------*/

char *newstring(const char *s)
{
  char *p;
  if (s == NULL)
    return NULL;
  p = malloc(strlen(s)+1);
  if (p)
    strcpy(p, s);
  else
    fputs("Error in routine \"newstring\": no memory available\n", stderr);
  return p;
}
