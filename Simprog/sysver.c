/*
C*SYSVER -- return name and version of operating system [Sun-UNIX]
C+
      SUBROUTINE SYSVER (STRING)
      CHARACTER*(*) STRING
C
C Return the name and version number of the operating system.
C
C Argument:
C  STRING  (output) : receives the information, or string '<unknown>' if
C                    it is not available.
C-----------------------------------------------------------------------
*/

#include <sys/utsname.h>
#include <stdio.h>
#include <string.h>

void sysver_(char string[], int len);

void sysver_(char string[], int len)
{
  struct utsname p;
  int i, junk, l;
  char t[128];
  junk = uname(&p);
  sprintf(t, "%s %s", p.sysname, p.release);
  l = strlen(t);
  if (l <= len)
    {
      for (i=0; i < l; i++)
	string[i] = t[i];
      for (i=l; i < len; i++)
	string[i] = ' ';
    }
  else
    {
      for (i=0; i < l; i++)
	string[i] = t[i];
    }
  return;
}
