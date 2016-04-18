C*PUTOUT -- write line on standard output [Convex-UNIX]
C+
      SUBROUTINE PUTOUT(TEXT)
      CHARACTER*(*) TEXT
C
C This subroutine writes one line on the standard output; the text to be
C written is supplied as a character-string argument.
C
C Argument:
C  TEXT   (input)  : character string for output.
C
C Subroutines required:
C  Fortran formatted I/O.
C
C History:
C  1987 Nov 11 - TJP
C-----------------------------------------------------------------------
      WRITE (6,'(A)') TEXT
      END
