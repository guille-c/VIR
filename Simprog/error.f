C*ERROR -- print error message on stdout, and abort [Convex-UNIX]
C+
      SUBROUTINE ERROR(STRING)
      CHARACTER*(*) STRING
C
C This subroutine prints its argument (a character string) on the 
C standard output, prefixed with '++ERROR++', and then terminates 
C execution of the program, returning an error status to its parent
C (shell or command interpreter).
C
C Argument:
C  STRING (input)  : text of message
C
C Subroutines required:
C  EXIT
C  PUTOUT
C
C History:
C  1987 Apr 8 - TJP.
C  2000 Dec 6 - TJP: mods for g77.
C-----------------------------------------------------------------------
      CHARACTER*128 TXT
      INTEGER LEN1
C
      TXT = '+++ERROR+++ '//STRING
      CALL PUTOUT(TXT(1:LEN1(TXT)))
      CALL EXIT(1)
C
      END
