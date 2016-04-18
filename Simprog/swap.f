      SUBROUTINE SWAP (A,B)
C-----------------------------------------------------------------------
C SWAP: interchange values of two quantities.
C
C Arguments:
C  A,B (real, input/output): variables whose values are to be swapped.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL A,B,C
C-----------------------------------------------------------------------
      C = A
      A = B
      B = C
      RETURN
      END
