C*SSPACE -- remove excess spaces from character string
C+
      SUBROUTINE SSPACE (S, L)
      CHARACTER*(*) S
      INTEGER L
C
C Remove excess spaces from a character string: each string of
C consecutive spaces is replaced by a single space, and all leading
C spaces are removed.
C
C Standard Fortran 77, with all variables declared.
C
C Arguments:
C  S (in/out) : the character string to be modified.
C  L (output) : the length of the string, excluding trailing spaces.
C
C History:
C   1-Nov-1986: T. J. Pearson.
C-----------------------------------------------------------------------
      CHARACTER*1 LAST
      INTEGER I
C
      LAST = ' '
      L = 0
      I = 1
   10 IF (I.GT.LEN(S)) RETURN
      IF (LAST.EQ.' ' .AND. S(I:I).EQ.' ') THEN
          CONTINUE
      ELSE
          L = L+1
          S(L:L) = S(I:I)
          LAST = S(L:L)
      END IF
      I = I+1
      GOTO 10
C-----------------------------------------------------------------------
      END
