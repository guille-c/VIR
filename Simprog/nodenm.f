C*NODENM -- return name of current host (node) [Convex-UNIX]
C+
      SUBROUTINE NODENM (NODE)
      CHARACTER*(*) NODE
C
C Return the name of the current host (node).
C
C Argument:
C  NODE   (output) : receives the node name, or string '<unknown>' if
C                    information is not available.
C
C Subroutines required:
C  HOSTNM (Convex)
C
C History:
C  1987 Nov 12 - TJP.
C-----------------------------------------------------------------------
      INTEGER HOSTNM
C
      IF (HOSTNM(NODE) .NE. 0) NODE = '<unknown>'
C
      END
