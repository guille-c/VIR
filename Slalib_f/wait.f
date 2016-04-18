      SUBROUTINE sla_WAIT (DELAY)
*+
*     - - - - -
*      W A I T
*     - - - - -
*
*  Interval wait
*
*  !!! Version for Sun SPARC and DECstation !!!
*
*  Given:
*     DELAY     real      delay in seconds
*
*  Called:  SLEEP (an RTL subroutine in DEC and Sun Fortran)
*
*  P.T.Wallace   Starlink   16 November 1993
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      REAL DELAY


      CALL SLEEP(NINT(DELAY))

      END
