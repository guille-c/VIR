      SUBROUTINE sla_MAP (RM, DM, PR, PD, PX, RV, EQ, DATE, RA, DA)
*+
*     - - - -
*      M A P
*     - - - -
*
*  Transform star RA,Dec from mean place to geocentric apparent
*
*  The reference frames and timescales used are post IAU 1976.
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Given:
*     RM,DM    dp     mean RA,Dec (rad)
*     PR,PD    dp     proper motions:  RA,Dec changes per Julian year
*     PX       dp     parallax (arcsec)
*     RV       dp     radial velocity (km/sec, +ve if receding)
*     EQ       dp     epoch and equinox of star data (Julian)
*     DATE     dp     TDB for apparent place (JD-2400000.5)
*
*  Returned:
*     RA,DA    dp     apparent RA,Dec (rad)
*
*  Called:
*     sla_MAPPA       star-independent parameters
*     sla_MAPQK       quick mean to apparent
*
*  Notes:
*
*  1)  EQ is the Julian epoch specifying both the reference frame and
*      the epoch of the position - usually 2000.  For positions where
*      the epoch and equinox are different, use the routine sla_PM to
*      apply proper motion corrections before using this routine.
*
*  2)  The distinction between the required TDB and TT is always
*      negligible.  Moreover, for all but the most critical
*      applications UTC is adequate.
*
*  3)  The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt.
*
*  4)  This routine may be wasteful for some applications because it
*      recomputes the Earth position/velocity and the precession-
*      nutation matrix each time, and because it allows for parallax
*      and proper motion.  Where multiple transformations are to be
*      carried out for one epoch, a faster method is to call the
*      sla_MAPPA routine once and then either the sla_MAPQK routine
*      (which includes parallax and proper motion) or sla_MAPQKZ (which
*      assumes zero parallax and proper motion).
*
*  5)  The accuracy is limited by imperfections in the IAU 1976/1980
*      models for precession and nutation.  Corrections are tabulated
*      in IERS Bulletin B and at the present epoch are of order 50 mas.
*      An improved precession-nutation model can be introduced by
*      using sla_MAPPA and sla_MAPQK (see the previous note) and
*      replacing the precession-nutation matrix into the parameter
*      array directly.
*
*  6)  The accuracy is further limited by the routine sla_EVP, called
*      by sla_MAPPA, which computes the Earth position and velocity
*      using the methods of Stumpff.  The maximum error is about
*      0.3 mas.
*
*  P.T.Wallace   Starlink   8 May 2000
*
*  Copyright (C) 2000 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      DOUBLE PRECISION RM,DM,PR,PD,PX,RV,EQ,DATE,RA,DA

      DOUBLE PRECISION AMPRMS(21)



*  Star-independent parameters
      CALL sla_MAPPA(EQ,DATE,AMPRMS)

*  Mean to apparent
      CALL sla_MAPQK(RM,DM,PR,PD,PX,RV,AMPRMS,RA,DA)

      END
