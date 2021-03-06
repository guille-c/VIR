      SUBROUTINE sla_UNPCD (DISCO,X,Y)
*+
*     - - - - - -
*      U N P C D
*     - - - - - -
*
*  Remove pincushion/barrel distortion from a distorted [x,y]
*  to give tangent-plane [x,y].
*
*  Given:
*     DISCO    d      pincushion/barrel distortion coefficient
*     X,Y      d      distorted coordinates
*
*  Returned:
*     X,Y      d      tangent-plane coordinates
*
*  Notes:
*
*  1)  The distortion is of the form RP = R*(1 + C*R**2), where R is
*      the radial distance from the tangent point, C is the DISCO
*      argument, and RP is the radial distance in the presence of
*      the distortion.
*
*  2)  For pincushion distortion, C is +ve;  for barrel distortion,
*      C is -ve.
*
*  3)  For X,Y in "radians" - units of one projection radius,
*      which in the case of a photograph is the focal length of
*      the camera - the following DISCO values apply:
*
*          Geometry          DISCO
*
*          astrograph         0.0
*          Schmidt           -0.3333
*          AAT PF doublet  +147.069
*          AAT PF triplet  +178.585
*          AAT f/8          +21.20
*          JKT f/8          +13.32
*
*  4)  The present routine is an approximate inverse to the
*      companion routine sla_PCD, obtained from two iterations
*      of Newton's method.  The mismatch between the sla_PCD and
*      sla_UNPCD routines is negligible for astrometric applications;
*      to reach 1 milliarcsec at the edge of the AAT triplet or
*      Schmidt field would require field diameters of 2.4 degrees
*      and 42 degrees respectively.
*
*  P.T.Wallace   Starlink   1 August 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      DOUBLE PRECISION DISCO,X,Y

      DOUBLE PRECISION CR2,A,CR2A2,F



      CR2=DISCO*(X*X+Y*Y)
      A=(2D0*CR2+1D0)/(3D0*CR2+1D0)
      CR2A2=CR2*A*A
      F=(2D0*CR2A2*A+1D0)/(3D0*CR2A2+1D0)
      X=X*F
      Y=Y*F

      END
