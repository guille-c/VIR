C-----------------------------------------------------------------------
C VLBINC:MAPCOM.INC -- include file for defining the common block for
C map header information. Corresponding FITS parameter names are given
C in parentheses.
C
C MPNDM       number of dimensions (1..7) (NAXIS)
C MPSIZ(7)    number of pixels in each dim (NAXISn)
C MPNIT       number of CLEAN iterations (NITER)
C
C MPOBJ       source name (OBJECT)
C MPTEL       telescope name, eg 'VLA', (TELESCOP)
C MPINS       instrument name, eg 'VLBI' (INSTRUME)
C MPOBS       observer name (OBSERVER)
C MPDOB       date of observation (DATE-OBS)
C MPDMP       date of creation of map (DATE-MAP)
C MPBUN       map units, eg 'Jy/beam' (BUNIT)
C MPCTY(7)    coordinate, eg 'RA', 'Dec' (CTYPEn)
C MPTYP       type of map (--)
C MPCRT       creating program (--)
C
C MPRV(7)     ref point value (CRVALn)
C MPRP(7)     ref point pixel locn (CRPIXn)
C MPIN(7)     coord value increment (CDELTn)
C MPRO(7)     rotation angle (CROTAn)
C MPBMJ       beam major axis (BMAJ)
C MPBMN       beam minor axis (BMIN)
C MPBPA       beam position angle (BPA)
C MPEPO       epoch of coordinates (EPOCH)
C MPFRQ       frequency, MHz (--)
C MPBWD       bandwidth, MHz (--)
C-----------------------------------------------------------------------
      INTEGER MPNDM, MPSIZ(7), MPNIT
C
      CHARACTER*32 MPTEL
      CHARACTER*16 MPOBJ, MPINS, MPOBS, MPDOB, MPDMP, MPBUN
      CHARACTER*16 MPCTY(7), MPTYP, MPCRT
C
      DOUBLE PRECISION MPRV(7), MPRP(7), MPIN(7), MPRO(7)
      DOUBLE PRECISION MPBMJ, MPBMN, MPBPA
      DOUBLE PRECISION MPEPO, MPFRQ, MPBWD
C
      COMMON /MAP_IHD/ MPNDM,MPSIZ,MPNIT
      COMMON /MAP_CHD/ MPOBJ,MPTEL,MPINS,MPOBS,MPDOB,MPDMP,MPBUN,
     1                 MPCTY,MPTYP,MPCRT
      COMMON /MAP_RHD/ MPRV, MPRP, MPIN, MPRO,
     2                 MPBMJ,MPBMN,MPBPA,MPEPO,MPFRQ,MPBWD
C-----------------------------------------------------------------------
