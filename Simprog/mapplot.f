      PROGRAM MAPPLT
C-----------------------------------------------------------------------
C MAPPLOT: contour a FITS-format image (on disk).
C
C Author:            T.J. Pearson (TJP)
C Language:          Fortran-77
C
C Modifications:
C Version 1.0:  ??????????? - (TJP).
C Version 1.1:  1983 May 24 - add FONT, LINEWIDTH (TJP).
C Version 2.0:  1984 Jan 23 - change map format to FITS, and remove
C                             special Grinnell plotting (TJP).
C Version 2.1:  1984 Jan 30 - handle default contours better; label
C                             MEM maps correctly; add LTYPE and BEAM
C                             parameters.
C Version 2.2:  1984 Jun 28 - fix to recognize AIPS-like coordinates
C                             RA---SIN and DEC--SIN (TJP).
C Version 2.3:  1985 Jun 26 - add more info to legend (TJP).
C Version 2.4:  1985 Jul 18 - add CHARSIZE parameter (TJP).
C Version 2.5:  1985 Nov 24 - add TOPLEFT and NOTITLE parameters (TJP).
C Version 2.6:  1986 Jul  7 - add BPOS parameter (TJP).
C Version 2.7:  1986 Jul 28 - add annotation (TJP).
C Version 2.8:  1986 Nov  1 - more annotation: TOPRIGHT, ANNOT, CONTLW,
C                             NPXY (Tasso Tzioumis).
C Version 2.9:  1987 Apr  4 - improved annotation of grey scales; add
C                             COLORS (TJP).
C Version 3.0:  1988 Apr 15 - changes for Convex (TJP).
C Version 3.1:  1988 Jul 20 - add additional labels, enable LTYPE=6
C                             (TJP).
C Version 3.2:  1989 Jun 30 - correct bug (beam annotation wrong if
C                             beam not plotted) (TJP).
C Version 3.3:  1989 Aug 23 - increase maximum array size (TJP).
C Version 3.4:  1989 Sep 10 - add WEDGE parameter (TJP).
C Version 3.5:  1990 Jan 10 - remove FORMAT=MEM option (MEM maps are
C                             now in FITS format).
C Version 3.6:  1991 Mar 25 - better placement of beam, change
C                             annotation (TJP).
C Version 3.7:  1991 Apr  1 - correct minor bug (TJP).
C Version 3.8:  1991 Apr 17 - handle rotated maps (TJP).
C Version 4.0:  1991 May 14 - add windows (TJP).
C Version 4.1:  1991 Dec 23 - add surface option; add more ANNOT,
C                             CHARSIZE, COLOR options (TJP).
C Version 4.2:  1992 Jul  6 - add FORMAT=INV option (Edward King); add
C                             statistics (TJP).
C Version 4.3:  1993 Mar 18 - add more LABELx parameters (TJP).
C Version 4.4:  1993 Apr 27 - allow longer telescope names (TJP).
C Version 4.5:  1993 May 14 - allow user to format legend (TJP).
C Version 4.6:  1993 Sep 20 - add linear scale (TJP).
C Version 4.7:  1993 Oct 16 - add OPAQUE option (TJP).
C Version 4.8:  1993 Nov  8 - add OVERLAY, SHIFT (TJP).
C Version 4.9:  1994 Feb  9 - add BOTLEFT, BOTRIGHT (TJP).
C Version 4.10: 1994 Mar 14 - handle negative maps better (TJP).
C Version 5.0:  1994 Aug  8 - add UNIT parameter to force scaling (TJP).
C Version 5.1:  1994 Nov 10 - fix bug (failing to call PGBEG) (TJP).
C Version 5.2:  1996 Feb  8 - label bar in kpc if necessary (TJP).
C Version 5.3:  1996 Dec 20 - change wedge display (TJP).
C Version 6.0:  2000 Sep 06 - Linux/PC port, add GRAYMAX (RLM), SLA int.
C Version 7.0:  2000 DEc 07 - use CFITSIO for input (TJP).
C-----------------------------------------------------------------------
C
C Constants
C
      CHARACTER*(*) VERSN
      INTEGER    INC,OUTC,NPARS,IDIM,MAXCON,NLAB
      REAL       VELC
      PARAMETER  (VERSN='7.0 - 2000 Dec 7')
C     -- array size
      PARAMETER  (IDIM=2048*2048)
C     -- unit numbers
      PARAMETER  (INC=5,OUTC=6)
C     -- number of parameters
      PARAMETER  (NPARS=303)
C     -- maximum number of contours
      PARAMETER  (MAXCON=20)
C     -- maximum number of labels
      PARAMETER  (NLAB=8)
C     -- speed of light (m/s)
      PARAMETER  (VELC=2.997924580E8)
C
C Variables.
C
      CHARACTER*128  MAPFIL,MAPOLD,PLTFIL,PLTOLD,TOPLEF,TOPRIG, TEXT
      CHARACTER*128  BOTLEF, BOTRIG
      CHARACTER*255  MAPNAM
      CHARACTER*40   XLABEL, YLABEL, XUNIT, YUNIT, FREQS, WAVES, EPOCHS
      CHARACTER*16   UNIT
      CHARACTER*128  TXTERR, TITLE1, TITLE2, STR
      CHARACTER*8    PSTYLE, MAPFMT, TEST
      CHARACTER*5    Z
      CHARACTER*48   LABELS(NLAB)
      CHARACTER*17   REVTIM
      DOUBLE PRECISION NAMES(NPARS),PARS(NPARS),ENDMRK
      DOUBLE PRECISION DJM, EPJ
      REAL    XLPOS(NLAB), YLPOS(NLAB)
      REAL    MAPARR(IDIM), MAPMAX, MAPMIN, ABSMAX
      REAL    TR(6),VCN(20),CLEV,PLEV,LEVS(20), CHSIZ1, CHSIZ2, BPOS(2)
      REAL    BLACK, BMAJ, BMIN, BPA, CTH, DLEV, DY, D
      REAL    PH, STH, VPX1, VPX2, VPY1, VPY2, WHITE, X1, X2
      REAL    XBC, XFACT, XM, Y1, Y2, YBC, YFACT, YM, Y
      REAL    LRTB(80), ANGLE, REDSHI, H0, Q0, SCL, DL, BAR
      REAL    XSHIFT, YSHIFT, XMARK, YMARK
      REAL    PGRND, GRAYMAX, GRAY1, GRAY2
      INTEGER ICON, IER, I, J, K, LENERR, LTYPE, L, MMAX, JUNK, IBAR
      INTEGER MMIN, MODE, MSHI, MSLO, M, NCONS, NMAX, NMIN
      INTEGER NSHI, NSLO, N, ID, IM, IY, LX, LY
      INTEGER LINEW, FONT, NLEVS, CONLW
      INTEGER NPX, NPY, NPXOLD, NPYOLD, COLOR(6)
      INTEGER LTEST, LL
      INTEGER LEN1, PGBEG, ANNOT, CTAB, NMARK, CMARK
      LOGICAL DOGRAY, DOCNTR, DOTITL, DOTICK, DOPIXL, DOBEAM
      LOGICAL WEDGE, SETW, DOSURF, OPAQUE, OVERLA
C
C Common blocks.
C
      INCLUDE 'mapcom.inc'
C
C Input parameters.
C
      DATA    ENDMRK/ 1H/ /
      DATA    NAMES / 5HSTyle, 6HFORmat,
     1                7HMAPfile, 9*1H ,
     2                20*4HLEvs, 5HBLack,5HWHite,4HPLEv,4HCLev,
     3                2*6HXrange, 2*6HYrange,
     4                8HPLOtfile, 9*1H , 4HEXit, 8HLInewidt, 4HFONt,
     5                1H*, 9*1H , 5HLType, 4HBEam, 5HANGle,
     6                7HTOPLeft,9*1H , 7HNOTItle, 2*8HBPositio,
     7                8HTOPRight, 9*1H , 5HANNot, 6HCONTLW, 2*4HNPxy,
     8                5HWEdge, 6*6HCOLors,
     9                3*6HLABEL1, 5*1H , 3*6HLABEL2, 5*1H ,
     A                3*6HLABEL3, 5*1H , 3*6HLABEL4, 5*1H ,
     B                3*6HLABEL5, 5*1H , 3*6HLABEL6, 5*1H ,
     C                3*6HLABEL7, 5*1H , 3*6HLABEL8, 5*1H ,
     D                8HSETWindo, 80*4HLRTB, 2*8HCHarsize,
     E                5HTItle, 9*1H , 5HTItle, 9*1H ,
     F                1HZ, 2HH0, 2HQ0, 8HBARLENGT, 6HOPaque, 7HOVerlay,
     G                2*5HSHift, 7HBOTLeft, 9*1H , 8HBOTRight, 9*1H ,
     H                5HUnits, 1H , 7HPalette, 4*4HMark, 7HGraymax /
      DATA    PARS  / 1HC, 4HFITS,
     1                3HMAP,9*1H ,
     2                20*0D0, 1.0D0, 0.0D0, 5D0, 0D0,
     3                2*0D0, 2*0D0 ,
     4                1H?, 9*1H ,-1D0 , 1D0, 1D0,
     5                10*1H , 4D0, -1D0, 25D0, 
     6                10*1H , -1D0, 2*0D0,
     7                10*1H , 1D0, 1D0, 2*1D0,
     8                -1D0, 9D0, 1D0, 2D0, 1D0, 5D0, 2D0,
     9                2*0D0, 6*1H , 2*0D0, 6*1H ,
     A                2*0D0, 6*1H , 2*0D0, 6*1H ,
     B                2*0D0, 6*1H , 2*0D0, 6*1H ,
     C                2*0D0, 6*1H , 2*0D0, 6*1H ,
     D                1D0, 80*0D0, 1D0, 0.6D0,
     E                10*1H , 7H%s  %t  ,6H%f  %d, 8*1H ,
     F                0D0, 100D0, 0.5D0, 0D0, -1D0, -1D0,
     G                2*0D0, 10*1H , 10*1H ,
     H                2*1H , 1D0, 3*0D0, 1D0, 0D0 /
C
C Statement function (from SLALIB)
C
        EPJ(DJM) = 2000D0 + (DJM-51544.5D0)/365.25D0

C
C Format statements.
C
  600 FORMAT (' MAPPLOT makes contour plots of FITS-format maps'/
     1        ' (Version ',A,')'//
     2        ' Control parameters (end with /):')
  660 FORMAT (/1X,'File: ',A)
  662 FORMAT (5X,'Source: ',A/
     2        5X,'Type: ',A,'   Created by: ',A/
     3        5X,'Telescope: ',2A/
     4        5X,'Observer: ',A/
     5        5X,'Size in pixels (x,y): ',I5,',',I5/
     6        5X,'Center at ',F7.1,',',F7.1/
     7        5X,'Maximum ',1PE12.3,'  in cell',2I5/
     8        5X,'Minimum ',1PE12.3,'  in cell',2I5)
  680 FORMAT (/' Contour levels (',A,'):'/(1X,6F12.6))
  681 FORMAT ( ' Contour levels (percent of peak):'/(1X,6F12.2))
  700 FORMAT ( ' +++ERROR+++ cannot send plot to ',A)
 709  FORMAT ( ' +++ERROR+++ UNIT must be "mas", "milliarcsec", ',
     :         '"arcsec", "arcmin", "degree",'/' "rad", or blank')
  710 FORMAT ( ' +++ERROR+++ STYLE must be "C", "G", "CG" or "S"')
  720 FORMAT ( ' +++ERROR+++ BLACK=WHITE is not allowed')
  730 FORMAT ( 'Maximum: ',1PG12.4,' ',A)
  740 FORMAT ( 'Contours (%): ',10F7.2)
  745 FORMAT ( 'Rotated by: ',F7.1,'\(718)')
  750 FORMAT ( 'Beam:  FWHM ',F7.2,' \\x',F7.2,1X,A,', p.a.',F7.1,
     1         '\(718)')
  760 FORMAT ( 'Beam:  FWHM ',F7.2,1X,A)
  770 FORMAT ( 'Grey scale: ',2F7.2)
  780 FORMAT (1X,'z = ',F6.4,',  H0 = ',F6.1,' km/s/Mpc,  q0 = ',F4.1)
  790 FORMAT (1X,'Linear scale = ',F12.2,' kpc/arcsec=pc/milliarcsec')
  791 FORMAT (1X,'Bar length = ',I6,' pc, =',F6.1,' mas')
 792  FORMAT (I4,' pc')
 793  FORMAT (I4,' kpc')
 794  FORMAT (1X,'Axis units are ',A,', ',A)
C-----------------------------------------------------------------------
C
C Start.
C
      WRITE (OUTC,600) VERSN
      MAPOLD = '*** UNSPECIFIED ***'
      PLTOLD = '*** UNSPECIFIED ***'
      NPXOLD = 1
      NPYOLD = 1
C
C Read control parameters.
C
 1000 MODE = 0
      PARS(165) = -1D0
      PARS(273) = -1D0
      CALL KEYIN(NAMES,PARS,NPARS,ENDMRK,MODE,INC,OUTC)
C     -- test for end-of-file or EXIT
      IF (MODE.EQ.1) GOTO 2000
      IF (PARS(51).NE.-1D0) GOTO 2000
C
C     -- MAPFILE
      WRITE (MAPFIL,'(10A8)') (PARS(I),I=3,12)
C     -- STYLE
      WRITE (PSTYLE,'(A8)') PARS(1)
      CALL UPCASE(PSTYLE)
      DOGRAY = INDEX(PSTYLE,'G').NE.0
      DOCNTR = INDEX(PSTYLE,'C').NE.0
      DOSURF = PSTYLE(1:1).EQ.'S'
C     -- FORMAT
      WRITE (MAPFMT,'(A8)') PARS(2)
      CALL UPCASE(MAPFMT)
C     -- BLACK and WHITE
      BLACK = PARS(33)
      WHITE = PARS(34)
      GRAYMAX = PARS(303)
C     -- LEVS
      DO 1010 ICON=1,20
          LEVS(ICON) = PARS(12+ICON)
 1010 CONTINUE
      NLEVS = 20
      DO WHILE (NLEVS.GT.0 .AND. LEVS(NLEVS).EQ.0.0)
          NLEVS = NLEVS-1
      END DO
      IF (NLEVS.EQ.0) THEN
C         -- Default contours +/1 1,2,4,8,...,1024
          NLEVS = 20
          J = 1
          DO 1030 ICON=1,10
            LEVS(10+ICON) = J
            LEVS(11-ICON) = -J
            J = J*2
 1030     CONTINUE
      END IF
C     -- PLEV and CLEV
      PLEV = PARS(35)
      CLEV = PARS(36)
C     -- XRANGE and YRANGE
      MSLO = PARS(37)
      MSHI = PARS(38)
      NSLO = PARS(39)
      NSHI = PARS(40)
      WRITE (PLTFIL,'(10A8)') (PARS(I),I=41,50)
C     -- LINEWIDTH
      LINEW = PARS(52)
C     -- FONT
      FONT  = PARS(53)
C     -- TITLE
      WRITE (TITLE1,'(10A8)') (PARS(J),J=248,257)
      WRITE (TITLE2,'(10A8)') (PARS(J),J=258,267)
C     -- LTYPE
      LTYPE = PARS(64)
      DOTITL = PARS(77).EQ.-1D0
      DOTICK = LTYPE.GT.2
      DOPIXL = LTYPE.EQ.6
C     -- BEAM
      DOBEAM = PARS(65).NE.-1D0
C     -- Angle
      ANGLE = PARS(66)
      WRITE (TOPLEF,'(10A8)') (PARS(I),I=67,76)
      BPOS(1) = PARS(78)
      BPOS(2) = PARS(79)
      WRITE (TOPRIG,'(10A8)') (PARS(I),I=80,89)
      WRITE (BOTLEF,'(10A8)') (PARS(I),I=276,285)
      WRITE (BOTRIG,'(10A8)') (PARS(I),I=286,295)
      ANNOT = NINT(PARS(90))
      CONLW = PARS(91)
      NPX   = PARS(92)
      NPY   = PARS(93)
C     -- WEDGE
      WEDGE = PARS(94).NE.-1D0
C     -- Color indices for annotation/+contours/-contours/grayscale
      COLOR(1) = PARS(95)
      COLOR(2) = PARS(96)
      COLOR(3) = PARS(97)
      COLOR(4) = PARS(98)
      COLOR(5) = PARS(99)
      COLOR(6) = PARS(100)
C     -- Additional labels
      DO 1040 I=1,NLAB
          XLPOS(I) = PARS(101 + (I-1)*8)
          YLPOS(I) = PARS(102 + (I-1)*8)
          WRITE(LABELS(I), '(6A8)') (PARS(J), J=103+(I-1)*8, 
     1      108+(I-1)*8)
 1040 CONTINUE
C     -- Windows
      SETW = PARS(165).GE.0D0
      DO 1050 I=1,80
          LRTB(I) = PARS(165+I)
 1050 CONTINUE
C     -- CHARSIZE
      CHSIZ1 = PARS(246)
      CHSIZ2 = PARS(247)
C     -- Cosmology
      REDSHI = PARS(268)
      WRITE (Z,'(F5.3)',IOSTAT=JUNK) REDSHI
      IF (REDSHI.LE.0.0) Z = ' '
      H0     = PARS(269)
      Q0     = PARS(270)
      BAR    = PARS(271)
      OPAQUE = PARS(272).GE.0D0
      OVERLA = PARS(273).GE.0D0
C     -- Shift: input in mas; applied in degrees
      XSHIFT = PARS(274)/3600000.0
      YSHIFT = PARS(275)/3600000.0
C     -- Angular units
      WRITE (UNIT,'(2A8)') (PARS(I),I=296,297)
      CALL UPCASE(UNIT)
      IF (UNIT.NE.'MAS'    .AND. UNIT.NE.'MILLIARCSEC' .AND.
     :    UNIT.NE.'ARCSEC' .AND. UNIT.NE.'ARCMIN'      .AND.
     :    UNIT.NE.'DEGREE' .AND. UNIT.NE.'RAD'         .AND.
     :    UNIT.NE.' '      ) THEN
          WRITE (OUTC,709)
          GOTO 1000
      END IF
C     -- Color palette
      CTAB = NINT(PARS(298))
      CTAB = MIN(6, MAX(1,CTAB))
C     -- Marker
      XMARK = PARS(299)
      YMARK = PARS(300)
      NMARK = PARS(301)
      CMARK = PARS(302)
C
C Check something is to be plotted.
C
      IF (.NOT.(DOCNTR.OR.DOGRAY.OR.DOSURF)) THEN
          WRITE (OUTC,710)
          GOTO 1000
      END IF
      IF (DOGRAY .AND. (BLACK.EQ.WHITE)) THEN
          WRITE (OUTC,720)
          GOTO 1000
      END IF
C
C Open plot file (if necessary).
C
      IF (PLTFIL.NE.PLTOLD .OR. NPX.NE.NPXOLD .OR. NPY.NE.NPYOLD) THEN
          IER = PGBEG(0,PLTFIL,NPX,NPY)
          IF (IER.NE.1) THEN
              WRITE (OUTC,700) PLTFIL(1:LEN1(PLTFIL))
              PLTOLD = '*** ILLEGAL ***'
              GOTO 1000
          END IF
          CALL PGASK(.FALSE.)
          PLTOLD = PLTFIL
          NPXOLD = NPX
          NPYOLD = NPY
          IF (OVERLA)
     :        WRITE (OUTC,'(A)') ' Starting new plot: OVERLAY ignored'
          OVERLA = .FALSE.
      END IF
C
C Read map from file MAPFIL (if necessary), locate maximum and minimum
C and type header.
C
      IF (MAPFIL.NE.MAPOLD) THEN
          MAPOLD = MAPFIL
          IER = 0
          CALL FTOPEN(2, MAPFIL, 0, JUNK, IER)
          IF (IER.NE.0) GOTO 3000
          MAPNAM = MAPFIL
          WRITE (OUTC,660) MAPNAM(1:LEN1(MAPNAM))
          CALL FITRD(2, MAPARR, IDIM, IER, TXTERR)
          IF (IER.NE.0) THEN
              LENERR = LEN1(TXTERR)
              CALL PUTOUT(TXTERR(1:LENERR))
              GOTO 3000
          END IF
          CLOSE (UNIT=2)
C
C Locate maximum, minimum.
C
          MAPMAX = -1E37
          MAPMIN =  1E37
          DO 1110 N=1,MPSIZ(2)
              K = (N-1)*MPSIZ(1)
              DO 1100 M=1,MPSIZ(1)
                  IF (MAPARR(K+M).GT.MAPMAX) THEN
                      MAPMAX = MAPARR(K+M)
                      NMAX = N
                      MMAX = M
                  END IF
                  IF (MAPARR(K+M).LT.MAPMIN) THEN
                      MAPMIN = MAPARR(K+M)
                      NMIN = N
                      MMIN = M
                  END IF
 1100         CONTINUE
 1110     CONTINUE
          WRITE (OUTC,662) MPOBJ,MPTYP,MPCRT,
     1                  MPTEL,MPINS,MPOBS,MPSIZ(1),MPSIZ(2),MPRP(1),
     2                  MPRP(2),MAPMAX,MMAX,NMAX, MAPMIN,MMIN,NMIN
          ABSMAX = MAPMAX
          IF (ABS(MAPMIN).GT.ABS(MAPMAX)) ABSMAX = MAPMIN
C         -- Format frequency string
          IF (MPFRQ.GT.1E3) THEN
              WRITE (FREQS,'(F8.3,'' GHz'')') MPFRQ*1E-3
          ELSE IF (MPFRQ.GT.0.0) THEN
              WRITE (FREQS,'(F8.3,'' MHz'')') MPFRQ
          ELSE
              FREQS = '*'
          END IF
          DO WHILE (FREQS(1:1).EQ.' ')
              FREQS = FREQS(2:)
          END DO
C         -- Format wavelength string
          IF (MPFRQ.LE.0.0) THEN
              WAVES = '*'
          ELSE IF (MPFRQ.GT.3E4) THEN
              WRITE (WAVES,'(F8.2,'' mm'')') 1E-3*VELC/MPFRQ
          ELSE
              WRITE (WAVES,'(F8.2,'' cm'')') 1E-4*VELC/MPFRQ
          END IF
          DO WHILE (WAVES(1:1).EQ.' ')
              WAVES = WAVES(2:)
          END DO
C         -- Format epoch string
          EPOCHS = MPDOB
C         -- Discern whether to decode old or new data format
          if (mpdob(3:3) .eq. '/') then
            READ (MPDOB, '(I2,1X,I2,1X,I2)', ERR=877) ID, IM, IY
          else
            read (MPDOB, '(2X,I2,1X,I2,1X,I2)', ERR=877) IY, IM, ID
          end if
          CALL CALDJ(IY, IM, ID, DJM, IER)
          WRITE (EPOCHS,'(F7.2)') EPJ(DJM)
  877     CONTINUE
      END IF
C
C Compute contour levels from PLEV, CLEV, LEVS and ABSMAX.
C Note that percentage contour levels are based on the
C maximum of the whole map, not just the region displayed.
C Type out contour levels in percent and in the units of
C the image. NCONS is the number of useful contour levels.
C
      IF (DOCNTR) THEN
          IF (PLEV.EQ.0.0) THEN
            DLEV = CLEV
          ELSE
            DLEV = PLEV*ABSMAX/100.0
          END IF
          NCONS = 0
          DO 1200 ICON=1,NLEVS
            IF (LEVS(ICON).NE.1E37) THEN
                NCONS = NCONS+1
                VCN(NCONS) = DLEV*LEVS(ICON)
                IF (VCN(NCONS).GT.MAPMAX .OR. VCN(NCONS).LT.MAPMIN)
     1                  NCONS = NCONS-1
            END IF
 1200     CONTINUE
          IF (NCONS.EQ.0) THEN
            WRITE (OUTC,'(A)')
     1          ' +++ERROR+++ requested contours outside range of map.'
            GOTO 1000
          END IF
          LL = MAX(1,LEN1(MPBUN))
          WRITE (OUTC,680) MPBUN(1:LL),
     1                  (VCN(ICON),ICON=1,NCONS)
          WRITE (OUTC,681) (VCN(ICON)*100.0/ABSMAX,ICON=1,NCONS)
      END IF
C
C Determine the ranges of the axes. The axis ranges are:
C MSLO to MSHI (horizontal) and NSLO to NSHI (vertical) in
C pixel numbers, and X1 to X2 (left to right) and Y1 to Y2
C (bottom to top) in real units.  We arrange for the 
C "real units" of the x-axis to increase from left to 
C right (except for RA-like quantities, which increase
C from right to left) and for the "real units" of the
C y-axis to increase from bottom to top.
C
      IF (MSLO.EQ.0 .AND. MSHI.EQ.0) THEN
          MSLO = 1
          MSHI = MPSIZ(1)
      ELSE IF (MSLO.LT.1 .OR. MSLO.GE.MSHI .OR. MSHI.GT.MPSIZ(1)) THEN
          WRITE (OUTC,'(A,I5)') ' +++ERROR+++ Requested XRANGE is '//
     1            'outside map pixel range of 1,',MPSIZ(1)
          GOTO 1000
      END IF
      IF (NSLO.EQ.0 .AND. NSHI.EQ.0) THEN
          NSLO = 1
          NSHI = MPSIZ(2)
      ELSE IF (NSLO.LT.1 .OR. NSLO.GE.NSHI .OR. NSHI.GT.MPSIZ(2)) THEN
          WRITE (OUTC,'(A,I5)') ' +++ERROR+++ Requested YRANGE is '//
     1            'outside map pixel range of 1,',MPSIZ(2)
          GOTO 1000
      END IF
      X1 = (MSLO-MPRP(1))*MPIN(1)
      X2 = (MSHI-MPRP(1))*MPIN(1)
      Y1 = (NSHI-MPRP(2))*MPIN(2)
      Y2 = (NSLO-MPRP(2))*MPIN(2)
      IF (MPCTY(1).EQ.'RA' .OR. MPCTY(1)(1:3).EQ.'RA-' .OR.
     1          MPCTY(1).EQ.'LL') THEN
          XLABEL = 'Relative R.A.'
          IF (MPRO(2).NE.0.0) XLABEL = 'x'
          CALL ANGLUN(UNIT,ABS(X2-X1),XUNIT,XFACT)
          MPRV(1) = XSHIFT
          IF (X1.LT.X2) CALL SWAP(X1,X2)
      ELSE
          XLABEL = MPCTY(1)
          XUNIT = ' '
          XFACT = 1.0
          IF (X1.GT.X2) CALL SWAP(X1,X2)
      END IF
      IF (MPCTY(2).EQ.'DEC' .OR. MPCTY(2)(1:4).EQ.'DEC-' .OR.
     1          MPCTY(2).EQ.'MM') THEN
          YLABEL = 'Relative Decl.'
          IF (MPRO(2).NE.0.0) YLABEL = 'y'
          CALL ANGLUN(UNIT,ABS(Y2-Y1),YUNIT,YFACT)
          MPRV(2) = YSHIFT
      ELSE
          YLABEL = MPCTY(2)
          YUNIT = ' '
          YFACT = 1.0
      END IF
      IF (Y1.GT.Y2) CALL SWAP(Y1,Y2)
      X1 = (X1+MPRV(1))*XFACT
      X2 = (X2+MPRV(1))*XFACT
      Y1 = (Y1+MPRV(2))*YFACT
      Y2 = (Y2+MPRV(2))*YFACT
C
C Determine the transformation between pixel numbers and real units.
C
      TR(1) = (MPRV(1)-MPRP(1)*MPIN(1))*XFACT
      TR(2) = MPIN(1)*XFACT
      TR(3) = 0.0
      TR(4) = (MPRV(2)-MPRP(2)*MPIN(2))*YFACT
      TR(5) = 0.0
      TR(6) = MPIN(2)*YFACT
C
C If LTYPE=6 requested, rescale everything in pixels.
C
      IF (DOPIXL) THEN
          X1 = MSLO
          X2 = MSHI
          Y1 = NSLO
          Y2 = NSHI
          XLABEL = MPCTY(1)
          XUNIT = 'pixels'
          YLABEL = MPCTY(2)
          YUNIT = XUNIT
          TR(1) = 0.0
          TR(2) = 1.0
          TR(3) = 0.0
          TR(4) = 0.0
          TR(5) = 0.0
          TR(6) = 1.0
      END IF
C
C Set window and viewport. The viewport is set to allow enough room
C for annotation. (This step omitted for an overlay plot.)
C
      IF (.NOT.OVERLA) THEN
          CALL PGPAGE
          CALL PGERAS
          CALL PGSVP(0.0, 1.0, 0.0, 1.0)
          CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
          D = MIN(VPX2-VPX1, VPY2-VPY1)/40.0
          VPX1 = VPX1 + 5.0*D
          VPX2 = VPX2 - 2.0*D
          VPY1 = VPY1 + 12.0*D
          VPY2 = VPY2 - 4.0*D
          CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
          IF (DOSURF) THEN
              CALL PGWNAD(0.0,1.0,0.0,1.0)
          ELSE
              CALL PGWNAD(X1,X2,Y1,Y2)
          END IF
      END IF
C     -- Omit the following statement if your version of PGPLOT
C     -- does not include this routine
      IF (OPAQUE) CALL PGSTBG(0)
C
C Draw a gray-scale plot if requested.
C
      IF (DOGRAY) THEN
          CALL PGBBUF
          CALL PGSLW(1)
          CALL PGSCI(COLOR(4))
          CALL PALETT(CTAB)
          IF (GRAYMAX .EQ. 0D0) THEN
                GRAY1 = BLACK*ABSMAX
                GRAY2 = WHITE*ABSMAX
          ELSE
                GRAY1 = BLACK*GRAYMAX
                GRAY2 = WHITE*GRAYMAX
          END IF
          CALL PGIMAG(MAPARR,MPSIZ(1),MPSIZ(2),MSLO,MSHI,NSLO,NSHI,
     1                GRAY1,GRAY2,TR)
          CALL PGEBUF
      END IF
C
C Draw a contour plot if requested.
C
      IF (DOCNTR) THEN
          CALL PGBBUF
          CALL PGSLW(CONLW)
          DO 1210 ICON=1,NCONS
            IF (VCN(ICON).LT.0.0) THEN
                CALL PGSCI(COLOR(3))
                CALL PGSLS(4)
            ELSE
                CALL PGSCI(COLOR(2))
                CALL PGSLS(1)
            END IF
              CALL PGCONT(MAPARR,MPSIZ(1),MPSIZ(2),MSLO,MSHI,NSLO,NSHI,
     1                      VCN(ICON),-1,TR)
 1210     CONTINUE
          CALL PGSLW(LINEW)
          CALL PGSLS(1)
          CALL PGEBUF
      END IF
C
C Draw surface plot if requested.
C
      IF (DOSURF) THEN
          CALL SURF(MAPARR,MPSIZ(1),MPSIZ(2),MSLO,MSHI,NSLO,NSHI,
     :              ANGLE, COLOR(5), COLOR(6), 14)
      END IF
C
C Draw windows.
C
      IF (DOGRAY.OR.DOCNTR) CALL MPLWIN(LRTB, 80, 2)
C
C Draw the beam FWHM contour if requested.
C
      IF (DOBEAM.AND.MPBMJ.GT.0.0.AND.XFACT.EQ.YFACT) THEN
          BMAJ = MPBMJ/2.0
          BMIN = MPBMN/2.0
          BPA = (90.-MPBPA)*3.1415926/180.0
C         -- Find bounding box of beam
          STH = SIN(BPA)
          CTH = COS(BPA)
          IF (STH.EQ.0.0) THEN
              XM = BMAJ
              YM = BMIN
          ELSE IF (CTH.EQ.0.0) THEN
              XM = BMIN
              YM = BMAJ
          ELSE
              PH = ATAN(-BMIN*STH/BMAJ/CTH)
              XM = ABS(BMAJ*COS(PH)*CTH - BMIN*SIN(PH)*STH)
              PH = ATAN(BMIN*CTH/BMAJ/STH)
              YM = ABS(BMIN*SIN(PH)*CTH + BMAJ*COS(PH)*STH)
          END IF
C         -- Choose location of beam
          IF (BPOS(1).EQ.0.0 .AND. BPOS(2).EQ.0.0) THEN
              XBC = X1 + SIGN(0.05*ABS(X2-X1)+XM*XFACT,X2-X1)
              YBC = Y1 + SIGN(0.05*ABS(Y2-Y1)+YM*YFACT,Y2-Y1)
          ELSE
              XBC = BPOS(1)
              YBC = BPOS(2)
          END IF
C         -- Plot it
          CALL PGBBUF
          CALL PGSCI(COLOR(1))
          CALL PGSLW(CONLW)
          CALL MPBEAM(XBC,YBC,XFACT*BMAJ,XFACT*BMIN,BPA)
          CALL PGSLW(LINEW)
          CALL PGEBUF
      END IF
C
C Draw the linear scale bar if requested.
C
      IF (REDSHI.GT.0.0 .AND. XFACT.EQ.YFACT) THEN
C         -- Compute linear scale SCL (pc/mas or kpc/arcsec)
          CALL COSMO(H0, Q0, REDSHI, SCL, DL)
          WRITE (OUTC,780) REDSHI, H0, Q0
          WRITE (OUTC,790) SCL
C         -- Get length of scale bar in pc (BAR) and mas (BAR/SCL)
          IF (BAR.LE.0.0) BAR = PGRND(0.05*SCL*ABS(X2-X1), JUNK)
          IBAR = NINT(BAR)
          IF (IBAR.LT.1) IBAR = 1
          IF (IBAR.GT.999) THEN
             IBAR = 1000*NINT(BAR/1000.0)
             WRITE (TEXT,793) IBAR/1000
          ELSE
             WRITE (TEXT,792) IBAR
          END IF
          BAR = IBAR
          WRITE (OUTC,791) IBAR, BAR/SCL
          BAR = BAR*XFACT/3600000.0
C         -- Plot scale bar
          CALL PGBBUF
          CALL PGSAVE
          CALL PGSCI(COLOR(1))
          CALL PGSLW(LINEW+2)
          CALL PGSCH(1.0)
          XBC = X2 - SIGN(0.05*ABS(X2-X1),X2-X1)
          YBC = Y1 + SIGN(0.05*ABS(Y2-Y1),Y2-Y1)
          CALL PGERRX(1, XBC, XBC-SIGN(BAR,X2-X1)/SCL, YBC, 1.0)
          CALL PGSLW(LINEW)
          CALL PGPTXT(XBC-0.5*SIGN(BAR,X2-X1)/SCL,
     :                Y1 + SIGN(0.07*ABS(Y2-Y1),Y2-Y1),
     :                0.0, 0.5, TEXT)
          CALL PGUNSA
          CALL PGEBUF
      END IF
C
C Annotation.
C
      CALL PGBBUF
      CALL PGSCI(COLOR(1))
      CALL PGSCH(CHSIZ1)
      CALL PGSLW(LINEW)
      CALL PGSCF(FONT)
      LX = MAX(1,LEN1(XUNIT))
      LY = MAX(1,LEN1(YUNIT))
      WRITE (OUTC,794) XUNIT(1:LX), YUNIT(1:LY)
      IF (DOSURF) THEN
      ELSE IF (DOTICK) THEN
          CALL PGBOX('BCNST',0.0,0,'BCNSTV',0.0,0)
          IF (DOTITL) CALL PGMTXT('B', 3.0, 0.5, 0.5, 
     1           XLABEL(1:LEN1(XLABEL))//' ('//XUNIT(1:LX)//')')
          IF (DOTITL) CALL PGMTXT('L', 3.0, 0.5, 0.5,
     1           YLABEL(1:LEN1(YLABEL))//' ('//YUNIT(1:LY)//')')
      ELSE
          CALL PGBOX('BC',0.0,0,'BC',0.0,0)
      END IF
      CALL LEGFMT(STR, TITLE1, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('T',2.5,0.0,0.0,STR)
      CALL LEGFMT(STR, TITLE2, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('T',0.8,0.0,0.0,STR)
      DO 1220 I=1,NLAB
          IF (LABELS(I).NE.' ') CALL
     1           PGPTXT(XLPOS(I), YLPOS(I), 0.0, 0.0, LABELS(I))
 1220 CONTINUE
      CALL PGSCH(CHSIZ1*1.2)
      CALL LEGFMT(STR, TOPLEF, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('T',-2.0,0.05,0.0,STR)
      CALL LEGFMT(STR, TOPRIG, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('T',-2.0,0.95,1.0,STR)
      CALL LEGFMT(STR, BOTLEF, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('B',-2.0,0.05,0.0,STR)
      CALL LEGFMT(STR, BOTRIG, MPOBJ, MPTEL, FREQS, MPDOB, WAVES,
     1                EPOCHS, Z)
      CALL PGMTXT('B',-2.0,0.95,1.0,STR)
      CALL PGSCH(CHSIZ1)
      IF (NMARK.NE.0) THEN
         CALL PGPT(1, XMARK, YMARK, NMARK)
      END IF
C
C Wedge display.
C
      Y = 4.0
      IF (DOGRAY.AND.WEDGE) THEN
         CALL PGWEDG('BI', Y, 3.0, GRAY1, GRAY2, MPBUN)
      END IF
      Y = Y+3.0
      Y = Y*(CHSIZ1/CHSIZ2)
C
C Additional annotation.
C
      IF (ANNOT.GE.0) THEN
          CALL PGSCH(CHSIZ2)
          CALL PGSCF(1)
          CALL PGSLW(1)
          DY = 1.5
C         -- Line 1: image maximum
          WRITE (TXTERR,730) ABSMAX, MPBUN
          CALL SSPACE(TXTERR,L)
          Y = Y+DY
          CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
C         -- Lines 2,3: contour levels
          IF (DOCNTR) THEN
              N = MIN(10,NCONS)
              WRITE (TXTERR,740) (VCN(ICON)*100.0/ABSMAX,ICON=1,N)
              CALL SSPACE(TXTERR,L)
              Y = Y+DY
              CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
              IF (NCONS.GT.10) THEN
                N = MIN(NCONS,20)
                WRITE (TXTERR,740) (VCN(ICON)*100.0/ABSMAX,ICON=11,N)
                CALL SSPACE(TXTERR,L)
                Y = Y+DY
                CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
              END IF
          END IF
C         -- Line 4: gray levels
          IF (DOGRAY) THEN
              WRITE (TXTERR,770) BLACK, WHITE
              CALL SSPACE(TXTERR,L)
              Y = Y+DY
              CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
          END IF
C         -- Line 5: rotation
          IF (MPRO(2).NE.0.0) THEN
              WRITE (TXTERR,745) MPRO(2)
              CALL SSPACE(TXTERR,L)
              Y = Y+DY
              CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
          END IF
C         -- Line 6: beam
          IF (MPBMJ.NE.0.0) THEN
              LL = MAX(1,LEN1(XUNIT))
              IF (MPBMJ.NE.MPBMN) THEN
                  WRITE (TXTERR,750) MPBMJ*XFACT, MPBMN*XFACT, 
     1                  XUNIT(1:LL), MPBPA-MPRO(2)
              ELSE
                  WRITE (TXTERR,760) MPBMJ*XFACT, XUNIT(1:LL)
              END IF
              CALL SSPACE(TXTERR,L)
              Y = Y+DY
              CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
          END IF
C         -- Line 7: file name
          IF (ANNOT.GE.1) THEN
              L = LEN1(MAPNAM)
              CALL FRDT(MAPNAM(1:L), REVTIM)
              Y = Y+DY
              CALL PGMTXT('B', Y,0.0,0.0,'File: '//MAPNAM(1:L)// ' ('//
     :             REVTIM//')')
C         -- Line 8: username and date
              CALL MKHIST('MAPPLOT', VERSN, TXTERR)
              L = LEN1(TXTERR)
              Y = Y+DY
              CALL PGMTXT('B',Y,0.0,0.0,TXTERR(1:L))
          END IF
C         -- end
          CALL PGSCH(CHSIZ1)
          CALL PGSCF(FONT)
          CALL PGSLW(LINEW)
      END IF
C      
C Flush graphics buffer. Do interactive window if requested.
C
      CALL PGEBUF
      CALL PGQINF('CURSOR', TEST, LTEST)
      IF (TEST.EQ.'NO')  SETW = .FALSE.
      IF (SETW) THEN
          CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
          CALL PGWNAD(X1,X2,Y1,Y2)
          CALL PGSCI(2)
          CALL MPEWIN(LRTB, 80, 2, IER)
          IF (IER.EQ.0) THEN
            DO 1230 I=1,80
              PARS(133+I) = LRTB(I)
 1230       CONTINUE
          END IF
      END IF
C
C Compute statistics.
C
      CALL MPSTAT(MAPARR, MPSIZ(1), MPSIZ(2), MSLO, MSHI, NSLO, NSHI,
     :          TR, LRTB, 80, OUTC)
C
C Get a new set of commands.
C
      CALL PGUPDT
      GOTO 1000
C
C Error exits.
C
 3000 CONTINUE
      CLOSE (UNIT=2,IOSTAT=IER)
      MAPOLD = '*** UNSPECIFIED ***'
      GOTO 1000
C
C Finish.
C
 2000 CALL PGEND
C-----------------------------------------------------------------------
      END

      SUBROUTINE MPBEAM (XC,YC,A,B,ANGLE)
C-----------------------------------------------------------------------
C Draw an ellipse in the current window.
C
C Arguments:
C
C XC (input, real) : world x-coordinate of the center of the ellipse.
C YC (input, real) : world y-coordinate of the center of the ellipse.
C A (input, real) : semi-major axis of ellipse (world coordinates).
C B (input, real) : semi-minor axis of ellipse (world coordinates).
C ANGLE (input, real) : the angle in radians that the major axis makes
C      with the x-axis of world coordinates; if x=0 or pi, the
C      major axis is parallel to the x-axis; if x = pi/2, the major
C      axis is parallel to the y-axis.
C-----------------------------------------------------------------------
      INTEGER  I, NP, CI, FS
      PARAMETER (NP=100)
      REAL     XP(NP), YP(NP)
      REAL     XC, YC, A, B, ANGLE, CA, SA, X, Y, DA
      REAL     PI
      PARAMETER (PI=3.141592653)
C
      CA = COS(ANGLE)
      SA = SIN(ANGLE)
      DA = 2*PI/NP
C
C Draw outline.
C
      DO 10 I=1,NP
          X = A*COS(I*DA)
          Y = B*SIN(I*DA)
          XP(I) = XC + X*CA - Y*SA
          YP(I) = YC + Y*CA + X*SA
   10 CONTINUE
      CALL PGBBUF
      CALL PGQCI(CI)
      CALL PGQFS(FS)
      CALL PGSCI(14)
      CALL PGSFS(1)
      CALL PGPOLY(NP, XP, YP)
      CALL PGSFS(2)
      CALL PGSCI(CI)
      CALL PGPOLY(NP, XP, YP)
      CALL PGSFS(FS)
      CALL PGEBUF
C
      END

      SUBROUTINE MPLWIN(LRTB, NDIM, CI)
      INTEGER NDIM, CI
      REAL LRTB(NDIM)
C-----------------------------------------------------------------------
      INTEGER I, NWIN
      REAL X1, X2, Y1, Y2
C
      CALL PGBBUF
      CALL PGSCI(CI)
      CALL PGSFS(2)
      NWIN = NDIM/4
      DO 10 I=1,NWIN
          X1 = -LRTB(4*I-3)
          X2 = -LRTB(4*I-2)
          Y1 = LRTB(4*I-1)
          Y2 = LRTB(4*I)
          IF (X1.NE.X2 .AND. Y1.NE.Y2) CALL PGRECT(X1, X2, Y1, Y2)
   10 CONTINUE
      CALL PGEBUF
      END

      SUBROUTINE MPEWIN(LRTB, NDIM, CI, IER)
      INTEGER NDIM, CI, IER
      REAL LRTB(NDIM)
C
C Routine for editing windows.
C-----------------------------------------------------------------------
      INTEGER I, J, NWIN, MAXWIN, STATE, IS, GETIN, JUNK, LEN1
      REAL X1, X2, Y1, Y2, XC, YC
      CHARACTER CH
      CHARACTER*255 FILE
C
      MAXWIN = NDIM/4
      NWIN = 0
      DO 10 I=1,NDIM,4
          X1 = -LRTB(I)
          X2 = -LRTB(I+1)
          Y1 = LRTB(I+2)
          Y2 = LRTB(I+3)
          IF (X1.EQ.X2 .OR. Y1.EQ.Y2) GOTO 15
          NWIN = NWIN+1
   10 CONTINUE
C
   15 CALL PUTOUT(' ')
      CALL PUTOUT('Editing windows:')
      CALL PUTOUT('  To DELETE a window, position the cursor in '//
     1            'the window and type D.')
      CALL PUTOUT('  To ADD a window, position the cursor at one '//
     1            'corner and type A,')
      CALL PUTOUT('    then position it at the opposite corner '//
     1            'and type A.')
      CALL PUTOUT('  Type Q to quit (exit without changing windows).')
      CALL PUTOUT('  Type X to exit and write lrtb.par.')
      CALL PUTOUT('  Type S to exit and save in named file.')
      CALL PUTOUT('[Mouse button 1=A, 2=D, 3=X]')
C
      CALL PGSCI(CI)
      XC = 0.0
      YC = 0.0
      STATE = 1
  20  CALL PGCURS(XC, YC, CH)
      CALL UPCASE(CH)
      IF (CH.EQ.'X') THEN
          IER = 0
          FILE = 'lrtb.par'
          GOTO 50
      ELSE IF (CH.EQ.'S') THEN
          IER = 0
          IS = GETIN(FILE, 'File name: ', JUNK)
          GOTO 50
      ELSE IF (CH.EQ.'Q') THEN
          IER = 1
          CALL PUTOUT('Window changes will not be saved!')
          RETURN
      ELSE IF (CH.EQ.'D') THEN
          IF (STATE.EQ.1 .AND. NWIN.LT.1) THEN
              CALL PUTOUT('No windows to delete')
          ELSE IF (STATE.EQ.2) THEN
              CALL PUTOUT('Current point cancelled')
              STATE = 1
              CALL PGSCI(0)
              CALL PGPT(1,X1,Y1,-1)
              CALL PGSCI(CI)
          ELSE
C             -- Find and delete window
              DO 30 I=1,NDIM,4
                  X1 = -LRTB(I)
                  X2 = -LRTB(I+1)
                  Y1 = LRTB(I+2)
                  Y2 = LRTB(I+3)
                  IF (XC.LE.X1 .AND. XC.GE.X2 .AND.
     1                YC.LE.Y1 .AND. YC.GE.Y2) THEN
                      CALL PGSCI(0)
                      CALL PGRECT(X1, X2, Y1, Y2)
                      CALL PGSCI(CI)
                      DO 25 J=I+4,NDIM
                          LRTB(J-4) = LRTB(J)
   25                 CONTINUE
                      LRTB(NDIM-4) = 0.0
                      LRTB(NDIM-3) = 0.0
                      LRTB(NDIM-2) = 0.0
                      LRTB(NDIM-1) = 0.0
                      NWIN = NWIN-1
                      GOTO 20
                  END IF
   30         CONTINUE
          END IF
      ELSE IF (CH.EQ.'A') THEN
          IF (NWIN.GE.MAXWIN) THEN
              CALL PUTOUT('Maximum number of windows exceeded')
          ELSE IF (STATE.EQ.1) THEN
C             -- first point
              X1 = XC
              Y1 = YC
              STATE = 2
              CALL PGPT(1,X1,Y1,-1)
          ELSE
C             -- second point
              X2 = XC
              Y2 = YC
              STATE = 1
              CALL PGRECT(X1,X2,Y1,Y2)
              NWIN = NWIN+1
              LRTB(4*NWIN-3) = MIN(-X1,-X2)
              LRTB(4*NWIN-2) = MAX(-X1,-X2)
              LRTB(4*NWIN-1) = MAX(Y1,Y2)
              LRTB(4*NWIN)   = MIN(Y1,Y2)
          END IF
      ELSE
          CALL PUTOUT('Type A, D, or X')
      END IF
      GOTO 20
C
C Done: save the current window list
C
   50 CALL PUTOUT('Saving windows in '//FILE(1:LEN1(FILE)))
      OPEN (UNIT=21, STATUS='UNKNOWN', FILE=FILE, IOSTAT=IS)
      IF (IS.NE.0) THEN
          CALL PUTOUT('Can''t open file!')
          IS = GETIN(FILE, 'File name: ', JUNK)
          GOTO 50
      END IF
      WRITE (21, '(A)') ' LRTB ='
      DO 60 I=1,NWIN
          IF (I.EQ.NWIN) THEN
              WRITE(21,'(3(F10.3,'',''),F10.3)') (LRTB(J),J=4*I-3,4*I)
          ELSE
              WRITE(21,'(4(F10.3,'',''))') (LRTB(J),J=4*I-3,4*I)
          END IF
   60 CONTINUE
      DO 70 I=4*NWIN+1,NDIM
          LRTB(I) = 0.0
   70 CONTINUE
      CLOSE (UNIT=21)
      END

      SUBROUTINE SURF(A, M, N, I1, I2, J1, J2, ANGLE, 
     :                CITOP, CIBOT, CIFRM)
      INTEGER M, N, I1, I2, J1, J2, CITOP, CIBOT, CIFRM
      REAL A(M,N), ANGLE
C-----------------------------------------------------------------------
C Surface mesh plot.
C
C A(M,N): array to be plotted
C I1,I2, J1,J2: range to be plotted (A(I,J): I1<=I<=I2, J1<=J<=J2)
C ANGLE: angle of view above horizontal (degrees)
C CITOP: color index for top surface
C CIBOT: color index for bottom surface
C CIFRM: color index for frame
C-----------------------------------------------------------------------
      REAL TMPARR(256*256)
      INTEGER NPIX, I, J, K, CI
C
      NPIX = (I2-I1+1)*(J2-J1+1)
      IF (NPIX.GT.256*256) THEN
          CALL PUTOUT('Image too big for 3D plot')
          RETURN
      END IF
C
C Copy subsection of array into temporary storage
C
      K = 0
      DO 20 J=J1,J2
          DO 10 I=I1,I2
              K = K+1
              TMPARR(K) = A(I,J)
   10     CONTINUE
   20 CONTINUE
C
C Plot lower surface
C
      CALL PGQCI(CI)
      CALL PGSCI(CIBOT)
      CALL FREDDY(TMPARR, I2-I1+1, J2-J1+1, 1.0, ANGLE, .TRUE., CIFRM)
C
C Copy subsection of array into temporary storage
C
      K = 0
      DO 40 J=J1,J2
          DO 30 I=I1,I2
              K = K+1
              TMPARR(K) = A(I,J)
   30     CONTINUE
   40 CONTINUE
C
C Plot upper surface
C
      CALL PGSCI(CITOP)
      CALL FREDDY(TMPARR, I2-I1+1, J2-J1+1, 1.0, ANGLE, .FALSE., CIFRM)
      CALL PGSCI(CI)
C
      END

      SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE,SURF,CIFRM)
      INTEGER KX, NY
      REAL ARRAY(KX,NY), SIZE, ANGLE
      LOGICAL SURF
      INTEGER CIFRM
C-----------------------------------------------------------------------
C Draws isometric plot of array
C
C ARRAY(KX,NY) : array to be plotted
C SIZE (real): width and height of plot window, assumed square, in
C PGPLOT world coordinates.
C SURF: .TRUE. for lower surface, .FALSE. for upper surface.
C-----------------------------------------------------------------------
      REAL FMAX,FMIN,DELTAX,DELTAY,DELTAV,SINE,PEAK,X,DX,HEIGHT
      REAL XBASE, YBASE
      INTEGER I,J,KI,KJ,NX,MX,MY,STEP,LEFT,RIGHT,IT,MN,INCX,IC
      LOGICAL VISBLE
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Statement functions define X and Y coordinates of the base of
C each pixel.
C
      XBASE(I,J) = DELTAX*(I-1+NY-J)
      YBASE(I,J) = DELTAY*(I+J-2)
C
      MN = KX*NY
      NX = KX
C     -- Check array size.
      IF(NX.LT.2 .OR. NY.LT.2) RETURN
C     -- Find extrema.
      FMAX = ARRAY(1,1)
      FMIN = FMAX
      DO 20 J=1,NY
          DO 10 I=1,NX
              FMIN = MIN(ARRAY(I,J),FMIN)
              FMAX = MAX(ARRAY(I,J),FMAX)
   10     CONTINUE
   20 CONTINUE
C     -- Plot scales
      DELTAX = SIZE/(NX+NY-2)
      SINE   = SIN(ANGLE/57.3)
      DELTAY = DELTAX*SINE
      HEIGHT = SIZE*(1.-ABS(SINE))
      DELTAV = HEIGHT
      FMAX   = FMAX-FMIN
      IF (FMAX.LT.0.0001) FMAX = DELTAV
      DELTAV = DELTAV/FMAX
      MX = NX+1
      MY = NY+1
      STEP = MX
C
C Start PGPLOT buffering.
C
      CALL PGBBUF
C
C Work our way down the Y axis, then up the X axis, calculating the Y
C plotter coordinates for each column of the plot, doing the hidden-line
C suppression at the same time. The array value is replaced by the
C y-coordinate of that pixel.
C
      DO 50 J=1,NY
          KJ = MY-J
          KI = 1
C               ( KI,KJ are coordinates of bottom of column)
          ARRAY(KI,KJ) = YBASE(KI,KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   30     PEAK = ARRAY(KI,KJ)
   40     KI = KI+1
          KJ = KJ+1
          IF (KI.GT.NX .OR. KJ.GT.NY) GOTO 50
          ARRAY(KI,KJ) = YBASE(KI,KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(SURF.XOR.(ARRAY(KI,KJ).GT.PEAK)) THEN
              GOTO 30
          ELSE
              ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          END IF
          GOTO 40
   50 CONTINUE
C
C Now to work our way up the X axis
C
      DO 80 I=2,NX
          KI = I
          KJ = 1
          ARRAY(KI,KJ) = YBASE(KI,KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   60     PEAK = ARRAY(KI,KJ)
   70     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 80
          ARRAY(KI,KJ) = YBASE(KI,KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(SURF.XOR.(ARRAY(KI,KJ).GT.PEAK)) THEN
              GOTO 60
          ELSE
              ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          END IF
          GOTO 70
   80 CONTINUE
C
C Draw the outline frame.
C
      IF (SURF) THEN
          CALL PGQCI(IC)
          CALL PGSCI(CIFRM)
          CALL PGMOVE(XBASE(NX,1),  YBASE(NX,1) )
          CALL PGDRAW(XBASE(1,1),   YBASE(1,1)  )
          CALL PGDRAW(XBASE(1,NY),  YBASE(1,NY) )
          CALL PGDRAW(XBASE(NX,NY), YBASE(NX,NY))
          CALL PGDRAW(XBASE(NX,1),  YBASE(NX,1) )
          CALL PGMOVE(XBASE(NX,1),  YBASE(NX,1) )
          CALL PGDRAW(XBASE(NX,1),  ABS(ARRAY(NX,1)) )
          CALL PGMOVE(XBASE(1,1),   YBASE(1,1)  )
          CALL PGDRAW(XBASE(1,1),   ABS(ARRAY(1,1))  )
          CALL PGMOVE(XBASE(1,NY),  YBASE(1,NY) )
          CALL PGDRAW(XBASE(1,NY),  ABS(ARRAY(1,NY)) )
          CALL PGMOVE(XBASE(NX,NY), YBASE(NX,NY))
          CALL PGDRAW(XBASE(NX,NY), ABS(ARRAY(NX,NY)))
          CALL PGSCI(IC)
      END IF
C
C Array is now ready for plotting.  If a point is
C positive, then it is to be plotted at that Y
C coordinate; if it is negative, then it is
C invisible, but at minus that Y coordinate (the point
C where the line heading towards it disappears has to
C be determined by finding the intersection of it and
C the cresting line).
C
C Plot rows:
C
      DO 110 J=1,NY,2
          KJ = MY-J
          DX = DELTAX*(J-2)
          X = DX+DELTAX
          CALL PGMOVE(X,ARRAY(1,KJ))
          VISBLE = .TRUE.
          DO 90 I=2,NX
              RIGHT = I+NX*(KJ-1)
              LEFT = RIGHT-1
              IT = RIGHT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
   90     CONTINUE
C
C Now at far end of row so come back
C
          KJ = KJ-1
          IF(KJ.LE.0) GOTO 170
          VISBLE = ARRAY(NX,KJ).GE.0.0
          DX = DELTAX*(NX+J)
          IF(VISBLE) CALL PGMOVE(DX-DELTAX,ARRAY(NX,KJ))
          DELTAX = -DELTAX
          DO 100 I=2,NX
              KI = MX-I
              LEFT = KI+NX*(KJ-1)
              RIGHT = LEFT+1
              IT = LEFT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
  100     CONTINUE
C
          X = DX+DELTAX*NX
C               (set DELTAX positive for return trip)
          DELTAX = -DELTAX
  110 CONTINUE
C
C Now do the columns:
C as we fell out of the last DO-loop we do the
C columns in ascending-X order
C
      INCX = 1
      KI = 1
C               (set DELTAX -ve since scanning R to L)
  120 DX = DELTAX*(KI+NY-1)
      DELTAX = -DELTAX
      X = DX+DELTAX
      CALL PGMOVE(X,ARRAY(1,1))
  130 VISBLE = .TRUE.
      DO 140 J=2,NY
          LEFT = KI+NX*(J-1)
          RIGHT = LEFT-NX
          IT = LEFT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
  140 CONTINUE
C
C At far end, increment X and check still inside array
C
      KI = KI+INCX
      IF(KI.LE.0 .OR. KI.GT.NX) GOTO 180
      VISBLE = ARRAY(KI,NY).GE.0.0
      DELTAX = -DELTAX
      DX = DELTAX*(KI-2)
      X = DX+DELTAX
      IF(VISBLE) CALL PGMOVE(X,ARRAY(KI,NY))
      DO 150 J=2,NY
          KJ = MY-J
          RIGHT = KI+NX*(KJ-1)
          LEFT = RIGHT+NX
          IT = RIGHT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
  150 CONTINUE
C
      X = DX+DELTAX*NY
      IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(KI,1))
      IF(KI.EQ.1) GOTO 180
      KI = KI+INCX
      IF(KI.GT.NX) GOTO 180
      IF(KI.EQ.1) GOTO 120
  160 DELTAX = -DELTAX
      DX = DELTAX*(1-KI-NY)
      X = DX+DELTAX
      CALL PGMOVE(X,ARRAY(KI,1))
      GOTO 130
C
C Do columns backwards because ended rows at far end of X
C
  170 KI = NX
      INCX = -1
      DX = DELTAX*(KI+NY)
      GOTO 160
C
C
  180 CALL PGEBUF
      END
C-----------------------------------------------------------------------
      SUBROUTINE FREDGO(ARRAY,MN)
      INTEGER MN
      REAL ARRAY(MN)
C
      INTEGER STEP,LEFT,RIGHT,IT,NX
      LOGICAL VISBLE
      REAL AL,AR,BL,EM,XX,X,Y,DELTAX
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Test visibility
C
      IF(ARRAY(IT).LT.0.0) GOTO 80
C
C This point is visible - was last?
C
      IF(VISBLE) GOTO 50
C
C No: calculate point where this line vanishes
C
   10 IF(LEFT.LE.NX .OR. MOD(LEFT-1,NX).EQ.0 .OR.
     1     RIGHT.LE.NX .OR. MOD(RIGHT-1,NX).EQ.0) GOTO 100
      AL = ABS(ARRAY(LEFT))
      AR = ABS(ARRAY(RIGHT))
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C               Right-hand point is crested
   20 RIGHT = RIGHT-STEP
      IF(ARRAY(RIGHT).LT.0.0) GOTO 20
C               Left-hand end of cresting line is either
C               RIGHT+NX or RIGHT-1
      LEFT = RIGHT+NX
      IF(ARRAY(LEFT).LT.0.0) LEFT = RIGHT-1
C
C               RIGHT and LEFT index into the endpoints of the
C               cresting line
   30 BL = ABS(ARRAY(LEFT))
      EM = ABS(ARRAY(RIGHT))-BL
      XX = EM-AR+AL
      IF(ABS(XX).LT.0.0001) GOTO 60
      XX = (AL-BL)/XX
   40 Y = EM*XX+BL
      IF(DELTAX.GT.0.0) XX = 1.0-XX
      XX = X-XX*DELTAX
      IF(VISBLE) GOTO 90
C               Drawing a line from an invisible point
C               to a visible one
      CALL PGMOVE(XX,Y)
      VISBLE = .TRUE.
   50 CALL PGDRAW(X,ARRAY(IT))
      RETURN
C
   60 XX = 0.5
      GOTO 40
C
C Left-hand point crested
C
   70 LEFT = LEFT-STEP
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C
C Right-hand end of cresting line is either LEFT+1 or LEFT-NX
C
      RIGHT = LEFT+1
      IF(ARRAY(RIGHT).LT.0.0) RIGHT = LEFT-NX
      GOTO 30
C
C This point is invisible; if last one was too, then forget it;
C else draw a line towards it
C
   80 IF(.NOT.VISBLE) RETURN
      GOTO 10
C
   90 CALL PGDRAW(XX,Y)
  100 VISBLE = .FALSE.
      RETURN
      END

      SUBROUTINE MPSTAT(A, IDIM, JDIM, I1, I2, J1, J2, TR, LRTB, 
     1                  LDIM, PR)
      INTEGER IDIM, JDIM, I1, I2, J1, J2, LDIM, PR
      REAL A(IDIM, JDIM), TR(6), LRTB(LDIM)
C
C Compute statistics of map within windows.
C-----------------------------------------------------------------------
      INTEGER NPIX, I, J, L, NWIN
      LOGICAL INSIDE, RANGE
      REAL X, Y, AMAX, AMIN, XMAX, YMAX, XMIN, YMIN, SUM, SUMSQ
      REAL SUMX, SUMY, MEAN, RMS
      REAL T, T1, T2
C
  800 FORMAT('  Number of windows:            ',I10/
     A       '  Number of cells in windows:   ',I10/
     1       '  Maximum in windows: ',1P,G13.4,'  at (x,y) ',0P,2F10.3/
     2       '  Minimum in windows: ',1P,G13.4,'  at (x,y) ',0P,2F10.3/
     3       '  Mean in windows:    ',0P,G13.4/
     4       '  RMS in windows:     ',0P,G13.4/
     5       '  Centroid:           ',0P,2F10.3)
C
      RANGE(T,T1,T2) = (T1.LE.T.AND.T.LT.T2) .OR. (T2.LE.T.AND.T.LT.T1)
C-----------------------------------------------------------------------
C
C Initialize counters.
C
      AMAX = -1.E37
      AMIN =  1.E37
      SUM  =  0.0
      SUMSQ = 0.0
      SUMX = 0.0
      SUMY = 0.0
      NPIX = 0
C
C How many windows?
C
      NWIN = 0
      DO 10 L=1,LDIM,4
          IF (LRTB(L).EQ.LRTB(L+1) .OR. LRTB(L+2).EQ.LRTB(L+3)) GOTO 20
          NWIN = NWIN+1
   10 CONTINUE
   20 CONTINUE
      IF (NWIN.LT.1) RETURN
C
C Search map for pixels in windows.
C
      DO 60 J=J1,J2
          DO 50 I=I1,I2
C             -- Find world coordinates of this pixel
              X = TR(1) + TR(2)*I + TR(3)*J
              Y = TR(4) + TR(5)*I + TR(6)*J
C             -- Is it within the windows?
              INSIDE = .FALSE.
              DO 30 L=1,4*NWIN,4
                  IF (RANGE(X,-LRTB(L),-LRTB(L+1)) .AND.
     :                RANGE(Y,LRTB(L+2),LRTB(L+3))) THEN
                      INSIDE = .TRUE.
                      GOTO 40
                  END IF
   30         CONTINUE
   40         CONTINUE
C             -- Accumulate statistics
              IF (INSIDE) THEN
                  NPIX = NPIX+1
                  IF (A(I,J).GT.AMAX) THEN
                     AMAX = A(I,J)
                     XMAX = X
                     YMAX = Y
                  END IF
                  IF (A(I,J).LT.AMIN) THEN
                     AMIN = A(I,J)
                     XMIN = X
                     YMIN = Y
                  END IF
                  SUM = SUM + A(I,J)
                  SUMSQ = SUMSQ + A(I,J)**2
                  SUMX = SUMX + X*A(I,J)
                  SUMY = SUMY + Y*A(I,J)
              END IF
   50     CONTINUE
   60 CONTINUE
C
C     Calculate and write statistics
C
      MEAN = SUM / REAL(NPIX)
      RMS = SQRT ( (SUMSQ/REAL(NPIX)) - MEAN**2)
      SUMX = SUMX/SUM
      SUMY = SUMY/SUM
      WRITE (PR, 800)  NWIN, NPIX, AMAX, XMAX, YMAX, AMIN, XMIN, YMIN,
     1                 MEAN, RMS, SUMX, SUMY
C
      END

      SUBROUTINE LEGFMT(OUT, FMT, S, T, F, D, W, E, Z)
      CHARACTER*(*) OUT, FMT, S, T, F, D, W, E, Z
C
C Make a legend: copy FMT string to OUT, replacing special strings
C as follows:
C    %S replace by string S (source name)
C    %T replace by string T (telescope name)
C    %F replace by string F (frequency)
C    %D replace by string D (date)
C    %W replace by string W (wavelength)
C    %E replace by string E (epoch)
C    %Z replace by string Z (redshift)
C    %% replace by '%'
C    %x ignore (where x is any character except the above)
C-----------------------------------------------------------------------
      INTEGER LEN1
      INTEGER LS, LT, LF, LD, LW, LE, LZ, I, J
      CHARACTER*1 TEST
C
      LS = LEN1(S)
      LT = LEN1(T)
      LF = LEN1(F)
      LD = LEN1(D)
      LW = LEN1(W)
      LE = LEN1(E)
      LZ = LEN1(Z)
C
      I = 0
      J = 0
   10 CONTINUE
          J = J+1
          IF (J.GT.LEN(FMT)) GOTO 100
          IF (FMT(J:J).EQ.'%') THEN
              J = J+1
              TEST = FMT(J:J)
              CALL UPCASE(TEST)
              IF (TEST.EQ.'S') THEN
                  IF (I+LS.GT.LEN(OUT)) GOTO 100
                  IF (LS.GT.0) OUT(I+1:I+LS) = S(1:LS)
                  I = I+LS
              ELSE IF (TEST.EQ.'T') THEN
                  IF (I+LT.GT.LEN(OUT)) GOTO 100
                  IF (LT.GT.0) OUT(I+1:I+LT) = T(1:LT)
                  I = I+LT
              ELSE IF (TEST.EQ.'F') THEN
                  IF (I+LF.GT.LEN(OUT)) GOTO 100
                  IF (LF.GT.0) OUT(I+1:I+LF) = F(1:LF)
                  I = I+LF
              ELSE IF (TEST.EQ.'D') THEN
                  IF (I+LD.GT.LEN(OUT)) GOTO 100
                  IF (LD.GT.0) OUT(I+1:I+LD) = D(1:LD)
                  I = I+LD
              ELSE IF (TEST.EQ.'W') THEN
                  IF (I+LW.GT.LEN(OUT)) GOTO 100
                  IF (LW.GT.0) OUT(I+1:I+LW) = W(1:LW)
                  I = I+LW
              ELSE IF (TEST.EQ.'E') THEN
                  IF (I+LE.GT.LEN(OUT)) GOTO 100
                  IF (LE.GT.0) OUT(I+1:I+LE) = E(1:LE)
                  I = I+LE
              ELSE IF (TEST.EQ.'Z') THEN
                  IF (I+LZ.GT.LEN(OUT)) GOTO 100
                  IF (LZ.GT.0) OUT(I+1:I+LZ) = Z(1:LZ)
                  I = I+LZ
              ELSE IF (TEST.EQ.'%') THEN
                  IF (I+1.GT.LEN(OUT)) GOTO 100
                  OUT(I+1:I+1) = '%'
                  I = I+1
              END IF
          ELSE
              IF (I.GE.LEN(OUT)) GOTO 100
              I = I+1
              OUT(I:I) = FMT(J:J)
          END IF
      GOTO 10
C
  100 DO 110 J=I+1,LEN(OUT)
          OUT(J:J) = ' '
  110 CONTINUE
C
      END

      SUBROUTINE COSMO(H, Q, Z,  SCL, DL)
      REAL H, Q, Z,  SCL, DL
C-----------------------------------------------------------------------
C Input parameters:
C     H    Hubble constant, km/s/Mpc [e.g., 100.0]
C     Q    deceleration parameter    [e.g., 0.0]
C     Z    redshift                  [e.g., 0.518]
C Output parameters:
C     SCL  linear scale, kpc/arcsec = pc/mas
C     DL   luminosity distance, Mpc
C-----------------------------------------------------------------------
      DOUBLE PRECISION RT, D, DF, X, XKM, XLY, VOVERC, ZZ
      DOUBLE PRECISION C, PI, THRAD, PC
C
C     Constants:
C        C = speed of light (km/sec)
C        PI = pi
C        THRAD = radians/milliarcsec
C        PC = km/pc
      PARAMETER (C=2.997924580D5)
      PARAMETER (PI=3.14159265359D0)
      PARAMETER (THRAD=PI/(1000D0*3600D0*180D0))
      PARAMETER (PC=30857D9)
C-----------------------------------------------------------------------
C
C DL is the luminosity distance (in units of c/H); formula from
C Peacock.
C
      ZZ = DBLE(Z)
      RT = SQRT(1D0 + 2D0*Q*ZZ)
      D = ZZ * (1D0 + RT + ZZ) / (1D0 + RT + Q*ZZ)
C
C DF is the angular-diameter distance (in units of c/H).
C
      DF = D/((1D0+ZZ)*(1D0+ZZ))
C
      X   = THRAD*C*DF/H*1D6
      XKM = X*PC*1D3
      XLY = XKM/C/(365.25*24*3600.)
      D   = D*C/H
      VOVERC = XLY*(1D0+ZZ)/1D3
C
      SCL = REAL(X)
      DL  = REAL(D)
      END

      SUBROUTINE ANGLUN (PREF,RANGE,UNIT,FACTOR)
C-----------------------------------------------------------------------
C Arguments:
C  PREF   (input, character): preferred unit, or blank for default.
C  RANGE  (input, real): coordinate range
C  UNIT   (output, character): units to use
C  FACTOR (output, real): scale factor to multiply coordinates by to
C      convert from degrees to UNIT.
C-----------------------------------------------------------------------
      CHARACTER*(*) PREF, UNIT
      REAL RANGE, FACTOR
C-----------------------------------------------------------------------
      IF (PREF.EQ.'MAS') THEN
          UNIT = 'mas'
          FACTOR = 3 600 000.
      ELSE IF (PREF.EQ.'MILLIARCSEC') THEN
          UNIT = 'milliarcsec'
          FACTOR = 3 600 000.
      ELSE IF (PREF.EQ.'ARCSEC') THEN
          UNIT = 'arcsec'
          FACTOR = 3600.
      ELSE IF (PREF.EQ.'ARCMIN') THEN
          UNIT = 'arcmin'
          FACTOR = 60.0
      ELSE IF (PREF.EQ.'DEGREE') THEN
          UNIT = 'degrees'
          FACTOR = 1.0
      ELSE IF (PREF.EQ.'RAD') THEN
          UNIT = 'rad'
          FACTOR = 0.01745329252
      ELSE IF (RANGE.LT.1/3600.) THEN
          UNIT = 'milliarcsec'
          FACTOR = 3 600 000.
      ELSE IF (RANGE.LT.1/60.) THEN
          UNIT = 'arcsec'
          FACTOR = 3600.
      ELSE IF (RANGE.LT.1.0) THEN
          UNIT = 'arcmin'
          FACTOR = 60.0
      ELSE
          UNIT = 'degrees'
          FACTOR = 1.0
      END IF
      RETURN
      END

      SUBROUTINE PALETT(TYPE)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE, I1, I2
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
      REAL TL(4), TR(4), TG(4), TB(4)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      DATA TL /0.0, 0.5, 0.5, 1.0/
      DATA TR /0.2, 0.6, 0.6, 1.0/
      DATA TG /0.0, 0.0, 0.5, 1.0/
      DATA TB /1.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, 1.0, 0.5)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, 1.0, 0.5)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, 1.0, 0.5)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, 1.0, 0.5)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, 1.0, 0.5)
      ELSE IF (TYPE.EQ.6) THEN
C        -- TJP
         CALL PGCTAB(TL, TR, TG, TB, 4, 1.0, 0.5)
      END IF
      CALL PGQCIR(I1,I2)
C      CALL PGSCR(I1, 0.5,0.6,0.6)
      CALL PGSCR(I1, 1.0,1.0,1.0)
      END
C
C     CALDJ, CLDJ subroutines from SLALIB
C
      SUBROUTINE CALDJ (IY, IM, ID, DJM, J)
*+
*     - - - - - -
*      C A L D J
*     - - - - - -
*
*  Gregorian Calendar to Modified Julian Date
*
*  (Includes century default feature:  use CLDJ for years
*   before 100AD.)
*
*  Given:
*     IY,IM,ID     int    year, month, day in Gregorian calendar
*
*  Returned:
*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
*     J            int    status:
*                           0 = OK
*                           1 = bad year   (MJD not computed)
*                           2 = bad month  (MJD not computed)
*                           3 = bad day    (MJD computed)
*
*  Acceptable years are 00-49, interpreted as 2000-2049,
*                       50-99,     "       "  1950-1999,
*                       100 upwards, interpreted literally.
*
*  Called:  CLDJ
*
*  P.T.Wallace   Starlink   November 1985
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

      INTEGER NY

*  Default century if appropriate
      IF (IY.GE.0.AND.IY.LE.49) THEN
         NY=IY+2000
      ELSE IF (IY.GE.50.AND.IY.LE.99) THEN
         NY=IY+1900
      ELSE
         NY=IY
      END IF

*  Modified Julian Date
      CALL CLDJ(NY,IM,ID,DJM,J)

      END



      SUBROUTINE CLDJ (IY, IM, ID, DJM, J)
*+
*     - - - - -
*      C L D J
*     - - - - -
*
*  Gregorian Calendar to Modified Julian Date
*
*  Given:
*     IY,IM,ID     int    year, month, day in Gregorian calendar
*
*  Returned:
*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
*     J            int    status:
*                           0 = OK
*                           1 = bad year   (MJD not computed)
*                           2 = bad month  (MJD not computed)
*                           3 = bad day    (MJD computed)
*
*  The year must be -4699 (i.e. 4700BC) or later.
*
*  The algorithm is derived from that of Hatcher 1984
*  (QJRAS 25, 53-55).
*
*  P.T.Wallace   Starlink   11 March 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31,28,31,30,31,30,31,31,30,31,30,31 /



*  Preset status
      J=0

*  Validate year
      IF (IY.LT.-4699) THEN
         J=1
      ELSE

*     Validate month
         IF (IM.GE.1.AND.IM.LE.12) THEN

*        Allow for leap year
            IF (MOD(IY,4).EQ.0) THEN
               MTAB(2)=29
            ELSE
               MTAB(2)=28
            END IF
            IF (MOD(IY,100).EQ.0.AND.MOD(IY,400).NE.0)
     :         MTAB(2)=28

*        Validate day
            IF (ID.LT.1.OR.ID.GT.MTAB(IM)) J=3

*        Modified Julian Date
            DJM=DBLE((1461*(IY-(12-IM)/10+4712))/4
     :               +(306*MOD(IM+9,12)+5)/10
     :               -(3*((IY-(12-IM)/10+4900)/100))/4
     :               +ID-2399904)

*        Bad month
         ELSE
            J=2
         END IF

      END IF

      END

      SUBROUTINE FITRD (INUNIT,AMAP,MAXPIX,STATUS,ERROR)
C-----------------------------------------------------------------------
C FITRD: read a FITS-format disk file.
C
C     Parameters -   (">" input, "<" output )
C
C     (>) INUNIT  (Integer) unit on which file has been opened with
C                 FTOPEN.
C     (<) AMAP    (Real array) Array to contain the map pixels.
C     (>) MAXPIX  (Integer) Dimension of AMAP.
C     (<) STATUS  (Integer) A status code for the read.
C                 0 => OK, 1 => error
C     (<) ERROR   (Character) Returns with a message describing
C                 the error.
C Revised: 12-May-2000
C-----------------------------------------------------------------------
C     IMPLICIT NONE
C
      INTEGER INUNIT,MAXPIX,STATUS
      REAL AMAP(MAXPIX)
      CHARACTER*(*) ERROR
      INCLUDE 'mapcom.inc'
      DOUBLE PRECISION  RVAL
      INTEGER I, BITPIX, PIXELS
      CHARACTER KEYW*8, CVAL*80, COMMEN*80
      LOGICAL NV, EXTEND, SIMPLE
      INTEGER MAXDIM
      PARAMETER (MAXDIM=7)
      INTEGER NAXIS, NAXES(MAXDIM), PCOUNT, GCOUNT
C-----------------------------------------------------------------------
C
C            Initial values
C
      MPNDM = 0
      MPNIT = 0
      MPOBJ = ' '
      MPTEL = ' '
      MPINS = ' '
      MPOBS = ' '
      MPDOB = ' '
      MPDMP = ' '
      MPBUN = ' '
      MPTYP = 'MAP'
      MPCRT = ' '
      MPBMJ = 0D0
      MPBMN = 0D0
      MPBPA = 0D0
      MPEPO = 0D0
      MPFRQ = 0D0
      MPBWD = 0D0
      DO 100 I=1,7
          MPSIZ(I) = 1
          MPCTY(I) = ' '
          MPRV(I) = 0D0
          MPRP(I) = 0D0
          MPIN(I) = 1D0
          MPRO(I) = 0D0
  100 CONTINUE
      MPCTY(1) = 'x'
      MPCTY(2) = 'y'
      MPCTY(3) = 'z'

C
C            Read the standard header keywords.
C
      STATUS = 0
C     WRITE (*,*) '--- Reading header'
      CALL FTGHPR(INUNIT, MAXDIM, SIMPLE, BITPIX, NAXIS, NAXES, 
     :     PCOUNT, GCOUNT, EXTEND, STATUS)
      IF (STATUS.NE.0 .OR. .NOT.SIMPLE) THEN
         ERROR = 'Bad FITS file'
         GOTO 550
      END IF
      IF (PCOUNT.NE.0. OR. GCOUNT.NE.1 .OR. NAXIS.LT.2) THEN
         ERROR = 'FITS file is not an image'
         GOTO 550
      END IF
      MPNDM = NAXIS
C
C Read additional keywords
C
C     WRITE (*,*) '--- Reading additional keywords'
      STATUS = 0
      CALL FTGKYS(INUNIT, 'OBJECT', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPOBJ = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'TELESCOP', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPTEL= CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'INSTRUME', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPINS = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'INSTRUME', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPINS = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'OBSERVER', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPOBS = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'DATE-OBS', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPDOB = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'DATE-MAP', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPDMP = CVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'BUNIT', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPBUN = CVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'BMAJ', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPBMJ = RVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'BMIN', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPBMN = RVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'BPA', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPBPA = RVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'BPA', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPBPA = RVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'EPOCH', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPEPO = RVAL
      STATUS = 0
      CALL FTGKYD(INUNIT, 'NITER', RVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPNIT = RVAL
      STATUS = 0
      CALL FTGKYS(INUNIT, 'ORIGIN', CVAL, COMMEN, STATUS)
      IF (STATUS.EQ.0) MPCRT = CVAL
C
C Read axis keyords
C
C     WRITE (*,*) '--- Reading axis keywords'
      DO I=1,NAXIS
         MPSIZ(I) = NAXES(I)
         WRITE (KEYW, '(A5,I1)') 'CRVAL', I
         STATUS = 0
         CALL FTGKYD(INUNIT, KEYW, RVAL, COMMEN, STATUS)
         IF (STATUS.EQ.0) MPRV(I) = RVAL
         WRITE (KEYW, '(A5,I1)') 'CRPIX', I
         STATUS = 0
         CALL FTGKYD(INUNIT, KEYW, RVAL, COMMEN, STATUS)
         IF (STATUS.EQ.0) MPRP(I) = RVAL
         WRITE (KEYW, '(A5,I1)') 'CDELT', I
         STATUS = 0
         CALL FTGKYD(INUNIT, KEYW, RVAL, COMMEN, STATUS)
         IF (STATUS.EQ.0) MPIN(I) = RVAL
         WRITE (KEYW, '(A5,I1)') 'CTYPE', I
         STATUS = 0
         CALL FTGKYS(INUNIT, KEYW, CVAL, COMMEN, STATUS)
         IF (STATUS.EQ.0) MPCTY(I) = CVAL
         WRITE (KEYW, '(A5,I1)') 'CROTA', I
         STATUS = 0
         CALL FTGKYD(INUNIT, KEYW, RVAL, COMMEN, STATUS)
         IF (STATUS.EQ.0) MPRO(I) = RVAL
      END DO
C
C            Look for frequency on an axis.
C
      DO I=1,MPNDM
         IF (MPCTY(I).EQ.'FREQ') MPFRQ = MPRV(I)*1E-6
      END DO
C
C            Compute size of map.
C
      PIXELS = 0
      IF (MPNDM.GT.0) THEN
         PIXELS = 1
         DO I=1,MPNDM
            PIXELS = PIXELS*MPSIZ(I)
         END DO
      END IF
      IF (PIXELS.LE.0) THEN
         ERROR = 'Map contains no pixels'
         GOTO 550
      ELSE IF (PIXELS.GT.MAXPIX) THEN
         ERROR = 'Map too big'
         GOTO 550
      END IF
C     WRITE (*,*) '--- Number of pixels = ', PIXELS
C
C            Read in PIXELS pixels.
C
C     WRITE (*,*) '--- Reading image pixels'
      STATUS = 0
      NV = .FALSE.
      CALL FTGPVE(INUNIT, 1, 1, PIXELS, 1E38, AMAP, NV, STATUS)
      IF (STATUS.NE.0) THEN
         ERROR = 'Error reading image array'
         GOTO 550
      END IF
      IF (NV) THEN
         CALL PUTOUT('Caution: image contains blanked pixels')
      END IF
C
C            All done: return.  STATUS=0 implies success,
C            STATUS=1 implies failure.
C
      STATUS = 0
      ERROR = ' '
  550 RETURN
      END
