      PROGRAM IMAGE
      REAL ARRAY(1024,1024)
      CHARACTER*120 INFILE, IMFILE
      REAL XINT, YINT, FREQ
      INTEGER NX, NY
C
      WRITE (*,'(A,$)') 'Input file name: '
      READ (*,'(A)', ERR=100) INFILE
      WRITE (*,'(A,$)') 'Output [FITS] file name: '
      READ (*,'(A)', ERR=100) IMFILE
C
      OPEN (UNIT=1, STATUS='OLD', FILE=INFILE)
      READ (1,*,ERR=100) NX, NY, XINT, YINT, FREQ
      WRITE (*,'(1X,2I6,3F10.3)') NX, NY, XINT, YINT, FREQ
      IF (NX.LT.4 .OR. NY.LT.4 .OR. NX*NY.GT.1024*1024
     :    .OR. XINT.LE.0.0 .OR. YINT.LE.0.0) THEN
         WRITE (*,*) 'Invalid image parameters in input file'
      ELSE
         CALL MKIMG(ARRAY, NX, NY, XINT, YINT, 1, IMFILE, FREQ*1E9)
      END IF
 100  CONTINUE
      END

      SUBROUTINE MKIMG(SKY, NX, NY, XINT, YINT, UNIT, FMAP, FREQ)
      INTEGER NX, NY
      REAL SKY(NX,NY)
      REAL XINT, YINT
      REAL FREQ
      INTEGER UNIT
      CHARACTER*(*) FMAP
C-----------------------------------------------------------------------
C Program to create a FITS image containing a temperature map with
C specified patterns.
C-----------------------------------------------------------------------
      CHARACTER*(*) PROGNM, VERSN
      PARAMETER (PROGNM='IMAGE')
      PARAMETER (VERSN='1998-02-10')
      REAL       PI,RPDEG,RPSA
      PARAMETER  (PI=3.1415926536)
      PARAMETER  (RPDEG=PI/180.0)
      PARAMETER  (RPSA=PI/648000.0)
C
      CHARACTER*120 MSG
      CHARACTER*12  OBSERV
      CHARACTER*32  STRING
      CHARACTER*10  DATOBS
      DOUBLE PRECISION RA0, DEC0
      INTEGER BITPIX
      INTEGER IX, IY, TYPE, K, IER, ID, IMON, IYR
      INTEGER NAXIS, NAXISN(4)
      INTEGER LEN1
      LOGICAL SIMPLE, EXTEND
      REAL DEL1, DEL2, PIX1, PIX2
      REAL FACT, AMP
      REAL SKYMIN, SKYMAX, SKYTOT
      REAL TEMP, X, Y, A, B, R, PHI, COSPHI, SINPHI, X0, Y0
      REAL XC, YC
      REAL XP, YP, RSQ
      REAL BETA, EXPON
C
      DATA RA0/0D0/, DEC0/0D0/
      FACT = 4.0*LOG(2.0)

C Image cell size (arcmin), center

      XC = NX/2
      YC = NY/2

C Set the image to zero

      DO IY=1,NY
         DO IX=1,NX
            SKY(IX,IY) = 0.0
         END DO
      END DO

C Read parameters of model component
C     TYPE = 1 (gaussian)
C            2 (elliptical disk)
C            3 (rectangle)
C            4 (SZE decrement)
C     TEMP = peak brightness (K)
C     X0, Y0 = location of center (arcmin)
C     A = major axis, or cluster core radius for SZE (arcmin)
C     R = axial ratio (<1)
C     PHI = position angle (deg)

      DO K=1,1000
C         READ (UNIT,*,ERR=100,END=100) TYPE, TEMP, X0, Y0, A, R, PHI
         X0 = 0.0
         Y0 = 0.0
         A = 0.0
         R = 0.0
         PHI = 0.0
         BETA = 0.0
         CALL GTMODL(UNIT, TYPE, TEMP, X0, Y0, A, R, PHI, BETA, 6)
         IF (TYPE.LT.0) GOTO 100
         IF (BETA.EQ.0.0) BETA = 2.0/3.0
         IF (TYPE.EQ.4) THEN
            WRITE (*,'(1X,I3,E12.4,7F10.3)')
     :           TYPE, TEMP, X0, Y0, A, R, PHI, BETA
         ELSE
            WRITE (*,'(1X,I3,E12.4,7F10.3)')
     :           TYPE, TEMP, X0, Y0, A, R, PHI
         END IF
         COSPHI = COS(PHI*RPDEG)
         SINPHI = SIN(PHI*RPDEG)
         B = A*R

C Add in this component

         DO IY=1,NY
            DO IX=1,NX
               X = (XC-IX)*XINT
               Y = (YC-IY)*YINT
               XP =  (X-X0)*SINPHI + (Y-Y0)*COSPHI
               YP = -(X-X0)*COSPHI + (Y-Y0)*SINPHI
               RSQ = (XP/A)**2 + (YP/B)**2
               IF (TYPE.EQ.1) THEN
C                 -- gaussian
                  IF ((ABS(XP) .LE. A*1.6557) .AND.
     :                 ABS(YP) .LE. B*1.6557) THEN
                     AMP = TEMP*EXP(-FACT*RSQ)
                     SKY(IX,IY) = SKY(IX,IY)+AMP
                  END IF
               ELSE IF (TYPE.EQ.2) THEN
C                 -- ellipse
                  IF (RSQ.LT.0.25)
     :                 SKY(IX,IY) = SKY(IX,IY) + TEMP
               ELSE IF (TYPE.EQ.3) THEN
C                 -- rectangle
                  IF ((ABS(XP).LE.A/2.0) .AND. ABS(YP).LE.B/2.0) 
     :                 SKY(IX,IY) = SKY(IX,IY) + TEMP
               ELSE IF (TYPE.EQ.4) THEN
C                 -- SZE isothermal beta model
                  BETA = 2.0/3.0
                  EXPON = (1.0 - 3.0*BETA)/2.0
                  AMP = TEMP*(1.0 + RSQ)**EXPON
                  SKY(IX,IY) = SKY(IX,IY) + AMP
               END IF
            END DO
         END DO
 90      CONTINUE
      END DO
 100  CLOSE (UNIT=UNIT)

C Statistics of image

      SKYTOT = 0.0
      SKYMIN = SKY(1,1)
      SKYMAX = SKY(1,1)
      DO IY=1,NY
         DO IX=1,NX
            SKYMIN = MIN(SKYMIN,SKY(IX,IY))
            SKYMAX = MAX(SKYMAX,SKY(IX,IY))
            SKYTOT = SKYTOT + SKY(IX,IY)
         END DO
      END DO
      WRITE (*,'(1X,A,2I6)')   'Image size:                ', NX,NY
      WRITE (*,'(1X,A,1PE12.4)') 'Minimum pixel in image (K):', SKYMIN
      WRITE (*,'(1X,A,1PE12.4)') 'Maximum pixel in image (K):', SKYMAX
C
C-----------------------------------------------------------------------
C Write the image to disk in FITS format.
C-----------------------------------------------------------------------
      WRITE (*,*) 'Writing FITS file'
      CALL USERNM(OBSERV)
C
      IER = 0
      CALL FTINIT(37, FMAP, 1, IER)
      MSG = 'Can''t open file MAPFILE: '//FMAP
      IF (IER.GT.0) CALL ERROR(MSG(1:LEN1(MSG)))
      SIMPLE = .TRUE.
      BITPIX = -32
      NAXIS = 4
      NAXISN(1) = NX
      NAXISN(2) = NY
      NAXISN(3) = 1
      NAXISN(4) = 1
      EXTEND = .TRUE.
C
C     Write required keywords.
C
      CALL FTPHPR(37, SIMPLE, BITPIX, NAXIS, NAXISN, 0, 1, EXTEND, IER)
      CALL FTPSCL(37, 1.0D0, 0.0D0, IER)
C
C     Write additional keywords.
C
      CALL FTPKYS(37, 'OBJECT',  'Model', 'Source name', IER)
      CALL FTPKYS(37, 'TELESCOP','Simulation',' ',IER)
      CALL FTPKYS(37, 'INSTRUME','Simulation',' ',IER)
      CALL FTPKYS(37, 'OBSERVER', OBSERV,' ',IER)
      CALL FTGSDT(ID, IMON, IYR, IER)
      WRITE (DATOBS,'(I4.4,''-'',I2.2,''-'',I2.2)') IYR,IMON,ID
      CALL FTPKYS(37, 'DATE-OBS', DATOBS,' ',IER)
      CALL FTPKYS(37, 'BUNIT',    'K','Brightness units',
     :            IER)
      CALL FTPKYS(37, 'RADECSYS', 'FK5','Mean place', IER)
      CALL FTPKYF(37, 'EQUINOX', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      CALL FTPKYF(37, 'EPOCH', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      CALL FTPKYS(37, 'ORIGIN', PROGNM//' (Caltech)',
     :                'Version '//VERSN, IER)
      CALL FTPKYE(37, 'DATAMAX', SKYMAX, 6, ' ', IER)
      CALL FTPKYE(37, 'DATAMIN', SKYMIN, 6, ' ', IER)
C
C     Axis 1 - RA [shift and rotation ignored!]
C
      PIX1 = XC+1
      DEL1 = -YINT/60.0
      STRING = 'RA  = '
      CALL PANGLE('H', RA0, STRING(6:23), 4)
      CALL FTPKYS(37, 'CTYPE1',   'RA---SIN',
     1     '--- Axis 1: Right Ascension ---', IER)
      CALL FTPKYG(37, 'CRVAL1',   RA0/RPDEG, 6, STRING, IER)
      CALL FTPKYF(37, 'CRPIX1',   PIX1,      6, ' ', IER)
      CALL FTPKYF(37, 'CDELT1',   DEL1,      6, ' ', IER)
      CALL FTPKYF(37, 'CROTA1',   0.0,       6, ' ', IER)
C     
C     Axis 2 - DEC
C     
      PIX2 = YC+1
      DEL2 = -YINT/60.0
      STRING = 'Dec = '
      CALL PANGLE('D', DEC0, STRING(6:22), 3)
      CALL FTPKYS(37, 'CTYPE2',   'DEC--SIN',
     1     '--- Axis 2: Declination ---', IER)
      CALL FTPKYG(37, 'CRVAL2',   DEC0/RPDEG, 6, STRING, IER)
      CALL FTPKYF(37, 'CRPIX2',   PIX2,       6, ' ', IER)
      CALL FTPKYF(37, 'CDELT2',   DEL2,       6, ' ', IER)
      CALL FTPKYF(37, 'CROTA2',   0.0,        6, ' ', IER)
C     
C     Axis 3 - Frequency (degenerate).
C     
      CALL FTPKYS(37, 'CTYPE3',   'FREQ',
     1     '--- Axis 3: Frequency ---', IER)
      CALL FTPKYE(37, 'CRVAL3',   FREQ,  6,    ' ', IER)
      CALL FTPKYF(37, 'CRPIX3',   1.0,   6,    ' ', IER)
      CALL FTPKYF(37, 'CDELT3',   0.0,   6,    ' ', IER)
      CALL FTPKYF(37, 'CROTA3',   0.0,   6,    ' ', IER)
C     
C     Axis 4 - Stokes parameter (degenerate).
C     
      CALL FTPKYS(37, 'CTYPE4',   'STOKES',
     1     '--- Axis 4: Stokes parameter ---', IER)
      CALL FTPKYF(37, 'CRVAL4',   1.0,   6, 'Dirty map', IER)
      CALL FTPKYF(37, 'CRPIX4',   1.0,   6,    ' ', IER)
      CALL FTPKYF(37, 'CDELT4',   1.0,   6,    ' ', IER)
      CALL FTPKYF(37, 'CROTA4',   0.0,   6,    ' ', IER)
C     
C     History.
C     
      CALL HISS(37, IER)
C     
C     End FITS header; write FITS data.
C     
      CALL FTPPRE(37, 1, 1, NX*NY, SKY, IER)
      CALL FTCLOS(37, IER)
      IF (IER.GT.0) CALL ERROR(
     1     'Something went wrong writing FITS file.')

      END

      SUBROUTINE HISS(UNIT, IER)
      INTEGER  UNIT, IER
C
C HISS: generate FITS history records giving username, node name, 
C operating system etc.
C-----------------------------------------------------------------------
      CHARACTER*20 USER, NODE, VER
      CHARACTER*20 TIME
      INTEGER      LEN1
C
C Find user name, etc.
C
      CALL USERNM(USER)
      CALL NODENM(NODE)
      CALL SYSVER(VER)
      CALL DATTIM(TIME)
C
C Write output.
C
      CALL FTPHIS(UNIT, 
     1      'Run by:          '//USER(1:LEN1(USER))//' on '//
     2      NODE(1:LEN1(NODE))//' ('//VER(1:LEN1(VER))//')', IER)
      CALL FTPHIS(UNIT, 
     1      'Date:            '//TIME, IER)
      END

C*GTMODL -- read one model component from file
C+
      SUBROUTINE GTMODL(INMOD,ITYPE,P1,P2,P3,P4,P5,P6,P7,OUTMOD)
      INTEGER INMOD, ITYPE, OUTMOD
      REAL    P1, P2, P3, P4, P5, P6, P7
C
C T.J.Pearson   1979 FEB 18
C Modified for VAX version  1979 JULY 5
C Modified to allow comments 1981 February 9
C Modified to copy input comments to an output file 1982 September 22
C
C Read one model component (card) from unit 'INMOD'.
C Each card carries up to 7 numbers:
C 6 floating-point parameters P1-P6, and an optional integer ITYPE.
C Copy comments following a '!' character to unit 'OUTMOD'
C 	Ignore zero-flux components (or blank lines)
C-----------------------------------------------------------------------
      REAL PARS(8)
      DOUBLE PRECISION KCTOR
      CHARACTER*128 CARD
      CHARACTER*120 MSG
      CHARACTER BLANK, EXCL
      INTEGER I, J, L
      INTEGER LEN1
      CHARACTER*64 FNAME
      DATA BLANK/' '/
      DATA EXCL/'!'/
C
   5  READ(INMOD, '(A)', END=40) CARD
      L = 128
   10 IF (CARD(L:L).EQ.BLANK) THEN
            L = L-1
            IF (L.LT.1) GOTO 5
        GOTO 10
      END IF
C
      I = 1
      DO 30 J=1,8
          CALL KSKIPB(CARD,I)
          PARS(J) = KCTOR(CARD,I)
   30 CONTINUE
      IF (I.GT.L) GOTO 50
      IF (CARD(I:I).EQ.BLANK) GOTO 50
      IF (CARD(I:I).EQ.EXCL) THEN
            IF (OUTMOD.GT.0) WRITE (OUTMOD,'(A)') CARD(I:L)
            GOTO 50
      END IF
      INQUIRE(UNIT=INMOD,NAME=FNAME)
      MSG = 'Invalid model format: file '//FNAME
      CALL ERROR(MSG(1:LEN1(MSG)))
   50 CONTINUE
      ITYPE = PARS(1)
      P1 = PARS(2)
      IF (P1.EQ.0.0) GOTO 5
      P2 = PARS(3)
      P3 = PARS(4)
      P4 = PARS(5)
      P5 = PARS(6)
      P6 = PARS(7)
      P7 = PARS(8)
      RETURN
C
   40 ITYPE = -1
      RETURN
      END
