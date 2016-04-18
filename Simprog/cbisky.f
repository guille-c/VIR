      PROGRAM CBISKY
C-----------------------------------------------------------------------
C CBISKY: main program.
C 
C Author:            T.J. Pearson
C
C Version 1.0: 1996 Jan 28 - adaptation of INVERT (TJP).
C Version 1.0: 1996 Oct  9 - use FITSIO.
C Version 1.0: 1997 Jun  2 - add plotting.
C Version 1.0: 1997 Jul 24 - correct scaling.
C Version 2.0: 2000 Apr 27 - for San Pedro; omit plotting.
C Version 2.1: 2000 May 15 - add FACTOR parameter to allow the
C                            program to read both Cmbfast and Sugiyama
C                            files.
C Version 3.0: 2000 Dec  6 - adapted for g77.
C Version 3.1: 2001 Jan 22 - fix BUG in conjugation; 
C                            add SHIFT option (no noise).
C Version 3.2: 2001 Jan 24 - save subset of image.
C Version 4.2: 2001 Jan 29 - better integration of C_l over cell;
C                            put map right way up.
C-----------------------------------------------------------------------
      CHARACTER*(*) PROGNM,VERSN
      PARAMETER (PROGNM='CBISKY')
      PARAMETER (VERSN='4.0 - 2001 Jan 29')
      INTEGER   INC, OUTC, NPARS
      PARAMETER (INC=5, OUTC=6, NPARS=45)
      INTEGER MAX
      PARAMETER (MAX=2*2048*1025)
C
      REAL      WORK(MAX)
      INTEGER   I, IDIM, MODE, ICENT
      INTEGER   NWORDS
      INTEGER   VALID(8)
C 
C Input control parameter definitions.
C 
      DOUBLE PRECISION PARS(NPARS),VALS(NPARS),ENDMRK
      DATA  ENDMRK/ 8H/        /
      DATA  PARS  / 8HSEED    , 8HFREQUENC, 8HPBCOR   , 8HSOURCE  ,
     :              8H        , 8HUVMAP   , 8HMAPSIZE , 8HXYINT   ,
     :              8HMAPFile , 9*8H        ,
     :              8HCMBfile , 9*8H        ,
     :              8HPlot    , 9*8H        ,
     :              8HColumn  , 8HT0      , 8HFactor  , 8HInterp  ,
     :              2*8HShift   , 8HSize    /
      DATA  VALS  / 1459865D0,  31D0,       -1D0,       8HCMB     ,
     :              8H        , -1D0,       512D0,      60D0,
     :              8Hcbisky.m, 8Hap      , 8*8H        ,
     :              8Hcl.txt  , 9*8H        ,
     :              10*8H        ,
     :              1D0,        2.728D0,    1D0,        1D0,
     :              2*0D0,      0D0 /
      DATA VALID/16,32,64,128,256,512,1024,2048/
C-----------------------------------------------------------------------
C 
C Introduction.
C 
      CALL PUTOUT('Program '//PROGNM//
     :     ' -- makes maps with specified power spectrum')
      CALL PUTOUT('(Version '//VERSN//')')
      CALL PUTOUT(' ')
C 
C Read the control parameters.
C 
      CALL PUTOUT('Control parameters:')
      MODE = 0
      CALL KEYIN(PARS,VALS,NPARS,ENDMRK,MODE,INC,OUTC)
C 
C Evaluate the parameters used for memory allocation.
C 
      IDIM = VALS(7)
      DO 10 I=1,8
          IF (IDIM.EQ.VALID(I)) GOTO 20
   10 CONTINUE
      CALL ERROR('MAPSIZE must be a power of 2 between 16 and 2048')
   20 CONTINUE
      ICENT = IDIM/2+1
      NWORDS = ICENT*IDIM*2
      IF (NWORDS.GT.MAX) CALL ERROR('Insufficient memory')
C 
C Call a subroutine to do the work.
C 
      CALL INVER1(IDIM, ICENT, WORK, WORK, VALS, NPARS, OUTC,
     :     PROGNM, VERSN)
      CALL PUTOUT(PROGNM//' completed')
C-----------------------------------------------------------------------
      END
  
      SUBROUTINE INVER1(IDIM,ICENT,UVARR,MAPARR,
     1                  VALS, NPARS, LIST, PROGNM, VERSN)
      INTEGER NPARS
      DOUBLE PRECISION VALS(NPARS)
      INTEGER LIST
      CHARACTER*(*) PROGNM, VERSN
C-----------------------------------------------------------------------
      REAL       PI,RPDEG,RPSA
      PARAMETER  (PI=3.1415926536)
      PARAMETER  (RPDEG=PI/180.0)
      PARAMETER  (RPSA=PI/648000.0)
C 
C               uv/map array
C               NB UVARR and MAPARR are equivalent arrays
C 
      INTEGER IDIM            ! array size for FFT
      INTEGER ICENT           ! position of center = IDIM/2+1
      COMPLEX UVARR(ICENT,IDIM)
      REAL    MAPARR(IDIM,IDIM), MAPMAX,MAPMIN
      INTEGER NDIM(2)
C 
C               Other variables
C 
      CHARACTER*16 SNAME
      INTEGER NAXISN(4), BITPIX
      CHARACTER*80  FMAP,FUV,STRING
      CHARACTER*80  MESSAGE
      CHARACTER*10  DATOBS
      CHARACTER*12  OBSERV
      CHARACTER*80  CLFILE,TEXT1, TEXT2
      REAL*8        RA0, DEC0
      REAL*8        FACTOR
C
      LOGICAL SIMPLE, EXTEND, PBCOR, UVMAP
      LOGICAL INTERP, SHIFT
      REAL    SHIFTX, SHIFTY, PHASE, CPH, SPH
      INTEGER I, J, K, ID, IER, IMON, ISEED, ISEED0
      INTEGER IU, IV, IX, IY
      INTEGER M, MMAX, MMIN, N, NMAX, NMIN, NAXIS
      INTEGER LEN1
      INTEGER L, LMIN, LMAX, ICOL
      INTEGER KEEPX, KEEPY, I1, J1, P
      REAL    DEL1, DEL2, FREQ, PIX1
      REAL    PIX2, UVINT, XYINT, RVIS, IVIS, RSQR
      REAL    GASDEV
      REAL    EL(5000), CL(5000), X, Y(10), ELR, SIGMA
      REAL    PIX, SUM, SUMSQ, TOTVAR, SMOVAR, WL, T0, R, VAR
      REAL    SUMWT
      INTEGER MC, NPT
      REAL    VOFF(5), UOFF(5), WT(5)
      DATA    VOFF/ 0.0, -0.5, -0.5,  0.5, 0.5 /
      DATA    UOFF/ 0.0, -0.5,  0.5, -0.5, 0.5 /
      DATA    WT  / 8.0,  1.0,  1.0,  1.0, 1.0 /
C-----------------------------------------------------------------------
      NDIM(1) = IDIM
      NDIM(2) = IDIM
C-----------------------------------------------------------------------
C Decode control parameters.
C-----------------------------------------------------------------------
C     Random number seed; negative integer required to start
C     the Numerical Recipes ran1 generator.
C
      ISEED0 = NINT(VALS(1))
      IF (ISEED0.EQ.0) ISEED0 = -1
      IF (ISEED0.GT.0) ISEED0 = -ISEED0
      ISEED = ISEED0
C
C     Frequency, for header of output image and primary beam.
C
      FREQ  = VALS(2)*1E9
      IF (FREQ.LE.0D0) CALL ERROR('Bad value for FREQUENCY')
C
C     Output image file name (and name of uv image)
C
      WRITE (FMAP,610)  (VALS(I),I=9,18)
      FUV ='uv-'//FMAP
C
C     "Source" name, for header of output image.
C
      WRITE (SNAME,'(2A8)') VALS(4),VALS(5)
C 
C     Image pixel size XYINT (arcsec)
C 
      XYINT = VALS(8)
      IF (XYINT.LE.0D0) CALL ERROR('Bad value for XYINT')
C
C     Do Primary Beam multiplication?
C
      PBCOR = VALS(3).GE.0.0
C
C     Create sd image of uv plane?
C
      UVMAP = VALS(6).GE.0.0
C
C     C_l spectrum file; column number; scale factor.
C
      WRITE (CLFILE,'(10A8)') (VALS(I),I=19,28)
      ICOL = VALS(39)
      IF (ICOL.LT.1 .OR. ICOL.GT.10) ICOL = 1
      FACTOR = VALS(41)
C
C     Temperature normalization
C
      T0 = VALS(40)
C
C     Interpolation mode
C
      INTERP = VALS(42).GT.0D0
C
C     Shift parameter
C
      SHIFTX = VALS(43)
      SHIFTY = VALS(44)
      SHIFT = (SHIFTX .NE. 0D0 .OR. SHIFTY .NE. 0D0)
      IF (SHIFT)
     :   WRITE (*,'(1X,A,2F8.3)') 'Shift in pixels: ', SHIFTX, SHIFTY
      SHIFTX = SHIFTX/IDIM
      SHIFTY = SHIFTY/IDIM
C
C Size of image to be kept.
C
      IF (VALS(45).GT.1D0) THEN
         KEEPX = VALS(45)
         KEEPY = VALS(45)
      ELSE
         KEEPX = IDIM
         KEEPY = IDIM
      END IF
C
C-----------------------------------------------------------------------
C Read the power spectrum from  the file specified and convert to C_l.
C Sugiyama files list 1E10 l(l+1) C_l / (2 pi), with a different
C model in each column. Factor should be 1E10.
C CMBFAST files list l(l+1)C_l/2 pi for T, E, [B,] TE. Factor should be
C 1.
C-----------------------------------------------------------------------
      WRITE (*,'(/2A)') ' Reading file ',CLFILE(1:LEN1(CLFILE))
      WRITE (*,'(A,I5)') ' Column number ', ICOL
      I = 0
      TOTVAR = 0.0
      SMOVAR = 0.0
      OPEN (UNIT=1, NAME=CLFILE, STATUS='OLD', IOSTAT=IER)
      IF (IER.NE.0) THEN
         WRITE (*,*) 'Unable to read file: ', CLFILE
         STOP
      END IF
      DO K=1,5000
         READ (1,*,ERR=20,END=20) L,(Y(J),J=1,ICOL)
         I = I+1
         X = L
         EL(I) = X
         IF (I.EQ.1) WRITE (*,*) Y(ICOL)
         CL(I)  = (Y(ICOL)/FACTOR)*2.0*PI/(X*(X+1))
         TOTVAR = TOTVAR + CL(I)*(2*X+1)
         R = -5.491E-3*(X+0.5)**2
         IF (R.GT.-20.0) THEN
            WL = EXP(R)
            SMOVAR = SMOVAR + WL*CL(I)*(2*X+1)
         END IF
      END DO
      CLOSE (UNIT=1)
 20   CONTINUE
      TOTVAR = TOTVAR/(4.0*PI)
      SMOVAR = SMOVAR/(4.0*PI)
      LMIN = EL(1)
      LMAX = EL(I)
      WRITE (*,*) ' l-range: ',               LMIN, LMAX
      WRITE (*,*) ' Total variance:        ', TOTVAR
      WRITE (*,*) ' Expected sky rms (µK): ', SQRT(TOTVAR)*1E6*T0
      WRITE (*,*) ' RMS in 10° beam (µK):  ', SQRT(SMOVAR)*1E6*T0
C-----------------------------------------------------------------------
      CALL USERNM(OBSERV)
C-----------------------------------------------------------------------
C 
C               Choose uv interval etc.
C               xy dimensions in milliarcsec, uv in wavelengths
C 
      UVINT = 1.0/(IDIM*XYINT*RPSA)
      WRITE (LIST,620) XYINT,UVINT,IDIM
      WRITE (LIST,621) 2.0*PI*IDIM/2*UVINT
C-----------------------------------------------------------------------
C Fill sd array (for output as image)
C-----------------------------------------------------------------------
      IF (UVMAP) THEN
      WRITE (*,'(/A)') 'Computing uv sd image'
      VAR = 0.0
      DO IV=1,IDIM
         DO IU=1,IDIM
            RSQR = REAL((IV-ICENT))**2 + REAL((IU-ICENT))**2
            ELR = 2.0*PI*SQRT(RSQR)*UVINT - 0.5
            L = NINT(ELR)
            IF (L.GE.LMIN .AND. L.LT.LMAX) THEN
               SIGMA = SQRT(CL(L))*UVINT/SQRT(2.0)
               VAR = VAR + SIGMA*SIGMA
            ELSE 
               SIGMA = 0.0
            END IF
            MAPARR(IU,IV) = SIGMA
         END DO
      END DO
      WRITE (*,*) 'Total variance in uv plane: ', VAR
      WRITE (*,*) 'Variance per pixel:         ', VAR/(IDIM*IDIM)
C
C              Plot sd image
C
C      CALL SHOMAP(IDIM, IDIM, MAPARR, GDEV, 'UV standard deviation')
C
C              Write sd image
C
      write (*,*) 'Writing sd image as FITS file'
      IER = 0
      CALL FTINIT(37, FUV, 1, IER)
      IF (IER.GT.0) CALL ERROR('Can''t open file: '//FUV)
      SIMPLE = .TRUE.
      BITPIX = -32
      NAXIS = 2
      NAXISN(1) = IDIM
      NAXISN(2) = IDIM
      EXTEND = .FALSE.
C
C     Write required keywords.
C
      CALL FTPHPR(37, SIMPLE, BITPIX, NAXIS, NAXISN, 0, 1, EXTEND, IER)
      CALL FTPSCL(37, 1.0D0, 0.0D0, IER)
      IF (IER.GT.0) write (*,*) 'error 1'
C
C     Write additional keywords.
C
      CALL FTPKYS(37, 'OBJECT',   SNAME, 'Source name', IER)
      CALL FTPKYS(37, 'TELESCOP','Simulation',' ',IER)
      CALL FTPKYS(37, 'INSTRUME','Simulation',' ',IER)
      CALL FTPKYS(37, 'OBSERVER', OBSERV,' ',IER)
      CALL FTGSDT(ID, IMON, IY, IER)
      WRITE (DATOBS,'(I4.4,''-'',I2.2,''-'',I2.2)') IY,IMON,ID
      CALL FTPKYS(37, 'DATE-OBS', DATOBS,' ',IER)
      CALL FTPKYS(37, 'BUNIT',    'JY/BEAM','Brightness units',
     :            IER)
      CALL FTPKYS(37, 'RADECSYS', 'FK5','Mean place', IER)
      CALL FTPKYF(37, 'EQUINOX', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      CALL FTPKYF(37, 'EPOCH', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      TEXT1 = PROGNM//' (Caltech)'
      TEXT2 = 'Version '//VERSN
      CALL FTPKYS(37, 'ORIGIN', TEXT1, TEXT2, IER)
      CALL FTPDAT(37, IER)
      IF (IER.GT.0) write (*,*) 'error 2'
C
C     Axis 1 - U
C
      PIX1 = ICENT
      DEL1 = -UVINT
      CALL FTPKYS(37, 'CTYPE1',   'UU---SIN',
     1     '--- Axis 1: U coordinate ---', IER)
      CALL FTPKYF(37, 'CRVAL1',   0.0,       6, ' ', IER)
      CALL FTPKYF(37, 'CRPIX1',   PIX1,      6, ' ', IER)
      CALL FTPKYF(37, 'CDELT1',   DEL1,      6, ' ', IER)
      CALL FTPKYF(37, 'CROTA1',   0.0,       6, ' ', IER)
      IF (IER.GT.0) write (*,*) 'error 3'
C     
C     Axis 2 - V
C     
      PIX2 = ICENT
      DEL2 = UVINT
      CALL FTPKYS(37, 'CTYPE2',   'VV---SIN',
     1     '--- Axis 2: V coordinate ---', IER)
      CALL FTPKYF(37, 'CRVAL2',   0.0,        6, ' ', IER)
      CALL FTPKYF(37, 'CRPIX2',   PIX2,       6, ' ', IER)
      CALL FTPKYF(37, 'CDELT2',   DEL2,       6, ' ', IER)
      CALL FTPKYF(37, 'CROTA2',   0.0,        6, ' ', IER)
      IF (IER.GT.0) write (*,*) 'error 4'
C     
C     History.
C     
      DO I=1,70
         MESSAGE(I:I) = '-'
      END DO
      TEXT1  = 'Spectrum file: '//CLFILE(1:LEN1(CLFILE))
      CALL FTPHIS(37, TEXT1, IER)
      WRITE (MESSAGE, '(''Seed: '',I12)') ISEED0
      CALL FTPHIS(37,  MESSAGE(1:70), IER)
      TEXT1 = 'Program:  '//PROGNM
      CALL FTPHIS(37, TEXT1, IER)
      TEXT1 = 'Version:  '//VERSN
      CALL FTPHIS(37, TEXT1, IER)
      CALL HISS(37, IER)
      IF (IER.GT.0) write (*,*) 'error 5'
C     
C     End FITS header; write FITS data.
C     
      CALL FTPPRE(37, 1, 1, IDIM*IDIM, MAPARR, IER)
      IF (IER.GT.0) write (*,*) 'error 6'
      CALL FTCLOS(37, IER)
      IF (IER.GT.0) CALL ERROR(
     1     'Something went wrong writing FITS file.')
      END IF ! (UVMAP)
C 
C               Fill visibility array
C
      WRITE (*,*) 'Computing visibilities according ',
     :     'to supplied power spectrum'
      IF (INTERP) THEN
         NPT = 5
         WRITE (*,*) 'Five-point cell average'
      ELSE
         NPT = 1
         WRITE (*,*) 'No cell averaging'
      END IF
      SUMWT = 0.0
      DO MC=1,NPT
         SUMWT = SUMWT+WT(MC)
      END DO
      DO 40 IV=1,IDIM
      DO 40 IU=1,ICENT
C        -- evaluate C_l at center and corners of cell and make
C        -- weighted average
         TOTVAR = 0.0
         DO MC=1,NPT
            RSQR = REAL((IV-ICENT+VOFF(MC)))**2 +
     :             REAL((IU-ICENT+UOFF(MC)))**2
            ELR = 2.0*PI*SQRT(RSQR)*UVINT - 0.5
            L = NINT(ELR)
            IF (L.GE.LMIN .AND. L.LT.LMAX) THEN
               VAR = CL(L)
            ELSE
               VAR = 0.0
            END IF
            TOTVAR = TOTVAR + WT(MC)*VAR
         END DO
C        -- Use the variance to generate random numbers
         IF (TOTVAR.GT.0.0) THEN
            SIGMA = SQRT(TOTVAR/SUMWT)*UVINT/SQRT(2.0)
            RVIS = SIGMA*GASDEV(ISEED)
            IVIS = SIGMA*GASDEV(ISEED)
         ELSE
            RVIS = 0.0
            IVIS = 0.0
         END IF
         IF (SHIFT) THEN
            PHASE = (2.0*PI)*(SHIFTX*(IU-ICENT) + SHIFTY*(IV-ICENT))
            CPH = COS(PHASE)
            SPH = SIN(PHASE)
            UVARR(IU,IV) = CMPLX(RVIS*CPH + IVIS*SPH,
     :                           IVIS*CPH - RVIS*SPH)
         ELSE
            UVARR(IU,IV) = CMPLX(RVIS, IVIS)
         END IF
 40   CONTINUE
C
C     Fix up central row for conjugate symmetry
C
      UVARR(ICENT, ICENT) = (0.0,0.0)
      DO 45 IV=ICENT+1,IDIM
C        write (*,*) iv, uvarr(icent, iv), UVARR(ICENT,ICENT-(IV-ICENT))
         UVARR(ICENT,IV) = CONJG(UVARR(ICENT,ICENT-(IV-ICENT)))
 45   CONTINUE
C     WRITE (*,*) 'Central uv pixel: ', UVARR(ICENT,ICENT)
C 
C               Do the Fourier Transform.
C               FOUR2 is Norman Brenner's routine, here used for
C               transforming conjugate symmetric complex data to real.
C               The sign changes shift the phase center from point
C               (1,1) to (ICENT,ICENT).
C 
      WRITE (*,*) 'Fourier transforming to sky plane'
      DO IV=1,IDIM,2
         DO IU=1,ICENT,2
            UVARR(IU,IV) = -UVARR(IU,IV)
         END DO
         DO IU=2,ICENT,2
            UVARR(IU,IV+1) = -UVARR(IU,IV+1)
         END DO
      END DO
      CALL FOUR2(UVARR,NDIM,2,-1,-1)
      DO IY=1,IDIM,2
         DO IX=1,IDIM,2
            MAPARR(IX,IY) = -MAPARR(IX,IY)
         END DO
         DO IX=2,IDIM,2
            MAPARR(IX,IY+1) = -MAPARR(IX,IY+1)
         END DO
      END DO
C
C     Multiply by Primary Beam.
C
      IF (PBCOR) THEN
         WRITE (*,*) 'Applying primary beam'
         CALL PBEAM(MAPARR,IDIM, IDIM, REAL(ICENT), REAL(ICENT),
     :        XYINT/60., XYINT/60.)
      END IF
C
C     Normalize with T0, and locate maximum and minimum
C 
      MAPMAX = -1E37
      MAPMIN =  1E37
      SUM = 0.0
      SUMSQ = 0.0
      DO N=1,IDIM
         DO M=1,IDIM
            MAPARR(M,N) = MAPARR(M,N)*T0
            SUM = SUM+MAPARR(M,N)
            SUMSQ = SUMSQ+MAPARR(M,N)**2
            IF (MAPARR(M,N).GT.MAPMAX) THEN
               MAPMAX = MAPARR(M,N)
               NMAX = N
               MMAX = M
            END IF
            IF (MAPARR(M,N).LT.MAPMIN) THEN
               MAPMIN = MAPARR(M,N)
               NMIN = N
               MMIN = M
            END IF
         END DO
      END DO
      WRITE (LIST,660) ICENT,ICENT,
     1     MAPMAX,MMAX,NMAX, MAPMIN,MMIN,NMIN
      PIX = REAL(IDIM)*REAL(IDIM)
      WRITE (LIST,661) IDIM*IDIM, SUM/PIX, SQRT(SUMSQ/PIX)
C
C Write clean map to disk in FITS format.
C
      WRITE (*,*) 'Writing FITS file'
      IER = 0
      CALL FTINIT(37, FMAP, 1, IER)
      IF (IER.GT.0) CALL ERROR('Can''t open file MAPFILE: '//FMAP)
      SIMPLE = .TRUE.
      BITPIX = -32
      NAXIS = 4
      NAXISN(1) = KEEPX
      NAXISN(2) = KEEPY
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
      CALL FTPKYS(37, 'OBJECT',   SNAME, 'Source name', IER)
      CALL FTPKYS(37, 'TELESCOP','Simulation',' ',IER)
      CALL FTPKYS(37, 'INSTRUME','Simulation',' ',IER)
      CALL FTPKYS(37, 'OBSERVER', OBSERV,' ',IER)
      CALL FTGSDT(ID, IMON, IY, IER)
      WRITE (DATOBS,'(I4.4,''-'',I2.2,''-'',I2.2)') IY,IMON,ID
      CALL FTPKYS(37, 'DATE-OBS', DATOBS,' ',IER)
      CALL FTPKYS(37, 'BUNIT',    'K','Brightness units',
     :            IER)
      CALL FTPKYS(37, 'RADECSYS', 'FK5','Mean place', IER)
      CALL FTPKYF(37, 'EQUINOX', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      CALL FTPKYF(37, 'EPOCH', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', IER)
      TEXT1 = PROGNM//' (Caltech)'
      TEXT2 = 'Version '//VERSN
      CALL FTPKYS(37, 'ORIGIN', TEXT1, TEXT2, IER)
      CALL FTPDAT(37, IER)
      CALL FTPKYE(37, 'DATAMAX', MAPMAX, 6, ' ', IER)
      CALL FTPKYE(37, 'DATAMIN', MAPMIN, 6, ' ', IER)
C
C     Axis 1 - RA [shift and rotation ignored!]
C
      PIX1 = KEEPX/2+1
      DEL1 = -XYINT/3600.0
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
      PIX2 = KEEPY/2+1
      DEL2 = XYINT/3600.0
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
      DO I=1,70
         MESSAGE(I:I) = '-'
      END DO
      TEXT1 = 'Spectrum file: '//CLFILE(1:LEN1(CLFILE))
      CALL FTPHIS(37, TEXT1, IER)
      WRITE (MESSAGE, '(''Column: '',I12)') ICOL
      CALL FTPHIS(37, MESSAGE(1:70), IER)
      WRITE (MESSAGE, '(''Seed: '',I12)') ISEED0
      CALL FTPHIS(37, MESSAGE(1:70), IER)
      TEXT1 = 'Program:  '//PROGNM
      CALL FTPHIS(37, TEXT1, IER)
      TEXT1 = 'Version:  '//VERSN
      CALL FTPHIS(37, TEXT1, IER)
      CALL HISS(37, IER)
C     
C     End FITS header; write FITS data.
C     
C      CALL FTPPRE(37, 1, 1, IDIM*IDIM, MAPARR, IER)
      I1 = IDIM/2 - KEEPX/2
      J1 = IDIM/2 - KEEPY/2
      P = 1
      DO J=J1+1, J1+KEEPY
         CALL FTPPRE(37, 1, P, KEEPX, MAPARR(I1,J), IER)
         P = P+KEEPX
      END DO
      CALL FTCLOS(37, IER)
      IF (IER.GT.0) CALL ERROR(
     1     'Something went wrong writing FITS file.')
C 
C Display map.
C
C      CALL SHOMAP(IDIM, IDIM, MAPARR, GDEV, 'Map')
      RETURN
C-----------------------------------------------------------------------
C       Format statements
C-----------------------------------------------------------------------
 610  FORMAT(10A8)
 620  FORMAT(/' Cellsize (XYINT):  ',F10.3,' arcsec'/
     1     ' u,v interval (UVINT):  ',F16.2,' wavelengths'/
     2     ' Array size for Fourier transform: ',I6)
 621  FORMAT(' Maximum sampled l: ',F8.1)
 660  FORMAT(4X,'Center at ',I5,',',I5/
     :       4X,'Maximum ',1PE12.3,'  in cell',2I5/
     :       4X,'Minimum ',1PE12.3,'  in cell',2I5)
 661  FORMAT(4X,'Number of pixels:', I8/
     :       4X,'Mean:',G12.4/
     :       4X,'rms: ',G12.4)
      END

      SUBROUTINE PBEAM(A, NX, NY, CX, CY, DELX, DELY)
      INTEGER NX, NY
      REAL A(NX,NY)
      REAL CX, CY, DELX, DELY
C
C Multiply map A(1.NX,1..NY) by primary beam. The pointing center is
C at pixel (CX,CY), possibly fractional; the grid spacing is (DELX,DELY)
C in arcmin.
C-----------------------------------------------------------------------
      INTEGER IX, IY
      REAL DX, DY, F
      REAL RS, RP, RM
      INTEGER N, M
      PARAMETER (N=241, M=201)
      REAL X(N), Y(N), XMIN, XMAX, YNORM
      REAL XPB(M), YPB(M), XX, R
      REAL BESSJ0
      INTEGER I, J, K
C
C Compute primary beam.
C
      RS = 15.5/2.0
      RP = 85.0/2.0
      RM = 90.0/2.0
      XMIN =   0.0
      XMAX =  50.0
C
C Compute aperture grading.
C
      DO I=1,N
         X(I) = XMIN + REAL(I-1)*(XMAX-XMIN)/REAL(N-1)
         IF (ABS(X(I)).GT.RM .OR. ABS(X(I)).LT.RS) THEN
            Y(I) = 0.0
         ELSE IF (RP.EQ.0.0) THEN
            Y(I) = 1.0
         ELSE
            Y(I) = EXP(-(1.15*X(I)/RP)**2)
         END IF
      END DO
C
C Do Hankel transform to compute primary beam.
C
      DO J=1,M
C        -- radius in arcmin
         XPB(J) = REAL(J-1)*100.0/REAL(M-1)
C        -- and radians
         XX = XPB(J)*3.141592653/(180.0*60.0)
         YPB(J) = 0.0
         DO I=1,N
            YPB(J) = YPB(J) + 
     :           Y(I)*X(I)*BESSJ0(2.0*3.141592653*XX*X(I))
         END DO
      END DO
C
C Normalize and square it.
C
      YNORM = YPB(1)
      DO J=1,M
         YPB(J) = (YPB(J)/YNORM)**2
      END DO
      WRITE (*,*) (YPB(J),J=1,M)
C
C Apply correction.
C
      DO IY=1,NY
         DY = ((IY-CY)*DELY)**2
         DO IX=1,NX
            DX = ((IX-CX)*DELX)**2
            R = SQRT(DX + DY)
            K = 1 + INT(R*REAL(M-1)/100.0)
            IF (K.GT.M-1) THEN
               F = 0.0
            ELSE
               CALL INTERP(XPB(K), XPB(K+1), YPB(K), YPB(K+1), R, F)
            END IF
C           IF (F.GT.0.0) WRITE (*,*) R, K, F
            A(IX,IY) = F*A(IX,IY)
         END DO
      END DO
      RETURN
      END

      SUBROUTINE INTERP(X1, X2, Y1, Y2, X, Y)
      REAL X1, X2, Y1, Y2, X, Y
C
C Linear interpolation between (X1,Y1) and (X2,Y2): return Y(X).
C
      REAL A, B
      A = (X-X1)/(X2-X1)
      B = (X-X2)/(X2-X1)
      Y = A*Y2 - B*Y1
      IF (ABS(A).GT.1.0 .OR. ABS(B).GT.1.0) 
     :     WRITE (*,*) 'INTERP failure'
      RETURN
      END

      SUBROUTINE HISS(UNIT, IER)
      INTEGER  UNIT, IER
C
C HISS: generate FITS history records giving username, node name, 
C operating system etc.
C-----------------------------------------------------------------------
      CHARACTER*20 USER, NODE, VER
      CHARACTER*20 TIME
      CHARACTER*80 TEXT
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
      TEXT = 'Run by: '//USER(1:LEN1(USER))//' on '//
     :      NODE(1:LEN1(NODE))//' ('//VER(1:LEN1(VER))//')' 
      CALL FTPHIS(UNIT, TEXT, IER)
      TEXT = 'Date: '//TIME 
      CALL FTPHIS(UNIT, TEXT, IER)
      END
