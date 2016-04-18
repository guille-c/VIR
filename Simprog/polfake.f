      PROGRAM POLFAK
C-----------------------------------------------------------------------
      CHARACTER*(*)  PROGNM, VERSN
      PARAMETER  (PROGNM='POLFAKE')
      PARAMETER  (VERSN='2002 Apr 8')
C
C 2000 Apr 28: make sure random number seed is negative.
C 2000 Apr 29: add integration time parameter in output FITS file.
C 2000 Apr 29: check at least 2 antennas.
C 2000 Apr 29: make Stokes parameter a user input (default LL).
C 2000 May 05: add NIGHTS parameter for simulating averaged data.
C 2000 May 05: create log file.
C 2000 May 06: allow random spectral indices.
C 2000 May 12: use fitsio to read sky image.
C 2000 May 15: shift UT start on each day for mosaic.
C 2000 Jul 15: add FGSI parameter.
C 2000 Nov 13: fix bug in Start_PA (wrong value used).
C 2000 Nov 13: add PALIST option.
C 2000 Dec  7: mods for g77.
C 2001 Jan 30: use bilinear interpolation.
C 2001 Feb 16: change flux densities from Jy to mJy in log file.
C 2001 Feb 20: fix check on size of input image.
C 2001 Feb 22: changed to use (faster) real->conjsym FFT instead of
C              complex FFT; max image size 2048.
C 2001 Feb 26: fix bug in GTMODL: ignored sources with RA=0.
C 2001 Feb 28: add option of foreground image (FOREGROUND, FGSI)
C 2001 Feb 28: change antenna names to RX0, RX1...
C 2001 Mar  5: correct a precision problem.
C 2001 Mar 14: allow negative point sources; add MINFLUX.
C 2001 Jun 22: increase MAXFLD from 100 to 200.
C 2001 Jul 16: add CBI beam option.
C 2001 Oct 25: add polarization (cbifake->polfake).
C 2001 Nov 08: save all Stokes combinations.
C 2002 Apr 08: fix bug in polarization assignment.
C
      INTEGER    MAXS, MAXF, NPARS, NXYZ, MAXDS
C     ! maximum number of stations
      PARAMETER  (MAXS=40)
C     ! maximum number of frequency channels
      PARAMETER  (MAXF=20)
C     ! maximum number of discrete sources
      PARAMETER  (MAXDS=5000)
C     ! number of input parameters
      PARAMETER  (NPARS=112+4*MAXS)
      PARAMETER  (NXYZ=4*MAXS)
C     ! maximum number of uv points
      INTEGER     KMAX
      PARAMETER  (KMAX=200000)
C     ! units for input, output, etc
      INTEGER    INC, OUTC, LOGUN, UNIT1, UNIT2, UNIT3, UNIT4
      PARAMETER  (INC=5)
      PARAMETER  (OUTC=6)
      PARAMETER  (LOGUN=7)
      PARAMETER  (UNIT1=17)
      PARAMETER  (UNIT2=18)
      PARAMETER  (UNIT3=47)
      PARAMETER  (UNIT4=18)
      INTEGER MAXFLD
      PARAMETER  (MAXFLD=200)
C     ! Maximum number of entries in PALIST
      INTEGER MAXPAL
      PARAMETER (MAXPAL=100)
      DOUBLE PRECISION    PI,RPDEG,RPHR,RPAMIN
      REAL       TWOPI
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (RPDEG=PI/180D0)
      PARAMETER (RPHR=PI/12D0)
      PARAMETER (TWOPI=2.0*PI)
      PARAMETER (RPAMIN=RPDEG/60D0)
C
      REAL CLIGHT
      PARAMETER (CLIGHT=2.997924580E8)
C
      INCLUDE 'polfake.inc'
C
      CHARACTER*16 FIELD, TPAVAL
      CHARACTER*10 RDSTR
      CHARACTER*8  STNAME(MAXS)
      CHARACTER*8  TIMEC
      CHARACTER*80 PREFIX, OUTDSN, DSDSN, CMBDSN, FGDSN, LOGFIL, PALFIL
      CHARACTER*80 CMBUDS, CMBQDS
      CHARACTER*8  BEAM
      DOUBLE PRECISION AZ, ZA, SQSZ, CQSZ, PA
      DOUBLE PRECISION FREQ(MAXF), BW(MAXF), REFDAT, RMJD, MJD, FD
      DOUBLE PRECISION LAST, HA, NIGHTS
      DOUBLE PRECISION OBSLAT, OBSLNG, ZALIM
      DOUBLE PRECISION PARS(NPARS), VALS(NPARS), ENDMRK, RPOL, LPOL
      DOUBLE PRECISION RA0D, DEC0D, PTRA(MAXFLD), PTDEC(MAXFLD)
      REAL             PTX(MAXFLD), PTY(MAXFLD)
      DOUBLE PRECISION SDECD, SRAD, GAST, APPRA, APPDEC, SLAPRM(21)
      DOUBLE PRECISION SH, CH, SD, CD, SP, CP, XX, YY, ZZ, R
      DOUBLE PRECISION SLA_GMST, EQEQX, SLA_EQEQX, SLA_AIRMAS
      DOUBLE PRECISION STX(MAXS), STY(MAXS), STZ(MAXS)
      DOUBLE PRECISION P1, P2, P3, P4, SPAC
      INTEGER BL(KMAX), POLZN(KMAX), UT(KMAX), T, MDATE
      INTEGER I, J, K, N, MODE, IF, NF, NSAMP, STATUS, IER, NSTART
      INTEGER RYEAR, RMONTH, RDAY, UTSTAR, UTSTOP, INTEG, NINTEG
      INTEGER PAC2, PAC3, IS, ISEED, SSEED, NDS
      INTEGER NPX, NPY, IFIELD, NFIELD
      INTEGER GTMODL, LEN1
      LOGICAL TRACK, IDLE, PALIST
      INTEGER ANTPOL(MAXS)
      INTEGER PAL, NPAL, PALCNT(MAXPAL)
      INTEGER IPOL, NPOL, KV
      REAL    PHASE, DELAY, PAC1, RUV
      REAL    THETA, U(KMAX), V(KMAX), PDATA(7), RDATA(3*4*MAXF)
      REAL    PASTAR, DELTPA, FWHM, SIGMA, ATTEN, RELF
      REAL    VISR(KMAX,MAXF,4), VISI(KMAX,MAXF,4), AIRMAS(KMAX), AM
      REAL    VISWT(KMAX,4)
      REAL    GASDEV
      REAL    TSYS, TATM, ANTARE, ETA, CHANBW, NFAC, RMS, SRCFAC
      REAL    SRCCUT, FL(MAXF), MAXFL, FLXLIM
      REAL    DSFLUX(MAXDS), DSSI(MAXDS), DSL(MAXDS), DSM(MAXDS)
      REAL    RAD(MAXDS)
      REAL    UU, VV
      COMPLEX VIS, QVIS, UVIS
      DOUBLE PRECISION DSRA(MAXDS), DSDEC(MAXDS)
      REAL    SIMEAN, SISD
      REAL    FGSI
      REAL    PALANG(MAXPAL)
      LOGICAL CBIBM
      REAL    CBIBEAM
      INTEGER COUNT(4)
      CHARACTER*4 FLDNAM(MAXFLD)
      LOGICAL ALLPOL
C
      DATA ENDMRK/8H/       /
      DATA PARS  /8HFitsfile, 9*8H        , 3*8HCycle_PA,
     :            8HN_ant   ,   8HTrack_PA,   8H        ,
     :          3*8HFrequenc,
     :            8HLatitude,   8HLongitud,   8HZAlimit ,
     :          2*8HAllPolar,   8HDeclinat,   8HRA      ,
     :            8HField_na,   8Hme      ,
     :            8HDate    , 2*8HUT_range,   8HIntegrat,
     :            8HStart_PA,   8HTsys    ,   8HTatm    ,
     :            8HArea    ,   8HEfficien,   8HBandwidt,
     :            8HSources , 9*8H        ,   8HFWHM    ,
     :          3*8HMosaic  ,   8HSFactor ,   8HSCutoff ,
     :            8HSeed    ,   8HCmbfile , 9*8H        ,
     :            8Hjunk    ,   8HNights  , 2*8HSpecInd ,
     :            8HFGSI    ,   8HPAlist  , 9*8H        ,
     :            8HForegrou, 9*8H        ,   8HMinflux ,
     :            8HBeam    ,
     :            8HQmap    , 9*8H        ,
     :            8HUmap    , 9*8H        ,
     :       NXYZ*8HXYZ_ant /
      DATA VALS / 8H?       , 9*8H        ,
     1     0D0, 1D0, 0D0,
     2     0D0,
     3     8HNo      , 8H        ,
     4     1D0, 26.5D0, 1D0,
     5     -23D0, 67.5D0, 45.0D0, 2*-1D0,
     6     -45D0, 0D0, 8HField-1 , 8H        ,
     7     0D0, 0D0, 60D0, 10D0,
     8     0D0,
     9     16D0, 3D0, 1.8D0, 0.95D0, 1.0D0,
     A     10*8H        , 43.93D0, 1D0, 1D0, 20D0, 1D0,
     B     1D5, -7654387D0, 10*8H        , -2D0, 1D0,
     C     -10D0, 0D0, 2D0, 10*8H        , 10*8H        ,
     :     1D-6,  8HGAUSS   , 20*8H        , NXYZ*0D0/
      DATA RPOL/8HR       /, LPOL/8HL       /
C-----------------------------------------------------------------------
C Read control parameters
C-----------------------------------------------------------------------
C
      CALL PUTOUT('Program '//PROGNM//' (Version '//VERSN//')')
    5 MODE = 0
      CALL PUTOUT(' ')
      CALL PUTOUT('Control parameters (end with /, EOF to finish):')
      VALS(29) = 3600*1998+60*4+18
      CALL KEYIN(PARS,VALS,NPARS,ENDMRK,MODE,INC,OUTC)
      IF (MODE.EQ.1) STOP
C
C Output file name
      WRITE (PREFIX, '(10A8)') (VALS(I),I=1,10)
C Logfile
      LOGFIL = PREFIX(1:LEN1(PREFIX))//'.log'
      OPEN (UNIT=LOGUN, STATUS='UNKNOWN', NAME=LOGFIL)
      WRITE (LOGUN,'(5A)') 'Program ',PROGNM,' (Version ',VERSN,')'
C Channel frequencies (Hz)
      NF      = MIN(MAXF,NINT(VALS(17)))
      WRITE (LOGUN,'(/1X,''Channel'',2X,''Center Frequency (GHz)'')')
      DO IF=1,NF
         FREQ(IF) = (VALS(18) + (IF-1)*VALS(19))*1E9
         BW(IF) = VALS(19)*1E9
         WRITE (LOGUN,'(1X,I7,2X,F12.4)') IF, FREQ(IF)*1E-9
      END DO
C Observatory location and limits
      OBSLAT = VALS(20)*RPDEG
      OBSLNG = VALS(21)*RPDEG
      ZALIM  = VALS(22)*RPDEG
C Do all polarizations?
      ALLPOL = VALS(23).NE.-1
C Field center (J2000) and name
      DEC0D  = VALS(25)/3600.0
      RA0D   = VALS(26)*15.0/3600.0
      WRITE (FIELD, '(2A8)') VALS(27), VALS(28)
C Reference date
      MDATE  = VALS(29)
      RYEAR  = MDATE/3600
      RMONTH = MOD(MDATE,3600)/60
      RDAY   = MOD(MDATE,60)
      IF (RYEAR.LE.50) THEN
         RYEAR = RYEAR + 2000
      ELSE IF (RYEAR.LT.100) THEN
         RYEAR = RYEAR + 1900
      END IF
      CALL SLA_CLDJ(RYEAR, RMONTH, RDAY, RMJD, IER)
      IF (IER.NE.0) CALL ERROR ('Invalid reference date')
C UT range (seconds from 0h UT on ref date)
      UTSTAR = VALS(30)
      UTSTOP = VALS(31)
C Integration time
      INTEG  = VALS(32)
C Nights: factor to increase effective integration time
      NIGHTS = VALS(67)
C Number of integrations
      NINTEG = (UTSTOP-UTSTAR)/INTEG
C Track PA? Start_PA
      WRITE (TPAVAL, '(2A8)') VALS(15), VALS(16)
      TRACK = TPAVAL(1:1).EQ.'y' .OR. TPAVAL(1:1).EQ.'Y'
      PASTAR = VALS(33)*RPDEG
      PAC1 = VALS(11)*RPDEG
      PAC2 = VALS(12)
      PAC3 = VALS(13)
C Noise characteristics
      TSYS =    VALS(34)
      TATM =    VALS(35)
      ANTARE =  VALS(36)
      ETA =     VALS(37)
      CHANBW =  VALS(38)*1E9
      NFAC = SQRT(2.0)*1.38E3/(ANTARE*ETA*SQRT(CHANBW*INTEG*NIGHTS))
      WRITE (LOGUN, *) 'Noise per integration at zenith (Jy): ', 
     :    NFAC*(TSYS+TATM)
C Mean and sd of the source spectral index distribution
      SIMEAN = VALS(68)
      SISD   = VALS(69)
C Spectral index of foreground image
      FGSI  = VALS(70)
C PA list file name
      WRITE (PALFIL , '(10A8)') (VALS(I),I=71,80)
C Antenna locations and polarizations
      N       = MIN(MAXS,NINT(VALS(14)))
      DO I=1,N
         IF (I.LE.10) THEN
            WRITE (STNAME(I), '(''RX'',I1)') I-1
         ELSE
            WRITE (STNAME(I), '(''RX'',I2.2)') I-1
         END IF
         STX(I) = VALS(112+4*I-3)
         STY(I) = VALS(112+4*I-2)
         STZ(I) = VALS(112+4*I-1)
         IF (VALS(112+4*I).EQ.LPOL) THEN
            ANTPOL(I) = 0
         ELSE IF (VALS(112+4*I).EQ.RPOL) THEN
            ANTPOL(I) = 1
         ELSE
            WRITE (*,*) 'Invalid polarization code'
            STOP
         END IF
         write (*,*) i, stx(i), sty(i), stz(i), antpol(i)
      END DO
C File name for source list, scale factor, upper cutoff, lower cutoff.
      WRITE (DSDSN, '(10A8)') (VALS(I),I=39,48)
      SRCFAC = VALS(53)
      SRCCUT = VALS(54)
      FLXLIM = VALS(91)
C Primary beam
      WRITE (BEAM, '(A8)') VALS(92)
      IF (BEAM.EQ.'CBI') THEN
         WRITE (*, '(1X,A)') 'Using CBI primary beam'
         CBIBM = .TRUE.
      ELSE IF (BEAM.EQ.'GAUSS') THEN
         FWHM = VALS(49)
         SIGMA = FWHM/1.66511
         WRITE (*, '(1X,A,F5.1,A)')
     :        'Using gaussian primary beam of FWHM',  FWHM, 'arcmin'
         CBIBM = .FALSE.
      ELSE
         CALL ERROR('Invalid option for "BEAM" parameter')
      END IF
C Mosaic parameters
      NPX = ABS(NINT(VALS(50)))
      NPY = ABS(NINT(VALS(51)))
      SPAC = VALS(52)
C Random number seed (must be negative to initialize Numerical Recipes
C ran1 generator)
      ISEED = NINT(VALS(55))
      IF (ISEED.GT.0) ISEED = -ISEED
      IF (ISEED.EQ.0) ISEED = -1
C File names for CMB sky maps
      WRITE (CMBDSN, '(10A8)') (VALS(I),I=56,65)
      WRITE (CMBQDS, '(10A8)') (VALS(I),I=93,102)
      WRITE (CMBUDS, '(10A8)') (VALS(I),I=103,112)
      DOCMB = .FALSE.
C File name for Foreground sky map
      WRITE (FGDSN, '(10A8)') (VALS(I),I=81,90)
      DOFG = .FALSE.
C-----------------------------------------------------------------------
C Check the number of antennas.
C-----------------------------------------------------------------------
      IF (N.LT.2) CALL ERROR('There should be at least 2 antennas.')
C-----------------------------------------------------------------------
C Read the PA list file.
C-----------------------------------------------------------------------
      PALIST = PALFIL .NE. ' '
      IF (PALIST) THEN
         WRITE (*,*) 'Reading ', PALFIL
         OPEN (UNIT=UNIT4, STATUS='OLD', ERR=511, FILE=PALFIL)
         DO I=1,MAXPAL
            READ (UNIT4,*,ERR=511, END=510) PALCNT(I), PALANG(I)
            NPAL = NPAL+1
         END DO
 510     CLOSE (UNIT=UNIT4)
         GOTO 512
 511     CALL ERROR('Problem reading PALIST file')
 512     CONTINUE
         IF (NPAL.LT.1) THEN
            PALIST = .FALSE.
         ELSE
            WRITE (*,*) NPAL, ' entries in PALIST file'
            DO I=1,NPAL
               WRITE (*,'(I8,F10.1)') PALCNT(I), PALANG(I)
            END DO
         END IF
      END IF
C-----------------------------------------------------------------------
C Read in the CMB images (if any).
C-----------------------------------------------------------------------
      IF (CMBDSN.NE.' ') THEN
         WRITE (*,*) 'Reading CMB "I" image from ',
     :        CMBDSN(1:LEN1(CMBDSN))
         WRITE (LOGUN,*) 'Reading CMB "I" image from ',
     :        CMBDSN(1:LEN1(CMBDSN))
         CALL ANREAD(1, CMBDSN(1:LEN1(CMBDSN)), IER)
         IF (IER.EQ.0) STOP
      END IF
      IF (CMBQDS.NE.' ') THEN
         WRITE (*,*) 'Reading CMB "Q" image from ',
     :        CMBQDS(1:LEN1(CMBQDS))
         WRITE (LOGUN,*) 'Reading CMB "Q" image from ',
     :        CMBQDS(1:LEN1(CMBQDS))
         CALL ANREAD(3, CMBQDS(1:LEN1(CMBQDS)), IER)
         IF (IER.EQ.0) STOP
      END IF
      IF (CMBUDS.NE.' ') THEN
         WRITE (*,*) 'Reading CMB "U" image from ',
     :        CMBUDS(1:LEN1(CMBUDS))
         WRITE (LOGUN,*) 'Reading CMB "U" image from ',
     :        CMBUDS(1:LEN1(CMBUDS))
         CALL ANREAD(4, CMBUDS(1:LEN1(CMBUDS)), IER)
         IF (IER.EQ.0) STOP
      END IF
C-----------------------------------------------------------------------
C Read in the Foreground image (if any).
C-----------------------------------------------------------------------
      IF (FGDSN.NE.' ') THEN
         WRITE (*,*) 'Reading foreground image from ',
     :        FGDSN(1:LEN1(FGDSN))
         WRITE (LOGUN,*) 'Reading foreground image from ',
     :        FGDSN(1:LEN1(FGDSN))
         CALL ANREAD(2, FGDSN(1:LEN1(FGDSN)), IER)
         IF (IER.EQ.0) STOP
      END IF
C-----------------------------------------------------------------------
C Loop for pointings.
C  For each pointing, we adjust the RA, Dec, date and file name.
C-----------------------------------------------------------------------
      NFIELD = MIN(MAXFLD,NPX*NPY)
      WRITE (LOGUN,*) ' '
      IF (NFIELD.LT.1) STOP
      CALL MOSAIC(NPX, NPY, SPAC, RA0D*RPDEG, DEC0D*RPDEG, MAXFLD,
     :  PTRA, PTDEC, PTX, PTY, FLDNAM)
      DO 1000 IFIELD=1,NFIELD
         WRITE (*,*) 'Doing field ', IFIELD, ' of ', NFIELD, ' ',
     :      FLDNAM(IFIELD)
C        -- Pointing center
         SRAD = PTRA(IFIELD)/RPDEG
         SDECD = PTDEC(IFIELD)/RPDEG
C        -- File name
         IF (NFIELD.GT.1) THEN
            WRITE (OUTDSN, '(A,A,A)')
     :           PREFIX(1:LEN1(PREFIX)), FLDNAM(IFIELD), '.uvf'
         ELSE
            WRITE (OUTDSN, '(A,A)')
     :           PREFIX(1:LEN1(PREFIX)), '.uvf'
         END IF
         WRITE (LOGUN,'(//''File: '',A)') OUTDSN(1:LEN1(OUTDSN))
C-----------------------------------------------------------------------
C Open FITS file.
C-----------------------------------------------------------------------
         STATUS = 0
         CALL FTINIT(UNIT3, OUTDSN, 1, STATUS)
         IF (STATUS .GT. 0) CALL ERROR('Cannot open output FITS file.')
C-----------------------------------------------------------------------
C Time and angle calculations.
C-----------------------------------------------------------------------
C Reference date: 0h UT on first day of observation.
C
      REFDAT = 2400000.5D0 + RMJD
      CALL SLA_DJCL(RMJD, RYEAR, RMONTH, RDAY, FD, IER)
      WRITE (RDSTR, '(I4.4,''-'',I2.2,''-'',I2.2)')
     :       RYEAR, RMONTH, RDAY
C
C Equation of equinoxes at reference date.
C
      EQEQX = SLA_EQEQX(RMJD)
C
C Precess source position to reference date; no proper motion or
C parallax.
C
      CALL SLA_MAPPA(2000D0, RMJD, SLAPRM)
      CALL SLA_MAPQKZ(SRAD*RPDEG, SDECD*RPDEG, SLAPRM, APPRA, APPDEC)
C-----------------------------------------------------------------------
C Read the list of point sources, and work out distance of each from
C pointing center.
C-----------------------------------------------------------------------
      IF (DSDSN.EQ.' ') THEN
         NDS = 0
         GOTO 77
      END IF
      WRITE (*,*) 'Reading point sources'
      IF (SIMEAN.LT.-2.0 .OR. SIMEAN.GT.2.0) THEN
         WRITE (*,*) 'Using spectral indices from input file'
      ELSE IF (SISD.LE.0.0) THEN
         WRITE (*,*) 'Using the same spectral index for all sources: ',
     :        SIMEAN
      ELSE
         WRITE (*,*) 'Using a gaussian distribution of spectral',
     :        ' indices with mean', SIMEAN, ' and sd ',SISD
         SSEED = -8756431
      END IF
      OPEN (UNIT=UNIT1, FILE=DSDSN, STATUS='OLD', IOSTAT=IER)
      IF (IER.NE.0) THEN
         WRITE (*,*) 'Cannot open source list'
         NDS = 0
      ELSE
         WRITE (LOGUN,*) 'Reading list of point sources from ',
     :        DSDSN(1:LEN1(DSDSN))
         WRITE (LOGUN,*) 'Scaling flux densities by', SRCFAC
         WRITE (LOGUN,*) 'Ignoring source stronger than', SRCCUT, ' mJy'
         NDS = 0
 22      IER = GTMODL(UNIT1, P1, P2, P3, P4)
         IF (IER.EQ.0) THEN
C           -- ignore zero-flux entries
            IF (P3.EQ.0.0) GOTO 22 
C           -- ignore sources above cutoff
            IF (ABS(P3*SRCFAC).GT.SRCCUT) GOTO 22
            NDS = NDS+1
            IF (NDS.GT.MAXDS) CALL ERROR('Too many point sources')
C           -- RA and Dec: precess to apparent
            P1 = 15D0*P1/3600D0
            P2 = P2/3600D0
            CALL SLA_MAPQKZ(P1*RPDEG, P2*RPDEG, SLAPRM, DSRA(NDS),
     :                      DSDEC(NDS))
C           -- compute direction cosines [ASP-Conf-6 p263]
            DSL(NDS) = COS(DSDEC(NDS))*SIN(DSRA(NDS)-APPRA)
            DSM(NDS) = SIN(DSDEC(NDS))*COS(APPDEC) -
     :                 COS(DSDEC(NDS))*SIN(APPDEC)*COS(DSRA(NDS)-APPRA)
C           -- distance from phase center (arcmin) and primary attenuation
            RAD(NDS) = ASIN(SQRT(DSL(NDS)**2+DSM(NDS)**2))
     :                 *60.0/RPDEG
C           -- modified flux density
            DSFLUX(NDS) = SRCFAC*
     :        0.001*P3/SQRT(1.0 - (DSL(NDS)**2 +DSM(NDS)**2))
C           -- spectral index
            IF (SIMEAN.LT.-2.0 .OR. SIMEAN.GT.2.0) THEN
               DSSI(NDS) = P4
            ELSE IF (SISD.LE.0.0) THEN
               DSSI(NDS) = SIMEAN
            ELSE
               DSSI(NDS) = SIMEAN + SISD*GASDEV(SSEED)
            END IF
C           -- summarize to log file
            WRITE (LOGUN,'(I4,F8.1,F7.3,F10.2)')
     :           NDS, P3, DSSI(NDS), RAD(NDS)
            GOTO 22
         END IF
         WRITE (LOGUN,*) 'Number of point sources:', NDS
         WRITE (*,*) 'Number of point sources:', NDS
         CLOSE(UNIT=UNIT1)
      END IF
C
C Write list of point sources in format acceptable to difmap.
C
      IF (NDS.GT.0) THEN
         OPEN(UNIT=UNIT2, FILE='sources.lm', STATUS='UNKNOWN')
         DO K=1,NDS
            WRITE (UNIT2, '(''pgpt '',g12.4,'', '',g12.4,'', 15'')') 
     :           DSL(K),DSM(K)
         END DO
         CLOSE(UNIT=UNIT2)
      END IF
 77   CONTINUE
C-----------------------------------------------------------------------
C Step through the observation integration by integration, accumulating
C uv points (but no data yet).
C-----------------------------------------------------------------------
      WRITE (*,*) 'Calculating uv coverage'
      WRITE (LOGUN,
     :  '(/1X,''Sample   UT       AZ      ZA   Airmass    PA'')')
      NSAMP = 0
      DELTPA = PASTAR
      IS = 0
      PAL = 0
      DO I=1,NINTEG
         IF (PALIST) THEN
            IF (IS.EQ.0) THEN
               PAL = PAL+1
               IF (PAL.GT.NPAL) PAL = 1
               IS = PALCNT(PAL)
            END IF
            IS = IS-1
            DELTPA = PALANG(PAL)*RPDEG
            IDLE = DELTPA .LT. 0.0
C           write (*,*) i, pal, deltpa/rpdeg, idle
         ELSE
            IS = IS+1
            IDLE = IS .GT. PAC2
         END IF
C
C        Time calculations. T is the number of seconds since REFDAT
C        for the center of the integration, and MJD the corresponding
C        modified Julian date (UT, strictly UT1). GAST is the
C        corresponding Greenwich apparent sidereal time (radians),
C        and HA the hour angle of the field center.
C
         T    = UTSTAR + (I-0.5)*INTEG
         MJD  = RMJD + DBLE(T)/86400D0
         GAST = SLA_GMST(MJD) + EQEQX
         LAST = GAST + OBSLNG
         HA   = LAST - APPRA
         IF (I.EQ.1) THEN
            WRITE(*, '(''At start LST='',F8.3,'' HA='',F8.3)')
     :           LAST/RPHR, HA/RPHR
         ELSE IF (I.EQ.NINTEG) THEN
            WRITE(*, '(''At end   LST='',F8.3,'' HA='',F8.3)')
     :           LAST/RPHR, HA/RPHR
         END IF
C
C        Azimuth (AZ), elevation (EL), parallactic angle (PA).
C
         SH = SIN(HA)
         CH = COS(HA)
         SD = SIN(APPDEC)
         CD = COS(APPDEC)
         SP = SIN(OBSLAT)
         CP = COS(OBSLAT)
         XX = -CH*CD*SP+SD*CP
         YY = -SH*CD
         ZZ = CH*CD*CP+SD*SP
         R  = SQRT(XX*XX+YY*YY)
         IF (R.EQ.0.0) THEN
            AZ = 0.0
         ELSE
            AZ = ATAN2(YY,XX)
         END IF
         ZA = PI/2D0 - ATAN2(ZZ,R)
         SQSZ = CP*SH
         CQSZ = SP*CD - CP*SD*CH
         IF (SQSZ.EQ.0D0.AND.CQSZ.EQ.0D0) CQSZ=1D0
         PA = ATAN2(SQSZ,CQSZ)
         IF (TRACK) THEN
            THETA = DELTPA
         ELSE
            THETA = MOD(DELTPA + PA, TWOPI)
         END IF
         IF (ZA.LE.ZALIM .AND. (.NOT.IDLE)) THEN
            AM = SLA_AIRMAS(ZA)
            NSTART = NSAMP+1
C
C           Add the uv points for this time.
C
            CALL UVCOV(N, STX, STY, ANTPOL, THETA, T, KMAX,
     :                 NSAMP, UT, BL, U, V, POLZN)
C
            WRITE (LOGUN,'(1X,I4,2X,A8,5F8.1, I8)') I, TIMEC(T),
     :           AZ/RPDEG, ZA/RPDEG, AM, PA/RPDEG, THETA/RPDEG, IS
            DO J=NSTART,NSAMP
               AIRMAS(J) = AM
            END DO
         END IF
         IF (.NOT.PALIST) THEN
            IF (IS.EQ.PAC2) DELTPA = MOD(DELTPA + PAC1, TWOPI)
            IF (IS.EQ.PAC2+PAC3) IS = 0
         END IF
      END DO
C-----------------------------------------------------------------------
C Check we have samples
C-----------------------------------------------------------------------
      IF (NSAMP.LT.1)
     :  CALL ERROR('Source is not visible in specified UT range')
C-----------------------------------------------------------------------
C Start with zero visibility.
C-----------------------------------------------------------------------
      DO I=1,NSAMP
         DO IPOL=1,4
            DO IF=1,NF
               VISR(I,IF,IPOL) = 0.0
               VISI(I,IF,IPOL) = 0.0
            END DO
            VISWT(I,IPOL) = 0.0
         END DO
      END DO
C-----------------------------------------------------------------------
C Add CMB contributions.
C-----------------------------------------------------------------------
      IF (CMBDSN.NE.' ' .OR. FGDSN.NE.' ') THEN
         WRITE (*,*) 'Computing CMB and foreground contribution'
C     -- Loop for frequency channels
      DO IF=1,NF
C        -- frequency = FREQ(IF) [Hz]
C        -- make a copy of the sky image plus the foreground
C        -- image multiplied by the primary
C        -- beam appropriate for this pointing. Note that width of 
C        -- primary beam scales with frequency. We also need to
C        -- scale from K to Jy/pixel, which depends on frequency.
         WRITE (LOGUN,*) SIGMA, FREQ(IF)
         CALL ANPBM(CBIBM, SIGMA*(30E9/REAL(FREQ(IF))), 
     :        PTX(IFIELD), PTY(IFIELD), REAL(FREQ(IF)), FGSI)
C        -- do FFT to put the image into the uv-plane
         CALL ANFFT
C        -- estimate visibility contribution for each uv point
         DO I=1,NSAMP
C           -- U(I), V(I) [in seconds]; convert to UU, VV in wavelengths
            UU = U(I)*FREQ(IF)
            VV = V(I)*FREQ(IF)
C           -- interpolation
            CALL GETVI2(UU, VV, VIS, MAPFT)
            CALL GETVI2(UU, VV, QVIS, MAPFTQ)
            CALL GETVI2(UU, VV, UVIS, MAPFTU)
C           -- RR
            VISR(I,IF,1) = VISR(I,IF,1) + REAL(VIS)
            VISI(I,IF,1) = VISI(I,IF,1) + AIMAG(VIS)
C           -- LL
            VISR(I,IF,2) = VISR(I,IF,2) + REAL(VIS)
            VISI(I,IF,2) = VISI(I,IF,2) + AIMAG(VIS)
C              -- RL = Q+iU
            VISR(I,IF,3) = VISR(I,IF,3) + REAL(QVIS) - AIMAG(UVIS)
            VISI(I,IF,3) = VISI(I,IF,3) + AIMAG(QVIS) + REAL(UVIS)
C              -- LR = Q-iU
            VISR(I,IF,4) = VISR(I,IF,4) + REAL(QVIS) + AIMAG(UVIS)
            VISI(I,IF,4) = VISI(I,IF,4) + AIMAG(QVIS) - REAL(UVIS)
         END DO
      END DO
      END IF
C-----------------------------------------------------------------------
C Add point sources.
C-----------------------------------------------------------------------
      IF (NDS.GT.0) THEN
         WRITE (*,*) 'Adding point sources (assumed unpolarized)'
         WRITE (LOGUN,'(/1X,2A)') 'Adding point sources ',
     :        '(effective fluxes in each channel in mJy)'
         WRITE (LOGUN,'(1X,A,G10.3,A)') 'Ignoring responses less than ',
     :        FLXLIM, ' Jy'
      END IF
      DO K=1,NDS
         MAXFL = 0.0
         DO IF=1,NF
            RELF = FREQ(IF)/30E9
            IF (CBIBM) THEN
               IF (RAD(K).GT.90.0) THEN
                  ATTEN = 0.0
               ELSE
                  ATTEN = CBIBEAM(REAL(RAD(K)*RPAMIN),
     :                            REAL(FREQ(IF)*1E-9))
               END IF
            ELSE
               ATTEN = EXP(-(RAD(K)*RELF/SIGMA)**2)
            END IF
            FL(IF) = DSFLUX(K)*ATTEN*RELF**DSSI(K)
C           write (*,*) k, if, freq(if)*1e-9, fl(if)
            MAXFL = MAX(MAXFL, ABS(FL(IF)))
         END DO
         IF (MAXFL.GT.FLXLIM) THEN
            WRITE(LOGUN,'(I4,2X,3P10F7.2)') K, (FL(IF),IF=1,NF)
            DO I=1,NSAMP
               DO IPOL=1,2 
C                 -- LL and RR only
                  DELAY = (U(I)*DSL(K) + V(I)*DSM(K))
                  DO IF=1,NF
                     PHASE = DELAY*TWOPI*FREQ(IF)
                     VISR(I,IF,IPOL) = VISR(I,IF,IPOL) +
     :                                 FL(IF)*COS(PHASE)
                     VISI(I,IF,IPOL) = VISI(I,IF,IPOL) +
     :                                 FL(IF)*SIN(PHASE)
                  END DO
               END DO
            END DO
         END IF
      END DO
C-----------------------------------------------------------------------
C Add gaussian noise.
C-----------------------------------------------------------------------
      WRITE (LOGUN,*) 'Adding noise (all polarizations)'
      RMS = NFAC*(TSYS+TATM*AIRMAS(I))
      DO IPOL=1,4
         DO I=1,NSAMP
            VISWT(I,IPOL) = 1.0/(RMS*RMS)
            DO IF=1,NF
               VISR(I,IF,IPOL) = VISR(I,IF,IPOL) + RMS*GASDEV(ISEED)
               VISI(I,IF,IPOL) = VISI(I,IF,IPOL) + RMS*GASDEV(ISEED)
            END DO
         END DO
      END DO
C-----------------------------------------------------------------------
C Calculate some statistics.
C-----------------------------------------------------------------------
      RUV = 0.0
      DO I=1,NSAMP
         RUV = MAX(RUV, U(I)*U(I) + V(I)*V(I))
      END DO
      WRITE (LOGUN,*) 'Max UV radius', SQRT(RUV), SQRT(RUV)*CLIGHT
C-----------------------------------------------------------------------
C Write FITS header.
C-----------------------------------------------------------------------
      write (*,*) 'writing FITS file'
      CALL FITSHD(UNIT3, NSAMP, NF, FIELD, RDSTR, SRAD, SDECD,
     :            FREQ(1), BW(1), REFDAT, -2, STATUS)
C     Additional history
      CALL FTPHIS(UNIT3, PROGNM//' version '//VERSN, STATUS)
      IF (CMBDSN.NE.' ') THEN
         CALL FTPHIS(UNIT3, 'CMB image: '//CMBDSN, STATUS)
      ELSE
         CALL FTPHIS(UNIT3, 'No CMB image', STATUS)
      END IF
      IF (FGDSN.NE.' ') THEN
         CALL FTPHIS(UNIT3, 'Foreground image: '//CMBDSN, STATUS)
      ELSE
         CALL FTPHIS(UNIT3, 'No foreground image', STATUS)
      END IF
      IF (STATUS.GT.0) CALL ERROR('Error while writing FITS header.')
C-----------------------------------------------------------------------
C Write FITS data. 4 polarizations at each uv point, only one filled.
C-----------------------------------------------------------------------
      NPOL = 4
      COUNT(1) = 0
      COUNT(2) = 0
      COUNT(3) = 0
      COUNT(4) = 0
      DO I=1,NSAMP
C         -- UVW
          PDATA(1) = U(I)
          PDATA(2) = V(I)
          PDATA(3) = 0.0
C         -- Baseline number: 256*NB1(IB) + NB2(IB)
          PDATA(4) = BL(I)
C         -- Days of time
          PDATA(5) = 0.0
C         -- Seconds of time
          PDATA(6) = UT(I)
C         -- Integration time (s)
          PDATA(7) = INTEG*NIGHTS
          CALL FTPGPE(UNIT3, I, 1, 7, PDATA, STATUS)
          KV = 0
          DO IF=1,NF
             DO IPOL=1,NPOL
C               -- order RR, LL, RL, LR
                RDATA(KV+1) = VISR(I,IF,IPOL)
                RDATA(KV+2) = VISI(I,IF,IPOL)
                IF (ALLPOL .OR. IPOL.EQ.POLZN(I)) THEN
                   RDATA(KV+3) = VISWT(I,IPOL)
                   COUNT(POLZN(I)) = COUNT(POLZN(I))+1
                ELSE
                   RDATA(KV+3) = -ABS(VISWT(I,IPOL))
                END IF
                KV = KV+3
             END DO
          END DO
          CALL FTPPRE(UNIT3, I, 1, 3*NPOL*NF, RDATA, STATUS)
       END DO
       IF (ALLPOL) WRITE (*,*) 'Writing all 4 polarizations'
       WRITE (*,*) COUNT(1)/NF, ' RR samples'
       WRITE (*,*) COUNT(2)/NF, ' LL samples'
       WRITE (*,*) COUNT(3)/NF, ' RL samples'
       WRITE (*,*) COUNT(4)/NF, ' LR samples'
C
C Write FITS-AIPS antenna table.
C
      CALL ANTAB(UNIT3, N,STNAME,STX,STY,STZ, RDSTR, FREQ, STATUS)
      IF (STATUS.GT.0) CALL ERROR('Error while writing AN table.')
C
C Write FITS-AIPS frequency table.
C
      CALL FQTAB(UNIT3, NF, FREQ, BW, STATUS)
      IF (STATUS.GT.0) CALL ERROR('Error while writing FQ table.')
C
C Close output file.
C
      CALL FTCLOS(UNIT3, STATUS)
      IF (STATUS.GT.0) CALL ERROR('Error while writing FITS file.')
C-----------------------------------------------------------------------
C End mosaic loop: increment date by one day; decrement start time by 4m
C-----------------------------------------------------------------------
         RMJD = RMJD+1
         UTSTAR = UTSTAR - 4*60
 1000 CONTINUE
      END

      SUBROUTINE UVCOV(N, X, Y, ANTPOL, THETA, T, KMAX, K, UT, BL,
     :                 U, V, POLZN)
      INTEGER N, KMAX, K
      DOUBLE PRECISION X(N), Y(N)
      INTEGER ANTPOL(N)
      REAL THETA
      INTEGER T
      INTEGER UT(KMAX), BL(KMAX)
      REAL U(KMAX), V(KMAX)
      INTEGER POLZN(KMAX)
C
C Compute uv coverage.
C Input:
C     N    = number of antennas
C     X,Y  = antenna coordinates on baseplate (meter)
C     ANTPOL = 0 for 'L' or 1 for 'R' polarization
C     THETA=plate rotation angle (rad)
C     T    = time of this sample (seconds from reference)
C     KMAX = dimension of arrays
C In/out:
C     K    = accumulated number of samples.
C Output:
C     UT   = time (seconds from reference)
C     BL   = baseline index
C     U, V = uv coordinates in seconds
C     POLZN = 1 for 'RR', 2 for 'LL', 3 for 'RL', 4 for 'LR'
C-----------------------------------------------------------------------
      INTEGER MAXS, MAXB
      PARAMETER (MAXS=40)
      PARAMETER (MAXB=MAXS*(MAXS-1)/2)
      REAL CLIGHT
      PARAMETER (CLIGHT=2.997924580E8)
      INTEGER I, J
      REAL UP, VP, CTH, STH
C
C Compute uv-coordinates for all antenna pairs; convert from meters to
C seconds using speed of light.
C
      CTH = COS(THETA)/CLIGHT
      STH = SIN(THETA)/CLIGHT
      DO 110 I=2,N
         DO 100 J=1,I-1
            IF (K.LT.KMAX) THEN
               K = K+1
               UT(K) = T
               BL(K) = 256*I + J
               UP = (X(I)-X(J))
               VP = (Y(I)-Y(J))
               U(K) =  UP*CTH + VP*STH
               V(K) = -UP*STH + VP*CTH
               IF (ANTPOL(I).EQ.0) THEN
                  IF (ANTPOL(J).EQ.0) THEN
                     POLZN(K) = 2 ! LL
                  ELSE
                     POLZN(K) = 3 ! RL
                  END IF
               ELSE ! ANTPOL(I).EQ.1
                  IF (ANTPOL(J).EQ.0) THEN
                     POLZN(K) = 4 ! LR
                  ELSE
                     POLZN(K) = 1 ! RR
                  END IF
               END IF
            END IF
 100     CONTINUE
 110  CONTINUE

      END

      SUBROUTINE FITSHD(UNIT, NSAMP, NIF, SOURCE, DATOBS, RA, DEC,
     :                  FREQ, BW, REFDAT, JUNK, STATUS)
      INTEGER UNIT, NSAMP, NIF, JUNK, STATUS
      CHARACTER*(*) SOURCE, DATOBS
      DOUBLE PRECISION RA, DEC, FREQ, BW, REFDAT
C
C Write FITS header. Version for 4 polarizations.
C
C UNIT:   I/O unit number
C NSAMP:  number of uv samples
C NIF:    number of IFs in each sample
C SOURCE: source/field name
C DATOBS: date of observation
C RA,DEC: field/pointing center (J2000, degrees)
C FREQ, BW: frequency and bandwidth for frequency axis
C REFDAT: reference date
C-----------------------------------------------------------------------
      DOUBLE PRECISION    PI,RPDEG
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (RPDEG=PI/180D0)
      INTEGER NAXIS
      PARAMETER (NAXIS=7)
      INTEGER PCOUNT
      PARAMETER (PCOUNT=7)
C-----------------------------------------------------------------------
      CHARACTER*12 OBSERV
      CHARACTER*23 STRRA
      CHARACTER*22 STRDEC
      INTEGER      NAXES(NAXIS)
      DATA NAXES/0, 3, 4, 1, 1, 1, 1/
C-----------------------------------------------------------------------
      STRRA = 'RA  = '
      CALL PANGLE('H', RA*RPDEG, STRRA(6:23), 4)
      STRDEC = 'Dec = '
      CALL PANGLE('D', DEC*RPDEG, STRDEC(6:22), 3)
      CALL USERNM(OBSERV)
C
C Write required keywords.
C
      NAXES(5) = NIF
      CALL FTPHPR(UNIT, .TRUE., -32, NAXIS, NAXES, 
     :            PCOUNT, NSAMP, .TRUE., STATUS)
      CALL FTPSCL(UNIT, 1.0D0, 0.0D0, STATUS)
      CALL FTPKYF(UNIT, 'BSCALE', 1.0, 1, ' ', STATUS)
      CALL FTPKYF(UNIT, 'BZERO',  0.0, 1, ' ', STATUS)
C
C Additional keywords.
C
      CALL FTPKYS(UNIT, 'OBJECT',   SOURCE,'Source name',STATUS)
      CALL FTPKYS(UNIT, 'TELESCOP','CBI',' ',STATUS)
      CALL FTPKYS(UNIT, 'INSTRUME','Simulation',' ',STATUS)
      CALL FTPKYS(UNIT, 'OBSERVER', OBSERV,' ',STATUS)
      CALL FTPKYS(UNIT, 'DATE-OBS', DATOBS,' ',STATUS)
      CALL FTPKYS(UNIT, 'BUNIT',    'JY','Correlated flux density',
     :            STATUS)
      CALL FTPKYS(UNIT, 'RADECSYS', 'FK5','Mean place', STATUS)
      CALL FTPKYF(UNIT, 'EQUINOX', 2000.0, 2,
     :            'Equinox of RA/Dec (J2000.0)', STATUS)
      CALL FTPKYF(UNIT, 'EPOCH', 2000.0, 2,
     :            'Alternate name for EQUINOX', STATUS)
      CALL FTPKYG(UNIT, 'OBSRA',RA, 12, 'Antenna pointing RA',STATUS)
      CALL FTPKYG(UNIT, 'OBSDEC',DEC, 12, 'Antenna pointing Dec',STATUS)
C
C FITS coordinate parameters.
C
      CALL FTPKYS(UNIT, 'CTYPE2','COMPLEX','1=real, 2=imag, 3=weight',
     1                STATUS)
      CALL FTPKYF(UNIT, 'CRVAL2', 1.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CDELT2', 1.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX2', 1.0,   1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'CTYPE3', 'STOKES',
     1            'Correlator: -1=RR, -2=LL, -3=RL, -4=LR',STATUS)
      CALL FTPKYF(UNIT, 'CRVAL3', -1.0,  1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CDELT3', -1.0,  1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX3',  1.0,  1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'CTYPE4',   'FREQ','Frequency, Hz',STATUS)
      CALL FTPKYG(UNIT, 'CRVAL4',     FREQ, 2, ' ',STATUS)
      CALL FTPKYG(UNIT, 'CDELT4',       BW, 2, ' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX4', 1.0,      1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'CTYPE5',     'IF','IF number',STATUS)
      CALL FTPKYF(UNIT, 'CRVAL5', 1.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CDELT5', 1.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX5', 1.0,   1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'CTYPE6','RA','Right ascension, degrees',STATUS)
      CALL FTPKYG(UNIT, 'CRVAL6', RA, 12, STRRA,STATUS)
      CALL FTPKYF(UNIT, 'CDELT6', 0.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX6', 1.0,   1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'CTYPE7','DEC','Declination, degrees',STATUS)
      CALL FTPKYG(UNIT, 'CRVAL7',DEC, 12, STRDEC,STATUS)
      CALL FTPKYF(UNIT, 'CDELT7', 0.0,   1,' ',STATUS)
      CALL FTPKYF(UNIT, 'CRPIX7', 1.0,   1,' ',STATUS)
C
C FITS random parameters.
C
      CALL FTPKYS(UNIT, 'PTYPE1','UU---SIN',
     1            'baseline u projection, seconds',STATUS)
      CALL FTPKYF(UNIT, 'PSCAL1',1.0,1,' ',STATUS)
      CALL FTPKYF(UNIT, 'PZERO1',0.0,1,' ',STATUS)
      CALL FTPKYS(UNIT, 'PTYPE2','VV---SIN',
     1            'baseline v projection, seconds',STATUS)
      CALL FTPKYF(UNIT, 'PSCAL2',1.0,1,' ',STATUS)
      CALL FTPKYF(UNIT, 'PZERO2',0.0,1,' ',STATUS)
      CALL FTPKYS(UNIT, 'PTYPE3','WW---SIN',
     1            'baseline w projection, seconds',STATUS)
      CALL FTPKYF(UNIT, 'PSCAL3',1.0,1,' ',STATUS)
      CALL FTPKYF(UNIT, 'PZERO3',0.0,1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'PTYPE4','BASELINE','256*ANT1 + ANT2',STATUS)
      CALL FTPKYF(UNIT, 'PSCAL4',1.0,1,' ',STATUS)
      CALL FTPKYF(UNIT, 'PZERO4',0.0,1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'PTYPE5','DATE','TAI Julian Date part 1',STATUS)
      CALL FTPKYD(UNIT, 'PSCAL5',1D0, 8, 'Days',STATUS)
      CALL FTPKYD(UNIT, 'PZERO5',REFDAT, 8, ' ',STATUS)
C
      CALL FTPKYS(UNIT, 'PTYPE6','DATE','TAI Julian Date part 2',STATUS)
      CALL FTPKYD(UNIT, 'PSCAL6',1D0/86400D0, 8,
     :            'Days/86400 (sec)', STATUS)
      CALL FTPKYF(UNIT, 'PZERO6',0.0,1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'PTYPE7','INTTIM','Integration time (sec)',
     :            STATUS)
      CALL FTPKYF(UNIT, 'PSCAL7',1.0,1,' ',STATUS)
      CALL FTPKYF(UNIT, 'PZERO7',0.0,1,' ',STATUS)
C
      CALL FTPKYS(UNIT, 'ORIGIN','POLFAKE Caltech',' ', STATUS)
      CALL FTPDAT(UNIT, STATUS)
      RETURN
C
      END

      SUBROUTINE ANTAB(UNIT,NS,STNAME,STX,STY,STZ, RDSTR, FREQ, STATUS)
C-----------------------------------------------------------------------
C Append AIPS antenna table as a FITS table extension to the FITS file.
C
C Arguments:
C  UNIT    (input, I*4)       I/O unit number
C  NS      (input, I*4)       Number of antennas.
C  STNAME  (input, C*8 array) Names of antennas.
C  STX     (input, R*8 array) Antenna x-coordinates (m).
C  STY     (input, R*8 array) Antenna y-coordinates (m).
C  STZ     (input, R*8 array) Antenna z-coordinates (m).
C  RDSTR   (input, C*)        Reference date (string).
C  FREQ    (input, R*8)       Reference frequency.
C  STATUS  (output, I*4)      Error flag.
C-----------------------------------------------------------------------
      INTEGER      STATUS, NS, UNIT
      CHARACTER*(*) STNAME(NS), RDSTR
      DOUBLE PRECISION       STX(NS), STY(NS), STZ(NS), FREQ
C-----------------------------------------------------------------------
      INTEGER MAXS
      PARAMETER (MAXS=40)
      DOUBLE PRECISION STXYZ(3,MAXS)
      INTEGER STNUM(MAXS), MNTSTA(MAXS)
      REAL STAXOF(MAXS), POLA(MAXS), POLCAL(4,MAXS)
      INTEGER I, TFIELD
      CHARACTER*8  EXTNAM
      CHARACTER*1  POLTY(MAXS)
      CHARACTER*16 TTYPE(12), TFORM(12), TUNIT(12)
      DATA TTYPE/ 'ANNAME', 'STABXYZ', 'ORBPARM', 'NOSTA',
     :            'MNTSTA', 'STAXOF',  'POLTYA',  'POLAA',
     :            'POLCALA','POLTYB',  'POLAB',   'POLCALB'/
      DATA TFORM/ '8A',     '3D',      '0D',      '1J',
     :            '1J',     '1E',      '1A',      '1E',
     :            '4E',     '1A',      '1E',      '4E'/
      DATA TUNIT/ ' ',      'METERS',  ' ',       ' ',
     :            ' ',      'METERS',  ' ',       'DEGREES',
     :            ' ',      ' ',       'DEGREES', ' '/
C
      DO I=1,NS
         STNUM(I) = I
         MNTSTA(I) = 4
C        -- MNTSTA=4 signifies 'bizarre'
         STAXOF(I) = 0.0
         POLTY(I) = ' '
         POLA(I) = 0.0
         POLCAL(1,I) = 0.0
         POLCAL(2,I) = 0.0
         POLCAL(3,I) = 0.0
         POLCAL(4,I) = 0.0
         STXYZ(1,I) = STX(I)
         STXYZ(2,I) = STY(I)
         STXYZ(3,I) = STZ(I)
      END DO
C
C Append an empty extension.
C
      CALL FTCRHD(UNIT, STATUS)
C
C Define parameters for table.
C
      TFIELD = 12
      EXTNAM = 'AIPS AN'
C
C Write required header parameters.
C
      CALL FTPHBN(UNIT, NS, TFIELD, TTYPE, TFORM, TUNIT,
     :            EXTNAM, 0, STATUS)
      CALL FTPKYJ(UNIT, 'EXTVER', 1,' ',STATUS)
C
C Write additional parameters.
C
      CALL FTPKYD(UNIT, 'ARRAYX', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'ARRAYY', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'ARRAYZ', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'GSTIA0', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'DEGPDY', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'FREQ',  FREQ, 18, ' ', STATUS)
      CALL FTPKYS(UNIT, 'RDATE',   RDSTR,  ' ', STATUS)
      CALL FTPKYD(UNIT, 'POLARX', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'POLARY', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'UT1UTC', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'DATUTC', 0D0, 18, ' ', STATUS)
      CALL FTPKYD(UNIT, 'IATUTC', 0D0, 18, ' ', STATUS)
      CALL FTPKYS(UNIT, 'TIMSYS', 'UTC',   ' ', STATUS)
      CALL FTPKYS(UNIT, 'ARRNAM', 'CBI',   ' ', STATUS)
      CALL FTPKYJ(UNIT, 'NUMORB', 0,       ' ', STATUS)
      CALL FTPKYJ(UNIT, 'NOPCAL', 4,       ' ', STATUS)
      CALL FTPKYJ(UNIT, 'FREQID', -1,      ' ', STATUS)
C
C Write the table data.
C
C ANNAME:
      CALL FTPCLS(UNIT, 1, 1, 1, NS, STNAME, STATUS)
C STABXYZ:
      CALL FTPCLD(UNIT, 2, 1, 1, NS*3, STXYZ, STATUS)
C ORBPARM (omitted):
C NOSTA:
      CALL FTPCLJ(UNIT, 4, 1, 1, NS, STNUM, STATUS)
C MNTSTA:
      CALL FTPCLJ(UNIT, 5, 1, 1, NS, MNTSTA, STATUS)
C STAXOF:
      CALL FTPCLE(UNIT, 6, 1, 1, NS, STAXOF, STATUS)
C POLTYA:
      CALL FTPCLS(UNIT, 7, 1, 1, NS, POLTY, STATUS)
C POLAA:
      CALL FTPCLE(UNIT, 8, 1, 1, NS, POLA, STATUS)
C POLCALA:
      CALL FTPCLE(UNIT, 9, 1, 1, NS*4, POLCAL, STATUS)
C POLTYB:
      CALL FTPCLS(UNIT,10, 1, 1, NS, POLTY, STATUS)
C POLAB:
      CALL FTPCLE(UNIT,11, 1, 1, NS, POLA, STATUS)
C POLCALB:
      CALL FTPCLE(UNIT,12, 1, 1, NS*4, POLCAL, STATUS)
C-----------------------------------------------------------------------
      END

      SUBROUTINE FQTAB(UNIT, NF, FREQ, BW, STATUS)
C-----------------------------------------------------------------------
C Append AIPS FQ table as a FITS table extension to the FITS file.
C
C Arguments:
C  UNIT    (input, I*4)       I/O unit number
C  NF      (input, I*4)       Number of IF's.
C  FREQ(*) (input, R*8)       Frequency
C  BW(*)   (input, R*8)       Bandwidth
C  STATUS  (in/out,I*4)
C-----------------------------------------------------------------------
      INTEGER      STATUS, NF, UNIT
      DOUBLE PRECISION FREQ(NF), BW(NF)
C-----------------------------------------------------------------------
      REAL FOFF(99), CHWID(99)
      INTEGER SB(99)
      INTEGER TFIELD, IF
      CHARACTER*2 STR
      CHARACTER*16 TTYPE(5), TFORM(5), TUNIT(5), EXTNAM
      DATA TTYPE/ 'FRQSEL', 'IF FREQ', 'CH WIDTH', 'TOTAL BANDWIDTH',
     :            'SIDEBAND'/
      DATA TFORM/ '1J    ', '*D     ', '*E      ', '*E             ',
     :            '*J    '/
      DATA TUNIT/ '      ', 'HZ     ', 'HZ      ', 'HZ             ',
     :            '      '/
C-----------------------------------------------------------------------
C     Convert number of IFs to character string and adjust format codes.
C
      IF (NF.LT.1 .OR. NF.GT.99) CALL ERROR('Too many IFs')
      WRITE (STR,'(I2)') NF
      IF (NF.GT.9) THEN
         TFORM(2) = STR(1:2)//'D'
         TFORM(3) = STR(1:2)//'E'
         TFORM(4) = STR(1:2)//'E'
         TFORM(5) = STR(1:2)//'J'
      ELSE
         TFORM(2) = STR(2:2)//'D'
         TFORM(3) = STR(2:2)//'E'
         TFORM(4) = STR(2:2)//'E'
         TFORM(5) = STR(2:2)//'J'
      END IF
C
C Append an empty extension.
C
      CALL FTCRHD(UNIT, STATUS)
C
C Define parameters for table.
C
      TFIELD = 5
      EXTNAM = 'AIPS FQ'
C
C Write required header parameters.
C
      CALL FTPHBN(UNIT, 1, TFIELD, TTYPE, TFORM, TUNIT,
     :            EXTNAM, 0, STATUS)
      CALL FTPKYJ(UNIT, 'EXTVER', 1,' ',STATUS)
      CALL FTPKYJ(UNIT, 'NO_IF', NF,' ',STATUS)
C
C Write the table data.
C
      DO 10 IF=1,NF
         FOFF(IF) = FREQ(IF) - FREQ(1)
         CHWID(IF) = BW(IF)
         SB(IF) = 1
 10   CONTINUE
         
C FRQSEL:
      CALL FTPCLJ(UNIT, 1, 1, 1, 1, 1, STATUS)
C IF FREQ:
      CALL FTPCLE(UNIT, 2, 1, 1, NF,  FOFF, STATUS)
C CH WIDTH:
      CALL FTPCLE(UNIT, 3, 1, 1, NF, CHWID, STATUS)
C TOTAL BANDWIDTH4
      CALL FTPCLE(UNIT, 4, 1, 1, NF, CHWID, STATUS)
C SIDEBAND:
      CALL FTPCLJ(UNIT, 5, 1, 1, NF, SB,    STATUS)
C-----------------------------------------------------------------------
      END

C*GTMODL -- read one model component from file
C+
      INTEGER FUNCTION GTMODL(INMOD,P1,P2,P3,P4)
      INTEGER INMOD
      DOUBLE PRECISION P1, P2, P3, P4
C
C Read one model component (card) from unit 'INMOD'.
C Each card carries up to 7 numbers:
C 6 floating-point parameters P1-P6, and an optional integer ITYPE.
C Copy comments following a '!' character to unit 'OUTMOD'
C-----------------------------------------------------------------------
      INTEGER OUTMOD
      PARAMETER (OUTMOD=6)
      DOUBLE PRECISION PARS(4)
      DOUBLE PRECISION KCTOR
      CHARACTER*80 CARD
      CHARACTER BLANK, EXCL
      INTEGER I, J, L
      CHARACTER*64 FNAME
      DATA BLANK/' '/
      DATA EXCL/'!'/
C
   5  READ(INMOD, '(A)', END=40) CARD
      L = 80
   10 IF (CARD(L:L).EQ.BLANK) THEN
            L = L-1
            IF (L.LT.1) GOTO 5
        GOTO 10
      END IF
C
      I = 1
      DO 30 J=1,4
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
      CALL ERROR('Invalid model format: file '//FNAME)
   50 CONTINUE
      P1 = PARS(1)
      P2 = PARS(2)
      P3 = PARS(3)
      P4 = PARS(4)
      GTMODL = 0
      RETURN 
C
   40 GTMODL = 1
      RETURN
      END

      SUBROUTINE MOSAIC(NX, NY, SPAC, RA0, DEC0, MAXFLD,
     :                  RA, DEC, PTX, PTY, FLDNAM)
      INTEGER NX, NY, MAXFLD
      DOUBLE PRECISION SPAC, RA0, DEC0, RA(MAXFLD), DEC(MAXFLD)
      REAL PTX(MAXFLD), PTY(MAXFLD)
      CHARACTER*(*) FLDNAM(MAXFLD)
C
      DOUBLE PRECISION    PI,RPDEG,RPAMIN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (RPDEG=PI/180D0)
      PARAMETER (RPAMIN=RPDEG/60D0)
C
      DOUBLE PRECISION COSDEC, SINDEC, XR, YR, ZR, X, Y
      INTEGER I, IX, IY
C-----------------------------------------------------------------------
C
      COSDEC = COS(DEC0)
      SINDEC = SIN(DEC0)
C     WRITE (*,'(20X,2F10.4)') RA0/RPDEG, DEC0/RPDEG
      WRITE (*,'(A)') ' List of fields in mosaic'
      WRITE (*,'(A5,3X,A2,6X,A2,6X,A2,8X,A3)')
     : 'Field', 'dX', 'dY', 'RA', 'Dec'
C
C Loop for pointings
C
      I = 0
      DO IX=1,NX
         X = SPAC*(IX-0.5D0*(NX+1))
         XR = X*RPAMIN
         DO IY=1,NY
            Y = SPAC*(IY-0.5D0*(NY+1))
            YR = Y*RPAMIN
            ZR = SQRT(1D0 - XR*XR - YR*YR)
            I = I+1
            IF (I.GT.MAXFLD) RETURN
            DEC(I) = ASIN(YR*COSDEC + ZR*SINDEC)
            RA(I)  = RA0 + ATAN2(XR, ZR*COSDEC - YR*SINDEC)
            PTX(I) = X
            PTY(I) = Y
            WRITE (FLDNAM(I), '(2I2.2)') IX, IY
            WRITE (*,'(I4,2F8.2, 2F10.4, 2X, A)') 
     :        I, PTX(I), PTY(I), RA(I)/RPDEG, DEC(I)/RPDEG, FLDNAM(I)
         END DO
      END DO
C
      END

      SUBROUTINE ANREAD(WHICH, FILE, STAT)
      INTEGER WHICH
      CHARACTER*(*) FILE
      INTEGER STAT
C
C Read image from FITS file FILE. 
C WHICH=1 for the CMB I image (MAPARR)
C       2 for the foreground image (FGARR)
C       3 for the CMB Q image (QMAP)
C       4 for the CMB U image (IMAP)
C Returns STAT=0 on error, =1 otherwise.
C-----------------------------------------------------------------------
      INCLUDE 'polfake.inc'
      INCLUDE 'mapcom.inc'
C
      REAL TWOPI
      PARAMETER (TWOPI=2.0*3.141592654)
C
      INTEGER IER, I
      INTEGER LEN1, JUNK
      CHARACTER*120 TXTERR
C
      REAL MAPMAX, MAPMIN
      REAL MEAN, RMS
      REAL DU, DV
      INTEGER NMAX, MMAX, NMIN, MMIN, NORM
C
 600  FORMAT(5X,'Number of pixels: ', I10/
     :       5X,'Maximum:          ', 1PG12.3,'  in cell',2I5/
     :       5X,'Minimum:          ', 1PG12.3,'  in cell',2I5/
     :       5X,'Mean:             ', 1PG12.3/
     :       5X,'rms:              ', 1PG12.3/
     :       5X,'Pixel units:      ',A)
C
      IER = 0
C     WRITE (*,*) '--- Reading FITS file'
      CALL FTOPEN(2, FILE, 0, JUNK, IER)
      IF (IER.NE.0) THEN
         WRITE (*,*) 'Error opening file'
         STAT = 0
         RETURN
      END IF
C     WRITE (*,*) '--- Reading FITS header and data'
      IF (WHICH.EQ.1) THEN
         CALL FITRD(2, MAPARR, IDIM, IER, TXTERR)
      ELSE IF (WHICH.EQ.2) THEN
         CALL FITRD(2, FGARR, IDIM, IER, TXTERR)
      ELSE IF (WHICH.EQ.3) THEN
         CALL FITRD(2, QMAP, IDIM, IER, TXTERR)
      ELSE IF (WHICH.EQ.4) THEN
         CALL FITRD(2, UMAP, IDIM, IER, TXTERR)
      END IF
      IF (IER.NE.0) THEN
         CALL PUTOUT(TXTERR(1:LEN1(TXTERR)))
         STAT = 0
         RETURN
      END IF
      IER = 0
C     WRITE (*,*) '--- Closing FITS file'
      CALL FTCLOS(2, IER)
      IF (IER.NE.0) THEN
         CALL PUTOUT('Error closing FITS image')
         STAT = 0
         RETURN
      END IF
C     -- Report information about image
      IF (MPBUN.NE.'K') CALL PUTOUT(
     :     ' ++WARNING++ units are not K so results will be incorrect')
      DU = ABS(360.0/(MPSIZ(1)*MPIN(1)))/TWOPI
      DV = ABS(360.0/(MPSIZ(2)*MPIN(2)))/TWOPI
      WRITE (*,'(5X,A,2I7)') 'Image size:   ', MPSIZ(1), MPSIZ(2)
      WRITE (*,'(5X,A,2F8.2)') 'Pixel size (arcmin): ', 
     :     MPIN(1)*60.0, MPIN(2)*60.0
      WRITE (*,'(5X,A,2F10.3)') 'UVpixel size (wavelengths): ',
     :     DU, DV
C     -- Save parameters and calculate statistics
      STAT = 1
      DO I=1,IDIM
         MAPWT(I) = 1.0
      END DO
      IF (WHICH.EQ.1) THEN
C        - we have a CMB image
         DOCMB = .TRUE.
C        - save its parameters
         MDIM(1) = MPSIZ(1)
         MDIM(2) = MPSIZ(2)
         PIXL(1) = MPIN(1)
         PIXL(2) = MPIN(2)
         REFP(1) = MPRP(1)
         REFP(2) = MPRP(2)
C        - calculate statistics
         CALL STATIS(MAPARR, MPSIZ(1), MPSIZ(2),
     :        MAPMAX, MAPMIN, NMAX, MMAX, NMIN, MMIN,
     :        MEAN, RMS)
      ELSE
C        - we have a foreground or Q/U image
         DOFG = .TRUE.
         IF (DOCMB) THEN
C           - we already have a CMB image, so check for compatibility
            IF (MDIM(1) .NE. MPSIZ(1) .OR. MDIM(2) .NE. MPSIZ(2) .OR.
     :          PIXL(1) .NE. MPIN(1)  .OR. PIXL(2) .NE. MPIN(2)  .OR.
     :          REFP(1) .NE. MPRP(1)  .OR. REFP(2) .NE. MPRP(2)) THEN
               WRITE (*,*) MDIM(1), MPSIZ(1), MDIM(2),  MPSIZ(2),
     :              PIXL(1),  MPIN(1), PIXL(2), MPIN(2),
     :              REFP(1),  MPRP(1), REFP(2), MPRP(2)
               CALL PUTOUT('Image does not match CMB image')
               STAT = 0
               RETURN
            END IF
         ELSE
C           - save its parameters
            MDIM(1) = MPSIZ(1)
            MDIM(2) = MPSIZ(2)
            PIXL(1) = MPIN(1)
            PIXL(2) = MPIN(2)
            REFP(1) = MPRP(1)
            REFP(2) = MPRP(2)
         END IF
         If (WHICH.EQ.2) THEN
            CALL STATIS(FGARR, MPSIZ(1), MPSIZ(2),
     :           MAPMAX, MAPMIN, NMAX, MMAX, NMIN, MMIN,
     :           MEAN, RMS)
         ELSE IF (WHICH.EQ.3) THEN
            CALL STATIS(QMAP, MPSIZ(1), MPSIZ(2),
     :           MAPMAX, MAPMIN, NMAX, MMAX, NMIN, MMIN,
     :           MEAN, RMS)
         ELSE IF (WHICH.EQ.4) THEN
            CALL STATIS(UMAP, MPSIZ(1), MPSIZ(2),
     :           MAPMAX, MAPMIN, NMAX, MMAX, NMIN, MMIN,
     :           MEAN, RMS)
         END IF
      END IF
C
      NORM = MPSIZ(1)*MPSIZ(2)
      WRITE (*,600) NORM, MAPMAX, MMAX, NMAX, MAPMIN, MMIN, NMIN,
     :     MEAN, RMS, MPBUN
      RETURN
      END

      SUBROUTINE STATIS(MAPARR, NX, NY,
     :                  MAPMAX, MAPMIN, NMAX, MMAX, NMIN, MMIN,
     :                  MEAN, RMS)
C     -- input
      INTEGER NX, NY
      REAL MAPARR(NX, NY)
C     -- output
      REAL MAPMAX, MAPMIN
      INTEGER NMAX, MMAX, NMIN, MMIN
      REAL MEAN, RMS
C
      INTEGER NORM, N, M
      REAL SUM, SUMSQ
C
C Locate maximum, minimum.
C
      MAPMAX = -1E37
      MAPMIN =  1E37
      NORM = 0
      SUM = 0.0
      SUMSQ = 0.0
      DO 220 N=1,NY
         DO 210 M=1,NX
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
            SUM = SUM + MAPARR(M,N)
            SUMSQ = SUMSQ + MAPARR(M,N)**2
            NORM = NORM + 1
 210     CONTINUE
 220  CONTINUE
C
C     Calculate and write statistics
C
      MEAN = SUM/NORM
      RMS = SQRT((SUMSQ/NORM) - MEAN**2)
      END

      SUBROUTINE ANPBM(CBIBM, SIGMA, XP, YP, FREQ, FGSI)
      LOGICAL CBIBM
      REAL SIGMA, XP, YP, FREQ, FGSI
C-----------------------------------------------------------------------
C Make a copy of the CMB image plus foreground image with primary beam
C correction applied, and scaling from K to Jy/pixel.
C input CMB image = MAPARR, input foreground image = FGARR,
C output image = MAPWT
C Make copies of QMAP and UMAP with primary beam correction in
C MAPWTQ, MAPWTU.
C
C Parameters:
C     CBIBM: TRUE for CBI theoretical beam, FALSE for gaussian
C     SIGMA: sigma of gaussian primary beam [arcmin]
C     XP: x-position of pointing center [arcmin]
C     YP: y-position of pointing center [arcmin]
C     FREQ: frequency in Hz
C     FGSI: spectral index of the foreground emission.
C-----------------------------------------------------------------------
      INCLUDE 'polfake.inc'
C
      INTEGER I, NX, NY, IX, IY
      INTEGER IXSHFT, IYSHFT, K
      INTEGER IDEL, JDEL
      REAL R, RSQ, ATTEN, DX, DY, X, Y, SCALE1, SCALE2, G, XF
      REAL CBIBEAM
      DOUBLE PRECISION    PI,RPDEG,RPAMIN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (RPDEG=PI/180D0)
      PARAMETER (RPAMIN=RPDEG/60D0)
C
C     -- image size (both images should be the same)
      NX = MDIM(1)
      NY = MDIM(2)
C     -- pixel size (arcmin)
      DX = PIXL(1)*60.0
      DY = PIXL(2)*60.0
C     -- scale factor (assuming T0 = 2.726K)
      XF = 0.5281*(FREQ/30E9)
C     -- normal blackbody
      G = XF*XF*EXP(XF)/(EXP(XF)-1.0)**2
      SCALE1 = (2.340*DX*DY)*((FREQ/30E9)**2)*G
C     -- power-law: omit g-factor
      SCALE2 = (2.340*DX*DY)*((FREQ/30E9)**FGSI)
      WRITE (*,*) 'ANPBM', SIGMA, XP, YP, XF, G, FGSI
      WRITE (*,*) 'Scalefactor: 1 K = ', SCALE1, SCALE2,
     :            ' Jy/pixel at', FREQ*1E-9, ' GHz'
C      WRITE (*,*) 'ANPBM', NX, NY, DX, DY, REFP(1), REFP(2)
C
C Shift array to pointing center
C
      IDEL = NINT(XP/DX)
      JDEL = NINT(YP/DY)
      WRITE (*,*) 'Shifting sky image by ', IDEL, JDEL, ' pixels'
C
C Add images and apply primary beam
C
      I = 0
      DO IY=1,NY
C        Position in input array, taking into account the shift         
         IYSHFT = IY + JDEL
         DO IX= 1,NX
            IXSHFT = IX + IDEL
C           -- I = pixel number in 1D output array 
            I = I+1
C           -- R = radius from pointing center in sigma
            X = (IX-REFP(1))*DX
            Y = (IY-REFP(2))*DY
            MAPWT(I)  = 0.0
            MAPWTQ(I) = 0.0
            MAPWTU(I) = 0.0
            IF (IXSHFT.GE.1 .AND. IXSHFT.LE.NX .AND.
     :          IYSHFT.GE.1 .AND. IYSHFT.LE.NY ) THEN
C              -- ATTEN = attenuation of primary beam
               IF (CBIBM) THEN
C                 -- theoretical beam                  
                  R = SQRT(X**2 + Y**2)
                  IF (R.GT.90.0) THEN
                     ATTEN = 0.0
                  ELSE
                     ATTEN = CBIBEAM(REAL(R*RPAMIN), REAL(FREQ/1E9))
                  END IF
               ELSE
C                 -- gaussian beam
                  RSQ = (X**2 + Y**2)/(SIGMA**2)
                  IF (RSQ.LE.14.0) THEN
                     ATTEN = EXP(-RSQ)
                  ELSE
                     ATTEN = 0.0
                  END IF
               END IF
C              -- K = pixel number in the input array
               K = IXSHFT + (IYSHFT-1)*NX
C              -- add the two images with appropriate scale factors
               IF (ATTEN.GT.0.0) THEN
                  IF (DOCMB) THEN
                     MAPWT(I) = MAPWT(I)+ATTEN*MAPARR(K)*SCALE1
                     MAPWTQ(I) = MAPWTQ(I)+ATTEN*QMAP(K)*SCALE1
                     MAPWTU(I) = MAPWTU(I)+ATTEN*UMAP(K)*SCALE1
                  END IF
                  IF (DOFG)  MAPWT(I) = MAPWT(I)+ATTEN*FGARR(K)*SCALE2
               END IF
            END IF
         END DO
      END DO
      RETURN
      END

      SUBROUTINE ANFFT
C-----------------------------------------------------------------------
C Do complex FFT of sky CMB images. 
C Input: MAPWT (real); output: MAPFT (complex)
C Input: MAPWTQ (real); output: MAPFTQ (complex)
C Input: MAPWTU (real); output: MAPFTU (complex)
C
      INCLUDE 'polfake.inc'
C
      INTEGER I, J, K, DIM(2), NX, NY
C
C Array dimensions (had better be powers of 2!).
C
      NX = MDIM(1)
      NY = MDIM(2)
C
C FFT of real 2D array. NX and NY must be powers of two. This
C routine phase-shifts so that the zero-frequency element is at
C (NX/2+1, NY/2+1).
C Input array is real of size MAPWT(NX,NY).
C Output array MAPFT is conjugate symmetric complex of
C size MAPFT(NX/2+1,NY).
C
      DIM(1) = NX
      DIM(2) = NY
      DO J=1,NY
         DO I=1,NX
            K = I + (J-1)*NX
            IF (MOD(I+J,2).EQ.1) THEN
               MAPWT(K) = -MAPWT(K)
               MAPWTQ(K) = -MAPWTQ(K)
               MAPWTU(K) = -MAPWTU(K)
            END IF
         END DO
      END DO
      CALL FOUR2(MAPFT, DIM, 2, 1, 0)
      CALL FOUR2(MAPFTQ, DIM, 2, 1, 0)
      CALL FOUR2(MAPFTU, DIM, 2, 1, 0)
      DO J=1,NY
         DO I=1,NX/2+1
            K = I + (J-1)*(NX/2+1)
C           -- the following is fudged to make the sign right
            IF (MOD(I+J,2).NE.1) THEN
               MAPFT(K) = -MAPFT(K)
               MAPFTQ(K) = -MAPFTQ(K)
               MAPFTU(K) = -MAPFTU(K)
            END IF
         END DO
      END DO
C
      END

      SUBROUTINE GETVI2(UU, VV, VIS, ARRAY)
      REAL UU, VV
      COMPLEX VIS
      COMPLEX ARRAY(*)
C-----------------------------------------------------------------------
C Return complex visibility VIS [Jy] for baseline (UU, VV) [wavelengths]
C By interpolation in array ARRAY (complex) (bilinear interpolation).
C-----------------------------------------------------------------------
      INCLUDE 'polfake.inc'
C     INCLUDE 'mapcom.inc'
C
      REAL TWOPI
      PARAMETER (TWOPI=2.0*3.141592654)
      INTEGER I, J, NX, NY
      COMPLEX V1, V2, V3, V4
      INTEGER ICENT, JCENT
      REAL    DU, DV
      REAL    X, Y
      COMPLEX GETVAL
C
C Array dimensions (had better be powers of 2!).
C
      NX = MDIM(1)
      NY = MDIM(2)
      ICENT = (NX/2)+1
      JCENT = (NY/2)+1
C
C Pixel dimensions.
C
      DU = (360.0/(MDIM(1)*PIXL(1)))/TWOPI
      DV = (360.0/(MDIM(2)*PIXL(2)))/TWOPI
C
C ARRAY has dimensions (NX/2+1, NY) with zero spatial frequency at
C (NX/2+1,  NY/2+1).
C
      X = UU/DU + ICENT
      Y = VV/DV + JCENT
      I = X
      J = Y
      X = X-I
      Y = Y-J
      IF (I.LT.1 .OR. I.GT.NX-1 .OR. J.LT.1 .OR. J.GT.NY-1 .OR.
     :    X.LT.0.0 .OR. Y.LT.0.0 .OR. X.GE.1.0 .OR. Y.GE.1.0) 
     :     CALL ERROR('in GETVI2')
      V1 = GETVAL(ARRAY, ICENT, NX, NY, I,   J  )
      V2 = GETVAL(ARRAY, ICENT, NX, NY, I+1, J  )
      V3 = GETVAL(ARRAY, ICENT, NX, NY, I+1, J+1)
      V4 = GETVAL(ARRAY, ICENT, NX, NY, I,   J+1)
      VIS = (1.0-X)*(1.0-Y)*V1 +
     :           X *(1.0-Y)*V2 +
     :           X *     Y *V3 +
     :      (1.0-X)*     Y *V4
C
      END


      COMPLEX FUNCTION GETVAL(MAPFT, ICENT, NX, NY, I, J)
C-----------------------------------------------------------------------
C Get value of MAPFT(I,J), taking into account conjugate symmetry.
C-----------------------------------------------------------------------
      INTEGER ICENT, NX, NY, I, J
      COMPLEX MAPFT(ICENT, NY)
      INTEGER JCENT, I1, J1
      JCENT = (NY/2) + 1
      IF (I.GT.ICENT) THEN
         I1 = 2*ICENT-I
         J1 = 2*JCENT-J
         GETVAL = CONJG(MAPFT(I1, J1))
      ELSE
         GETVAL = MAPFT(I, J)
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
         WRITE (*,*) '--- Number of pixels = ', PIXELS
         WRITE (*,*) '--- Maximum number of pixels = ', MAXPIX
         ERROR = 'Map too big'
         GOTO 550
      END IF
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
      RETURN
  550 STATUS = 1      
      RETURN
      END




