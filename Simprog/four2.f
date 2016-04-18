      SUBROUTINE FOUR2 (DATA, N, NDIM, ISIGN, IFORM)
C-----------------------------------------------------------------------
C Cooley-Tukey Fast Fourier Transform.
C Multi-dimensional transform, each dimension a power of two,
C complex or real data.
C
C Transform(K1,K2,...) = Sum(Data(J1,J2,...)*exp(ISIGN*2*pi*sqrt(-1)
C     *((J1-1)*(K1-1)/N(1) + (J2-1)*(K2-1)/N(2) + ...))), 
C summed for all J1 and K1 from 1 to N(1), J2 and K2 from 1 to N(2),
C etc. for all NDIM subscripts.
C
C NDIM must be positive and each N(IDIM) must be a power of two.
C
C ISIGN is +1 or -1.
C
C Let NTOT = N(1)*N(2)*...*N(NDIM).  Then a -1 transform followed by
C a +1 one (or vice versa) returns NTOT times the original data.
C
C IFORM = 1, 0 or -1, as data is complex, real or the first half of
C a complex array.  Transform values are returned to array DATA.
C They are complex, real or the first half of a complex array, as
C IFORM = 1, -1 or 0.  The transform of a real array (IFORM = 0)
C dimensioned N(1) by N(2) by ... will be returned in the same
C array, now considered to be complex of dimensions N(1)/2+1 by N(2)
C by ....  Note that if IFORM = 0 or -1, N(1) must be even, and
C enough room must be reserved.  The missing values may be obtained
C by complex conjugation.  The reverse transformation, of a half
C complex array dimensioned N(1)/2+1 by N(2) by ..., is
C accomplished by setting iform to -1.  In the N array, N(1) must be
C the true N(1), not N(1)/2+1.  The transform will be real and
C returned to the input array.
C
C Running time is proportional to NTOT*log2(NTOT), rather than the
C naive NTOT**2.  Furthermore, less error is built up.
C
C Written by Norman Brenner of MIT Lincoln Laboratory, June 1968.
C Edited for Fortran-77 by Tim Pearson, February 2001.
C
C See-- IEEE Audio Transactions (June 1967), Special Issue On FFT.
C-----------------------------------------------------------------------
      REAL DATA(*)
      INTEGER N(*), NDIM, ISIGN, IFORM, NTOT, IDIM, NREM, NPREV
      INTEGER NCURR, JDIM

      NTOT = 1
      DO IDIM=1,NDIM
         NTOT = NTOT*N(IDIM)
      END DO
      IF (IFORM.GE.0) THEN
         NREM = NTOT
         DO IDIM=1,NDIM
            NREM = NREM/N(IDIM)
            NPREV = NTOT/(N(IDIM)*NREM)
            NCURR = N(IDIM)
            IF ((IDIM-1+IFORM).LE.0) THEN
               NCURR = NCURR/2
            END IF
            CALL BITRV (DATA,NPREV,NCURR,NREM)
            CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)
            IF ((IDIM-1+IFORM).LE.0) THEN
               CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)
               NTOT = (NTOT/N(1))*(N(1)/2+1)
            END IF
         END DO
      ELSE
         NTOT = (NTOT/N(1))*(N(1)/2+1)
         NREM = 1
         DO JDIM=1,NDIM
            IDIM = NDIM+1-JDIM
            NCURR = N(IDIM)
            IF ((IDIM-1).LE.0) THEN
               NCURR = NCURR/2
               CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)
               NTOT = NTOT/(N(1)/2+1)*N(1)
            END IF
            NPREV = NTOT/(N(IDIM)*NREM)
            CALL BITRV (DATA,NPREV,NCURR,NREM)
            CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)
            NREM = NREM*N(IDIM)
         END DO
      END IF
      RETURN
      END

      SUBROUTINE FIXRL (DATA,N,NREM,ISIGN,IFORM)
C-----------------------------------------------------------------------
C For IFORM = 0, convert the transform of a doubled-up real array,
C considered complex, into its true transform.  Supply only the
C first half of the complex transform, as the second half has
C conjugate symmetry.  For IFORM = -1, convert the first half
C of the true transform into the transform of a doubled-up real
C array.  N must be even.
C
C Using complex notation and subscripts starting at zero, the
C transformation is--
C     DIMENSION DATA(N,NREM)
C     ZSTP = EXP(ISIGN*2*PI*I/N)
C     DO 10 I2=0,NREM-1
C       DATA(0,I2) = CONJ(DATA(0,I2))*(1+I)
C       DO 10 I1=1,N/4
C         Z = (1+(2*IFORM+1)*I*ZSTP**I1)/2
C         I1CNJ = N/2-I1
C         DIF = DATA(I1,I2)-CONJ(DATA(I1CNJ,I2))
C         TEMP = Z*DIF
C         DATA(I1,I2) = (DATA(I1,I2)-TEMP)*(1-IFORM)
C 10      DATA(I1CNJ,I2) = (DATA(I1CNJ,I2)+CONJ(TEMP))*(1-IFORM)
C If I1=I1CNJ, the calculation for that value collapses into
C a simple conjugation of DATA(I1,I2).
C-----------------------------------------------------------------------
      REAL DATA(*), TWOPI,  TEMPR, THETA, SINTH, ZSTPR, ZSTPI, ZR, ZI
      REAL DIFR, DIFI, TEMPI
      INTEGER N, NREM, ISIGN, IFORM, IP0, IP1, IP2, J1, I2MIN, I2
      INTEGER I1MIN, I1, I2CNJ, I1MAX

      TWOPI=6.283185307*FLOAT(ISIGN)
      IP0=2
      IP1=IP0*(N/2)
      IP2=IP1*NREM
      IF (IFORM) 10,70,70
C     PACK THE REAL INPUT VALUES (TWO PER COLUMN)
 10   J1=IP1+1
      DATA(2)=DATA(J1)
      IF (NREM-1) 70,70,20
 20   J1=J1+IP0
      I2MIN=IP1+1
      DO 60 I2=I2MIN,IP2,IP1
      DATA(I2)=DATA(J1)
      J1=J1+IP0
      IF (N-2) 50,50,30
 30   I1MIN=I2+IP0
      I1MAX=I2+IP1-IP0
      DO 40 I1=I1MIN,I1MAX,IP0
      DATA(I1)=DATA(J1)
      DATA(I1+1)=DATA(J1+1)
 40   J1=J1+IP0
 50   DATA(I2+1)=DATA(J1)
 60   J1=J1+IP0
 70   DO 80 I2=1,IP2,IP1
      TEMPR=DATA(I2)
      DATA(I2)=DATA(I2)+DATA(I2+1)
 80   DATA(I2+1)=TEMPR-DATA(I2+1)
      IF (N-2) 200,200,90
 90   THETA=TWOPI/FLOAT(N)
      SINTH=SIN(THETA/2.)
      ZSTPR=-2.*SINTH*SINTH
      ZSTPI=SIN(THETA)
      ZR=(1.-ZSTPI)/2.
      ZI=(1.+ZSTPR)/2.
      IF (IFORM) 100,110,110
 100  ZR=1.-ZR
      ZI=-ZI
 110  I1MIN=IP0+1
      I1MAX=IP0*(N/4)+1
      DO 190 I1=I1MIN,I1MAX,IP0
         DO 180 I2=I1,IP2,IP1
            I2CNJ=IP0*(N/2+1)-2*I1+I2
            IF (I2-I2CNJ) 150,120,120
 120        IF (ISIGN*(2*IFORM+1)) 130,140,140
 130        DATA(I2+1)=-DATA(I2+1)
 140        IF (IFORM) 170,180,180
 150        DIFR=DATA(I2)-DATA(I2CNJ)
            DIFI=DATA(I2+1)+DATA(I2CNJ+1)
            TEMPR=DIFR*ZR-DIFI*ZI
            TEMPI=DIFR*ZI+DIFI*ZR
            DATA(I2)=DATA(I2)-TEMPR
            DATA(I2+1)=DATA(I2+1)-TEMPI
            DATA(I2CNJ)=DATA(I2CNJ)+TEMPR
            DATA(I2CNJ+1)=DATA(I2CNJ+1)-TEMPI
            IF (IFORM) 160,180,180
 160        DATA(I2CNJ)=DATA(I2CNJ)+DATA(I2CNJ)
            DATA(I2CNJ+1)=DATA(I2CNJ+1)+DATA(I2CNJ+1)
 170        DATA(I2)=DATA(I2)+DATA(I2)
            DATA(I2+1)=DATA(I2+1)+DATA(I2+1)
 180     CONTINUE
         TEMPR=ZR-.5
         ZR=ZSTPR*TEMPR-ZSTPI*ZI+ZR
 190     ZI=ZSTPR*ZI+ZSTPI*TEMPR+ZI
C     RECURSION SAVES TIME, AT A SLIGHT LOSS IN ACCURACY.  IF AVAILABLE,
C     USE DOUBLE PRECISION TO COMPUTE ZR AND ZI.
 200  IF (IFORM .GE. 0) THEN
C     UNPACK THE REAL TRANSFORM VALUES (TWO PER COLUMN)
         I2=IP2+1
         I1=I2
         J1=IP0*(N/2+1)*NREM+1
         GO TO 250
 220     DATA(J1)=DATA(I1)
         DATA(J1+1)=DATA(I1+1)
         I1=I1-IP0
         J1=J1-IP0
 230     IF (I2-I1) 220,240,240
 240     DATA(J1)=DATA(I1)
         DATA(J1+1)=0.
 250     I2=I2-IP1
         J1=J1-IP0
         DATA(J1)=DATA(I2+1)
         DATA(J1+1)=0.
         I1=I1-IP0
         J1=J1-IP0
         IF (I2-1) 260,260,230
 260     DATA(2)=0.
      END IF
      RETURN
      END

      SUBROUTINE BITRV(DATA,NPREV,N,NREM)
C-----------------------------------------------------------------------
C        SHUFFLE THE DATA BY 'BIT REVERSAL'.
C        DIMENSION DATA(NPREV,N,NREM)
C        DATA(I1,I2REV,I3) = DATA(I1,I2,I3), FOR ALL I1 FROM 1 TO
C        NPREV, ALL I2 FROM 1 TO N (WHICH MUST BE A POWER OF TWO), AND
C        ALL I3 FROM 1 TO NREM, WHERE I2REV-1 IS THE BITWISE REVERSAL
C        OF I2-1.  FOR EXAMPLE, N = 32, I2-1 = 10011 AND I2REV-1 =
C        11001.
C-----------------------------------------------------------------------
      COMPLEX DATA(*), TEMP
      INTEGER NPREV, N, NREM, IP0, IP1, IP4, IP5, I4REV, I4MAX, I4
      INTEGER I1MAX, I1, I5, I5REV, IP2

      IP0=1
      IP1=IP0*NPREV
      IP4=IP1*N
      IP5=IP4*NREM
      I4REV=1
      I4MAX=IP4
      DO 60 I4=1,I4MAX,IP1
      IF(I4-I4REV)10,30,30
   10 I1MAX=I4+IP1-IP0
      DO 20 I1=I4,I1MAX,IP0
      DO 20 I5=I1,IP5,IP4
      I5REV=I4REV+I5-I4
      TEMP=DATA(I5)
      DATA(I5)=DATA(I5REV)
   20 DATA(I5REV)=TEMP
   30 IP2=IP4/2
   40 IF(I4REV-IP2)60,60,50
   50 I4REV=I4REV-IP2
      IP2=IP2/2
      IF(IP2-IP1)60,40,40
   60 I4REV=I4REV+IP2
      RETURN
      END

      SUBROUTINE COOL2(DATA,NPREV,N,NREM,ISIGN)
C-----------------------------------------------------------------------
C        FOURIER TRANSFORM, LENGTH N, BY THE COOLEY-TUKEY ALGORITHM, IN
C        PLACE, BIT-REVERSED TO NORMAL ORDER, SANDE-TUKEY PHASE SHIFTS.
C        COMPLEX DATA
C        DIMENSION DATA(NPREV,N,NREM)
C        DATA(I1,J2,I3) = SUM(DATA(I1,I2,I3)CEXP(ISIGNC2CPICIC((I2-1)C
C        (J2-1)/N))), SUMMED OVER I2 = 1 TO N FOR ALL I1 FROM 1 TO
C        NPREV, J2 FROM 1 TO N AND I3 FROM 1 TO NREM.  N MUST BE A
C        POWER OF TWO.  FACTORING N BY FOUR'S DECREASES RUNNING TIME
C        BY ABOUT TWENTY FIVE PERCENT OVER FACTORING BY TWO'S.
C        GENERATING THE PHASE SHIFT FACTORS BY RECURSION SAVES ABOUT
C        TWENTY FIVE PERCENT OVER COMPUTING THEM, WHILE TABLING THEM
C        IS UNLIKELY TO SAVE MORE THAN ANOTHER TEN PERCENT.
C-----------------------------------------------------------------------
      COMPLEX DATA(*),TEMP,WSTP,W,W2,W3,T0,T1,T2,T3
      INTEGER NPREV, N, NREM, ISIGN, IP0, IP1, IP2, IP3, IP4, IP5
      INTEGER NPART, I1MAX, I1, I2, I5, I3A, I3B, I3C, I3D
      REAL TWOPI, THETA, SINTH

      TWOPI=6.28318530717958647692*FLOAT(ISIGN)
      IP0=1
      IP1=IP0*NPREV
      IP4=IP1*N
      IP5=IP4*NREM
      IP2=IP1
      NPART=N
   10 IF(NPART-2)50,30,20
   20 NPART=NPART/4
      GO TO 10
C
C        DO A FOURIER TRANSFORM OF LENGTH TWO
C
   30 IP3=IP2*2
      I1MAX=IP1
      DO 40 I1=1,I1MAX,IP0
      DO 40 I5=I1,IP5,IP3
      I3A=I5
      I3B=I3A+IP2
      TEMP=DATA(I3B)
      DATA(I3B)=DATA(I3A)-TEMP
   40 DATA(I3A)=DATA(I3A)+TEMP
      GO TO 140
C
C        DO A FOURIER TRANSFORM OF LENGTH FOUR (WITHOUT BIT REVERSAL)
C
   50 IP3=IP2*4
      THETA=TWOPI/FLOAT(IP3/IP1)
      SINTH=SIN(THETA/2.)
      WSTP=CMPLX(-2.*SINTH*SINTH,SIN(THETA))
C
C        COS(THETA)-1, FOR ACCURACY
C
      W=1.
      DO 130 I2=1,IP2,IP1
      IF(I2-1)70,70,60
   60 W2=W*W
      W3=W2*W
   70 I1MAX=I2+IP1-IP0
      DO 120 I1=I2,I1MAX,IP0
      DO 120 I5=I1,IP5,IP3
      I3A=I5
      I3B=I3A+IP2
      I3C=I3B+IP2
      I3D=I3C+IP2
      IF(I2-1)90,90,80
C
C        MULTIPLY BY THE PHASE SHIFT FACTORS
C
   80 DATA(I3B)=W2*DATA(I3B)
      DATA(I3C)=W*DATA(I3C)
      DATA(I3D)=W3*DATA(I3D)
   90 T0=DATA(I3A)+DATA(I3B)
      T1=DATA(I3A)-DATA(I3B)
      T2=DATA(I3C)+DATA(I3D)
      T3=DATA(I3C)-DATA(I3D)
      DATA(I3A)=T0+T2
      DATA(I3C)=T0-T2
      TEMP=(0.,1.)*T3
      IF(ISIGN)100,100,110
  100 TEMP=-TEMP
  110 DATA(I3B)=T1+TEMP
  120 DATA(I3D)=T1-TEMP
  130 W=W*WSTP+W
  140 IP2=IP3
      IF(IP3-IP4)50,150,150
  150 RETURN
      END
