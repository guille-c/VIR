C*TIMEC  -- convert time in seconds to 'hh:mm:ss'
C+
      CHARACTER*8 FUNCTION TIMEC(TIME)
      INTEGER TIME
C
C Converts TIME in seconds to 'hh:mm:ss' character string. Times
C outside the range 0 to 86400 (24 hours) give a blank string.
C
C Argument:
C  TIME   (input) : count of seconds since midnight (i.e., integer
C                   in range 0-86400).
C
C Returns:
C  TIMEC          : character string in form 'hh:mm:ss'; e.g. TIME=0
C                   gives '00:00:00'; TIME=86399 gives '23:59:59'.
C
C History:
C  1987 Nov 12 - TJP.
C-----------------------------------------------------------------------
      INTEGER HOURS,MINS,SECS
      CHARACTER*1 DIGIT(0:9),SEP
      CHARACTER*8 OUT
      DATA DIGIT/'0','1','2','3','4','5','6','7','8','9'/
      DATA SEP/':'/
C
      IF (TIME.LT.0 .OR. TIME.GT.86400) THEN
          OUT = ' '
      ELSE
          HOURS = TIME/3600
          MINS = MOD(TIME,3600)/60
          SECS = MOD(TIME,60)
          OUT(1:1) = DIGIT(HOURS/10)
          OUT(2:2) = DIGIT(MOD(HOURS,10))
          OUT(3:3) = SEP
          OUT(4:4) = DIGIT(MINS/10)
          OUT(5:5) = DIGIT(MOD(MINS,10))
          OUT(6:6) = SEP
          OUT(7:7) = DIGIT(SECS/10)
          OUT(8:8) = DIGIT(MOD(SECS,10))
      END IF
      TIMEC = OUT
C
      END
