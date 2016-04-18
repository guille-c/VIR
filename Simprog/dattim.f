C*DATTIM -- return current date and time [UNIX]
C+
      SUBROUTINE DATTIM (STRING)
      CHARACTER*(*) STRING
C
C Return the current date and time in a character string in standard
C format: 'dd-MMM-yyyy hh:mm:ss'. The caller can request the full
C string, by supplying a string of length 20 characters, or just the
C date, by supplying a string of length 11 characters.
C
C Argument:
C  STRING  (output) : receives the date and time, or string '<unknown>'
C                    if the information is not available; should be
C                    either 11 characters or 20 characters.
C
C Subroutines required:
C  UDATE (UNIX -- in udate.c)
C
C History:
C  1988 Apr 2 - TJP.
C  1989 Aug 10 - DLM modified for Sun (would be better to use standard
C                UNIX calls in C).
C  1989 Aug 14 - DLM modified for general UNIX system (using time()/
C                localtime() calls in udate.c)
C-----------------------------------------------------------------------
      INTEGER INOW(6), I
      CHARACTER DATE*11, TIME*8
      CHARACTER*3 MLIST(12)
      SAVE MLIST
      DATA MLIST/ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     1             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
C
      CALL UDATE(INOW)
C
C     Write integer day of month in DATE string
      WRITE (DATE(1:2),'(I2)')  INOW(3)
      DATE(3:3) = '-'
C
C     Write month character string in DATE string
      DATE(4:6) = MLIST(INOW(2))
      DATE(7:7) = '-'
C
C     Write integer year in DATE string
      WRITE (DATE(8:11),'(I4)') INOW(1)
C
C     Write hour, minute, second into TIME string; pad with zeros
      WRITE (TIME(1:2),'(I2)')   INOW(4)
      TIME(3:3) = ':'
      WRITE (TIME(4:5),'(I2)')  INOW(5)
      TIME(6:6) = ':'
      WRITE (TIME(7:8),'(I2)')  INOW(6)
      DO I=1,8
	IF (TIME(I:I) .EQ. ' ') TIME(I:I) = '0'
      END DO
C
      STRING = DATE(1:11)
      IF (LEN(STRING) .GE. 20) STRING = DATE(1:11)//' '//TIME(1:8)
C
      RETURN
      END
