C*FRDT   -- revision date of file [Convex-UNIX]
C+
      SUBROUTINE FRDT (FNAM, RDATE)
      CHARACTER*(*) FNAM, RDATE
C
C Determine the "last modified" (revision) date of a named file.
C (Specific to Convex-UNIX).
C
C Arguments:
C   FNAM (input, character): file name.
C   RDATE (output, character): ASCII string representing the date/time
C         (up to 20 characters). Result is string 'Undefined' if date
C         cannot be determined (e.g., file does not exist).
C
C Author:
C   Tim Pearson
C
C Subroutines required:
C   STAT (3F) [Convex-UNIX]
C   CTIME (3F) [Convex-UNIX]
C
C Modifications:
C   1988 Jul 14 - new subroutine.
C   2000 Sep 06 - RLM - change dim of STATB(13)
C 
C-----------------------------------------------------------------------
      INTEGER STATB(13), IER
      INTEGER STAT
      CHARACTER*24 CTIME, UTIME
      CHARACTER*20 VTIME
C
      IER = STAT(FNAM, STATB)
      IF (IER.NE.0) THEN
          RDATE = 'Undefined'
      ELSE
          UTIME = CTIME(STATB(10))
          VTIME(1:2) = UTIME(9:10)
          VTIME(3:3) = '-'
          VTIME(4:6) = UTIME(5:7)
          VTIME(7:7) = '-'
          VTIME(8:11) = UTIME(21:24)
          VTIME(12:12) = ' '
          VTIME(13:20) = UTIME(12:19)
          RDATE = VTIME
      END IF
C
      END
