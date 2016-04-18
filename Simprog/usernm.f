C*USERNM -- return name of current user [Convex-UNIX]
C+
      SUBROUTINE USERNM (USER)
      CHARACTER*(*) USER
C
C Return the username currently associated with the calling process.
C
C Argument:
C  USER   (output) : receives the node name, or string '<unknown>' if
C                    information is not available; should be at least
C                    12 characters.
C
C Subroutines required:
C  GETLOG (Convex)
C  GETENV (Convex)
C
C  1988 Apr 14 - TJP.
C-----------------------------------------------------------------------
C
C First see what GETLOG will return:
C
      CALL GETLOG(USER)
      IF (USER.NE.' ') RETURN
C
C If this fails, use GETENV to get the environment variable USER:
C
      CALL GETENV('USER', USER)
      IF (USER.NE.' ') RETURN
C
C Otherwise give up:
C
      USER = '<unknown>'
C
      END
