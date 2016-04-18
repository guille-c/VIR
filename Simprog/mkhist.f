C*MKHIST -- create a standard history record
C+
      SUBROUTINE MKHIST(PROGNM, VERSN, HIST)
      CHARACTER*(*) PROGNM
      CHARACTER*(*) VERSN
      CHARACTER*(*) HIST
C
C Create a history record for a Merge file in standard format. This
C format is:
C "<Program> (v<version>) run by <user> [<node>], <date> <time>"
C
C Arguments:
C  PROGNM  (input) : program name.
C  VERSN   (input) : program version.
C  HIST   (output) : the history record; N.B. the length of HIST must
C                    be 80 characters.
C-----------------------------------------------------------------------
      CHARACTER    UDATE*20, USER*12, NODE*12
      INTEGER      LEN1
C
      CALL DATTIM(UDATE)
      CALL USERNM(USER)
      CALL NODENM(NODE)
      HIST = PROGNM//' (v'//VERSN//') run by '//
     1       USER(1:LEN1(USER))//' ['//NODE(1:LEN1(NODE))//'], '//UDATE
      END
