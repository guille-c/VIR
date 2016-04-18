C*GETIN -- get one line from standard input, with prompt [Convex-UNIX]
C+
      INTEGER FUNCTION GETIN(TEXT,PROMPT,LENG)
      CHARACTER*(*) TEXT
      CHARACTER*(*) PROMPT
      INTEGER LENG
C
C Read a character string from standard input with optional prompt.
C
C Arguments:
C  TEXT   (output) : text read from standard input.
C  PROMPT (input)  : prompt string; if this is blank, or if the
C                    standard input is not a terminal, no prompt is
C                    issued.
C  LENG   (output) : number of characters readand returned; GETIN
C                    reads a complete line, but only returns as many
C                    characters as will fit in TEXT.
C
C Returns:
C  GETIN           : 1 if successful;
C                    0 if end-of-file or an error occurs.
C
C Non-standard extensions to Fortran-77:
C  Uses Q and $ format descriptors.
C
C Subroutines required:
C  Fortran I/O.
C
C History:
C  1987 Nov 11 - TJP.
C  1989 Aug 23 - DLM. Remove Q descriptor for Sun-3/Sun-4 operation.
C-----------------------------------------------------------------------
      INTEGER IER
C
      IF (PROMPT.NE.' ') WRITE (6, '(A,$)') PROMPT
      READ (5, '(A)', IOSTAT=IER) TEXT
      LENG = LEN1(TEXT)
      IF (IER.EQ.0) THEN
          GETIN = 1
          IF (LENG.GT.LEN(TEXT)) LENG = LEN(TEXT)
      ELSE
          GETIN = 0
          LENG = 0
      END IF
      END
