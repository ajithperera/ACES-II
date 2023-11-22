      SUBROUTINE AROUND(HEAD)
      CHARACTER HEAD*(*)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      LNG = LEN(HEAD) + 2
      IND = (72 - LNG)/2 + 1
      WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '! ', HEAD, ' !'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
C     WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '| ', HEAD, ' |'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
      WRITE (LUPRI, '()')
      RETURN
      END
