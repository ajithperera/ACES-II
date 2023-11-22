      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
      PARAMETER ( LUPRI = 6 )
C
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(/,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      RETURN
      END
