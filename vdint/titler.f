      SUBROUTINE TITLER(HEAD,A,IN)
      CHARACTER HEAD*(*), A
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
C
      LENGTH = LEN(HEAD)
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = IN - 100
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//,80A)')
     *      (' ',I=1,INDENT), '.', ('-',I=1,LENGTH), '.'
         WRITE (LUPRI,'(80A)') (' ',I=1,INDENT),'| ', HEAD, ' |'
         WRITE (LUPRI, '(80A)')
     *      (' ',I=1,INDENT), '`', ('-',I=1,LENGTH), ''''
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD, ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      RETURN
      END
