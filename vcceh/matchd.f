C      
      FUNCTION MATCHD(I, J, A, B, IIRREP, JIRREP, ASPIN, BSPIN, 
     &                INDEX, IRREP)
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION INDEX(2)
      LOGICAL MATCHD
C
      MATCHD = .FALSE.
C
      IF ((IIRREP .EQ. IRREP .AND. I .EQ. INDEX(ASPIN)) .OR.
     &    (JIRREP .EQ. IRREP .AND. J .EQ. INDEX(BSPIN)) ) THEN
         MATCHD = .TRUE.
      ENDIF
C     
      RETURN
      END





