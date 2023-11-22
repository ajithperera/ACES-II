      FUNCTION MATCH2(I,J,A,B,IIRREP,JIRREP,AIRREP,BIRREP,
     &   ASPIN,BSPIN,INDEX,IRREP,PHTYPE)
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 PHTYPE
      DIMENSION INDEX(2)
      LOGICAL MATCH2
C
      MATCH2 = .FALSE.
      IF (PHTYPE .EQ. 'HOLE') THEN
        IF ( (IIRREP. EQ. IRREP .AND.  I .EQ. INDEX(ASPIN)) .OR.
     &     (JIRREP. EQ. IRREP . AND. J .EQ. INDEX(BSPIN)) ) THEN
          MATCH2 = .TRUE.
        ENDIF
      ELSEIF (PHTYPE .EQ. 'PART' ) THEN
        IF ( (AIRREP. EQ. IRREP .AND.  A .EQ. INDEX(ASPIN)) .OR.
     &     (BIRREP. EQ. IRREP . AND. B .EQ. INDEX(BSPIN)) ) THEN
          MATCH2 = .TRUE.
        ENDIF
      ENDIF
C
      RETURN
      END





