      FUNCTION ORDERED(I,J,A,B,IIRREP,JIRREP,AIRREP,BIRREP,ICASE)
C
      INTEGER I,J,A,B, IIRREP, JIRREP, AIRREP, BIRREP
      LOGICAL ORDERED
C
      ORDERED = .FALSE.
C
      IF (ICASE .EQ. 1) THEN
      ORDERED = .TRUE.
      RETURN
C  AB SPIN CASE
        IF (IIRREP.LT.JIRREP) THEN
          ORDERED = .TRUE.
        ELSEIF (IIRREP.EQ.JIRREP) THEN
          IF (I .LT. J) THEN
            ORDERED = .TRUE.
          ELSEIF (I. EQ. J) THEN
            IF (AIRREP .LT. BIRREP) THEN
              ORDERED = .TRUE.
            ELSEIF (AIRREP .EQ. BIRREP) THEN
              ORDERED = A. LE. B
            ENDIF
          ENDIF
        ENDIF
      ELSEIF (ICASE .EQ. 2 .OR. ICASE .EQ. 3) THEN
C AA AND BB SPINCASE
        IF (IIRREP.LT.JIRREP) THEN
          IF (AIRREP .LT. BIRREP) THEN
            ORDERED = .TRUE.
          ELSEIF (AIRREP.EQ.BIRREP) THEN
            IF (A .LT. B) ORDERED = .TRUE.
          ENDIF
        ELSEIF (IIRREP.EQ.JIRREP) THEN
          IF (I .LT. J) THEN
            IF (AIRREP .LT. BIRREP) THEN
              ORDERED = .TRUE.
            ELSEIF (AIRREP.EQ.BIRREP) THEN
              IF (A .LT. B) ORDERED = .TRUE.
            ENDIF
          ENDIF
        ENDIF
C
      ENDIF
      RETURN
      END
