C
      SUBROUTINE FIXDBLCONT(SCR, MAXCOR, IRREPX, NSIZEC)
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(MAXCOR), ZERO
C
      DATA ONE, ZERO /1.0D0, 0.0D0/
C    
      I000 = 1
      I010 = I000 + NSIZEC
      I020 = I010 + NSIZEC
C
C Let's get the core, valence vectors to the memory,
C
      CALL GETLST(SCR(I000), 1, 1, 1, IRREPX, 375)
      CALL GETLST(SCR(I010), 2, 1, 1, IRREPX, 375)
C 
      DO 10 ICOUNT = 0, NSIZEC - 1 
C         
         IF (SCR(I000 + ICOUNT) .EQ. SCR(I010 + ICOUNT)) THEN
            IF (SCR(I010 + ICOUNT) .EQ. ONE) SCR(I010 + ICOUNT) = ZERO
         ENDIF
C
 10   CONTINUE
C
      CALL PUTLST(SCR(I010), 2, 1, 1, IRREPX, 375)
C     
C Now get the inner-valence  vector to the memory.
C
      CALL GETLST(SCR(I010), 3, 1, 1, IRREPX, 375)

      DO 20 ICOUNT = 0, NSIZEC - 1 
C         
         IF (SCR(I000 + ICOUNT) .EQ. SCR(I010 + ICOUNT)) THEN
            IF (SCR(I010 + ICOUNT) .EQ. ONE) SCR(I010 + ICOUNT) = ZERO
         ENDIF
C
 20   CONTINUE
C
C Now get the modified valence vectors to the memory
C
      CALL GETLST(SCR(I000), 2, 1, 1, IRREPX, 375)
C
      DO 30 ICOUNT = 0, NSIZEC - 1 
C         
         IF (SCR(I000 + ICOUNT) .EQ. SCR(I010 + ICOUNT)) THEN
            IF (SCR(I010 + ICOUNT) .EQ. ONE) SCR(I010 + ICOUNT) = ZERO
         ENDIF
C
 30   CONTINUE
C
      CALL PUTLST(SCR(I010), 3, 1, 1, IRREPX, 375)
C
      RETURN
      END

