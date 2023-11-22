C     
      SUBROUTINE PROJECTCHECK(SCR, MAXCOR, IRREPX, NSIZEC)
C
C This subroutine puts calcualted excitation pattern on the appropriate
C columns of LISTPTRN (375).
C
      IMPLICIT INTEGER (A - Z)
      DOUBLE PRECISION SCR(MAXCOR), ONE
C
      DATA ONE /1.0D0/
C Put projection vectors to appropriate lists
C    
      I000 = 1
      I010 = I000 + NSIZEC
      I020 = I010 + NSIZEC
C
      CALL GETLST(SCR(I000), 1, 1, 1, IRREPX, 375)

C      DO 10 I = 0, NSIZEC-1
C         WRITE(6, *) SCR(I000 + I)
C
C 10   CONTINUE

      CALL GETLST(SCR(I010), 2, 1, 1, IRREPX, 375)
C
C      DO 20 I = 0, NSIZEC-1
C         WRITE(6, *) SCR(I010 + I)
C
C 20   CONTINUE

      CALL SAXPY(NSIZEC, ONE, SCR(I000), 1, SCR(I010), 1)
C
C      DO 30 I = 0, NSIZEC-1
C         WRITE(6, *) SCR(I010 + I)
C
C 30   CONTINUE
C
      CALL GETLST(SCR(I000), 3, 1, 1, IRREPX, 375)
C
      CALL SAXPY(NSIZEC, ONE, SCR(I000), 1, SCR(I010), 1)

      DO 40 I = 0, NSIZEC - 1
         WRITE(6, *) SCR(I010 + I)

 40   CONTINUE
C
      RETURN
      END

