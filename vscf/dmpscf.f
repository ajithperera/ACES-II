      SUBROUTINE DMPSCF(FNEW,FSCR,LENF,ETOT,ITER,IUHF,DAMP,DMPFLG,
     &                  DMPTOL,DEAVG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CVGING,DMPFLG
C
      DIMENSION FNEW((IUHF+1)*LENF),FSCR((IUHF+1)*LENF)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /DMPCOM/ ITDLST
C
      DATA CVGTOL /0.05D+00/
C
C     Highly daring and experimental damping on F.
C
      PT2    = 0.2D+00
      TWOPT2 = 2.2D+00
C
      E0 = ETOT
      IF(ITER.EQ.1)THEN
        CALL PUTREC(20,'JOBARC','E0      ',IINTFP,E0)
c        CALL PUTREC(20,'JOBARC','FA0     ',IINTFP*ITRILN(NIRREP+1),
c     &              ICORE(I040))
        CALL PUTREC(20,'JOBARC','FA0     ',IINTFP*LENF,FNEW)
       IF(IUHF.EQ.1)THEN
c        CALL PUTREC(20,'JOBARC','FB0     ',IINTFP*ITRILN(NIRREP+1),
c     &              ICORE(I040 + IINTFP*ITRILN(NIRREP+1)))
        CALL PUTREC(20,'JOBARC','FB0     ',IINTFP*LENF,FNEW(1+LENF))
       ENDIF
        DE    = 0.0D+00
        DEP   = 0.0D+00
        DEAVG = 0.0D+00
      ENDIF
C
C     ----- End of ITER = 1 block -----
C
      IF(ITER.EQ.2)THEN
        CALL GETREC(20,'JOBARC','E0      ',IINTFP,E1)
        CALL PUTREC(20,'JOBARC','E1      ',IINTFP,E1)
        CALL PUTREC(20,'JOBARC','E0      ',IINTFP,E0)
C
        CALL GETREC(20,'JOBARC','FA0     ',IINTFP*LENF,FSCR)
        CALL PUTREC(20,'JOBARC','FA1     ',IINTFP*LENF,FSCR)

        IF(IUHF.EQ.1)THEN
          CALL GETREC(20,'JOBARC','FB0     ',IINTFP*LENF,FSCR(1+LENF))
          CALL PUTREC(20,'JOBARC','FB1     ',IINTFP*LENF,FSCR(1+LENF))
        ENDIF
C
        CALL PUTREC(20,'JOBARC','FA0     ',IINTFP*LENF,FNEW)
C
        IF(IUHF.EQ.1)THEN
          CALL PUTREC(20,'JOBARC','FB0     ',IINTFP*LENF,FNEW(1+LENF))
        ENDIF
C
        DE    = E0 - E1
        DEP   = 0.0D+00
        DEAVG = DABS(DE)
      ENDIF
C
C     ----- End of ITER = 2 block -----
C
      IF(ITER.GT.2)THEN
        CALL GETREC(20,'JOBARC','E1      ',IINTFP,E2)
        CALL PUTREC(20,'JOBARC','E2      ',IINTFP,E2)
        CALL GETREC(20,'JOBARC','E0      ',IINTFP,E1)
        CALL PUTREC(20,'JOBARC','E1      ',IINTFP,E1)
        CALL PUTREC(20,'JOBARC','E0      ',IINTFP,E0)
C
        CALL GETREC(20,'JOBARC','FA0     ',IINTFP*LENF,FSCR)
        CALL PUTREC(20,'JOBARC','FA1     ',IINTFP*LENF,FSCR)
C
        IF(IUHF.EQ.1)THEN
        CALL GETREC(20,'JOBARC','FB0     ',IINTFP*LENF,FSCR(1+LENF))
        CALL PUTREC(20,'JOBARC','FB1     ',IINTFP*LENF,FSCR(1+LENF))
        ENDIF
C
CSSS      Write(*,*) E0
CSSS      Write(*,*) E1
CSSS      Write(*,*) E2

      DE  = E0 - E1
      DEP = E1 - E2

      DEAVG = (DABS(DE) + DABS(DEP) + PT2*DEAVG)/TWOPT2
C
      IF(IFLAGS2(109).EQ.1)THEN
        CALL DAMPD(DE,DEP,DEAVG,DAMP)
CSSS      WRITE(LUOUT,5350) DAMP
CSSS 5350 FORMAT(T8,'Dynamical damping parameter is ',F15.10)

      ENDIF
C check extra control here (cvging)
      IF(DAMP.GE.DMPTOL)THEN
      CALL SCFDMP(FNEW,FSCR,DAMP,LENF)
C
      IF(IUHF.EQ.1)THEN
      CALL SCFDMP(FNEW(1+LENF),FSCR(1+LENF),DAMP,LENF)
      ENDIF
C
      ENDIF
C
      CALL PUTREC(20,'JOBARC','FA0     ',IINTFP*LENF,FNEW)
C
      IF(IUHF.EQ.1)THEN
      CALL PUTREC(20,'JOBARC','FB0     ',IINTFP*LENF,FNEW(1+LENF))
      ENDIF
C
      ENDIF
C
C     ----- End of ITER > 2 block -----
C
      CVGING = DABS(DE).LT.CVGTOL
      DMPFLG = IFLAGS2(109).GT.0.AND.(DAMP.GE.DMPTOL .OR. ITER.LE.2 
     &                          . OR. .NOT.CVGING)
      IF(DMPFLG)THEN
        ITDLST     = ITER + IFLAGS(20)
      ENDIF
C
      RETURN
      END
