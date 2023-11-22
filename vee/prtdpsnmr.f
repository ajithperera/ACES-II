      SUBROUTINE PRTDPSNMR (DIPSTR, TMNMRL, TMNMRR, ROOT, NTPERT,
     &                      NATOM)
C
C This routine calculate dipole strength over NMR operators.
C S. Ajith Perera 04/96
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION TMNMRL, TMNMRR
      DIMENSION DIPSTR (NTPERT, NTPERT), TMNMRL(NTPERT),
     &           TMNMRR(NTPERT)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/IFLAGS(100)
C
      LENGTH = NTPERT*IINTFP
      CALL GETREC(20, 'JOBARC', 'NMRTMLHT', LENGTH, TMNMRL)
      CALL GETREC(20, 'JOBARC', 'NMRTMRHT', LENGTH, TMNMRR)
C
C Calculate dipole strength for each perturbation
C 
      DO 30 IPERT = 1, NTPERT
C
         DO 40 JPERT = 1, NTPERT
C     
            DIPSTR(IPERT, JPERT) = TMNMRL(IPERT)*TMNMRR(JPERT)/ROOT
C
 40      CONTINUE
 30   CONTINUE
C
C Multiply by appropriate factors, take appropriate tensor averages
C and print the results.
C
      IF (IFLAGS(18) .EQ. 8) THEN
         IFERMI = 0
         ISDIP  = 0
         IPSO   = 1
         IDSO   = 0
      ELSE IF (IFLAGS(18) .EQ. 9) THEN
         IFERMI = 1
         ISDIP  = 0
         IPSO   = 0
         IDSO   = 0
      ELSE IF (IFLAGS(18) .EQ. 10) THEN
         IFERMI = 0
         ISDIP  = 1
         IPSO   = 0
         IDSO   = 0
      ENDIF
C
      CALL FACTOR(DIPSTR, NTPERT, IFERMI, ISDIP, IPSO, IDSO)
C 
      WRITE (6,*)
C
      RETURN
      END
