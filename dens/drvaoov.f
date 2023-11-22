C DRIVER FOR AO-BASED I(ai) = G(aecd)*<ie||cd> contraction
C modified slightly from DRAOLAD
C
c Started by PR 12/2003, but left incomplete. Substantial
C amount of work was need to salvage this to a correctly functioning
C driver for the AO-based formation of I(ae).  Ajith Perera, 04/2005.
C
C
      SUBROUTINE DRVAOOV(AIOV,ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      LOGICAL TAU,MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &        QCISD,ROHF4,ITRFLG,LAMBDA,UCC,CCD,CCSD,M4DQ,M4SDQ,
     &        M4SDTQ
      DIMENSION AIOV(*), IT1OFF(2)
      INTEGER POP, VRT
      COMMON /SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
C
      CC   = CCD .OR. CCSD .OR. QCISD
      MBPT4= M4DQ .OR. M4SDQ .OR. M4SDTQ
      
      IF (MBPT3 .OR. UCC) THEN
         NTIMES = 1
      ELSE 
         NTIMES = 2
      ENDIF
C
      TAU = .FALSE. 
      IF (CC) THEN
         TAU = .TRUE. 
         IT1OFF(1) = 1
         IT1OFF(2) = IT1OFF(1) + IINTFP*NTAA*IUHF
         IT1END    = IT1OFF(2) + IINTFP*NTBB 
      ELSE
         IT1END    = 1 
      ENDIF 
C
      CALL Z2INIOV(AIOV,ICORE,MAXCOR,IUHF,TAU,IT1OFF,IT1END, 
     &             NTIMES)

      RETURN
      END
