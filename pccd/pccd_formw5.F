      SUBROUTINE PCCD_FORMW5(ICORE,MAXCOR,IUHF,TERM1,TERM2,TERM3,TERM4,
     &                       TERM5,TERM6,TAU,CCSD)
  
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL TERM1,TERM2,TERM3,TERM4,TERM5,TERM6,TAU,CCSDT2,CCSD
      LOGICAL bRedundant, CC2, HBAR_4LCCSD, NONHF
      DIMENSION ICORE(MAXCOR)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/IFLAGS2(500)
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      IOFFLIST=0
C
C ZERO OUT THESE LISTS BECAUSE ROUTINES USE FANCYPUT 
C
      CALL ACES_COM_SYMLOC
      CALL ZEROLIST(ICORE,MAXCOR,130+IOFFLIST)
      NONHF        = .FALSE.
      NONHF        = (IFLAGS(38) .NE. 0)
C
C REPLACE IAJK INTEGRALS WITH CORRESPONDING HBAR ELEMENTS
C FOR EVERYTHING EXCEPT ROHF-MBPT AND CCSDT-2
C
      CALL MODIAJK(ICORE,MAXCOR,IUHF,ONE)
C
C EVALUATE THE CONTRIBUTION FROM HELL
C
 
      IF(TERM1.OR.TERM5)THEN
         CALL W5INRHF(ICORE,MAXCOR,IUHF,TERM1,TERM5,IOFFLIST)
      ENDIF
      
      IF(TERM2.OR.TERM3.OR.TERM6)THEN
       IF(TERM6)CALL W5TAU(ICORE,MAXCOR,IUHF,TAU,IOFFLIST,CCSD)
       IF(TERM3)CALL W5F(ICORE,MAXCOR,IUHF,IOFFLIST,HBAR_4LCCSD,NONHF)

        CALL W5AB2RHF(ICORE,MAXCOR,IUHF,TERM3,TERM2,TERM6,
     &                         TAU,IOFFLIST,HBAR_4LCCSD,NONHF)
      ENDIF
      CALL MODIAJK(ICORE,MAXCOR,IUHF,ONEM)
C
      RETURN
      END
