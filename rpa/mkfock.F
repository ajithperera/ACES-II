      SUBROUTINE MKFOCK(WORK, MAXCOR, IUHF)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION WORK(MAXCOR)
      INTEGER POP, VRT, DIRPRD
      CHARACTER*1 SP(2)
C
#include "machsp.com"
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA SP /'A', 'B'/
C      
      CALL GETREC(20, 'JOBARC', 'NBASTOT ', 1, NBAS)
C
      IFOCK = 1
      ISCR1 = IFOCK + NBAS*NBAS
      ISCR2 = ISCR1 + NBAS*NBAS
      IEND  = ISCR2 + 2*NBAS*NBAS
C      
      DO ISPIN = 1, 1 + MIN(1,IUHF)
         CALL GETREC(20,'JOBARC','FOCK'//SP(ISPIN),NBAS*NBAS, 
     &               WORK(IFOCK))
         CALL AO2MO2(WORK(IFOCK), WORK(IFOCK), WORK(ISCR1), 
     &               WORK(ISCR2), NBAS, NBAS, ISPIN)

#ifdef _DEBUG_LVLM1
         Write(6,*) "The MO basis Fock matrix"
         call output(work(ifock), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
         CALL MKFOCK_BLCK(WORK(IFOCK), WORK(ISCR1), WORK(ISCR2),
     &                    NOCCO(ISPIN), NVRTO(ISPIN), NBAS, ISPIN)
      ENDDO
      
      RETURN
      END
 
