      SUBROUTINE H4INX2(ICORE,MAXCOR,IUHF)
C
C     SUBROUTINE CALLS H4X2ALL FOR THE THREE SPIN CASES
C     
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DEBUG
      DIMENSION ICORE(MAXCOR)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC

#include "parcc.com"
C
      IBOT=1
      IF(IUHF.EQ.0) IBOT=3
C
      IF (Ispar .AND. Coulomb) Then

      Write(6,"(a)") " Entering DCC H4X2ALL Blocks"

         DO ISPIN=IBOT,3
C            CALL H4X2ALL_4PDCC_1(ICORE,MAXCOR,ISPIN,IUHF,CCSD,
C     &                           .FALSE.)

            CALL H4X2ALL_4PDCC_11(ICORE,MAXCOR,ISPIN,IUHF,CCSD,
     &                            .FALSE.)
            CALL H4X2ALL_4PDCC_12(ICORE,MAXCOR,ISPIN,IUHF,CCSD,
     &                            .FALSE.)
            CALL H4X2ALL_4PDCC_2(ICORE,MAXCOR,ISPIN,IUHF,CCSD,
     &                           .TRUE.)
            CALL H4X2ALL_4PDCC_3(ICORE,MAXCOR,ISPIN,IUHF,CCSD,
     &                           .FALSE.)
         ENDDO 

      ELSE

      DO 100 ISPIN=IBOT,3
            CALL H4X2ALL(ICORE,MAXCOR,ISPIN,IUHF,CCSD)
 100  CONTINUE

      ENDIF 
C
CSSS      Call check_phph(Icore,Maxcor,Iuhf)
      RETURN
      END
