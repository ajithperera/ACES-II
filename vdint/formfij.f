      SUBROUTINE FORMFIJ(F,CMO,SCR,ISPIN)
C
C FORMS THE OCCUPIED-OCCUPIED BLOCK OF THE FOCK MATRICES
C REQUIRED TO FORM THE INTERMEDIATES I WHICH ARE CONTRACTED
C WITH THE OVERLAP MATRIX DERIVATIVES
C
CEND
C
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       INTEGER POP,VRT,DIRPRD
C
       DIMENSION F(1),CMO(1),SCR(1)
C
       COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
       COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
       COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
       DATA AZERO,ONE /0.D0,1.D0/
C
C   LOOP OVER ALL IRREPS ON THE RIGHT HAND SIDE OF DAO
C
        IOFFC=1
        IOFFAO=1
        IOFFMO=1
        DO 900 IRREP=1,NIRREP
C
C   PERFORM MULTIPLICATION
C
         CALL XGEMM('N','N',NBAS5(IRREP),POP(IRREP,ISPIN),
     &              NBAS5(IRREP),ONE,
     &              F(IOFFAO),NBAS5(IRREP),CMO(IOFFC),
     &              NBAS5(IRREP),AZERO,SCR,NBAS5(IRREP))
C
         CALL XGEMM('T','N',POP(IRREP,ISPIN),POP(IRREP,ISPIN),
     &              NBAS5(IRREP),ONE,
     &              CMO(IOFFC),NBAS5(IRREP),SCR,NBAS5(IRREP),
     &              AZERO,F(IOFFMO),POP(IRREP,ISPIN))
C
         IOFFMO=IOFFMO+POP(IRREP,ISPIN)*POP(IRREP,ISPIN)
         IOFFAO=IOFFAO+NBAS5(IRREP)*NBAS5(IRREP)
         IOFFC=IOFFC+NBAS5(IRREP)*NBAS5(IRREP)
C
900     CONTINUE 
C
       RETURN
       END