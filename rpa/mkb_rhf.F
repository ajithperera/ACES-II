      SUBROUTINE MKB_RHF(AMAT,SCRATCH,IRREP)
C
C   THIS ROUTINE MAKES THE FULL SPIN ORBITAL A MATRIX
C   TO TEST RHF->UHF STABILITY FOR A PARTICULAR DPD IRREP.
C
C   NOTE, THAT ONLY ONE BLOCK OF THE A MATRIX IS COMPUTED 
C
C   THE FORMULAS ARE GIVEN BY
C
C   A(AI,BJ) = <AB//IJ> + <AJ//IB>      AA-BLOCK
C
C   A(ai,bj) = <ab//ij> + <aj//ib>      BB-BLOCK
C
C   A(AI,bj) = 2  <Ab//Ij>              AB-BLOCK
C
C   A(ai,BJ) = 2  <aB//iJ>              BA-BLOCK (TRANPOSE OF AB-BLOCK)
C  
C AMAT HOLDS ALL FOUR PARTS OF THE A MATRIX 
C SCRATCH IS A SCRATCH ARRAY REQUIRED FOR CONTRUCTING A
C IRREP IS THE CORRESPONDING IRREP OF THE FORMED A MATRIX
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYW,DIRPRD,POP,VRT
      DIMENSION AMAT(1),SCRATCH(1)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NJUNK(6)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SHIFT/ ISHIFT,NDRGEO
C
      DATA ONE,ONEM,TWO /1.0D+0,-1.0D+0,2.0D+0/
C
C     FORM FIRST THE AA AND BB-BLOCK
C
C  DEAL WITH AA INTEGRALS
C
C
C FIRST READ IN THE INTEGRALS <AB//IJ> (ORDERING AI,BJ)
C THEY ARE STORED ON LIST 19 AND THEN ADD IN THE <AI//BJ> INTS
C
      LISTW=19 + ISHIFT 
      NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
C
      IOFF=1
C
      DO 1001 IDIS=1,NUMSYW
       CALL GETLST(AMAT(IOFF),IDIS,1,1,IRREP,LISTW)
#ifdef _DEBUG_LVLM1
       Write(6,"(a)") "The list 19"
       Write(6, "(F10.8)") (AMAT(IOFF-1+I),I=1, DISSYW)
#endif
       IOFF=IOFF+DISSYW
1001  CONTINUE
C
C   NOW FORM THE AB BLOCK (THE BA BLOCK IS THE TRANSPOSED AS LONG
C   AS WE NOT SCALE WITH THE ORBITAL ENERGY DENOMINATOR
C
    
C#ifdef _NOSKIP
      LISTW=18 + ISHIFT 
      NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
C
      IOFF  = NUMSYW*DISSYW+1
      IOFF1 = 1
      DO 1002 IDIS=1,NUMSYW
       CALL GETLST(AMAT(IOFF),IDIS,1,1,IRREP,LISTW)
#ifdef _DEBUG_LVLM1
       Write(6,"(a)") "The list 18"
       Write(6, "(F10.8)") (AMAT(IOFF-1+I),I=1, DISSYW)
#endif
       CALL DAXPY(DISSYW, ONE, AMAT(IOFF), 1, AMAT(IOFF1),1)
       IOFF  = IOFF  + DISSYW
       IOFF1 = IOFF1 + DISSYW 
1002  CONTINUE
C#endif
C
C  ALL DONE, RETURN
C
      RETURN
      END
