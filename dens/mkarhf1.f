      SUBROUTINE MKARHF1(AMAT,SCRATCH,ICORE,MAXCOR,bRedundant)
C
C   THIS ROUTINE CALCULATES THE A-MATRIX USED IN THE
C   SOLUTION OF THE CPHF EQUATION FOR RHF WAVE FUNCTIONS
C   NOTE, ONLY THE BLOCK OF THE A MATRIX CORRESPONDING TO
C   IRREP 1 IS FORMED. FURTHER, FOR RHF WE ARE USING
C   A SPIN-ADAPTED VERSION. THE FORMULA IS GIVEN BY
C
C   A(AI,BJ) = <AB//IJ> + <AJ//IB> + <Ab//Ij> + <Aj//Ib>
C
C            = <AB//IJ> - <JA//IB> + 2 <Ab//Ij>
C
C
C   AMAT : WILL RETURN THE CORRESPONDING PART OF THE A MATRIX
C   SCRATCH : SCRATCH VECTOR REQUIRED TO CONSTRUCT THE A MATRIX
C
CEND
C
C   CODED JULY/90   JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,POP,VRT,reflist
      LOGICAL bRedundant
      DIMENSION AMAT(1),SCRATCH(1),icore(1)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP2/IRPDPD(8,22)
      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM2/POP(8,2),VRT(8,2),NJUNK(6)
      common /dropgeo/ ndrgeo
      COMMON /SHIFT/ ISHIFT 
C
      DATA ONE,ONEM,TWO /1.0D+0,-1.0D+0,2.0D+0/
C
C   READ IN FIRST THE INTEGRALS <AB//IJ> (ORDERING AI,BJ) 
C   THEY ARE STORED ON LIST 19
C
      LISTW=19  + ISHIFT 

      NUMSYW=IRPDPD(1,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(1,ISYTYP(1,LISTW))
C
      IF(bRedundant) THEN
        CALL GETLST(AMAT,1,NUMSYW,2,1,LISTW)
      ELSE
        CALL GETLST_NR(AMAT,ICORE(1),MAXCOR,LISTW,1)
      ENDIF
      
C
C   NOW ADD THE INTEGRALS <AJ//BI> (ORDERING AI,BJ)
C   THEY ARE STORED ON LIST 23, HOWEVER TAKE CARE OF
C   THE SIGN, SINCE THEY ARE ACTUALLY <JA//BI>
C
      LISTW=23 + ISHIFT 
      NUMSYW=IRPDPD(1,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(1,ISYTYP(1,LISTW))
C
      CALL GETLST(SCRATCH,1,NUMSYW,2,1,LISTW)
C
      CALL SAXPY(NUMSYW*DISSYW,ONEM,SCRATCH,1,AMAT,1)
C
C   NOW ADD THE INTEGRALS <Ab//Ij> (ORDERING AI,bj)
C   WITH A FACTOR 2 TO A
C   THEY ARE STORED ON LIST 18
C
      LISTW=18 + ISHIFT 

      NUMSYW=IRPDPD(1,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(1,ISYTYP(1,LISTW))
C
      IF(bRedundant) THEN
       CALL GETLST(SCRATCH,1,NUMSYW,2,1,LISTW)
      ELSE
       CALL GETLST_NR(SCRATCH,ICORE(1),MAXCOR,LISTW,1)
      ENDIF

C
      CALL SAXPY(NUMSYW*DISSYW,TWO,SCRATCH,1,AMAT,1)    
C
C  ALL DONE
C
      RETURN
      END      
