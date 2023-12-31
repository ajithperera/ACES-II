
      SUBROUTINE REORD(A,SCR,NBAST,POP,VRT,NBAS)
C
C  THIS SUBROUTINE REORDERS A GIVEN MATRIX IN SUCH A WAY THAT
C  ABACUS ( CALCULATOR AND IN THE FUTURE CRAY) CAN HANDLE THE
C  BACK TRANSFORMATION OF THE ONE ELECTRON DENSITY MATRIX
C
C  A ..... INPUT AND OUTPUT MATRIX 
C  SCR ... SCRATCH ARRAY
C  NBAST . TOTAL NUMBER OF BASIS FUNCTION
C  POP ... NUMBER OF OCCUPIED ORBITALS PER IRREP
C  VRT ... NUMBER OF VIRTUAL ORBITALS PER IRREP
C  NBAS .. NUMBER OF BASIS FUNCTIONS PER IRREP
C
CEND
C
C  CODED OCT/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION A(NBAST,NBAST),SCR(NBAST,NBAST)
      DIMENSION POP(8),VRT(8),NBAS(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C  REORDER FIRST COLUMS
C
C
C  OFFSET FOR OCCUPIED AND VIRTUAL BLOCKS
C
      IOFFO=1
      IOFFV=1
      JOFFV=1
      DO 10 IRREP=1,NIRREP   
       JOFFV=JOFFV+POP(IRREP) 
       IOFFV=IOFFV+POP(IRREP)*NBAST
10    CONTINUE
C
      IOFFT=1
C
C  LOOP OVER ALL IRREPS
C
      DO 100 IRREP=1,NIRREP
C
C  FIRST FILL WITh OCCUPIED ORBITALS OF THIS BLOCK
C
c YAU : old
c      CALL ICOPY(IINTFP*POP(IRREP)*NBAST,A(IOFFO,1),1,SCR(IOFFT,1),1)
c YAU : new
       CALL DCOPY(POP(IRREP)*NBAST,A(IOFFO,1),1,SCR(IOFFT,1),1)
c YAU : end
C
       IOFFO=IOFFO+NBAST*POP(IRREP)
       IOFFT=IOFFT+NBAST*POP(IRREP)
C
C  FILL NOW WITH THE VIRTUAL BLOCK OF THIS IRREP
C
c YAU : old
c      CALL ICOPY(IINTFP*VRT(IRREP)*NBAST,A(IOFFV,1),1,SCR(IOFFT,1),1)
c YAU : new
       CALL DCOPY(VRT(IRREP)*NBAST,A(IOFFV,1),1,SCR(IOFFT,1),1)
c YAU : end
C       
       IOFFV=IOFFV+NBAST*VRT(IRREP)
       IOFFT=IOFFT+NBAST*VRT(IRREP)
C
100   CONTINUE
C
C  NOW DEAL WITH THE ROWS
C
      IOFFO=1
      IOFFV=JOFFV
C
      IOFFT=1
C
C  LOOP OVER ALL IRREPS
C
      DO 200 IRREP=1,NIRREP
C
C  FIRST FILL WITH THE OCCUPIED BLOCKS OF THIS IRREP
C
       DO 150 I=1,POP(IRREP)
C
        CALL SCOPY(NBAST,SCR(IOFFO,1),NBAST,A(IOFFT,1),NBAST)
C
        IOFFO=IOFFO+1
        IOFFT=IOFFT+1
C
150    CONTINUE
C
C  NOW FILL WITH THE VIRTUAL BLOCK OF THIS IRREP
C
       DO 160 IA=1,VRT(IRREP)
C
        CALL SCOPY(NBAST,SCR(IOFFV,1),NBAST,A(IOFFT,1),NBAST)
C
        IOFFV=IOFFV+1
        IOFFT=IOFFT+1
C
160    CONTINUE
C
200   CONTINUE 
C
C  ALL DONE, RETURN
C
      RETURN
      END
