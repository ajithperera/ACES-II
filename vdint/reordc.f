
      SUBROUTINE REORDC(A,SCR,NBAST,POP,VRT,NBAS)
C
C   THIS SUBROUTINE REORDERS A GIVEN MATRIX IN SUCH A WAY THAT
C   ABACUS (CALCULATOR AND IN THE FUTURE CRAY) CAN HANDLE THE
C   BACK TRANSFORMATION OF THE ONE-ELECTRON DENSITY MATRIX
C
C   A ... .... INPUT ARRAY
C   SCR ...... OUTPUT ARRAY
C   NBAST .... TOTAL NUMBER OF BASIS FUNCTIONS
C   POP ...... NUMBER OF OCCUPIED ORBITALS PER IRREP
C   VRT ...... NUMBER OF VIRTUAL ORBITALS PER IRREP
C   NBAS ..... NUMBER OF BASIS FUNCTIONS PER IRREP
C
CEND
C
C CODED OCT/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION A(NBAST,NBAST),SCR(NBAST,NBAST)
      DIMENSION POP(8),VRT(8),NBAS(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C COPY FIRT CFULL TO SCRATCH
C
c YAU : old
c     CALL ICOPY(NBAST*NBAST*IINTFP,A,1,SCR,1)
c YAU : new
      CALL DCOPY(NBAST*NBAST,A,1,SCR,1)
c YAU : end
C
C  OFFSET FOR OCCUPIED AND VIRTUAL BLOCK
C
      IOFFO=1
      IOFFV=1
      DO 10 IRREP=1,NIRREP
       IOFFV=IOFFV+POP(IRREP)*NBAST
10    CONTINUE
C
C  OFFSET FOR TARGET ARRAY
C
      IOFFT=1
C
C  LOOP OVER ALL IRREPS
C
      DO 100 IRREP=1,NIRREP
C
C FILL FIRST WITH OCCUPIED ORBITALS OF THIS BLOCK
C
c YAU : old
c      CALL ICOPY(IINTFP*POP(IRREP)*NBAST,SCR(IOFFO,1),1,A(IOFFT,1),1)
c YAU : new
       CALL DCOPY(POP(IRREP)*NBAST,SCR(IOFFO,1),1,A(IOFFT,1),1)
c YAU : end
C
       IOFFO=IOFFO+NBAST*POP(IRREP)
       IOFFT=IOFFT+NBAST*POP(IRREP)
C
C FILL NOW WITH VIRTUAL ORBITALS OF THIS BLOCK
C
c YAU : old
c      CALL ICOPY(IINTFP*VRT(IRREP)*NBAST,SCR(IOFFV,1),1,A(IOFFT,1),1)
c YAU : new
       CALL DCOPY(VRT(IRREP)*NBAST,SCR(IOFFV,1),1,A(IOFFT,1),1)
c YAU : end
C
       IOFFV=IOFFV+NBAST*VRT(IRREP)
       IOFFT=IOFFT+NBAST*VRT(IRREP)
C
100   CONTINUE
C
C  ALL DONE, RETURN
C
      RETURN
      END
