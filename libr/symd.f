      SUBROUTINE SYMD(DFULL,DSYM,NBAST,NBAS,ISPIN) 
C
C  THIS ROUTINE SYMMETRY PACKS THE DENSITY (AND INTERMEDIATE)
C  MATRIX SUCH THAT ABACUS ( AND LATER CALCULATOR AND CRAY)
C  ARE ABLE TO DEAL WITH IT
C
C     DFULL ...  UNPACKED FULL DENSITY (INTERMEDIATE) MATRIX
C     DSYM ....  SYMMETRY PACKED DENSITY (INTERMEDIATE) MATRIX
C     NBAST ...  TOTAL NUMBER OF BASIS FUNCTIONS
C     NBAS ....  NUMBER OF BASIS FUNCTIONS PER IRREP
C     ISPIN ...  SPIN CASE
C
CEND
C
C  CODED OCT/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION DFULL(NBAST,NBAST),DSYM(1),NBAS(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
C  REORDER FIRST ROWS AND COLUMNS
C
      CALL REORD(DFULL,DSYM,NBAST,POP(1,ISPIN),
     &           VRT(1,ISPIN),NBAS)
C
C  NOW DEAL WITH THE SYMMETRY PACKING !
C
C  SET OFFSETS
C
      IOFFR=1
      IOFFC=0
      IOFFT=1

C
C LOOP OVER IRREPS
C
      DO 100 IRREP=1,NIRREP
C
       N=NBAS(IRREP)
       IOFFP=IOFFT
       DO 10 IMO=1,N
        IOFFC=IOFFC+1
        CALL SCOPY(N,DFULL(IOFFR,IOFFC),1,DSYM(IOFFP),1)
        IOFFP=IOFFP+N
10     CONTINUE
C
       CALL SQUEEZ(DSYM(IOFFT),N,0)
C
       IOFFT=IOFFT+N*(N+1)/2
       IOFFR=IOFFR+N
C
100   CONTINUE
C
C  ALL DONE
C
    
      RETURN
      END