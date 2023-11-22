      SUBROUTINE SYMI(FFULL,FSYM,NBAST,NBAS,ISPIN) 
C
C  THIS ROUTINE SYMMETRY PACKS THE INTERMEDIATE
C  MATRIX SUCH THAT ABACUS ( AND LATER CALCULATOR AND CRAY)
C  ARE ABLE TO DEAL WIth IT
C
C     FFULL ...  UNPACKED FULL INTERMEDIATE MATRIX
C     FSYM ....  SYMMETRY PACKED INTERMEDIATE MATRIX
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
      DIMENSION FFULL(NBAST,NBAST),FSYM(1),NBAS(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
C  REORDER FIRST ROWS AND COLUMNS
C
      CALL REORD(FFULL,FSYM,NBAST,POP(1,ISPIN),
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
        CALL SCOPY(N,FFULL(IOFFR,IOFFC),1,FSYM(IOFFP),1)
        IOFFP=IOFFP+N
10     CONTINUE
C
       IOFFT=IOFFT+N*N
       IOFFR=IOFFR+N
C
100   CONTINUE
C
C  ALL DONE
C
      RETURN
      END
