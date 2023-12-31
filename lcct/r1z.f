      SUBROUTINE R1Z(G,R1A,R1B,Z1A,Z1B,DISSYG,NUMSYG,POP1,POP2,
     &               VRT1,VRT2,IRREP,IRREPX,ISPIN)
C
C THIS SUBROUTINE FORMS THE SYMMETRY PACKED G(AB,IJ) 
C DENSITY MATRIX CONTRIBUTION FOR EOM-CCSD GRADIENTS
C GIVEN THE T1(J,B) AND Z(I,A) INTERMEDIATE
C
C  G(AB,IJ) = R(A,I)*Z(B,J)- R(A,J)*Z(B,I)
C
C             + Z(A,I)*R(B,J)- Z(A,J)*R(B,I)
C
C FOR ISPIN =1 (AAAA CASE) AND ISPIN =2 (BBBB CASE)
C THE EQUATION GIVEN APLIES DIRECTLY. FOR ISPIN=3
C (ABAB CASE) IT REDUCES TO
C
C G(Ab,Ij) = R(A,I)*Z(b,j)
C
C            + Z(A,I)*R(b,j)
C
CEND 
C
C CODED OCTOBER/93 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP1,POP2,VRT1,VRT2,DIRPRD,DISSYG
      DIMENSION G(DISSYG,NUMSYG),R1A(1),R1B(1),Z1A(1),Z1B(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8)
      DIMENSION IOFFT2O(8),IOFFT2V(8),IOFFT1A(8),
     &          IOFFT1B(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      NNM1O2(I)=((I-1)*I)/2
C
      IF(ISPIN.LT.3) THEN
C
C AAAA OR BBBB SPIN CASE
C
C DETERMINE FIRST OFFSETS FOR T2 AND T1
C
       IOFFT2O(1)=0
       IOFFT2V(1)=0
       IOFFT1A(1)=0
       DO 1 IRREPJR=1,NIRREP-1
        IRREPJL=DIRPRD(IRREPX,IRREPJR)
        IOFFT1A(IRREPJR+1)=IOFFT1A(IRREPJR)+POP1(IRREPJR)*VRT1(IRREPJL)
1      CONTINUE
       IF(IRREP.EQ.1) THEN
        DO 2 IRREPJ=1,NIRREP-1
         IOFFT2O(IRREPJ+1)=IOFFT2O(IRREPJ)+NNM1O2(POP1(IRREPJ))
         IOFFT2V(IRREPJ+1)=IOFFT2V(IRREPJ)+NNM1O2(VRT1(IRREPJ))
2       CONTINUE
       ELSE
        DO 3 IRREPJ=1,NIRREP-1
         IRREPI=DIRPRD(IRREP,IRREPJ)
         IF(IRREPI.LT.IRREPJ) THEN
          IOFFT2O(IRREPJ+1)=IOFFT2O(IRREPJ)+POP1(IRREPI)*POP1(IRREPJ)
          IOFFT2V(IRREPJ+1)=IOFFT2V(IRREPJ)+VRT1(IRREPI)*VRT1(IRREPJ)
         ELSE
          IOFFT2O(IRREPJ+1)=IOFFT2O(IRREPJ)
          IOFFT2V(IRREPJ+1)=IOFFT2V(IRREPJ)
         ENDIF
3       CONTINUE
       ENDIF
C
       DO 100 IRREPJ=1,NIRREP
C
        NUMJ=POP1(IRREPJ)
        INDIJ=IOFFT2O(IRREPJ)
        INDJ0=IOFFT1A(IRREPJ)
        IRREPI=DIRPRD(IRREP,IRREPJ)
        IF(IRREPI.EQ.IRREPJ) THEN
C
         IF(NNM1O2(NUMJ).EQ.0) GO TO 100
C
          IRREPB=DIRPRD(IRREPJ,IRREPX)
          IRREPA=DIRPRD(IRREPB,IRREP)
C
          NUMB=VRT1(IRREPB)
C
          INDAB0=IOFFT2V(IRREPB)
C
         DO 10 J=2,NUMJ
         INDJ=INDJ0+(J-1)*NUMB
         DO 10 I=1,J-1
          INDI=INDJ0+(I-1)*NUMB
C
          INDIJ=INDIJ+1
C 
          INDAB=INDAB0
C
          DO 20 IB=2,NUMB
          INDBI=INDI+IB
          INDBJ=INDJ+IB
          DO 20 IA=1,IB-1
           INDAI=INDI+IA
           INDAJ=INDJ+IA
           INDAB=INDAB+1
C
           G(INDAB,INDIJ)=G(INDAB,INDIJ)+R1A(INDAI)*Z1A(INDBJ)
     &                                  -R1A(INDAJ)*Z1A(INDBI)
     &                                  +Z1A(INDAI)*R1A(INDBJ)
     &                                  -Z1A(INDAJ)*R1A(INDBI)
20        CONTINUE
10       CONTINUE
C
        ELSE IF(IRREPI.LT.IRREPJ) THEN
C
         NUMI=POP1(IRREPI)
C
         INDIJ=IOFFT2O(IRREPJ)
C
         IRREPB=DIRPRD(IRREPJ,IRREPX)
         IRREPA=DIRPRD(IRREPB,IRREP)
C
         NUMB=VRT1(IRREPB)
         NUMA=VRT1(IRREPA)
         IF(IRREPA.LT.IRREPB) THEN
          INDAB0=IOFFT2V(IRREPB)
         ELSE
          INDAB0=IOFFT2V(IRREPA)
         ENDIF
C
         DO 30 J=1,NUMJ
         INDJ=IOFFT1A(IRREPJ)+(J-1)*NUMB
         DO 30 I=1,NUMI
C
          INDIJ=INDIJ+1         
C
          INDI=IOFFT1A(IRREPI)+(I-1)*NUMA
C
          INDAB=INDAB0
C
          IF(IRREPA.LT.IRREPB) THEN
C
          DO 50 IB=1,NUMB
           INDBJ=INDJ+IB
           DO 60 IA=1,NUMA
            INDAI=INDI+IA
            INDAB=INDAB+1
            G(INDAB,INDIJ)=G(INDAB,INDIJ)+R1A(INDAI)*Z1A(INDBJ)
     &                                   +Z1A(INDAI)*R1A(INDBJ)
60         CONTINUE
50        CONTINUE
C
          ELSE
C
          DO 150 IA=1,NUMA
           INDAI=INDI+IA
           DO 160 IB=1,NUMB
            INDBJ=INDJ+IB
            INDAB=INDAB+1
            G(INDAB,INDIJ)=G(INDAB,INDIJ)-R1A(INDAI)*Z1A(INDBJ)
     &                                   -Z1A(INDAI)*R1A(INDBJ)
160         CONTINUE
150        CONTINUE
C
          ENDIF
30       CONTINUE
        ENDIF
100    CONTINUE
      ELSE
C
C  ABAB SPIN CASE
C
C  GET FIRST OFFSETS OF T1 AND T2
C
       IOFFT2O(1)=0
       IOFFT2V(1)=0
       IOFFT1A(1)=0
       IOFFT1B(1)=0
       DO 101 IRREPJ=1,NIRREP-1
        IRREPI=DIRPRD(IRREPJ,IRREP)
        IRREPJL=DIRPRD(IRREPX,IRREPJ)
        IOFFT2O(IRREPJ+1)=IOFFT2O(IRREPJ)+POP2(IRREPJ)*POP1(IRREPI)
        IOFFT2V(IRREPJ+1)=IOFFT2V(IRREPJ)+VRT2(IRREPJ)*VRT1(IRREPI)
        IOFFT1A(IRREPJ+1)=IOFFT1A(IRREPJ)+POP1(IRREPJ)*VRT1(IRREPJL)
        IOFFT1B(IRREPJ+1)=IOFFT1B(IRREPJ)+POP2(IRREPJ)*VRT2(IRREPJL)
101    CONTINUE
C
       DO 200 IRREPJ=1,NIRREP

        IRREPB=DIRPRD(IRREPX,IRREPJ)
        NUMJ=POP2(IRREPJ)
        INDIJ=IOFFT2O(IRREPJ)
        INDAB0=IOFFT2V(IRREPB)
        NUMB=VRT2(IRREPB)
C
        IRREPI=DIRPRD(IRREP,IRREPJ)
        IRREPA=DIRPRD(IRREP,IRREPB)
        NUMI=POP1(IRREPI)
        NUMA=VRT1(IRREPA)
C
        DO 210 J=1,NUMJ
        INDJ=IOFFT1B(IRREPJ)+(J-1)*NUMB
        DO 210 I=1,NUMI
         INDI=IOFFT1A(IRREPI)+(I-1)*NUMA
         INDIJ=INDIJ+1
         INDAB=INDAB0
         DO 220 IB=1,NUMB
         INDBJ=INDJ+IB
         DO 220 IA=1,NUMA
          INDAB=INDAB+1
          INDAI=INDI+IA
          G(INDAB,INDIJ)=G(INDAB,INDIJ)+R1A(INDAI)*Z1B(INDBJ)
     &                                 +Z1A(INDAI)*R1B(INDBJ)
220      CONTINUE
210     CONTINUE
200    CONTINUE
      ENDIF
      RETURN
      END
