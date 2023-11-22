      SUBROUTINE GTAU(T2,T1A1,T1A2,T1B1,T1B2,DISSYT,NUMSYT,POP1,POP2,
     &               VRT1,VRT2,IRREP,ISPIN,FACT)
C
C THIS SUBROUTINE FORMS THE SYMMETRY PACKED TAU(AB,IJ) 
C AMPLITUDES FOR CCSD GIVEN THE SYMMETRY PACKED T2(AB,IJ)
C AND T1(A,I) AMPLITUDES. TAU(AB,IJ) IS DEFINED AS
C
C  TAU(AB,IJ) = T2(AB,IJ) + T1(A,I)*T1(B,J)*FACT
C                         - T1(A,J)*T1(B,I)*FACT
C
C FOR ISPIN =1 (AAAA CASE) AND ISPIN =2 (BBBB CASE)
C THE EQUATION GIVEN APLIES DIRECTLY. FOR ISPIN=3
C (ABAB CASE) IT REDUCES TO
C
C TAU(Ab,Ij) = T2(Ab,Ij) + T1(A,I)*T1(b,j)*FACT
C
C NOTE THAT SYMMETRY PACKING IS USED AND THAT THE
C SYMMETRY INFORMATION IS ALSO USED IN ORDER TO
C DECIDE IF THERE ARE ANY SINGLE CONTRIBUTION OR NOT.
C FOR THE ABAB SPIN CASE, THERE ARE ONLY CONTRIBUTIONS
C WHEN THE IRREP OF A IS EQUAL TO THE IRREP OF I (AND
C IRREPB EQUAL TO IRREPJ, WHICH IS FORCED BY THE REQUIREMENT
C THAT THE T2 AMPLITUDES ARE TOTAL SYMMETRIC) IN THE AAAA
C AND BBBB SPIN CASES, THERE ARE CONTRIBUTIONS IF EITHERE
C IRREPA EQUALS IRREPI (FIRST TERM) OR IRREPA EQUALS IRREPJ 
C (SECOND TERM)
C
C
C  This routine is a more generalized version of FTAU in that one
C  can use T1 increments of different order.
CEND 
C
C CODED JULY/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LINCC,CICALC
      INTEGER POP1,POP2,VRT1,VRT2,DIRPRD,DISSYT
      DIMENSION T2(DISSYT,NUMSYT),T1A1(1),T1B1(1),T1A2(1),T1B2(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8)
      DIMENSION IOFFT2O(8),IOFFT2V(8),IOFFT1A(8),
     &          IOFFT1B(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /LINEAR/ LINCC,CICALC
C
      NNM1O2(I)=((I-1)*I)/2
C
C      IF(LINCC)RETURN
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
       DO 1 IRREPJ=1,NIRREP-1
        IOFFT1A(IRREPJ+1)=IOFFT1A(IRREPJ)+POP1(IRREPJ)*VRT1(IRREPJ)
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
          NUMB=VRT1(IRREPJ)
C
          INDAB0=IOFFT2V(IRREPJ)
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
           T2(INDAB,INDIJ)=T2(INDAB,INDIJ)+FACT*T1A1(INDAI)*T1A2(INDBJ)
     &                                    -FACT*T1A1(INDAJ)*T1A2(INDBI)
     &                                    +FACT*T1A2(INDAI)*T1A1(INDBJ)
     &                                    -FACT*T1A2(INDAJ)*T1A1(INDBI)
20        CONTINUE
10       CONTINUE
C
        ELSE IF(IRREPI.LT.IRREPJ) THEN
C
         NUMI=POP1(IRREPI)
C
         INDIJ=IOFFT2O(IRREPJ)
C
C
         NUMB=VRT1(IRREPJ)
         NUMA=VRT1(IRREPI)
         INDAB0=IOFFT2V(IRREPJ)
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
C
          DO 50 IB=1,NUMB
           INDBJ=INDJ+IB
           DO 60 IA=1,NUMA
            INDAI=INDI+IA
            INDAB=INDAB+1
            T2(INDAB,INDIJ)=T2(INDAB,INDIJ)+FACT*
     &                     (T1A1(INDAI)*T1A2(INDBJ)
     &                      +T1A2(INDAI)*T1A1(INDBJ))
60         CONTINUE
50        CONTINUE
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
        IOFFT2O(IRREPJ+1)=IOFFT2O(IRREPJ)+POP2(IRREPJ)*POP1(IRREPI)
        IOFFT2V(IRREPJ+1)=IOFFT2V(IRREPJ)+VRT2(IRREPJ)*VRT1(IRREPI)
        IOFFT1A(IRREPJ+1)=IOFFT1A(IRREPJ)+POP1(IRREPJ)*VRT1(IRREPJ)
        IOFFT1B(IRREPJ+1)=IOFFT1B(IRREPJ)+POP2(IRREPJ)*VRT2(IRREPJ)
101    CONTINUE
C
       DO 200 IRREPJ=1,NIRREP

        NUMJ=POP2(IRREPJ)
        INDIJ=IOFFT2O(IRREPJ)
        INDAB0=IOFFT2V(IRREPJ)
        NUMB=VRT2(IRREPJ)
C
        IRREPI=DIRPRD(IRREP,IRREPJ)
        NUMI=POP1(IRREPI)
        NUMA=VRT1(IRREPI)
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
          T2(INDAB,INDIJ)=T2(INDAB,INDIJ)+FACT*(T1A1(INDAI)*T1B2(INDBJ)
     &                                         +T1A2(INDAI)*T1B1(INDBJ))
220      CONTINUE
210     CONTINUE
200    CONTINUE
      ENDIF
      RETURN
      END
