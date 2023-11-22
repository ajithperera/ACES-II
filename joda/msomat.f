      SUBROUTINE MSOMAT(IAOM,IAOR,IAOC,ISOM,ISOR,ISOC,ICHAO,ICHSO,NAT,
     +      NGENG)
C
C     CONSTRUCT AO --> SO MATRIXES
C
      IMPLICIT INTEGER (A-Z)
C
      INTEGER IAOM(IAOR,IAOC),ISOM(ISOR,ISOC),ICHAO(3,IAOC),ICHSO
     1 (ISOR),ISCR(15,3),ISCS(3)
C     SCRATCH VECTOR ISCR SUFFICIENT UP TO G-AO
      COMMON/IB/IBING(3,3),ICHAT(8,8),ICHGEN(8,3),IMAGE(200,3),NEQA
      COMMON/JJJ/ JJ(8),JJPP(8)
C
C     CONSTRUCT CURRENT REDUCTION MATRIX (AO TO SO)
C     AS DIRECT PRODUCT OF MATRICES IAOM * ICHAT
C
      II=0
      DO 40 I=1,NEQA
         DO 30 J=1,IAOR
            II=II+1
            KK=0
            DO 20 K=1,NEQA
               DO 10 L=1,IAOC
                  KK=KK+1
                  ISOM(KK,II)=IAOM(J,L)*ICHAT(K,I)
 10            CONTINUE
 20         CONTINUE
 30      CONTINUE
 40   CONTINUE
C
C     DETERMINE CHARACTER OF SO'S
C
      IF (NGENG.GT.0) GO TO 60
C
C     TRIVIAL CASE OF NO SYMMETRY
C
      DO 50 I=1,ISOR
         ICHSO(I)=1
 50   CONTINUE
      RETURN
C
   60 CONTINUE
C
C     DETERMINE FIRST CHARACTERS OF AO W.R. TO GENERATORS
C
      DO 80 I=1,NGENG
         DO 70 J=1,IAOR
            ISCR(J,I)= ICHAO(1,J)**IBING(1,I)
     1         *ICHAO(2,J)**IBING(2,I)
     2         *ICHAO(3,J)**IBING(3,I)
 70      CONTINUE
 80   CONTINUE
C
C     NOW CHARACTERS OF SO'S
C
      IFLW=NAT-NEQA
      DO 140 L=1,ISOR
         DO 110 I=1,NGENG
            II=0
            ISUM=0
            DO 100 J=1,NEQA
               KK=(IMAGE(IFLW+J,I)-1)*IAOR
               DO 90 K=1,IAOR
                  II=II+1
                  ISUM=ISUM+ISOM(L,II)*ISOM(L,KK+K)*ISCR(K,I)
 90            CONTINUE
 100        CONTINUE
            ISCS(I)=1
            IF (ISUM.LT.0) ISCS(I)=-1
 110     CONTINUE
C
C     DETERMINE NUMBER OR IRREP BY COMPARISON WITH CHARACT. OF GENER.
C     OF IRREP AS GIVEN IN ICHGEN
C
         M=2**NGENG
         DO 130 I=1,M
            ISUM=0
            DO 120 J=1,NGENG
               JCHGEN=ICHGEN(I,1)**IBING(1,J)*
     1            ICHGEN(I,2)**IBING(2,J)*
     2            ICHGEN(I,3)**IBING(3,J)
               ISUM=ISUM+IABS(ISCS(J)-JCHGEN)
 120        CONTINUE
            IF (ISUM.EQ.0) ICHSO(L)=JJ(I)
 130     CONTINUE
 140  CONTINUE
C
      RETURN
      END