      SUBROUTINE DUMP9(AODINT,BUF,BUF1,BUF2,NSTART,NEND,IRREPX)
C
C THIS ROUTINE DUMPS TRANSFORMED INTEGRAL DERIVATIVES
C TO DISK. HERE, INTEGRALS FOR NON-SYMMETRIC PERTURBATIONS
C OF TYPE JKII ARE HANDLED
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,AOPOP,POP,VRT
C
      DIMENSION AODINT(1),BUF(1),BUF1(1),
     &          BUF2(1),NSTART(8),NEND(8)
      DIMENSION IOFFOO(8,8),IOFFVO(8,8),IOFFOV(8,8),
     &          IOFFVV(8,8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/AOPOPS/AOPOP(8)
      COMMON/INFOI/TIMEI,NWRITI
C
      INDXF(I,J,N)=I+N*(J-1)
C
C INITIALIZE OFFSETS
C
      IOFF=1
C
      DO 1 IRREP=1,NIRREP
       IOFF1=0
       IOFF2=0
       IOFF3=0
       IOFF4=0
       DO 2 IRREP1=1,NIRREP
        IRREP2=DIRPRD(IRREP1,IRREP)
        IOFFOO(IRREP2,IRREP1)=IOFF1
        IOFFVO(IRREP2,IRREP1)=IOFF2
        IOFFOV(IRREP2,IRREP1)=IOFF3
        IOFFVV(IRREP2,IRREP1)=IOFF4
        IOFF1=IOFF1+POP(IRREP1,1)*POP(IRREP2,1)
        IOFF2=IOFF2+VRT(IRREP2,1)*POP(IRREP1,1)
        IOFF3=IOFF3+VRT(IRREP1,1)*POP(IRREP2,1)
        IOFF4=IOFF4+VRT(IRREP1,1)*VRT(IRREP2,1)
2      CONTINUE
1     CONTINUE
C
C LOOP OVER ALL IRREPS
C
      DO 1000 IRREP1=1,NIRREP
C
       IRREP2=DIRPRD(IRREP1,IRREPX)
       NBAS1=AOPOP(IRREP1)
       NBAS2=AOPOP(IRREP2)
       NOCC1=POP(IRREP1,1)
       NSTART1=NSTART(IRREP1)
       NEND1=NEND(IRREP1)
       NVRT1=VRT(IRREP1,1) 
       NOCC2=POP(IRREP2,1)
       NVRT2=VRT(IRREP2,1) 
C
       LEN=0
       DO 1001 IRREP3=1,NIRREP
        IF(IRREP3.NE.IRREP1.AND.IRREP3.NE.IRREP2) THEN
         LEN=LEN+AOPOP(IRREP3)*AOPOP(IRREP3)
        ENDIF
1001   CONTINUE
C
       DO 100 IND=NSTART1,NEND1
C
        IOFFA=IOFF
        DO 2000 IRREP3=1,NIRREP
C
         IF(IRREP3.EQ.IRREP1.OR.IRREP3.EQ.IRREP2) GO TO 2000
         IRREP13=DIRPRD(IRREP1,IRREP3)
         IRREP23=DIRPRD(IRREP2,IRREP3)
         NBAS3=AOPOP(IRREP3)
         NOCC3=POP(IRREP3,1)
         NVRT3=VRT(IRREP3,1)
         NBAS33=NBAS3*NBAS3
         NBAS13=NBAS1*NBAS3
         NBAS23=NBAS2*NBAS3
         NBAS233=NBAS2*NBAS33
C
C GET INTEGRALS FOR GIVEN IND AS IRREP2,IRREP1,IRREP2
C
         IOFF2=IOFFA
         IOFF1=1
         DO 2001 IND2=1,NBAS2
          CALL SCOPY(NBAS33,AODINT(IOFF2),1,BUF(IOFF1),1)
          IOFF1=IOFF1+NBAS33
          IOFF2=IOFF2+LEN
2001     CONTINUE
C
C REORDER THEM AS IRREP1,IRREP3,IRREP2,IRREP3 (slowest to fastest)
C
        DO 3010 IP=1,NBAS3
         CALL SCOPY(NBAS23,BUF(IP),NBAS3,BUF2,1)
         CALL TRANSP(BUF2,BUF2(1+NBAS23),NBAS2,NBAS3)
         CALL SCOPY(NBAS23,BUF2(1+NBAS23),1,BUF(IP),NBAS3) 
3010    CONTINUE
C
C THE ORDERING IS NOW :  (FASTEST TO SLOWEST) IRREP3,IRREP2; IRREP3,
C                                             FIXED OCC(IRREP1)
C
        LENGTH=NOCC3*NOCC2
        NWRITI=NWRITI+LENGTH*NOCC3
        IF(LENGTH.NE.0) THEN
         DO 15 INDQ=1,NOCC3
          IOFFIJKL=INDXF(1,INDQ,NBAS23)
          IDISIJKL=IOFFOO(IRREP3,IRREP1)+INDXF(INDQ,IND,NOCC3)
          LOFFIJKL=IOFFOO(IRREP3,IRREP2)
          CALL BLKCPY2(BUF(IOFFIJKL),NBAS3,NBAS2,BUF1,NOCC3,NOCC2,1,1)
          CALL ULIST(BUF1,BUF2,LENGTH,IDISIJKL,LOFFIJKL,
     &               1,IRREP23,IRREP13,313,1)
15       CONTINUE
        ENDIF
C
        LENGTH=NOCC3*NOCC2
        NWRITI=NWRITI+LENGTH*NVRT3
        IF(LENGTH.NE.0) THEN
         DO 35 IQ=1,NVRT3
          INDQ=IQ+NOCC3
          IOFFIJKA=INDXF(1,INDQ,NBAS23)
          IDISIJKA=IOFFOV(IRREP1,IRREP3)+INDXF(IND,IQ,NOCC1)
          LOFFIJKA=IOFFOO(IRREP2,IRREP3)
          CALL BLKCPY2(BUF(IOFFIJKA),NBAS3,NBAS2,BUF1,NOCC3,NOCC2,1,1)
          CALL TRANSP(BUF1,BUF1(1+LENGTH),NOCC2,NOCC3)
          CALL ULIST(BUF1(1+LENGTH),BUF2,LENGTH,IDISIJKA,LOFFIJKA,
     &               1,IRREP23,IRREP13,310,1)
35       CONTINUE
        ENDIF
C
        LENGTH=NVRT3*NVRT2
        NWRITI=NWRITI+LENGTH*NOCC3
        IF(LENGTH.NE.0) THEN
         DO 45 INDQ=1,NOCC3
          IOFFABIJ=INDXF(1,INDQ,NBAS23)
          IDISABIJ=IOFFOO(IRREP3,IRREP1)+INDXF(INDQ,IND,NOCC3)
          LOFFABIJ=IOFFVV(IRREP3,IRREP2)
          CALL BLKCPY2(BUF(IOFFABIJ),NBAS3,NBAS2,BUF1,NVRT3,NVRT2,
     &                 NOCC3+1,NOCC2+1)
          CALL ULIST(BUF1,BUF2,LENGTH,IDISABIJ,LOFFABIJ,
     &               1,IRREP23,IRREP13,316,1)
45       CONTINUE
        ENDIF
C
        LENGTH=NVRT3*NVRT2
        NWRITI=NWRITI+LENGTH*NVRT3
        IF(LENGTH.NE.0) THEN
         DO 65 IQ=1,NVRT3
          INDQ=IQ+NOCC3
          IOFFABCI=INDXF(1,INDQ,NBAS23)
          IDISABCI=IOFFVO(IRREP3,IRREP1)+INDXF(IQ,IND,NVRT3)
          LOFFABCI=IOFFVV(IRREP3,IRREP2)
          CALL BLKCPY2(BUF(IOFFABCI),NBAS3,NBAS2,BUF1,NVRT3,NVRT2,
     &                 NOCC3+1,NOCC2+1)
          CALL ULIST(BUF1,BUF2,LENGTH,IDISABCI,LOFFABCI,
     &               1,IRREP23,IRREP13,330,1)
65       CONTINUE
        ENDIF
C
        LENGTH=NOCC2*NVRT3
        NWRITI=NWRITI+LENGTH*NVRT3
        IF(LENGTH.NE.0) THEN
         DO 75 IQ=1,NVRT3
          INDQ=IQ+NOCC3
          IOFFAIBJ=INDXF(1,INDQ,NBAS23)
          IDISAIBJ=IOFFVO(IRREP3,IRREP1)+INDXF(IQ,IND,NVRT3)
          LOFFAIBJ=IOFFVO(IRREP3,IRREP2)
          CALL BLKCPY2(BUF(IOFFAIBJ),NBAS3,NBAS2,BUF1,NVRT3,NOCC2,
     &                 NOCC3+1,1)
          CALL ULIST(BUF1,BUF2,LENGTH,IDISAIBJ,LOFFAIBJ,
     &               1,IRREP23,IRREP13,325,1)
75       CONTINUE
        ENDIF
C
        LENGTH=NOCC3*NVRT2
        NWRITI=NWRITI+LENGTH*NVRT3
        IF(LENGTH.NE.0) THEN
         DO 76 IQ=1,NVRT3
          INDQ=IQ+NOCC3
          IOFFIABJ=INDXF(1,INDQ,NBAS23)
          IDISIABJ=IOFFVO(IRREP3,IRREP1)+INDXF(IQ,IND,NVRT3)
          LOFFIABJ=IOFFVO(IRREP2,IRREP3)
          CALL BLKCPY2(BUF(IOFFIABJ),NBAS3,NBAS2,BUF1,NOCC3,NVRT2,
     &                 1,NOCC2+1)
          CALL TRANSP(BUF1,BUF1(1+LENGTH),NVRT2,NOCC3)
          CALL ULIST(BUF1(1+LENGTH),BUF2,LENGTH,IDISIABJ,LOFFIABJ,
     &               1,IRREP23,IRREP13,321,1)
76       CONTINUE
        ENDIF
C
        IOFFA=IOFFA+NBAS33
C
2000    CONTINUE
C
        IOFF=IOFF+LEN*NBAS2
C
100    CONTINUE
C
1000  CONTINUE
C
      RETURN
      END
