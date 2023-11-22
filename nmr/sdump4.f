      SUBROUTINE SDUMP4(AODINT,BUF,BUF1,NSTART,NEND,IRREPX)
C
C THIS ROUTINE DUMPS TRANSFORMED INTEGRAL DERIVATIVES
C TO DISK. HERE, INTEGRALS FOR SYMMETRIC PERTURBATIONS
C OF TYPE IJKL ARE HANDLED
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,AOPOP,POP,VRT
C
      DIMENSION AODINT(1),BUF(1),BUF1(1),NSTART(8),NEND(8)
      DIMENSION IOFFVV(8,8),IOFFVO(8,8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/AOPOPS/AOPOP(8)
C
      INDXF(I,J,N)=I+N*(J-1)
C
C INITIALIZE OFFSETS
C
      IOFF=1
C
C DETERMINE OFFSETS
C
      DO 1 IRREP=1,NIRREP
C
       IVV=0
       IVO=0
C
       DO 1 IRREPR=1,NIRREP 
C
        IRREPL=DIRPRD(IRREP,IRREPR)
        IOFFVV(IRREPL,IRREPR)=IVV
        IOFFVO(IRREPL,IRREPR)=IVO
        IVV=IVV+VRT(IRREPR,1)*VRT(IRREPL,1)
        IVO=IVO+VRT(IRREPL,1)*POP(IRREPR,1)
C
1     CONTINUE
C
C LOOP OVER ALL IRREPS
C
      DO 1000 IRREP1=1,NIRREP
C
       IRREP1R=DIRPRD(IRREP1,IRREPX)
       NSTARTR=NSTART(IRREP1R) 
       NENDR=NEND(IRREP1R) 
C
       DO 2000 IRREP2=1,NIRREP
C
        IF(IRREP1.EQ.IRREP2) GO TO 2000
C
        IRREP12=DIRPRD(IRREP1,IRREP2)
C
        LEN=0
        DO 2100 IRREP3=1,NIRREP
         IF(IRREP3.NE.IRREP1.AND.IRREP3.NE.IRREP2) THEN
          IRREP4=DIRPRD(IRREP12,IRREP3)
          IF(IRREP4.LT.IRREP3) THEN
           LEN=LEN+VRT(IRREP3,1)*VRT(IRREP4,1)
          ENDIF
         ENDIF
2100    CONTINUE 
         
        NVRT2=VRT(IRREP2,1)
C
        DO 100 IND=NSTARTR,NENDR
         IOFFA=IOFF
         DO 3000 IRREP3=1,NIRREP
          IF(IRREP3.NE.IRREP1.AND.IRREP3.NE.IRREP2) THEN
           IRREP4=DIRPRD(IRREP12,IRREP3)
           IF(IRREP4.LT.IRREP3) THEN
            IRREP13=DIRPRD(IRREP1,IRREP3)
            IRREP13R=DIRPRD(IRREP1R,IRREP3)
            IRREP14R=DIRPRD(IRREP1R,IRREP4)
            IRREP24=DIRPRD(IRREP2,IRREP4)
            IRREP23=DIRPRD(IRREP2,IRREP3)
            NVRT3=VRT(IRREP3,1)
            NVRT4=VRT(IRREP4,1)
C
            NVRT23=NVRT2*NVRT3
            NVRT24=NVRT2*NVRT4
            NVRT34=NVRT3*NVRT4
            NVRT234=NVRT2*NVRT3*NVRT4
C
C 
C GET INTEGRALS FOR GIVEN IND AS IRREP2,IRREP3,IRREP4
C
            IOFF1=1 
            IOFF2=IOFFA
            DO 200 IND2=1,NVRT2
             CALL SCOPY(NVRT34,AODINT(IOFF2),1,BUF(IOFF1),1)
             IOFF1=IOFF1+NVRT34
             IOFF2=IOFF2+LEN
200         CONTINUE
C
C REORDER THEM AS IRREP1,IRREP3,IRREP2,IRREP4
C
            DO 3010 IP=1,NVRT4
             CALL SCOPY(NVRT23,BUF(IP),NVRT4,BUF1,1)
             CALL TRANSP(BUF1,BUF1(1+NVRT23),NVRT2,NVRT3)
             CALL SCOPY(NVRT23,BUF1(1+NVRT23),1,BUF(IP),NVRT4) 
3010        CONTINUE
C
C THE ORDERING IS NOW :  (FASTEST TO SLOWEST) IRREP4,IRREP2; IRREP3,
C                                             FIXED OCC(IRREP1)
            LENGTH=NVRT2*NVRT4
            IF(LENGTH.NE.0) THEN
             DO 65 IQ=1,NVRT3
              IOFFABCI=INDXF(1,IQ,NVRT24)
              IDISABCI=IOFFVO(IRREP3,IRREP1R)+INDXF(IQ,IND,NVRT3)
              LOFFABCI=IOFFVV(IRREP4,IRREP2)
              CALL ULIST(BUF(IOFFABCI),BUF1,LENGTH,IDISABCI,LOFFABCI,
     &                   1,IRREP24,IRREP13R,330,2)
65           CONTINUE
            ENDIF
C
C REORDER THEM AS IRREP1,IRREP4,IRREP2,IRREP3
C
            IOFF1=1 
            IOFF2=IOFFA
            DO 201 IND2=1,NVRT2
             CALL SCOPY(NVRT34,AODINT(IOFF2),1,BUF(IOFF1),1)
             IOFF1=IOFF1+NVRT34
             IOFF2=IOFF2+LEN
201         CONTINUE
C
            IOFF1=1
            DO 3011 IP=1,NVRT2
             CALL TRANSP(BUF(IOFF1),BUF1,NVRT3,NVRT4)
             CALL SCOPY(NVRT34,BUF1,1,BUF(IOFF1),1)
             IOFF1=IOFF1+NVRT34
3011        CONTINUE
C
            DO 3012 IP=1,NVRT3
             CALL SCOPY(NVRT24,BUF(IP),NVRT3,BUF1,1)
             CALL TRANSP(BUF1,BUF1(1+NVRT24),NVRT2,NVRT4)
             CALL SCOPY(NVRT24,BUF1(1+NVRT24),1,BUF(IP),NVRT3) 
3012        CONTINUE
C
C THE ORDERING IS NOW :  (FASTEST TO SLOWEST) IRREP3,IRREP2; IRREP4,
C                                             FIXED OCC(IRREP1)
            LENGTH=NVRT2*NVRT3
            IF(LENGTH.NE.0) THEN
             DO 165 IQ=1,NVRT4
              IOFFABCI=INDXF(1,IQ,NVRT23)
              IDISABCI=IOFFVO(IRREP4,IRREP1R)+INDXF(IQ,IND,NVRT4)
              LOFFABCI=IOFFVV(IRREP3,IRREP2)
              CALL ULIST(BUF(IOFFABCI),BUF1,LENGTH,IDISABCI,LOFFABCI,
     &                   1,IRREP23,IRREP14R,330,2)
165          CONTINUE
            ENDIF
C
            IOFFA=IOFFA+NVRT34
C
           ENDIF
          ENDIF
3000     CONTINUE
C
         IOFF=IOFF+LEN*NVRT2
C
100     CONTINUE
2000   CONTINUE
1000  CONTINUE
C
      RETURN
      END