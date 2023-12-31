      SUBROUTINE LOADA3(CMO,W,W2,BUF,IBUF,ISYMAO,NBAS,
     &                  NSTART,NEND,ISIZE3,IOFFAO,IOFFI,
     &                  IOFF3,NSIZE,ILNBUF,NPERT,ISPIN)
C
C THIS ROUTINE LOADS AO INTEGRALS FROM
C THE INTEGRAL FILE IJIJ AND TRANSFORMS THE FIRST
C INDEX.
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,AND
      CHARACTER*8 NAMXYZ
      CHARACTER*80 FNAME
C
      DIMENSION IBUF(ILNBUF),BUF(ILNBUF),W(1),W2(1),CMO(1)
      DIMENSION NBAS(8),NSTART(8),NEND(8),IOFFAO(8),IOFFI(8)
      DIMENSION ISIZE3(8),IOFF3(8,8,2)
      DIMENSION ISYMAO(1) 
C
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/PERT4/IXYZ
      COMMON/BFILE/NAMXYZ(4,3)
      COMMON/INFOL/NREAD,NPASS,TIME
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      DATA HALF/0.5D0/
C
      IUPKI(IX)=AND(IX,IALONE)
      IUPKJ(IX)=AND(ISHFT(IX,-IBITWD),IALONE)
      IUPKK(IX)=AND(ISHFT(IX,-2*IBITWD),IALONE)
      IUPKL(IX)=AND(ISHFT(IX,-3*IBITWD),IALONE)
      INDXF(I,J,N)=I+N*(J-1)
C
      CALL TIMER(1)
C
C  OPEN FILE CONTAINING 2e DERIVATIVE INTEGRALS
C
      CALL GFNAME(NAMXYZ(3,IXYZ),FNAME,ILENGTH)
      OPEN(UNIT=30,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      REWIND(UNIT=30)
      NAOINT=0
C
C READ IN INTEGRALS
C
      CALL ZERO(W,NSIZE)
1     READ(30) BUF,IBUF,NUT
      IF(NUT.EQ.-1) GO TO 11
      DO 10 INT=1,NUT
       X=BUF(INT)
       ITMP=IBUF(INT)
C
C  X VALUE, IX,JX,KX,LX INDICES  
C
       IX=IUPKI(ITMP)
       JX=IUPKJ(ITMP)
       KX=IUPKK(ITMP)
       LX=IUPKL(ITMP)
C
C GET IRREP OF INTEGRAL
C
       IRREPX1=ISYMAO(IX)
       IRREPX2=ISYMAO(JX)
       IRREPX3=ISYMAO(KX)
       IRREPX4=ISYMAO(LX)
C
       IRREPX=DIRPRD(IRREPX1,IRREPX2)
C
C  OFFSET WITHIN BASIS FUNCTIONS
C
       IOFFX1=IOFFAO(IRREPX1)
       IOFFX2=IOFFAO(IRREPX2)
       IOFFX3=IOFFAO(IRREPX3)
       IOFFX4=IOFFAO(IRREPX4)
C
C  INDICES WITHIN IRREP
C
       I=IX-IOFFX1
       J=JX-IOFFX2
       K=KX-IOFFX3
       L=LX-IOFFX4
C
C  OFFSET FOR CMO WITHIN IRREPX1 AND IRREPX2
C
       IOFFC1=INDOCC(IRREPX1,ISPIN)
       IOFFC2=INDOCC(IRREPX2,ISPIN)
       IOFFC3=INDOCC(IRREPX3,ISPIN)
       IOFFC4=INDOCC(IRREPX4,ISPIN)
C
C  NUMBER OF OCCUPIED ORBITALS WITHIN IRREPX1 AND IRREPX2
C
       NSTART1=NSTART(IRREPX1)
       NEND1=NEND(IRREPX1)
       NSTART2=NSTART(IRREPX2)
       NEND2=NEND(IRREPX2)
       NSTART3=NSTART(IRREPX3)
       NEND3=NEND(IRREPX3)
       NSTART4=NSTART(IRREPX4)
       NEND4=NEND(IRREPX4)
C
C  NUMBER OF BASIS FUNCTIONS WITHIN IRREPX 
C
       NBASX1=NBAS(IRREPX1)
       NBASX2=NBAS(IRREPX2)
       NBASX3=NBAS(IRREPX3)
       NBASX4=NBAS(IRREPX4)
C
C  OFFSET WITHIN W
C
       IOFFW41=IOFFI(IRREPX1)
       IOFFW42=IOFFI(IRREPX2)
       IOFFW43=IOFFI(IRREPX3)
       IOFFW44=IOFFI(IRREPX4)
C
       ISIZE31=ISIZE3(IRREPX1)
       ISIZE32=ISIZE3(IRREPX2)
       ISIZE33=ISIZE3(IRREPX3)
       ISIZE34=ISIZE3(IRREPX4)
C
C DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
C
       IF(IRREPX1.EQ.IRREPX3) THEN
        ITYP=1
        IF(I.EQ.K.AND.J.EQ.L) X=X*HALF
       ELSE
        ITYP=2
        IF(I.EQ.L.AND.J.EQ.K) X=X*HALF
       ENDIF
C
C ADDITIONAL OFFSET WITHIN W
C
       IOFFW41=IOFFW41+IOFF3(IRREPX1,IRREPX2,ITYP)
       IOFFW42=IOFFW42+IOFF3(IRREPX2,IRREPX1,ITYP)
       IOFFW43=IOFFW43+IOFF3(IRREPX3,IRREPX4,ITYP)
       IOFFW44=IOFFW44+IOFF3(IRREPX4,IRREPX3,ITYP)
C
C THERE ARE A TOTAL OF EIGHT CONTRIBUTIONS
C
C   (IJ|KL) (JI|KL) (IJ|LK) (JI|LK)
C   (KL|IJ) (KL|JI) (LK|IJ) (LK|JI)
C
C HOWEVER, WE NEED ONLY FOUR SINCE WE STORE THE AS
C
C   NU, SIGMA >= RHO,  I   
C
C FIRST ELEMENT:     TRANSFORM INDEX L TO MO BASIS
C
        X1=X
        IND=INDXF(I,J,NBASX1)
        DO 9001 IOCC=NSTART4,NEND4
         W2(IOCC)=CMO(IOFFC4+L-1+(IOCC-1)*NBASX4)*X1 
9001    CONTINUE 
C
C INCREMENT L(KX,IX,JX;I)
C
        DO 9002 IOCC=NSTART4,NEND4
         IADR=IOFFW44+K+NBASX3*(IND-1)+ISIZE34*(IOCC-NSTART4)
         W(IADR)=W(IADR)+W2(IOCC)
9002    CONTINUE
C
C SECOND ELEMENT:     TRANSFORM INDEX K TO MO BASIS
C
        X1=-X
        IND=INDXF(J,I,NBASX2)
        DO 9003 IOCC=NSTART3,NEND3
         W2(IOCC)=CMO(IOFFC3+K-1+(IOCC-1)*NBASX3)*X1
9003    CONTINUE 
C
C INCREMENT L(LX,JX,IX;I)
C
        DO 9004 IOCC=NSTART3,NEND3
         IADR=IOFFW43+L+NBASX4*(IND-1)+ISIZE33*(IOCC-NSTART3)
         W(IADR)=W(IADR)+W2(IOCC)
9004    CONTINUE
C
C THIRD ELEMENT:     TRANSFORM INDEX J TO MO BASIS
C
        X2=X
        IND=INDXF(K,L,NBASX3)
        DO 9005 IOCC=NSTART2,NEND2
         W2(IOCC)=CMO(IOFFC2+J-1+(IOCC-1)*NBASX2)*X2
9005    CONTINUE 
C
C INCREMENT L(IX,KX,LX;I)
C
        DO 9006 IOCC=NSTART2,NEND2
         IADR=IOFFW42+I+NBASX1*(IND-1)+ISIZE32*(IOCC-NSTART2)
         W(IADR)=W(IADR)+W2(IOCC)
9006    CONTINUE
C
C FOURTH ELEMENT:     TRANSFORM INDEX I TO MO BASIS
C
        X2=-X
        IND=INDXF(L,K,NBASX4)
        DO 9007 IOCC=NSTART1,NEND1
         W2(IOCC)=CMO(IOFFC1+I-1+(IOCC-1)*NBASX1)*X2
9007    CONTINUE 
C
C INCREMENT L(JX,LX,KX;I)
C
        DO 9008 IOCC=NSTART1,NEND1
         IADR=IOFFW41+J+NBASX2*(IND-1)+ISIZE31*(IOCC-NSTART1)
         W(IADR)=W(IADR)+W2(IOCC)
9008    CONTINUE
10    CONTINUE
C
      NAOINT=NAOINT+NUT
      GO TO 1
C
11    CLOSE(UNIT=30,STATUS='KEEP')
C
      CALL TIMER(1)
      TIME=TIME+TIMENEW
      NREAD=NAOINT
      NPASS=NPASS+1
C
      RETURN
900   WRITE(6,1000)
1000  FORMAT(T3,'@LOADA3-F, Unexpected end-of-file on integral file.')
      CALL ERREX
      RETURN
901   WRITE(6,1000)
1001  FORMAT(T3,'@LOADA3-F, I/O error on unit ',I5,'.')
      CALL ERREX
      RETURN
      END
