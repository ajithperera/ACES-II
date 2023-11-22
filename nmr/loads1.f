      SUBROUTINE LOADS1(PCMO,W,W2,BUF,IBUF,ISYMAO,NBAS,
     &                  NSTART,NEND,ISIZE,IOFFAO,IOFFI,NSIZE,
     &                  ILNBUF,IPERT,ISPIN,IRREPX)
C
C THIS ROUTINE LOADS AO INTEGRALS FROM
C THE INTEGRAL FILE IIII AND TRANSFORMS THE FIRST
C INDEX.
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LAST
      INTEGER DIRPRD,AND,POP,VRT
      CHARACTER*80 FNAME
C
      DIMENSION IBUF(ILNBUF),BUF(ILNBUF),W(1),W2(1),PCMO(1)
      DIMENSION NBAS(8),NSTART(8),NEND(8),
     &          ISIZE(8),IOFFAO(8),IOFFI(8)
      DIMENSION ISYMAO(100),IOFFPC(8,2) 
C
      COMMON/SYM/POP(8,2),VRT(8,2),NDUM(6)
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFOT/NPASS1,NPASS2,NPASS3,NPASS4,
     &             NREAD1,NREAD2,NREAD3,NREAD4,
     &             TIME1L,TIME2L,TIME3L,TIME4L,
     &             TIME1S,TIME2S,TIME3S,TIME4S
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      DATA HALF/0.5D0/
C
      IUPKI(IX)=AND(IX,IALONE)
      IUPKJ(IX)=AND(ISHFT(IX,-IBITWD),IALONE)
      IUPKK(IX)=AND(ISHFT(IX,-2*IBITWD),IALONE)
      IUPKL(IX)=AND(ISHFT(IX,-3*IBITWD),IALONE)
      INDXT(I,J)=J+(I*(I-1))/2
      INDXF(I,J,N)=I+N*(J-1)
C
      CALL TIMER(1)
C
      IOFF=1
      DO 11 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       IOFFPC(IRREPR,1)=IOFF
       IOFF=IOFF+NBAS(IRREPL)*POP(IRREPR,1)
11    CONTINUE
C 
C  OPEN FILE CONTAINING 2e INTEGRALS
C
      LAST=.TRUE.
      LUINT=30
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL',STATUS='OLD')
      CALL LOCATE(LUINT,'TWOELSUP')
C
C READ IN INTEGRALS
C
      CALL ZERO(W,NSIZE)
      NAOBUF=0
      NAOINT=0
1     READ(LUINT) BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
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
       IRREPL=ISYMAO(IX)
       IRREPR=DIRPRD(IRREPL,IRREPX)
C
C  OFFSET WITHIN BASIS FUNCTIONS
C
       IOFF=IOFFAO(IRREPL)
C
C  INDICES WITHIN IRREP
C
       IX=IX-IOFF
       JX=JX-IOFF
       KX=KX-IOFF
       LX=LX-IOFF
C
       I=MAX(IX,JX)
       J=MIN(IX,JX)
       K=MAX(KX,LX)
       L=MIN(KX,LX)
C
       IF(I.EQ.K.AND.J.EQ.L) THEN
        X=X*HALF
       ENDIF
C
C  OFFSET FOR PCMO
C
       IOFFC=IOFFPC(IRREPR,ISPIN)
C
C  NUMBER OF OCCUPIED ORBITALS WITHIN IRREPX
C
C  NUMBER OF BASIS FUNCTIONS WITHIN IRREPX 
C
       NBASX=NBAS(IRREPL)
       NSTARTX=NSTART(IRREPR)
       NENDX=NEND(IRREPR)
C
C  OFFSET WITHIN W
C
       IOFFW=IOFFI(IRREPR)
C
       ISIZEX=ISIZE(IRREPR)
C
C DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
C
C THERE ARE A TOTAL OF EIGHT CONTRIBUTIONSa
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
        IND=INDXT(I,J)
        X1=X
        IF(K.EQ.L) X1=X1*HALF 
        DO 9001 IFIRST=NSTARTX,NENDX
         W2(IFIRST)=PCMO(IOFFC+LX-1+(IFIRST-1)*NBASX)*X1 
9001    CONTINUE 
C
C INCREMENT L(KX,IX,JX;I)
C
        DO 9002 IFIRST=NSTARTX,NENDX
         IADR=IOFFW+KX+NBASX*(IND-1)+ISIZEX*(IFIRST-NSTARTX)
         W(IADR)=W(IADR)+W2(IFIRST)
9002    CONTINUE
        
C
C SECOND ELEMENT:     TRANSFORM INDEX K TO MO BASIS
C
        IND=INDXT(I,J) 
        DO 9003 IFIRST=NSTARTX,NENDX
         W2(IFIRST)=PCMO(IOFFC+KX-1+(IFIRST-1)*NBASX)*X1
9003    CONTINUE 
C
C INCREMENT L(LX,JX,IX;I)
C
        DO 9004 IFIRST=NSTARTX,NENDX
         IADR=IOFFW+LX+NBASX*(IND-1)+ISIZEX*(IFIRST-NSTARTX)
         W(IADR)=W(IADR)+W2(IFIRST)
9004    CONTINUE
C
C THIRD ELEMENT:     TRANSFORM INDEX J TO MO BASIS
C
        IND=INDXT(K,L)
        X2=X
        IF(I.EQ.J) X2=X2*HALF
        DO 9005 IFIRST=NSTARTX,NENDX
         W2(IFIRST)=PCMO(IOFFC+JX-1+(IFIRST-1)*NBASX)*X2
9005    CONTINUE 
C
C INCREMENT L(IX,KX,LX;I)
C
        DO 9006 IFIRST=NSTARTX,NENDX
         IADR=IOFFW+IX+NBASX*(IND-1)+ISIZEX*(IFIRST-NSTARTX)
         W(IADR)=W(IADR)+W2(IFIRST)
9006    CONTINUE
C
C FOURTH ELEMENT:     TRANSFORM INDEX I TO MO BASIS
C
        IND=INDXT(K,L)
        DO 9007 IFIRST=NSTARTX,NENDX
         W2(IFIRST)=PCMO(IOFFC+IX-1+(IFIRST-1)*NBASX)*X2
9007    CONTINUE 
C
C INCREMENT L(JX,LX,KX;I)
C
        DO 9008 IFIRST=NSTARTX,NENDX
         IADR=IOFFW+JX+NBASX*(IND-1)+ISIZEX*(IFIRST-NSTARTX)
         W(IADR)=W(IADR)+W2(IFIRST)
9008    CONTINUE
10    CONTINUE
C
      
      IF(NUT.EQ.ILNBUF) THEN
       NAOINT=NAOINT+NUT
       GO TO 1
      ENDIF
      IF(NUT.EQ.-1) NUT=0
C
      NAOINT=NAOINT+NUT
      IF(LAST) THEN
       CLOSE(UNIT=LUINT,STATUS='KEEP')
      ELSE
       REWIND(LUINT)
       CALL LOCATE(LUINT,'TWOELSUP')
      ENDIF
C
      NREAD1=NAOINT
      NPASS1=NPASS1+1
C      
      CALL TIMER(1)
      TIME1L=TIME1L+TIMENEW
C
      RETURN
900   WRITE(6,1000)
1000  FORMAT(T3,'@LOADS1-F, Unexpected end-of-file on integral file.')
      CALL ERREX
      RETURN
901   WRITE(6,1000)
1001  FORMAT(T3,'@LOADS1-F, I/O error on unit ',I5,'.')
      CALL ERREX
      RETURN
      END
