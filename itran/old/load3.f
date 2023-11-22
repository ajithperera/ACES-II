      SUBROUTINE LOAD3(CMO,W,W2,BUF,IBUF,ISYMAO,NBAS,
     &                 NOCC,NSTART,NEND,IOFFAO,IOFF4,ISIZE3,
     &                 NSIZE,ILNBUF,ISPIN,LUINT,
     &                 LAST,IOFFSET)
C
C THIS ROUTINE LOADS AO INTEGRALS FROM
C THE INTEGRAL FILE IJIJ AND TRANSFORMS THE FIRST
C INDEX.
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LAST
      INTEGER DIRPRD,AND,OR
C
      DIMENSION IBUF(ILNBUF),BUF(ILNBUF),W(1),W2(1),CMO(1)
      DIMENSION NBAS(8),NOCC(8),IOFFAO(8),IOFF4(8,8)
      DIMENSION ISIZE3(8,8),NSTART(8),NEND(8)
      DIMENSION ISYMAO(100) 
C
      COMMON/SYMINF/NSTRAT,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/AOOFST/INDOCC(8,2)
      COMMON/VTINFO/NPASS1,NPASS2,NPASS3,NPASS4,
     &              NLOAD1,NLOAD2,NLOAD3,NLOAD4,
     &              NWRIT1,NWRIT2,NWRIT3,NWRIT4,
     &              NWRIT1A,NWRIT2A,NWRIT3A,NWRIT4A,
     &              NWRIT1B,NWRIT2B,NWRIT3B,NWRIT4B
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
      CALL ZERO(W,NSIZE)
C
      NPASS3=NPASS3+1
C
C READ IN INTEGRALS
C
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
       IRREPX1=ISYMAO(IX)
       IRREPX2=ISYMAO(JX)
       IRREPX=DIRPRD(IRREPX1,IRREPX2)
C
       IRREPX3=ISYMAO(KX)
       IRREPX4=ISYMAO(LX)
       IF(IRREPX3.NE.IRREPX1) THEN 
        KXX=LX
        LX=KX
        KX=KXX
       ENDIF
       IF(IRREPX1.GT.IRREPX2) THEN 
        ITMP=IX
        IX=JX
        JX=ITMP
        KTMP=KX
        KX=LX
        LX=KTMP
        IRREPX1T=IRREPX1
        IRREPX1=IRREPX2
        IRREPX2=IRREPX1T
       ENDIF
C
C  OFFSET WITHIN BASIS FUNCTIONS
C
       IOFF1=IOFFAO(IRREPX1)
       IOFF2=IOFFAO(IRREPX2)
C
C  INDICES WITHIN IRREP
C
       I=IX-IOFF1
       J=JX-IOFF2
       K=KX-IOFF1
       L=LX-IOFF2
C
C  OFFSET FOR CMO WITHIN IRREPX1 AND IRREPX2
C
       IOFFC1=INDOCC(IRREPX1,ISPIN)
       IOFFC2=INDOCC(IRREPX2,ISPIN)
C
C  NUMBER OF OCCUPIED ORBITALS WITHIN IRREPX1 AND IRREPX2
C
       NOCCX1=NOCC(IRREPX1)
       NOCCX2=NOCC(IRREPX2)
       NSTARTX1=NSTART(IRREPX1)
       NENDX1=NEND(IRREPX1)
       NSTARTX2=NSTART(IRREPX2)
       NENDX2=NEND(IRREPX2)
C
C  NUMBER OF BASIS FUNCTIONS WITHIN IRREPX 
C
       NBASX1=NBAS(IRREPX1)
       NBASX2=NBAS(IRREPX2)
C
C  OFFSET WITHIN W
C
       IOFFW41=IOFF4(IRREPX1,IRREPX2)
       IOFFW42=IOFF4(IRREPX2,IRREPX1)
C
       ISIZE31=ISIZE3(IRREPX1,IRREPX2)
       ISIZE32=ISIZE3(IRREPX2,IRREPX1)
C
C DETERMINE REDUNDANCY FACTOR FOR PLUGGING IN INTEGRALS
C
       IF(I.EQ.K.AND.J.EQ.L) X=X*HALF
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
        IND=INDXF(I,J,NBASX1)
        DO 9001 IOCC=NSTARTX2,NENDX2
         W2(IOCC)=CMO(IOFFC2+L-1+(IOCC-1)*NBASX2)*X 
9001    CONTINUE 
C
C INCREMENT L(KX,IX,JX;I)
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9002 IOCC=NSTARTX2,NENDX2
         IADR=IOFFW42+K+NBASX1*(IND-1)+ISIZE32*(IOCC-NSTARTX2)
         W(IADR)=W(IADR)+W2(IOCC)
9002    CONTINUE
C
C SECOND ELEMENT:     TRANSFORM INDEX K TO MO BASIS
C
        IND=INDXF(I,J,NBASX1)
        DO 9003 IOCC=NSTARTX1,NENDX1
         W2(IOCC)=CMO(IOFFC1+K-1+(IOCC-1)*NBASX1)*X
9003    CONTINUE 
C
C INCREMENT L(LX,JX,IX;I)
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9004 IOCC=NSTARTX1,NENDX1
         IADR=IOFFW41+L+NBASX2*(IND-1)+ISIZE31*(IOCC-NSTARTX1)
         W(IADR)=W(IADR)+W2(IOCC)
9004    CONTINUE
C
C THIRD ELEMENT:     TRANSFORM INDEX J TO MO BASIS
C
        IND=INDXF(K,L,NBASX1)
        DO 9005 IOCC=NSTARTX2,NENDX2
         W2(IOCC)=CMO(IOFFC2+J-1+(IOCC-1)*NBASX2)*X
9005    CONTINUE 
C
C INCREMENT L(IX,KX,LX;I)
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9006 IOCC=NSTARTX2,NENDX2
         IADR=IOFFW42+I+NBASX1*(IND-1)+ISIZE32*(IOCC-NSTARTX2)
         W(IADR)=W(IADR)+W2(IOCC)
9006    CONTINUE
C
C FOURTH ELEMENT:     TRANSFORM INDEX I TO MO BASIS
C
        IND=INDXF(K,L,NBASX1)
        DO 9007 IOCC=NSTARTX1,NENDX1
         W2(IOCC)=CMO(IOFFC1+I-1+(IOCC-1)*NBASX1)*X
9007    CONTINUE 
C
C INCREMENT L(JX,LX,KX;I)
C
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9008 IOCC=NSTARTX1,NENDX1
         IADR=IOFFW41+J+NBASX2*(IND-1)+ISIZE31*(IOCC-NSTARTX1)
         W(IADR)=W(IADR)+W2(IOCC)
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
C
CJDW 2/1/95. Delete integral file if possible.
C
       IF(IFLAGS(18).GE.3.OR.IFLAGS(26).EQ.1.OR.IFLAGS(93).EQ.2.OR.
     &   (ISPIN.EQ.1.AND.IFLAGS(11).GT.0)   .OR.IFLAGS2(103).EQ.1)THEN
        CLOSE(UNIT=LUINT,STATUS='KEEP')
       ELSE
        CLOSE(UNIT=LUINT,STATUS='DELETE')
       ENDIF
C
      ELSE
       REWIND(LUINT)
       CALL LOCATE(LUINT,'TWOELSUP')
      ENDIF
      NLOAD3=NAOINT
C      
      RETURN
900   WRITE(6,1000)
1000  FORMAT(T3,'@LOAD3-F, Unexpected end-of-file on integral file.')
      CALL ERREX
      RETURN
901   WRITE(6,1000)
1001  FORMAT(T3,'@LOAD3-F, I/O error on unit ',I5,'.')
      CALL ERREX
      RETURN
      END
