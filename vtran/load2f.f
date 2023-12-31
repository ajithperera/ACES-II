      SUBROUTINE LOAD2F(CMO,W,W2,BUF,IBUF,ISYMAO,NBAS,
     &                  NMO,NSTART,NEND,ISIZE,ISIZT,IOFFAO,
     &                  IOFFI,NSIZE,ILNBUF,ISPIN,LUINT,LAST)
C
C THIS ROUTINE LOADS AO INTEGRALS FROM
C THE INTEGRAL FILE IIJJ AND TRANSFORMS THE FIRST
C INDEX.
C
C FULL TRANSFORMATION 
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER AND,OR
      LOGICAL LAST
C
      DIMENSION IBUF(ILNBUF),BUF(ILNBUF),W(NSIZE),W2(1),
     &          CMO(1000)
      DIMENSION NBAS(8),NMO(8),ISIZE(8),IOFFAO(8),IOFFI(8)
      DIMENSION NSTART(8),NEND(8)
      DIMENSION ISIZT(8,8)
      DIMENSION ISYMAO(100) 
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
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
      NPASS2=NPASS2+1
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
       IRREPX2=ISYMAO(KX)
C
C ENSURE THAT IRREPX1 IS LOWER THAN IRREPX2
C
       IF(IRREPX2.LT.IRREPX1) THEN
C
        IRREPX1T=IRREPX1
        IRREPX1=IRREPX2
        IRREPX2=IRREPX1T
        IXT=IX
        IX=KX
        KX=IXT
        JXT=JX
        JX=LX
        LX=JXT
C
       ENDIF
C
C  OFFSET WITHIN BASIS FUNCTIONS
C
       IOFF1=IOFFAO(IRREPX1)
       IOFF2=IOFFAO(IRREPX2)
C
C  INDICES WITHIN IRREP
C
       IX=IX-IOFF1
       JX=JX-IOFF1
       KX=KX-IOFF2
       LX=LX-IOFF2
C
       K=MAX(KX,LX)
       L=MIN(KX,LX)
C
C  OFFSET FOR CMO WITHIN IRREPX1 
C
       IOFFC1=INDOCC(IRREPX1,ISPIN)-1
C
C  NUMBER OF  ORBITALS WITHIN IRREPX1 
C
       NMOX1=NMO(IRREPX1)
       NSTARTX1=NSTART(IRREPX1)
       NENDX1=NEND(IRREPX1)
C
C  NUMBER OF BASIS FUNCTIONS WITHIN IRREPX1
C
       NBASX1=NBAS(IRREPX1)
C
C  OFFSET WITHIN W
C
       IOFFW1=IOFFI(IRREPX1)
C
       ISIZEX1=ISIZE(IRREPX1)
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
C WE NEED ONLY TWO !
C
C
C FIRST ELEMENT:     TRANSFORM INDEX J TO MO BASIS
C
        IND=INDXT(K,L)+ISIZT(IRREPX1,IRREPX2)-1
        IF(IX.EQ.JX) X=X*HALF
        IOFFC1J=IOFFC1+JX
        DO 9005 IMO=NSTARTX1,NENDX1
         W2(IMO)=CMO(IOFFC1J+(IMO-1)*NBASX1)*X
9005    CONTINUE 
C
C INCREMENT L(IX,KX,LX;I)
C
        IADRR=IOFFW1+IX+NBASX1*IND
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9006 IMO=NSTARTX1,NENDX1
         IADR=IADRR+ISIZEX1*(IMO-NSTARTX1)
         W(IADR)=W(IADR)+W2(IMO)
9006    CONTINUE
C
C SECOND ELEMENT:     TRANSFORM INDEX I TO MO BASIS
C
        IOFFC1I=IOFFC1+IX
        DO 9007 IMO=NSTARTX1,NENDX1
         W2(IMO)=CMO(IOFFC1I+(IMO-1)*NBASX1)*X
9007    CONTINUE 
C
C INCREMENT L(JX,LX,KX;I)
C
        IADRR=IOFFW1+JX+NBASX1*IND
*VOCL LOOP,NOVREC
CDIR$ IVDEP
        DO 9008 IMO=NSTARTX1,NENDX1
         IADR=IADRR+ISIZEX1*(IMO-NSTARTX1)
         W(IADR)=W(IADR)+W2(IMO)
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
c      WRITE(6,1002)NAOINT
      NLOAD2=NAOINT
1002  FORMAT(T3,'@LOAD2F-I, ',I8,' integrals read in from ',
     &          'file IIJJ')
c      ind=0
c      do 234 j=1,noccx
c      do 234 i=1,nbasx
c      do 234 l=1,nbasx
c      do 234 k=1,nbasx
c       ind=ind+1
c       write(*,235) j,i,l,k,ind,w(ind)
c235    format('indices I,MU,RHO,SIGMA',4I3,I5,F10.7)
c234   continue
      
      RETURN
900   WRITE(6,1000)
1000  FORMAT(T3,'@LOAD-F, Unexpected end-of-file on integral file.')
      CALL ERREX
      RETURN
901   WRITE(6,1000)
1001  FORMAT(T3,'@LOAD-F, I/O error on unit ',I5,'.')
      CALL ERREX
      RETURN
      END
