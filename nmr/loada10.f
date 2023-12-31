      SUBROUTINE LOADA10(CMO,W,W2,BUF,IBUF,ISYMAO,NBAS,
     &                   NSTART,NEND,IOFFAO,IOFFW,
     &                   NSIZE,ILNBUF,IRREPX,NPERT,ISPIN)
C
C THIS ROUTINE LOADS GIAO AO INTEGRALS FROM
C THE FILE IJIK AND IJKI AND TRANSFORMS THE FIRST
C INDEX TO THE MO BASIS. 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W,W2,BUF,X,X1,X2,CMO,HALF
      DOUBLE PRECISION TCPU,TSYS,TIME, TIMEIN,
     &                 TIMENOW, TIMETOT, TIMENEW 
      CHARACTER*8 NAMXYZ
      CHARACTER*80 FNAME
C
      DIMENSION NSTART(8),NEND(8),NBAS(8),ISYMAO(1),IOFFAO(8),IOFFW(8,8)
      DIMENSION IBUF(ILNBUF),BUF(ILNBUF),W(1),W2(1),CMO(1)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NDUMMY,NIRREP,IRREPA(255,2),DIRPRD(8,8)
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
C CHAIN IN THE INTEGRALS
C
      CALL ZERO(W,NSIZE)
1     READ(30)BUF,IBUF,NUT
      IF(NUT.EQ.-1) GO TO 11 
      DO 10 INT=1,NUT
       X=BUF(INT)
       ITMP=IBUF(INT)
C
C  X VALUE OF GIAO INTEGRAL 
C
C  I,K INDICES OF COMPLEX CONJUGATE ORBITALS
C
       IXX=IUPKI(ITMP)
       JXX=IUPKJ(ITMP)
       KXX=IUPKK(ITMP)
       LXX=IUPKL(ITMP)
C
C DETERMINE IRREP OF ORBITALS
C
       IRREPX1=ISYMAO(IXX)
       IRREPX2=ISYMAO(JXX)
       IRREPX3=ISYMAO(KXX)
       IRREPX4=ISYMAO(LXX)
C
       IF(IRREPX1.EQ.IRREPX3) THEN
        IRREPXI=IRREPX1
        IRREPXJ=IRREPX2
        IRREPXK=IRREPX4 
        IX=IXX
        JX=JXX
        KX=KXX
        LX=LXX
        ITYPE=1
       ELSE IF(IRREPX1.EQ.IRREPX4) THEN
        IRREPXI=IRREPX1
        IRREPXJ=IRREPX2
        IRREPXK=IRREPX3
        IX=IXX
        JX=JXX
        KX=KXX
        LX=LXX
        ITYPE=2
       ELSE IF(IRREPX2.EQ.IRREPX3) THEN
        IRREPXI=IRREPX2
        IRREPXJ=IRREPX1
        IRREPXK=IRREPX4
        IX=JXX
        JX=IXX
        KX=LXX
        LX=KXX
        ITYPE=2
        X=-X
       ELSE IF(IRREPX2.EQ.IRREPX4) THEN
        IRREPXI=IRREPX2
        IRREPXJ=IRREPX1
        IRREPXK=IRREPX3
        IX=JXX
        JX=IXX
        KX=LXX
        LX=KXX
        ITYPE=1
        X=-X
       ENDIF
C
C OFFSET WITHIN BASIS FUNCTIONS
C
       IOFFI=IOFFAO(IRREPXI)
       IOFFJ=IOFFAO(IRREPXJ)
       IOFFK=IOFFAO(IRREPXK)
C
C OFFSET FOR CMO
C
       IOFFC=INDOCC(IRREPXI,ISPIN)-1
C
C OFFSET WITHIN W
C
       IOFFWJ=IOFFW(IRREPXI,IRREPXJ)
       IOFFWK=IOFFW(IRREPXI,IRREPXK)
C
C NUMBER OF OCCUPIED ORBITALS AND NUMBER OF BASIS FUNCTIONS WITHIN IRREPX
C
       NSTARTI=NSTART(IRREPXI)
       NENDI=NEND(IRREPXI)
       NBASXI=NBAS(IRREPXI)
       NBASXJ=NBAS(IRREPXJ)
       NBASXK=NBAS(IRREPXK)
C
C SIZE OF THREE DIMENSIONAL ARRAY WITHIN IRREPX
C 
       ISIZEX=2*NBASXI*NBASXJ*NBASXK
C
       IF(ITYPE.EQ.1) THEN
C
C INDICES WITHIN IRREP
C
        IX=IX-IOFFI
        JX=JX-IOFFJ
        KX=KX-IOFFI
        LX=LX-IOFFK
C
C THIRD ELEMENT:     TRANSFORM INDEX K TO MO BASIS
C
        IND1=INDXF(JX,IX,NBASXJ)
        X1=-X
        DO 9005 I=NSTARTI,NENDI
         W2(I)=CMO(IOFFC+KX+(I-1)*NBASXI)*X1
9005    CONTINUE 
C
C INCREMENT L(LX,JX,IX;I)
C
CDIR$ IVDEP
        DO 9006 I=NSTARTI,NENDI
         IADR=IOFFWK+LX+NBASXK*(IND1-1)+ISIZEX*(I-NSTARTI)
         W(IADR)=W(IADR)+W2(I)
9006    CONTINUE
C
C FOURTH ELEMENT:     TRANSFORM INDEX I TO MO BASIS
C
        IND1=INDXF(LX,KX,NBASXK)
        DO 9007 I=NSTARTI,NENDI
         W2(I)=CMO(IOFFC+IX+(I-1)*NBASXI)*X1
9007    CONTINUE 
C
C INCREMENT L(JX,LX,KX;I)
C
CDIR$ IVDEP
        DO 9008 I=NSTARTI,NENDI
         IADR=IOFFWJ+JX+NBASXJ*(IND1-1)+ISIZEX*(I-NSTARTI)
         W(IADR)=W(IADR)+W2(I)
9008    CONTINUE
C
       ELSE IF(ITYPE.EQ.2) THEN
C
C INDICES WITHIN IRREP
C
        IX=IX-IOFFI
        JX=JX-IOFFJ
        KX=KX-IOFFK
        LX=LX-IOFFI
C
        IND1=INDXF(IX,JX,NBASXI)+NBASXI*NBASXJ
        X1=X
        DO 9015 I=NSTARTI,NENDI
         W2(I)=CMO(IOFFC+LX+(I-1)*NBASXI)*X1
9015    CONTINUE
C
C INCREMENT L(    )
C
CDIR$ IVDEP
        DO 9016 I=NSTARTI,NENDI
         IADR=IOFFWK+KX+NBASXK*(IND1-1)+ISIZEX*(I-NSTARTI)
         W(IADR)=W(IADR)+W2(I)
9016    CONTINUE
C
        IND1=INDXF(LX,KX,NBASXI)+NBASXI*NBASXK
        X1=-X
        DO 9017 I=NSTARTI,NENDI
         W2(I)=CMO(IOFFC+IX+(I-1)*NBASXI)*X1
9017    CONTINUE
C
C INCREMENT L(   )
C
CDIR$ IVDEP
        DO 9018 I=NSTARTI,NENDI
         IADR=IOFFWJ+JX+NBASXJ*(IND1-1)+ISIZEX*(I-NSTARTI)
         W(IADR)=W(IADR)+W2(I)
9018    CONTINUE
C
       ENDIF 
10    CONTINUE
      NAOINT=NAOINT+NUT
      GO TO 1
11    CLOSE(UNIT=30,STATUS='KEEP')
C
      CALL TIMER(1)
      TIME=TIME+TIMENEW
      NREAD=NAOINT
      NPASS=NPASS+1
      RETURN
C
900   WRITE(6,1000)
1000  FORMAT(T3,'@LOADA10-F, Unexpected end-of-file on integral file.')
      CALL ERREX
      RETURN
901   WRITE(6,1000)
1001  FORMAT(T3,'@LOADA10-F, I/O error on unit ',I5,'.')
      CALL ERREX
      RETURN
      END
