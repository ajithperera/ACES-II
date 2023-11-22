

      SUBROUTINE STT16(W,W2,BUF,IBUF,NSIZA,NSIZB,NSIZ2,NORB,
     &                 ITYPE,ISPIN,IUHF)
C
C THIS ROUTINE DOES AN IN-CORE SORT AND LIST GENERATION FOR
C  <ij|kl> OR <ab|cd> INTEGRALS (ALPHA-ALPHA SPIN CASE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND
      CHARACTER*4 SPCASE(3)
      CHARACTER*4 INTYPE(6)
      CHARACTER*8 INAME
      CHARACTER*80 FNAME
      DIMENSION W(NSIZA,NSIZB),W2(NSIZ2*NSIZ2)
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
C
C IINTLN IS THE INTEGER LENGTH IN BYTES.
C IFLTLN IS THE FLOATING POINT LENGTH IN BYTES.
C IINTFP IS THE RATIO IFLTLN/IINTLN
C
      COMMON /BUFLEN/ ILNBUF
      DATA INTYPE /'HHHH','PHHH','PPHH','PHPH','PPPH','PPPP'/
      DATA SPCASE /'AA  ','BB  ','AB  '/
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      INDX(I,J)=I+(J*(J-1)/2)
      IF(IFLAGS(1).Gt.1)THEN
       WRITE(LUOUT,2000)INTYPE(ITYPE),SPCASE(ISPIN)
2000   FORMAT(T3,'@STT16-I, Processing integral type ',
     & A4,' spin case ',A2,'.')
      ENDIF
      IX=0
      CALL ZERO(W,NSIZA*NSIZB)
      CALL ZERO(W2,NSIZ2*NSIZ2)
      INAME=INTYPE(ITYPE)//SPCASE(ISPIN)
      CALL GFNAME(INAME,FNAME,ILENGTH)
      OPEN(UNIT=15,FILE=FNAME(1:ILENGTH),
     &FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD')
C
C PICK UP INTEGRALS ONE BY ONE AND PLACE THEM IN A TRIANGULAR ARRAY.
C  <ij|kl> (OR <ab|cd>) ARE STORED i<=k;j<=l (a<=c;b<=d) IN W.
C
1     READ(15)BUF,IBUF,NUT
      DO 100 INT=1,NUT
       I=IUPKI(IBUF(INT))
       J=IUPKJ(IBUF(INT))
       K=IUPKK(IBUF(INT))
       L=IUPKL(IBUF(INT)) 
       I1=INDX(I,K)
       I2=INDX(J,L)
       W(I1,I2)=BUF(INT)
100   CONTINUE
      IF(NUT.EQ.ILNBUF)GOTO 1
      IF(ISPIN.NE.3)THEN
       DO 101 I=1,NSIZA
        DO 102 J=1,I-1
         W(I,J)=W(J,I)
102     CONTINUE
101    CONTINUE
      ENDIF
C
C FORM ANTISYMMETRIZED INTEGRALS IF THIS IS AA OR BB SPIN CASE.  THESE
C   ARE STORED i<=j;k<=l (a<=b,c<=d).
C
      IF(ISPIN.LT.3)THEN
       DO 10 I=1,NORB
        DO 11 J=1,I-1
         DO 20 K=1,NORB
          DO 21 L=1,K-1
           IX=IX+1
           I1=INDX(MIN(I,K),MAX(I,K))
           I2=INDX(MIN(J,L),MAX(J,L))
           I3=INDX(MIN(J,K),MAX(J,K))
           I4=INDX(MIN(I,L),MAX(I,L))
           W2(IX)=W(I1,I2)-W(I3,I4)
21        CONTINUE
20       CONTINUE
11      CONTINUE
10     CONTINUE
C
C NOW WRITE THEM OUT
C
       ISTART=1
       CALL PUTLST(W2,1,NSIZ2,2,ISPIN,ITYPE)
C
C IF THIS IS RHF, THEN WRITE THE UNANTISYMMETRIZED INTEGRALS 
C  INTO MOIO(2,ITYPE).
C
       IF(IUHF.EQ.0)THEN
        CALL PUTLST(W,1,NSIZA,2,2,ITYPE)
       ENDIF
      ELSEIF(ISPIN.EQ.3)THEN
C
C A-B SPIN CASE.
C
       CALL PUTLST(W,1,NSIZB,2,ISPIN,ITYPE)
      ELSE
       WRITE(LUOUT,200)ISPIN
200    FORMAT(T3,'@STT1XX-F, Spin case ',I3,' unknown.')
       CALL ERREX
      ENDIF
      CLOSE(UNIT=15,STATUS='DELETE')
      RETURN
      END 