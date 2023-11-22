      SUBROUTINE SABCI2(BUCK,BUF,IBUCK,IBUF,NINBCK,ICHAIN,
     &                  NBUCK,NDBCK,NBKINT,NREC,IRECL,
     &                  ILNBUF,ISYM,BUF2,IBUF2,IWHERE,ILOOKUP,
     &                  IMASTER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISIL,DISKL,A,B,C,AND,OR
      CHARACTER*80 FNAME
      PARAMETER (LUINT=10)
      PARAMETER (LUSRT=15)
      DIMENSION BUCK(NBKINT,NBUCK),IBUCK(NBKINT,NBUCK),ISYM(1)
      DIMENSION IBUF(ILNBUF),NINBCK(NBUCK),ICHAIN(NBUCK),BUF(ILNBUF)
      DIMENSION BUF2(2*ILNBUF),IBUF2(2*ILNBUF),IWHERE(2*ILNBUF)
      DIMENSION ILOOKUP(2*ILNBUF),IMASTER(NDBCK*NBUCK)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
C
C SUBROUTINE CREATES A SORT FILE FOR <AB|IC> INTEGRALS, AB SPIN
C  CASE, I ALPHA SPIN.
C
C   BUCK   - SCRATCH ARRAY TO HOLD BUCKET BUFFER VALUES.
C
C   IBUCK  - SCRATCH ARRAY TO HOLD BUCKET BUFFER INDICES.
C
C   BUF    - HOLDS A BUFFER OF INTEGRAL VALUES.
C
C   IBUF   - HOLDS A BUFFER OF INTEGRAL INDICES.
C
C   NINBCK- BUCKET POPULATION VECTOR -- KEEPS TRACK OF THE NUMBER OF ELE
C           IN EACH BUCKET AT A GIVEN TIME.
C
C   ICHAIN- SORT CHAIN VECTOR.  ICHAIN(J) ALWAYS HOLDS THE VALUE OF THE
C           RECORD NUMBER WHICH CORRESPONDS TO THE LAST RECORD HOLDING
C           MEMBERS OF DISTRIBUTION J.
C
C   NBUCK - THE NUMBER OF BUCKETS USED IN THE SORT.
C
C   NDBCK - THE NUMBER OF DISTRIBUTIONS PER BUCKET.
C
C   NBKINT- NUMBER OF INTEGRALS PER BUCKET.
C
C   NREC  - NUMBER OF RECORDS IN SORTFILE (RETURNED).
C
C   IRECL - RECORD LENGTH FOR SORTFILE.
C
C   ILNBUF- THE SIZE OF AN INTEGRAL BUFFER. (600 FOR VMOL).
C
C
C STATEMENT FUNCTIONS FOR BIT UNPACKING.
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
C
C STATEMENT FUNCTIONS FOR SQUARE TRIANGULAR INDEXING.
C
      INDX(I,J,N)=I+(J-1)*N
C
C OPEN UP SORTFILE AND INTEGRAL FILE.
C
      NOCA=NOCCO(1)
      NOCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      CALL GFNAME('PPPH1H  ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD',ACCESS='SEQUENTIAL')
      CALL GFNAME('SRTFIL  ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
C
C INITIALIZE SOME CRAP.
C
      NREC=1
      NUMINT=0
      CALL IZERO(NINBCK,NBUCK)
      CALL IZERO(ICHAIN,NBUCK)
      IOFF=0
      DO 5 I=1,NBUCK
       DO 6 J=1,NDBCK
        IOFF=IOFF+1
        IMASTER(IOFF)=I
6      CONTINUE
5     CONTINUE
C
C READ A BUFFER OF INTEGRALS AND DEAL WITH IT.
C
1     READ(LUINT)BUF,IBUF,NUT
      ISTICK=0
      CALL IZERO(IWHERE,2*ILNBUF)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 10 INT=1,NUT
       I=IUPKI(IBUF(INT))
       B=IUPKJ(IBUF(INT))
       A=IUPKK(IBUF(INT))
       C=IUPKL(IBUF(INT))
C
C ASSIGN DISTRIBUTION NUMBERS.  USE EQUIVALENCE OF <AB|IC> AND <AC|BI>,
C  BUT DO <BC|IA> FIRST.
C
       DISKL=ISYM(INDX(I,C,NOCA))
C
C PUT THIS INTEGRAL INTO A BUCKET BASED ON THE KL INDEX.  DO YOU
C   UNDERSTAND THE LOGIC OF THE ASSIGNMENT?
C
       IBKET=1+(DISKL-1)/NDBCK
       IBUF(INT)=IPACK(A,B,I,C)
       IWHERE(ISTICK+1)=IMASTER(DISKL)
       IBUF2(ISTICK+1) =IBUF(INT)
       BUF2(ISTICK+1)  =BUF(INT)
C
C IF THE INTEGRAL IS OF TYPE (AB|IB), THEN IT WILL ONLY BE PLUNKED
C  INTO ONE BUCKET.  ANY OTHER TYPE GOES INTO TWO.  NOTE THAT THE
C  <CB|AI> INTEGRAL IS USED HERE.
C
       IF(B.NE.C)THEN
        DISIL=ISYM(INDX(I,B,NOCA))
        IBUF(INT)=IPACK(A,C,I,B)
        IBKET=1+(DISIL-1)/NDBCK
        IWHERE(ISTICK+2)=IMASTER(DISIL)
        IBUF2(ISTICK+2) =IBUF(INT)
        BUF2(ISTICK+2)  =BUF(INT)
       ENDIF
       ISTICK=ISTICK+2
10    CONTINUE
C
C DEAL WITH STUFF NOW
C
      DO 102 ITYPE=1,NBUCK
       CALL WHENEQ(2*ILNBUF,IWHERE,1,ITYPE,ILOOKUP,NVAL)
       ioff=1
       inbuck=ninbck(itype)
11     idump=min(nval,nbkint-inbuck)
       call gather(idump,buck(inbuck+1,itype),buf2,ilookup(ioff))
       call igather(idump,ibuck(inbuck+1,itype),ibuf2,ilookup(ioff))
       inbuck=inbuck+idump 
       if(inbuck.eq.nbkint)then
        call plunk(lusrt,buck(1,itype),ibuck(1,itype),ichain(itype),
     &              inbuck,nbkint,nrec)
        nval=nval-idump
        ioff=ioff+idump
        if(nval.ne.0)goto 11
       endif
       ninbck(itype)=inbuck
102   CONTINUE
      NUMINT=NUMINT+NUT
      IF(NUT.EQ.ILNBUF)GOTO 1
      WRITE(6,1400)NUMINT
1400  FORMAT(T3,'@SABCI2-I, ',I8,' PPPH integrals sorted.')
C
C FLUSH REMAINING BUFFERS.
C
      DO 20 I=1,NBUCK
       CALL PLUNK(LUSRT,BUCK(1,I),IBUCK(1,I),ICHAIN(I),NINBCK(I),
     &            NBKINT,NREC)
20    CONTINUE
      NREC=NREC-1
      CLOSE(UNIT=LUSRT,STATUS='KEEP')
      CLOSE(UNIT=LUINT,STATUS='DELETE')
      RETURN
      END
