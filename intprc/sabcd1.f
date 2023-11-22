

      SUBROUTINE SABCD1(BUCK,BUF,IBUCK,IBUF,NINBCK,ICHAIN,
     &                  NBUCK,NDBCK,NBKINT,NREC,IRECL,
     &                  ILNBUF,ISYM,ISPIN,BUF2,IBUF2,IWHERE,
     &                  ILOOKUP,IMASTER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*4 SPCASE(2)
      CHARACTER*8 INAME
      CHARACTER*80 FNAME
      INTEGER DISCD,DISAD,DISCB,DISAB,DISDA,DISBC,DISBA,DISDC
      INTEGER A,B,C,D,AND,OR
      PARAMETER (LUINT=10)
      PARAMETER (LUSRT=15)
      DIMENSION BUCK(NBKINT,NBUCK),IBUCK(NBKINT,NBUCK),ISYM(1)
      DIMENSION IBUF(ILNBUF),NINBCK(NBUCK),ICHAIN(NBUCK),BUF(ILNBUF)
      DIMENSION BUF2(4*ILNBUF),IBUF2(4*ILNBUF),IWHERE(4*ILNBUF)
      DIMENSION IMASTER(NDBCK*NBUCK),ILOOKUP(4*ILNBUF)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
C
      DATA SPCASE /'AA  ','BB  '/
C
C SUBROUTINE CREATES A SORT FILE FOR <AB|CD> INTEGRALS, AA SPIN CASE,
C   RHF REFERENCE.
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
c
C STATEMENT FUNCTIONS FOR BIT UNPACKING.
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
C
C STATEMENT FUNCTIONS FOR INDEXING IN A TRIANGULAR MATRIX WITH A<B
C
      INDX(I,J)=I+((J-1)*(J-2))/2
C
C OPEN UP SORTFILE AND INTEGRAL FILE.
C
      NOCC=NOCCO(ISPIN)
      NVRT=NVRTO(ISPIN)
      INAME='PPPP'//SPCASE(ISPIN)
      CALL GFNAME(INAME,FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD',ACCESS='SEQUENTIAL')
      CALL GFNAME('SRTFIL  ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
C
C INITIALIZE SOME CRAP AND SET UP MASTER INDEX.
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
      CALL IZERO(IWHERE,4*ILNBUF)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 10 INT=1,NUT
       A=IUPKI(IBUF(INT))
       B=IUPKJ(IBUF(INT))
       C=IUPKK(IBUF(INT))
       D=IUPKL(IBUF(INT))
C
C INTEGRAL WILL BE ONE OF THREE TYPES:
C   <AB|CD>, <AB|AD> OR <AB|CB>, <AB|AB> AND <AA|AA>.
C
C  FOR TYPE I: WE ALSO MUST WRITE OUT SEVEN OTHER PERMUTATIONS:
C                  <CB|AD>,<AD|CB>,<CD|AB>,<BC|DA>,<DA|BC>,<DC|BA> 
C                  AND <BA|DC>
C
C  FOR TYPE II: WE NEED TO WRITE OUT:
C                  <AB|AD> OR <BA|DA> AND <AD|AB> OR <DA|BA>.
C
C  FOR TYPE III: WE NEED TO WRITE OUT:
C                  <AB|AB> OR <BA|BA> WITH A<B
C
C  FOR TYPE IV: NO CONTRIBUTION
C
C ASSIGN THE DISTRIBUTION NUMBER AND DETERMINE WHICH CLASS IT BELONGS TO
C
       ITYPE=1
       IF(A.EQ.C.OR.B.EQ.D)ITYPE=2
       IF(A.EQ.C.AND.B.EQ.D)THEN
        ITYPE=3
       ENDIF
       ITYPE=1
C
C PUT THIS INTEGRAL INTO A BUCKET BASED ON ITS CD INDEX.
C
      IF(C.LT.D) THEN
       DISCD=ISYM(INDX(C,D))
       IBKET=IMASTER(DISCD)
       ITMP=IPACK(A,B,C,D)
       iwhere(istick+1)=ibket
       ibuf2(istick+1)=itmp
       buf2(istick+1)=buf(int)
C
      ELSE IF(D.LT.C) THEN
       DISDC=ISYM(INDX(D,C))
       ITMP=IPACK(B,A,D,C)
       IBKET=IMASTER(DISDC)
       iwhere(istick+1)=ibket
       ibuf2(istick+1)=itmp
       buf2(istick+1)=buf(int)
      ENDIF
C
C IF INTEGRAL IS OF TYPE 1 OR 2, STILL MORE TO GO.
C
       IF(ITYPE.LE.2)THEN
        IF(A.LT.B) THEN
         DISAB=ISYM(INDX(A,B))
         ITMP=IPACK(C,D,A,B)
         IBKET=IMASTER(DISAB)
         iwhere(istick+2)=ibket
         ibuf2(istick+2)=itmp
         buf2(istick+2)=buf(int)
        ELSE IF(B.LT.A) THEN
         DISBA=ISYM(INDX(B,A))
         ITMP=IPACK(D,C,B,A)
         IBKET=IMASTER(DISBA)
         iwhere(istick+2)=ibket
         ibuf2(istick+2)=itmp
         buf2(istick+2)=buf(int)
        ENDIF
C
C IF TYPE 1, STILL A FEW MORE.
C
        IF(ITYPE.EQ.1)THEN
         IF(A.LT.D) THEN
          DISAD=ISYM(INDX(A,D))
          ITMP=IPACK(C,B,A,D)
          IBKET=IMASTER(DISAD)
          iwhere(istick+3)=ibket
          ibuf2(istick+3)=itmp
          buf2(istick+3)=buf(int)
         ELSE IF(D.LT.A) THEN
          DISDA=ISYM(INDX(D,A))
          ITMP=IPACK(B,C,D,A)
          IBKET=IMASTER(DISDA)
          iwhere(istick+3)=ibket
          ibuf2(istick+3)=itmp
          buf2(istick+3)=buf(int)
         ENDIF
         IF(B.LT.C) THEN
          DISBC=ISYM(INDX(B,C))
          ITMP=IPACK(D,A,B,C)
          IBKET=IMASTER(DISBC)
          iwhere(istick+4)=ibket
          ibuf2(istick+4)=itmp
          buf2(istick+4)=buf(int)
         ELSE IF(C.LT.B) THEN
          DISCB=ISYM(INDX(C,B))
          ITMP=IPACK(A,D,C,B)
          IBKET=IMASTER(DISCB)
          iwhere(istick+4)=ibket
          ibuf2(istick+4)=itmp
          buf2(istick+4)=buf(int)
         ENDIF
        ENDIF
       ENDIF
       istick=istick+4
10    CONTINUE
      NUMINT=NUMINT+NUT
C
C DEAL WITH STUFF NOW
C
      DO 15 ITYPE=1,NBUCK
       CALL WHENEQ(4*ILNBUF,IWHERE,1,ITYPE,ILOOKUP,NVAL)
       ioff=1
       inbuck=ninbck(itype)
21     idump=min(nval,nbkint-inbuck)
       call gather(idump,buck(inbuck+1,itype),buf2,ilookup(ioff))
       call igather(idump,ibuck(inbuck+1,itype),ibuf2,ilookup(ioff))
       inbuck=inbuck+idump 
       if(inbuck.eq.nbkint)then
        call plunk(lusrt,buck(1,itype),ibuck(1,itype),ichain(itype),
     &              inbuck,nbkint,nrec)
        nval=nval-idump
        ioff=ioff+idump
        if(nval.ne.0)goto 21
       endif
       ninbck(itype)=inbuck
15    CONTINUE
      IF(NUT.EQ.ILNBUF)GOTO 1
      WRITE(6,1400)NUMINT
1400  FORMAT(T3,'@SABCD1-I, ',I8,' PPPP integrals sorted.')
C
C FLUSH REMAINING BUFFERS.
C
      DO 20 I=1,NBUCK
       CALL PLUNK(LUSRT,BUCK(1,I),IBUCK(1,I),ICHAIN(I),NINBCK(I),
     &            NBKINT,NREC)
20    CONTINUE
      NREC=NREC-1
      CLOSE(UNIT=LUINT,STATUS='DELETE')
      CLOSE(UNIT=LUSRT,STATUS='KEEP')
      RETURN
      END
