      SUBROUTINE SABCD0(BUCK,BUF,IBUCK,IBUF,NINBCK,ICHAIN,
     &                  NBUCK,NDBCK,NBKINT,NREC,IRECL,
     &                  ILNBUF,ISYM,IUHF,GRAD,BUF2,IBUF2,IWHERE,
     &                  ILOOKUP,IMASTER)
C
C SUBROUTINE CREATES A SORT FILE FOR <AB|CD> INTEGRALS, AB SPIN CASE,
C   RHF AND UHF REFERENCE.
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
C   ISYM -  INVERSES SYMMETRY VECTOR REQUIRE TO DETERMINE THE
C           NUMBER OF A DISTRIBUTION IN THE SYMMETRY PACKED
C           INTEGRAL LIST
C
C   IUHF -  RHF/UHF FLAG
C
CEND
C
C  CODED JULY/90   JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GRAD
      CHARACTER*8 INAME 
      CHARACTER*80 FNAME 
      INTEGER DISCD,DISAD,DISCB,DISAB,DISDA,DISBC,DISBA,DISDC
      INTEGER A,B,C,D,AND,OR
      PARAMETER (LUINT=10)
      PARAMETER (LUSRT=15)
      DIMENSION BUCK(NBKINT,NBUCK),IBUCK(NBKINT,NBUCK)
     &,ISYM(1),IWHERE(8*ILNBUF),BUF2(8*ILNBUF),IBUF2(8*ILNBUF)
      dimension ilookup(8*ilnbuf),imaster(ndbck*nbuck)
      DIMENSION IBUF(ILNBUF),NINBCK(NBUCK),ICHAIN(NBUCK),BUF(ILNBUF)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
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
C STATEMENT FUNCTIONS FOR SQUARE TRIANGULAR INDEXING.
C
      INDX(I,J,N)=I+(J-1)*N
C
C OPEN UP SORTFILE AND INTEGRAL FILE.
C
      NOCC=NOCCO(1)
      NVRT=NVRTO(1)
      IF(IUHF.EQ.0) THEN
       INAME='PPPPAA  '
      ELSE
       INAME='PPPPAB  '
      ENDIF
      CALL GFNAME(INAME,FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD',ACCESS='SEQUENTIAL')
      CALL GFNAME('SRTFL1  ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
      CALL GFNAME('SRTFL2  ',FNAME,ILENGTH)
      OPEN(UNIT=LUSRT+1,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)

C
C INITIALIZE SOME CRAP AND COMPUTE MASTER INDEX.
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
C DIFFER HERE BETWEEN RHF AND UHF CASES SINCE THE INTEGRALS
C ARE TREATED IN A DIFFERENT WAY
C
C HERE RHF :
C
      IF(IUHF.EQ.0.AND..NOT.GRAD) THEN
C
C READ A BUFFER OF INTEGRALS AND DEAL WITH IT.
C
1     READ(LUINT)BUF,IBUF,NUT
      istick=0
      call izero(iwhere,8*ilnbuf)
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
C                  <CB|AD>,<AD|CB>,<CD|AB>,<BC|DA>,<DA|BC>,<DC|BA> AND <
C
C  FOR TYPE II: WE NEED TO WRITE OUT:
C                  <AB|AD>,<AD|AB>,<DA|BA> AND <BA|DA>.
C
C  FOR TYPE III: WE NEED TO WRITE OUT:
C                  <AB|AB> AND <BA|BA>
C
C  FOR TYPE IV: WE ONLY NEED TO WRITE ONE OUT.
C
C ASSIGN THE DISTRIBUTION NUMBER AND DETERMINE WHICH CLASS IT BELONGS TO
C
C HOWEVER, IN THE NEW SPIN-ADAPTED SCHEME, WE STORE ONLY <AB//CD> WITH
C A.LE.B,C,D, SO CHECK ALWAYS THE FIRST TWO INDICES IF THEY MEET THIS
C REQUIREMENT. THIS PROCEDURE REDUCES THE SIZE OF THE MOINTS FILE
C  APPROXIMATELY BY A FACTOR OF TWO !
C
       ITYPE=1
       IF(A.EQ.C.OR.B.EQ.D)ITYPE=2
       IF(A.EQ.C.AND.B.EQ.D)THEN
        ITYPE=3
        IF(A.EQ.B)ITYPE=4
       ENDIF
C
C PUT THIS INTEGRAL INTO A BUCKET BASED ON ITS CD INDEX.
C  
       IF(A.LE.B) THEN
       DISCD=ISYM(INDX(C,D,NVRT))
       IBKET=IMASTER(DISCD)
       ITMP=IPACK(A,B,C,D)
       iwhere(istick+1)=ibket
       ibuf2(istick+1)=itmp
       buf2(istick+1)=buf(int)
       ENDIF
C
C IF THE INTEGRAL IS NOT OF TYPE IV, WE ARE NOT DONE.
C

       IF(ITYPE.NE.4)THEN
        IF(B.LE.A) THEN
        DISDC=ISYM(INDX(D,C,NVRT))
        ITMP=IPACK(B,A,D,C)
        IBKET=IMASTER(DISDC)
        iwhere(istick+2)=ibket
        ibuf2(istick+2)=itmp
        buf2(istick+2)=buf(int)
        ENDIF
       ENDIF
C
C IF INTEGRAL IS OF TYPE 1 OR 2, STILL MORE TO GO.
C
       IF(ITYPE.LE.2)THEN
        IF(C.LE.D) THEN 
        DISAB=ISYM(INDX(A,B,NVRT))
        ITMP=IPACK(C,D,A,B)
        IBKET=IMASTER(DISAB)
        iwhere(istick+3)=ibket
        ibuf2(istick+3)=itmp
        buf2(istick+3)=buf(int)
        ENDIF
  
        IF(D.LE.C) THEN
        DISBA=ISYM(INDX(B,A,NVRT))
        ITMP=IPACK(D,C,B,A)
        IBKET=IMASTER(DISBA)
        iwhere(istick+4)=ibket
        ibuf2(istick+4)=itmp
        buf2(istick+4)=buf(int)
        ENDIF 
       ENDIF
C
C IF TYPE 1, STILL A FEW MORE.
C
       IF(ITYPE.EQ.1)THEN
        IF(C.LE.B) THEN
        DISAD=ISYM(INDX(A,D,NVRT))
        ITMP=IPACK(C,B,A,D)
        IBKET=IMASTER(DISAD)
        iwhere(istick+5)=ibket
        ibuf2(istick+5)=itmp
        buf2(istick+5)=buf(int)
        ENDIF
 
        IF(B.LE.C) THEN
        DISDA=ISYM(INDX(D,A,NVRT))
        ITMP=IPACK(B,C,D,A)
        IBKET=IMASTER(DISDA)
        iwhere(istick+6)=ibket
        ibuf2(istick+6)=itmp
        buf2(istick+6)=buf(int)
        ENDIF
 
        IF(D.LE.A) THEN
        DISBC=ISYM(INDX(B,C,NVRT))
        ITMP=IPACK(D,A,B,C)
        IBKET=IMASTER(DISBC)
        iwhere(istick+7)=ibket
        ibuf2(istick+7)=itmp
        buf2(istick+7)=buf(int)
        ENDIF
 
        IF(A.LE.D) THEN
        DISCB=ISYM(INDX(C,B,NVRT))
        ITMP=IPACK(A,D,C,B)
        IBKET=IMASTER(DISCB)
        iwhere(istick+8)=ibket
        ibuf2(istick+8)=itmp
        buf2(istick+8)=buf(int)
        ENDIF
       ENDIF
        istick=istick+8
10    CONTINUE
      NUMINT=NUMINT+NUT
C
C DEAL WITH STUFF NOW
C
      DO 102 ITYPE=1,NBUCK
       CALL WHENEQ(8*ILNBUF,IWHERE,1,ITYPE,ILOOKUP,NVAL)
       ioff=1
       inbuck=ninbck(itype)
11     idump=min(nval,nbkint-inbuck)
       call gather(idump,buck(inbuck+1,itype),buf2,ilookup(ioff))
       call igather(idump,ibuck(inbuck+1,itype),ibuf2,ilookup(ioff))
       inbuck=inbuck+idump 
       if(inbuck.eq.nbkint)then
        call plunk0(lusrt,buck(1,itype),ibuck(1,itype),ichain(itype),
     &              inbuck,nbkint,nrec)
        nval=nval-idump
        ioff=ioff+idump
        if(nval.ne.0)goto 11
       endif
       ninbck(itype)=inbuck
102   CONTINUE
C
C NOW CLEAR OUT BUFFERS, WRITING OUT WHATEVER IS LEFT.
C
      IF(NUT.EQ.ILNBUF)GOTO 1
C
      ELSE IF(IUHF.EQ.0.AND.GRAD) THEN
C
C READ A BUFFER OF INTEGRALS AND DEAL WITH IT.
C
101   READ(LUINT)BUF,IBUF,NUT
      call izero(iwhere,8*ilnbuf)
      istick=0
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 110 INT=1,NUT
       A=IUPKI(IBUF(INT))
       B=IUPKJ(IBUF(INT))
       C=IUPKK(IBUF(INT))
       D=IUPKL(IBUF(INT))
C
C INTEGRAL WILL BE ONE OF THREE TYPES:
C   <AB|CD>, <AB|AD> OR <AB|CB>, <AB|AB> AND <AA|AA>.
C
C  FOR TYPE I: WE ALSO MUST WRITE OUT SEVEN OTHER PERMUTATIONS:
C                  <CB|AD>,<AD|CB>,<CD|AB>,<BC|DA>,<DA|BC>,<DC|BA> AND <
C
C  FOR TYPE II: WE NEED TO WRITE OUT:
C                  <AB|AD>,<AD|AB>,<DA|BA> AND <BA|DA>.
C
C  FOR TYPE III: WE NEED TO WRITE OUT:
C                  <AB|AB> AND <BA|BA>
C
C  FOR TYPE IV: WE ONLY NEED TO WRITE ONE OUT.
C
C ASSIGN THE DISTRIBUTION NUMBER AND DETERMINE WHICH CLASS IT BELONGS TO
C
       ITYPE=1
       IF(A.EQ.C.OR.B.EQ.D)ITYPE=2
       IF(A.EQ.C.AND.B.EQ.D)THEN
        ITYPE=3
        IF(A.EQ.B)ITYPE=4
       ENDIF
C
C PUT THIS INTEGRAL INTO A BUCKET BASED ON ITS CD INDEX.
C  
       DISCD=ISYM(INDX(C,D,NVRT))
       IBKET=IMASTER(DISCD)
       ITMP=IPACK(A,B,C,D)
       iwhere(istick+1)=ibket
       ibuf2(istick+1)=itmp
       buf2(istick+1)=buf(int)
C
C IF THE INTEGRAL IS NOT OF TYPE IV, WE ARE NOT DONE.
C

       IF(ITYPE.NE.4)THEN
        DISDC=ISYM(INDX(D,C,NVRT))
        ITMP=IPACK(B,A,D,C)
        IBKET=IMASTER(DISDC)
        iwhere(istick+2)=ibket
        ibuf2(istick+2)=itmp
        buf2(istick+2)=buf(int)
       ENDIF
C
C IF INTEGRAL IS OF TYPE 1 OR 2, STILL MORE TO GO.
C
       IF(ITYPE.LE.2)THEN
        DISAB=ISYM(INDX(A,B,NVRT))
        ITMP=IPACK(C,D,A,B)
        IBKET=IMASTER(DISAB)
        iwhere(istick+3)=ibket
        ibuf2(istick+3)=itmp
        buf2(istick+3)=buf(int)
  
        DISBA=ISYM(INDX(B,A,NVRT))
        ITMP=IPACK(D,C,B,A)
        IBKET=IMASTER(DISBA)
        iwhere(istick+4)=ibket
        ibuf2(istick+4)=itmp
        buf2(istick+4)=buf(int)
       ENDIF
C
C IF TYPE 1, STILL A FEW MORE.
C
       IF(ITYPE.EQ.1)THEN
        DISAD=ISYM(INDX(A,D,NVRT))
        ITMP=IPACK(C,B,A,D)
        IBKET=IMASTER(DISAD)
        iwhere(istick+5)=ibket
        ibuf2(istick+5)=itmp
        buf2(istick+5)=buf(int)
 
        DISDA=ISYM(INDX(D,A,NVRT))
        ITMP=IPACK(B,C,D,A)
        IBKET=IMASTER(DISDA)
        iwhere(istick+6)=ibket
        ibuf2(istick+6)=itmp
        buf2(istick+6)=buf(int)
 
        DISBC=ISYM(INDX(B,C,NVRT))
        ITMP=IPACK(D,A,B,C)
        IBKET=IMASTER(DISBC)
        iwhere(istick+7)=ibket
        ibuf2(istick+7)=itmp
        buf2(istick+7)=buf(int)

        DISCB=ISYM(INDX(C,B,NVRT))
        ITMP=IPACK(A,D,C,B)
        IBKET=IMASTER(DISCB)
        iwhere(istick+8)=ibket
        ibuf2(istick+8)=itmp
        buf2(istick+8)=buf(int)
       ENDIF
       istick=istick+8
110    CONTINUE
      NUMINT=NUMINT+NUT
C
C DEAL WITH STUFF NOW
C
      DO 103 ITYPE=1,NBUCK
       CALL WHENEQ(8*ILNBUF,IWHERE,1,ITYPE,ILOOKUP,NVAL)
       ioff=1
       inbuck=ninbck(itype)
12     idump=min(nval,nbkint-inbuck)
       call gather(idump,buck(inbuck+1,itype),buf2,ilookup(ioff))
       call igather(idump,ibuck(inbuck+1,itype),ibuf2,ilookup(ioff))
       inbuck=inbuck+idump 
       if(inbuck.eq.nbkint)then
        call plunk0(lusrt,buck(1,itype),ibuck(1,itype),ichain(itype),
     &              inbuck,nbkint,nrec)
        nval=nval-idump
        ioff=ioff+idump
        if(nval.ne.0)goto 12
       endif
       ninbck(itype)=inbuck
103   CONTINUE
C
      IF(NUT.EQ.ILNBUF)GOTO 101
C
C  HERE UHF :
C
C  THE DIFFERENCE IS THAT WE CANNOT SWITCH A AND B SPINS AS IN RHF
C
      ELSE IF(IUHF.EQ.1) THEN
C
C READ A BUFFER OF INTEGRALS AND DEAL WITH IT.
C
1001  READ(LUINT)BUF,IBUF,NUT
      istick=0
      call izero(iwhere,4*ilnbuf)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 1010 INT=1,NUT
       A=IUPKI(IBUF(INT))
       B=IUPKJ(IBUF(INT))
       C=IUPKK(IBUF(INT))
       D=IUPKL(IBUF(INT))
C
C INTEGRAL WILL BE ONE OF THREE TYPES:
C   <AB|CD>, <AB|AD> OR <AB|CB>, <AB|AB> (=<AA|AA>).
C
C  FOR TYPE I: WE ALSO MUST WRITE OUT THREE OTHER PERMUTATIONS:
C                  <CB|AD>,<AD|CB>,<CD|AB>.
C
C  FOR TYPE II: WE NEED TO WRITE OUT:
C                  <AB|AD>,<AD|AB>.
C
C  FOR TYPE III: WE NEED TO WRITE OUT:
C                  <AB|AB>.
C
C ASSIGN THE DISTRIBUTION NUMBER AND DETERMINE WHICH CLASS IT BELONGS TO
C
       ITYPE=1
       IF(A.EQ.C.OR.B.EQ.D)ITYPE=2
       IF(A.EQ.C.AND.B.EQ.D)THEN
        ITYPE=3
       ENDIF
C
C PUT THIS INTEGRAL INTO A BUCKET BASED ON ITS CD INDEX.
C
       DISCD=ISYM(INDX(C,D,NVRT))
       IBKET=IMASTER(DISCD)
       ITMP=IPACK(A,B,C,D)
       iwhere(istick+1)=ibket
       ibuf2(istick+1)=itmp
       buf2(istick+1)=buf(int)
C
C IF INTEGRAL IS OF TYPE 1 OR 2, STILL MORE TO GO.
C
       IF(ITYPE.LE.2)THEN
        DISAB=ISYM(INDX(A,B,NVRT))
        ITMP=IPACK(C,D,A,B)
        IBKET=IMASTER(DISAB)
        iwhere(istick+2)=ibket
        ibuf2(istick+2)=itmp
        buf2(istick+2)=buf(int)
       ENDIF
C
C IF TYPE 1, STILL A FEW MORE.
C
       IF(ITYPE.EQ.1)THEN
        DISAD=ISYM(INDX(A,D,NVRT))
        ITMP=IPACK(C,B,A,D)
        IBKET=IMASTER(DISAD)
        iwhere(istick+3)=ibket
        ibuf2(istick+3)=itmp
        buf2(istick+3)=buf(int)
 
        DISCB=ISYM(INDX(C,B,NVRT))
        ITMP=IPACK(A,D,C,B)
        IBKET=IMASTER(DISCB)
        iwhere(istick+4)=ibket
        ibuf2(istick+4)=itmp
        buf2(istick+4)=buf(int)
       ENDIF
       istick=istick+4
1010  CONTINUE
      NUMINT=NUMINT+NUT
C
C DEAL WITH STUFF NOW
C
      DO 104 ITYPE=1,NBUCK
       CALL WHENEQ(4*ILNBUF,IWHERE,1,ITYPE,ILOOKUP,NVAL)
       ioff=1
       inbuck=ninbck(itype)
13     idump=min(nval,nbkint-inbuck)
       call gather(idump,buck(inbuck+1,itype),buf2,ilookup(ioff))
       call igather(idump,ibuck(inbuck+1,itype),ibuf2,ilookup(ioff))
       inbuck=inbuck+idump 
       if(inbuck.eq.nbkint)then
        call plunk0(lusrt,buck(1,itype),ibuck(1,itype),ichain(itype),
     &              inbuck,nbkint,nrec)
        nval=nval-idump
        ioff=ioff+idump
        if(nval.ne.0)goto 13
       endif
       ninbck(itype)=inbuck
104   CONTINUE
      IF(NUT.EQ.ILNBUF)GOTO 1001
      ENDIF
      WRITE(6,1400)NUMINT
1400  FORMAT(T3,'@SABCD0-I, ',I8,' PPPP integrals sorted.')
C
C FLUSH REMAINING BUFFERS.
C
      DO 20 I=1,NBUCK
       CALL PLUNK0(LUSRT,BUCK(1,I),IBUCK(1,I),ICHAIN(I),NINBCK(I),
     &            NBKINT,NREC)
20    CONTINUE
      NREC=NREC-1
      CLOSE(UNIT=LUINT,STATUS='DELETE')
      CLOSE(UNIT=LUSRT,STATUS='KEEP')
      CLOSE(UNIT=LUSRT+1,STATUS='KEEP')
      RETURN
      END
