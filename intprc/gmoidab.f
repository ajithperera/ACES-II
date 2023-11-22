      SUBROUTINE GMOIDAB(IORDRA,IORDRB,BUCK,IBUCK,BUF,IBUF,
     &                  NODA,NODB,NBKINT,ISPIN,NIRREP,NORB,
     &          IWHERE,ILOOKUP,IPNTRA,IPNTRB,IPNTDA,IPNTDB,NBSALL)
C
C SUBROUTINE TO READ MO INTEGRALS FROM HF2 AND PUT THEM IN
C  ONE OF THE FOLLOWING CATEGORIES: PPPP,PPPH,PPHH,PHPH,PHHH,HHHH
C  AND THEN WRITE THEM OUT TO DISK.  NOTE: IN HF2 THE INDICES ARE
C  STORED IN MULLIKEN FORMAT.  SWITCH HERE TO DIRAC TO FACILITATE
C  CORRELATED CALCULATIONS.  THIS ROUTINE IS SPECIFIC FOR THE
C  ALPHA-BETA BLOCK, WHICH IS STORED SQUARE TRIANGULAR INSTEAD
C  OF CANONICAL.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,ABCDTYPE
      LOGICAL ISOPN,DOALL,ABIJ,NOABCD,YESNO,UNLOAD(9)
      CHARACTER*10 STRING
      CHARACTER*4 SPCASE(3)
      CHARACTER*80 FNAME
      PARAMETER (NFPPPP   = 16)
      PARAMETER (NFPPPH1H = 17)
      PARAMETER (NFPPPH2H = 18)
      PARAMETER (NFPPHH   = 19)
      PARAMETER (NFPHPH1P = 20)
      PARAMETER (NFPHPH2P = 21)
      PARAMETER (NFPHHH1P = 22)
      PARAMETER (NFPHHH2P = 23)
      PARAMETER (NFHHHH   = 24)
C
C IGRNBF IS THE LENGTH OF INTEGRAL RECORDS IN THE HF2 FILE WRITTEN BY
C  GRNFNC.  ILNBUF IS THE LENGTH OF RECORDS WRITTEN BY GMOIAB.
C
      PARAMETER (IGRNBF =600)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /ABCD/ ABCDTYPE
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      DIMENSION BUCK(NBKINT,9),IBUCK(NBKINT,9)
      DIMENSION BUF(IGRNBF),IBUF(IGRNBF)
      DIMENSION IORDRA(NORB),IORDRB(NORB)
      DIMENSION IPNTRA(NBSALL),IPNTRB(NBSALL)
      DIMENSION IPNTDA(NBSALL),IPNTDB(NBSALL) 
      DIMENSION NTOT(9),NINBCK(9),ILOOKUP(IGRNBF),IWHERE(IGRNBF)
      DATA SPCASE /'AA  ','BB  ','AB  '/
      IPCK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),
     &ISHFT(K,2*IBITWD)),ISHFT(L,3*IBITWD))
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      INDP(I)=(I*(I-1))/2
c-------------------------------------------------------------------
c make ipntd for drop-mo analytic gradient  ---  Mar. '95  --  KB --
c-------------------------------------------------------------------
      call izero(ipntda,nbsall) 
      call izero(ipntdb,nbsall) 
      call getrec(20,'JOBARC','NUMDROPA',1,ndrop) 
      call getrec(20,'JOBARC','NDROTPOP',1,ndrco) 
      call getrec(20,'JOBARC','NDROTPOP',1,ndrco) 
      call getrec(20,'JOBARC','IDROPA  ',NBSALL,IPNTRA) 
      call getrec(20,'JOBARC','IDROPB  ',NBSALL,IPNTRB) 
      NOA = NODA + ndrco
      NOB = NODB + ndrco
      iidra = 0
      iidrb = 0
      do 50 i=1,nbsall
       ipntda(i) = iidra 
       ipntdb(i) = iidrb 
       if (ipntra(i).ne.0) then
        iidra = iidra + 1
       endif
       if (ipntrb(i).ne.0) then
        iidrb = iidrb + 1
       endif
 50   continue
      if ((iidra.ne.ndrop) .or. (iidrb.ne.ndrop)) then
       write (6,55) iidra,iidrb,ndrop
 55    format (4x,' there must be some error,  iidra iidrb=',2i3,3x,
     &         'must be equal to ndrop=',i3)
       call errex
      endif 
c----------------------
c     write (6,57) (ipntra(i),i=1,nbsall)
c     write (6,57) (ipntrb(i),i=1,nbsall)
c57   format (2x,'ipntr=',30i2) 
c     write (6,58) (ipntda(j),j=1,nbsall)
c     write (6,58) (ipntdb(j),j=1,nbsall)
c58   format (2x,'ipntd=',30i2) 
c     write (6,59) iidr,ndrop
c59   format (2x,'iidr  ndrop =',2i5) 
c-------------------------------------------------------------------
      WRITE(LUOUT,1500)SPCASE(ISPIN)
1500  FORMAT(T3,'@GMOIAB-I, Processing MO integrals for spin ',
     &          'case ',A2,'.')
C
C LOCATE FIRST INTEGRAL RECORD FOR GIVEN SPIN CASE.
C
      CALL GFNAME('HF2AB   ',FNAME,ILENGTH)
      OPEN(UNIT=25,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
 1    READ(25)STRING
      IF(INDEX(STRING,'**').EQ.0)GOTO 1
      IF(.NOT.NOABCD.AND.ABCDTYPE.EQ.0)THEN
       CALL GFNAME('PPPPAB  ',FNAME,ILENGTH)
       OPEN(UNIT=16,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
      ENDIF
      IF(.NOT.ABIJ)THEN
       CALL GFNAME('HHHHAB  ',FNAME,ILENGTH)
       OPEN(UNIT=24,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PPPH1H  ',FNAME,ILENGTH)        
       OPEN(UNIT=17,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PPPH2H  ',FNAME,ILENGTH)        
       OPEN(UNIT=18,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PHPH1P  ',FNAME,ILENGTH)        
       OPEN(UNIT=20,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PHPH2P  ',FNAME,ILENGTH)        
       OPEN(UNIT=21,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PHHH1P  ',FNAME,ILENGTH)        
       OPEN(UNIT=22,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       CALL GFNAME('PHHH2P  ',FNAME,ILENGTH)        
       OPEN(UNIT=23,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
      ENDIF
      CALL GFNAME('PPHHAB  ',FNAME,ILENGTH)        
      OPEN(UNIT=19,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &STATUS='NEW')
C
C CODE TO PREVENT WRITES TO UNOPENED FILES
C
      DO 340 IUNIT=1,9
       INQUIRE(UNIT=15+IUNIT,OPENED=YESNO)
       UNLOAD(IUNIT)=YESNO
340   CONTINUE
C
      NINT=0
      CALL IZERO(NINBCK,9)
      CALL IZERO(NTOT,9)
C
      norec = 1 
 2    READ(25,END=999)BUF,IBUF,NUT
      IF(NUT.EQ.-1) GO TO 103
      DO 11 IK=1,NUT
C
C CRAPS PUTS THE AB INTEGRALS OUT AS BETA-ALPHA INSTEAD OF THE
C  LOGICALLY SENSIBLE ALPHA-BETA.  SWITCH SENSE HERE.
C
       K=IUPKI(IBUF(IK))
       L=IUPKJ(IBUF(IK))
       I=IUPKK(IBUF(IK))
       J=IUPKL(IBUF(IK))
c-------------------------------------------------
       IWHERE(ik)=0
       if (ipntrb(l).ne.0) go to 10 
       if (ipntrb(k).ne.0) go to 10 
       if (ipntra(j).ne.0) go to 10 
       if (ipntra(i).ne.0) go to 10 
c-------------------------------------------------
       IF(NIRREP.NE.1)THEN
C
C ASSURE THAT J>I, L>K.
C
        IF(I.GT.J)THEN
         IT=J
         J=I
         I=IT
        ENDIF
        IF(K.GT.L)THEN
         IT=L
         L=K
         K=IT
        ENDIF
       ENDIF
       I0=I - ipntda(i) 
       J0=J - ipntda(j) 
       K0=K - ipntdb(k) 
       L0=L - ipntdb(l) 
       IF(I.GT.NOA)I0=I-NODA - ipntda(i) 
       IF(J.GT.NOA)J0=J-NODA - ipntda(j) 
       IF(K.GT.NOB)K0=K-NODB - ipntdb(k) 
       IF(L.GT.NOB)L0=L-NODB - ipntdb(l) 
C
C Indices will be stored using Dirac convention.
C
C
C  Now classify integral according to Mulliken indices (IJ|KL).
C  (These are stored upper triangular on IJ and on KL).
C    If J and L are occupied                   it must be HHHHAB.
C    If I and K are unoccupied                 it must be PPPPAB.
C
C  After branching for these possibilities, only (PP|xH),(PH|xP) and
C                                                (PH|xH) are left.
C   Now, if the possibilities are tested in order selection is based
C    on the first satisfied criterion:
C
C    If J is occupied and K is unoccupied      it must be PHPH1P.
C    If L is occupied and I is unoccupied      it must be PHPH2P.
C    If J is occupied and L is unoccupied      it must be PHHH2P.
C    If L is occupied and I is unoccupied      it must be PHHH1P.
C    If J and K are unoccupied                 it must be PPPH1H.
C    If L and I are unoccupied                 it must be PPPH2H.
C    If J and L are occupied                   it must be PPHH.
C
C  Designations above refer to Dirac indices.
C
       IBUF(IK)=IPCK(I0,K0,J0,L0)
       IF(J.LE.NOA.AND.L.LE.NOB)THEN
        IWHERE(IK)=9
       ELSEIF(I.GT.NOA.AND.K.GT.NOB)THEN
        IWHERE(IK)=1
       ELSEIF(L.LE.NOB.AND.I.GT.NOA)THEN
        IWHERE(IK)=5
       ELSEIF(J.LE.NOA.AND.K.GT.NOB)THEN
        IWHERE(IK)=6
       ELSEIF(J.LE.NOA.AND.L.GT.NOB)THEN
        IWHERE(IK)=8
       ELSEIF(L.LE.NOB.AND.J.GT.NOA)THEN
        IWHERE(IK)=7
       ELSEIF(J.GT.NOA.AND.K.GT.NOB)THEN
        IWHERE(IK)=2
       ELSEIF(L.GT.NOB.AND.I.GT.NOA)THEN
        IWHERE(IK)=3
       ELSEIF(J.GT.NOA.AND.L.GT.NOB)THEN
        IWHERE(IK)=4
       ENDIF
C
 10   CONTINUE
 11   CONTINUE
C
C DEAL WITH THIS INFORMATION NOW
C
      DO 102 ITYPE=1,9
       IF(UNLOAD(ITYPE))THEN
        CALL WHENEQ(NUT,IWHERE,1,ITYPE,ILOOKUP,NVAL)
        INBUCK=NINBCK(ITYPE)
        IOFF=1
12      IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUF,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUF,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT)THEN
         CALL DROPIT(ITYPE+15,BUCK(1,ITYPE),IBUCK(1,ITYPE),NBKINT,
     &               NBKINT)
         NTOT(ITYPE)=NTOT(ITYPE)+NBKINT
         INBUCK=0
         NVAL=NVAL-IDUMP
         IOFF=IOFF+IDUMP
         IF(NVAL.NE.0)GOTO 12
        ENDIF
        NINBCK(ITYPE)=INBUCK
       ENDIF
102   CONTINUE
C  
      GO TO 2
103   CONTINUE 
C
C NOW CLEAR OUT BUFFERS, WRITING OUT WHATEVER IS LEFT.
C
      DO 104 ITYPE=1,9
       IF(UNLOAD(ITYPE))THEN
        CALL DROPIT(ITYPE+15,BUCK(1,ITYPE),IBUCK(1,ITYPE),
     &              NBKINT,NINBCK(ITYPE))
        NTOT(ITYPE)=NTOT(ITYPE)+NINBCK(ITYPE)
        NINT=NINT+NTOT(ITYPE)
       ENDIF
104   CONTINUE
C
C CLOSE FILES.
C
      DO 13 I=16,24,1
       INQUIRE(UNIT=I,OPENED=ISOPN)
       IF(ISOPN)CLOSE(UNIT=I,STATUS='KEEP')
13    CONTINUE
c     IF(IFLAGS(84).NE.2) THEN 
c      CLOSE(UNIT=25,STATUS='DELETE')
c     ELSE
       CLOSE(UNIT=25,STATUS='KEEP')
c     ENDIF
C
C WRITE SUMMARY.
C
      WRITE(LUOUT,1001)(NTOT(ITYPE),ITYPE=1,9),NINT
 1000 FORMAT(T3,'@GMOIAB-I, Generation of integral list completed.')
 1001 FORMAT(T23,'TYPE',T39,'NUMBER',/,T23,'----',T38,'--------',/,
     &T23,'PPPP',T37,I9,/,T23,
     &'PPPH1H',T37,I9,/,T23,'PPPH2H',T37,I9,/,T23,'PPHH',T37,I9,/,
     &T23,'PHPH1P',T37,I9,/,T23,'PHPH2P',T37,I9,/,
     &T23,'PHHH1P',T37,I9,/,T23,'PHHH2P',T37,I9,/,
     &T23,'HHHH  ',T37,I9,//,T22,'TOTAL',T37,
     &I9)
      RETURN
999   WRITE(LUOUT,1400)
1400  FORMAT(T3,'@GMOIAB-F, Unexpected end-of-file on HF2.')
      CALL ERREX
      END
