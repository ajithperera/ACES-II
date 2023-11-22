      SUBROUTINE GMOIDAA(IORDER,BUCK,IBUCK,BUF,IBUF,NO,NBKINT,
     &                  ISPIN,NIRREP,NSTO,IWHERE,ILOOKUP,
     &                  IPNTR,IPNTD,NBSALL,IUHF)
C
C SUBROUTINE TO READ MO INTEGRALS FROM HF2 AND PUT THEM IN
C  ONE OF THE FOLLOWING CATEGORIES: PPPP,PPPH,PPHH,PHPH,PHHH,HHHH
C  AND THEN WRITE THEM OUT TO DISK.  NOTE: IN HF2 THE INDICES ARE
C  STORED IN MULLIKEN FORMAT.  SWITCH HERE TO DIRAC TO FACILITATE
C  CORRELATED CALCULATIONS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,ABCDTYPE
      LOGICAL ISOPN,DOALL,ABIJ,NOABCD,YESNO,UNLOAD(7)
      CHARACTER*10 STRING
      CHARACTER*4 SPCASE(3)
      CHARACTER*1 SP(2)
      CHARACTER*8 INAME,HF2NAME
      CHARACTER*80 FNAME
      PARAMETER (NFPPPP = 16)
      PARAMETER (NFPPPH = 17)
      PARAMETER (NFPPHH = 18)
      PARAMETER (NFPHPH = 19)
      PARAMETER (NFPHHH = 20)
      PARAMETER (NFHHHH = 21)
C
C IGRNBF IS THE LENGTH OF INTEGRAL RECORDS IN THE HF2 FILE WRITTEN BY
C  GRNFNC.  NBKINT IS THE LENGTH OF RECORDS WRITTEN BY GMOIAA.
C
      PARAMETER (IGRNBF =600)
      DIMENSION BUCK(NBKINT,7),IBUCK(NBKINT,7),ILOOKUP(IGRNBF)
      DIMENSION BUF(IGRNBF),IBUF(IGRNBF),IORDER(NSTO),IWHERE(IGRNBF)
      DIMENSION IPNTR(NBSALL),IPNTD(NBSALL) 
      DIMENSION NFULL(7),IU(7),NINBCK(7)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      COMMON /ABCD/ ABCDTYPE
      DATA SPCASE /'AA  ','BB  ','AB  '/
      DATA SP /'A','B'/
      DATA IU /21,20,19,18,17,17,16/
      IPCK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),
     &ISHFT(K,2*IBITWD)),ISHFT(L,3*IBITWD))
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      indp(i)=(i*(I+1))/2
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
c-------------------------------------------------------------------
c make ipntd for drop-mo analytic gradient  ---  Mar. '95  --  KB --
c-------------------------------------------------------------------
      call izero(ipntd,nbsall)
      call getrec(20,'JOBARC','IDROP'//SP(ISPIN),NBSALL,IPNTR) 
      call getrec(20,'JOBARC','NUMDROPA',1,ndrop) 
      call getrec(20,'JOBARC','NDROTPOP',1,NDRCO)
      NOa = NO + NDRCO 
      iidr = 0
      do 50 i=1,nbsall
       ipntd(i) = iidr 
       if (ipntr(i).ne.0) then
        iidr = iidr + 1
       endif
 50   continue
      if (iidr.ne.ndrop) then
       write (6,55) iidr,ndrop
 55    format (4x,'==== there must be some error,  iidr =',i3,3x,
     &         'must be equal to ndrop=',i3)
       call errex
      endif 
c----------------------
c     write (6,57) (ipntr(i),i=1,nbsall)
c57   format (2x,'ipntr=',30i2) 
c     write (6,58) (ipntd(j),j=1,nbsall)
c58   format (2x,'ipntd=',30i2) 
c     write (6,59) iidr,ndrop
c59   format (2x,'iidr  ndrop =',2i5) 
c-------------------------------------------------------------------
      WRITE(LUOUT,1500)SPCASE(ISPIN)
1500  FORMAT(T3,'@GMOIAA-I, Processing MO integrals for spin ',
     &          'case ',A2,'.')
      CALL IZERO(NINBCK,7)
      CALL IZERO(NFULL,7)
C
C LOCATE FIRST INTEGRAL RECORD FOR GIVEN SPIN CASE.
C
      IF(ISPIN.EQ.1) THEN
       HF2NAME='HF2AA   '
      ELSE IF(ISPIN.EQ.2) THEN
       HF2NAME='HF2BB   '
      ENDIF    
      IF(IUHF.EQ.0) HF2NAME='HF2     '
      CALL GFNAME(HF2NAME,FNAME,ILENGTH)
      OPEN(UNIT=22,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='OLD')
      IF(IUHF.EQ.0) THEN
 1     READ(22)STRING
       IF(INDEX(STRING,'**').EQ.0)GOTO 1
      ENDIF
      IF(.NOT.NOABCD.AND.ABCDTYPE.EQ.0)THEN
       INAME='PPPP'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
       OPEN(UNIT=16,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
      ENDIF
      IF(.NOT.ABIJ)THEN
       INAME='HHHH'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
       OPEN(UNIT=21,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       INAME='PPPH'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
       OPEN(UNIT=17,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       INAME='PHPH'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
       OPEN(UNIT=19,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
       INAME='PHHH'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
       OPEN(UNIT=20,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     & STATUS='NEW')
      ENDIF
       INAME='PPHH'//SPCASE(ISPIN)
       CALL GFNAME(INAME,FNAME,ILENGTH)
      OPEN(UNIT=18,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &STATUS='NEW')
C
C NEW CODE TO PREVENT WRITES TO UNOPENED FILES
C
      DO 340 IUNIT=1,7
       INQUIRE(UNIT=IU(IUNIT),OPENED=YESNO)
       UNLOAD(IUNIT)=YESNO
340   CONTINUE
C      
 2    READ(22,END=999)BUF,IBUF,NUT
      IF(NUT.EQ.-1) GO TO 103
      DO 11 IK=1,NUT
       I=IUPKI(IBUF(IK))
       J=IUPKJ(IBUF(IK))
       K=IUPKK(IBUF(IK))
       L=IUPKL(IBUF(IK))
c-------------------------------------------------
       IWHERE(IK)=0
       if (ipntr(l).ne.0)  go to 10
       if (ipntr(k).ne.0)  go to 10
       if (ipntr(j).ne.0)  go to 10
       if (ipntr(i).ne.0)  go to 10
c-------------------------------------------------
c
       IF(NIRREP.NE.1)THEN
C
C INSURE CANONICAL ORDERING OF INDICES TO FACILITATE TESTING.
C
        IF(J.LT.I)THEN
         IT=I
         I=J
         J=IT
        ENDIF
        IF(L.LT.K)THEN
         IT=K
         K=L
         L=IT
        ENDIF
        IF(INDP(J)+I.GT.INDP(L)+K)THEN
         IT=I
         I=K
         K=IT
         IT=J
         J=L
         L=IT
        ENDIF
       ENDIF
C
       I0=I - ipntd(i)
       J0=J - ipntd(j)
       K0=K - ipntd(k)
       L0=L - ipntd(l)
       IF(I.GT.NOa)I0=I-NO -ipntd(i) 
       IF(J.GT.NOa)J0=J-NO -ipntd(j) 
       IF(K.GT.NOa)K0=K-NO -ipntd(k) 
       IF(L.GT.NOa)L0=L-NO -ipntd(l) 
C
       ibuf(IK)=ipck(i0,k0,j0,l0)
C
C  Designations above refer to Dirac indices.
C
       IF(I.GT.NOa.AND.K.GT.NOa)THEN
        IWHERE(IK)=7 
       ELSEIF(K.GT.NOa.AND.J.GT.NOa.AND.L.GT.NOa)THEN
        IWHERE(IK)=5
       ELSEIF(I.GT.NOa.AND.J.GT.NOa.AND.L.GT.NOa)THEN
        IWHERE(IK)=5
        ibuf(IK)=ipck(k0,i0,l0,j0)
       ELSEIF(L.LE.NOa)THEN
         IWHERE(IK)=1
       ELSEIF(K.LE.NOa.AND.J.LE.NOa)THEN
         IWHERE(IK)=2
       ELSEIF(K.GT.NOa.AND.J.LE.NOa)THEN
         IWHERE(IK)=3
       ELSEIF(K.LE.NOa.AND.J.GT.NOa.AND.I.LE.NOa)THEN
         IWHERE(IK)=4
       ENDIF
 10   CONTINUE
 11   CONTINUE
C
C DEAL WITH STUFF NOW
C
      DO 102 ITYPE=1,7
       IF(UNLOAD(ITYPE))THEN
        CALL WHENEQ(NUT,IWHERE,1,ITYPE,ILOOKUP,NVAL)
        INBUCK=NINBCK(ITYPE)
        IOFF=1
21      IDUMP=MIN(NVAL,NBKINT-INBUCK)
        CALL GATHER(IDUMP,BUCK(INBUCK+1,ITYPE),BUF,ILOOKUP(IOFF))
        CALL IGATHER(IDUMP,IBUCK(INBUCK+1,ITYPE),IBUF,ILOOKUP(IOFF))
        INBUCK=INBUCK+IDUMP
        IF(INBUCK.EQ.NBKINT)THEN
         CALL DROPIT(IU(ITYPE),BUCK(1,ITYPE),IBUCK(1,ITYPE),NBKINT,
     &               NBKINT)
         INBUCK=0
         NFULL(ITYPE)=NFULL(ITYPE)+1
         NVAL=NVAL-IDUMP
         IOFF=IOFF+IDUMP
         IF(NVAL.NE.0)GOTO 21
        ENDIF
        NINBCK(ITYPE)=INBUCK
       ENDIF
102   CONTINUE
c     write (6,9861) (nfull(i),i=1,7),(ninbck(j),j=1,7),nbkint
C        
      GOTO 2
C
103   CONTINUE
C
C NOW CLEAR OUT BUFFERS, WRITING OUT WHATEVER IS LEFT.
C
      DO 333 ITYPE=1,7
       if(itype.eq.6)goto 333
       IF(UNLOAD(ITYPE))THEN
        CALL DROPIT(IU(ITYPE),BUCK(1,ITYPE),IBUCK(1,ITYPE),NBKINT,
     &              NINBCK(ITYPE))
       ENDIF
333   CONTINUE
C
C CLOSE FILES.
C
      DO 12 I=16,21,1
       INQUIRE(UNIT=I,OPENED=ISOPN)
       IF(ISOPN)CLOSE(UNIT=I,STATUS='KEEP')
 12   CONTINUE
c     IF(IFLAGS(84).NE.2) THEN
c      CLOSE(UNIT=22,STATUS='DELETE')
c     ELSE 
       CLOSE(UNIT=22,STATUS='KEEP')
c     ENDIF
C
C WRITE SUMMARY.
C
c--------------------
c     write (6,9861) (nfull(i),i=1,7),(ninbck(j),j=1,7),nbkint
c9861 format (4x,'==== nfull =',7i7/
c    x        4x,'==== ninbck=',7i7/
c    x        4x,'==== nbkint=',i7)
c--------------------
      NIPPPP=NFULL(7)*NBKINT+NINBCK(7)
      NIPPPH=NFULL(6)*NBKINT+NINBCK(6)
      NIPPPH=NIPPPH+NFULL(5)*NBKINT+NINBCK(5)
      NIPPHH=NFULL(4)*NBKINT+NINBCK(4)
      NIPHPH=NFULL(3)*NBKINT+NINBCK(3)
      NIPHHH=NFULL(2)*NBKINT+NINBCK(2)
      NIHHHH=NFULL(1)*NBKINT+NINBCK(1)
      NINT=NIPPPP+NIPPPH+NIPPHH+NIPHPH+NIPHHH+NIHHHH
      WRITE(LUOUT,1000)
      WRITE(LUOUT,1001)NIPPPP,NIPPPH,NIPPHH,NIPHPH,NIPHHH,NIHHHH,
     &NINT
 1000 FORMAT(T3,'@GMOIAA-I, Generation of integral list completed.')
 1001 FORMAT(T23,'TYPE',T39,'NUMBER',/,T23,'----',T38,'--------',/,
     &T23,'PPPP',T37,I9,/,T23,
     &'PPPH',T37,I9,/,T23,'PPHH',T37,I9,/,T23,'PHPH',T37,I9,/,
     &T23,'PHHH',T37,I9,/,T23,'HHHH',T37,I9,//,T22,'TOTAL',T37,
     &I9)
      RETURN
999   WRITE(LUOUT,1400)
1400  FORMAT(T3,'@GMOIAA-F, Unexpected end-of-file on HF2.')
      CALL ERREX
      END
