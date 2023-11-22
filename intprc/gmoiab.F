      SUBROUTINE GMOIAB(IORDRA,IORDRB,BUCK,IBUCK,BUF,IBUF,
     &                  NOA,NOB,NBKINT,ISPIN,NIRREP,NORB,
     &                  IWHERE,ILOOKUP)
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
      DIMENSION NTOT(9),NINBCK(9),ILOOKUP(IGRNBF),IWHERE(IGRNBF)
      DATA SPCASE /'AA  ','BB  ','AB  '/
cjp
cjp
#include "bwcc.com"
cjp
cjp data for integral sorting  wrto general fermi vacuum
       integer iquad,indx,iperm,intype
       dimension iquad(4),iperm(4,0:15),intype(0:15)
c     lkji    type permutation to get aces ordering
c                  (must be element of symmetry group of index permutations
c                   of tei in mullikan conv.(1234,2134,1243,2143)),
c at the moment of sort still mulliken order is valid and kl indices are beta-spin
c
c reversed 
c mulliken
c     1234
c     lkji
cspin bbaa    code	perm	aces-notation-Dirac
c     HHhh	9	1234	HHHHAB
c     HHhp   	7	1243	->PHHH1P
c     HHph   	7	1234	PHHH1P
c     HHpp   	5	1234	PHPH1P
c     HPhh   	8	2134	->PHHH2P
c     HPhp   	4	2143	->PPHH
c     HPph   	4	2134	->PPHH
c     HPpp   	3	2134	->PPPH2H
c     PHhh   	8	1234	PHHH2P
c     PHhp  	4	1243	->PPHH
c     PHph  	4	1234 	PPHH
c     PHpp  	3	1234	PPPH2H
c     PPhh  	6	1234	PHPH2P
c     PPhp  	2 	1243	->PPPH1H
c     PPph  	2	1234	PPPH1H
c     PPpp 	1	1234	PPPPAB 
c
cjp corresponding codes
      data intype/9,7,7,5,8,4,4,3,8,4,4,3,6,2,2,1/
      data iperm/
     +1,2,3,4,
     +1,2,4,3,
     +1,2,3,4,
     +1,2,3,4,
     +2,1,3,4,
     +2,1,4,3,
     +2,1,3,4,
     +2,1,3,4,
     +1,2,3,4,
     +1,2,4,3,
     +1,2,3,4,
     +1,2,3,4,
     +1,2,3,4,
     +1,2,4,3,
     +1,2,3,4,
     +1,2,3,4/
cjp additional data for fock contribution construction
      integer fperm,fockind,nfperm
      parameter(nfperm=4)
      dimension fperm(4,nfperm)
      data fperm/
     +1,2,3,4,
     +2,1,3,4,
     +1,2,4,3,
     +2,1,4,3/
      dimension fockind(4,nfperm)
      logical b1,b2
      integer findex
      findex(i,j)=(i*(i-1))/2+j
cjp
      IPCK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),
     &ISHFT(K,2*IBITWD)),ISHFT(L,3*IBITWD))
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      INDP(I)=(I*(I-1))/2
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
 2    READ(25,END=999)BUF,IBUF,NUT
      IF(NUT.EQ.-1) GO TO 103
      DO 11 IK=1,NUT
C
C CRAPS PUTS THE AB INTEGRALS OUT AS BETA-ALPHA INSTEAD OF THE
C  LOGICALLY SENSIBLE ALPHA-BETA.  SWITCH SENSE HERE.
C
       IF(NIRREP.EQ.1)THEN
        K=IUPKI(IBUF(IK))
        L=IUPKJ(IBUF(IK))
        I=IUPKK(IBUF(IK))
        J=IUPKL(IBUF(IK))
       ELSE
        K=IUPKI(IBUF(IK))
        L=IUPKJ(IBUF(IK))
        I=IUPKK(IBUF(IK))
        J=IUPKL(IBUF(IK))
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
cjp
cjpccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjp fock contributions; here such integral can contrib to both beta and alpha
cjp fock matrices, but only 4 index permutations are allowed and only two
cjp possibilities for equal indices apply
cjp
       if(isbwcc) then
cjp check for all possibilities of equal at least two indices-only such integrals can contribute
       b1=i.eq.j
       b2=k.eq.l
       if(b1.or.b2) then
cjp generate all unique index permuted integrals
cjp  accumulate contributions to other fock matrices
       nperm=1
       iquad(1)=i
       iquad(2)=j
       iquad(3)=k
       iquad(4)=l
       do kperm=1,nfperm
       do ii=1,4
       fockind(ii,nperm)=iquad(fperm(ii,kperm))
       enddo
ctest if it is not identical with an older generated one
       do jperm=1,nperm-1
       isequal=1
       do ii=1,4
       if(fockind(ii,nperm).ne.fockind(ii,jperm)) goto 9910
       enddo
cequal found
       goto 9920
9910   continue
       enddo
c non equal to any previous
       nperm=nperm+1
9920   continue
       enddo
       nperm=nperm-1
c end kperm loop
cjp
cjp now we have all non-equivalent integrals (ij,kl)
cjp loop over such integrals
       do kperm=1,nperm
       myi=fockind(1,kperm)
       myj=fockind(2,kperm)
       myk=fockind(3,kperm)
       myl=fockind(4,kperm)
cjp add to one triangle of fock matrix difference
cjp here all such contributions are considered, addfock then itself selects
cjp which of them belong to the difference due to different occupation in
cjp reference of SCF and iref
c coulomb
cjp here only coulomb contributes, exchange vanishes due to different spins
cjp but coulomb must be considered for both alpha and beta blocks
c
       if(myi.ge.myj .and. myk.eq.myl)
     &     call addfock(myi,myj,myk,buf(ik),1,2,1)
c
       if(myi.eq.myj .and. myk.ge.myl)
     &     call addfock(myk,myl,myi,buf(ik),2,1,1)
c
       enddo
cjp        kperm loop
       endif
cjp      b1.or.b2
       endif
cjp         isbwcc
cjp end fock contribution
cjpccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cjp
cjp choose between old logic and general multireference
cjp
      if(isbwcc) then
cjp general algorithm for mr-bwcc case, in fact, it does not require canonical ordering at all
      iquad(4)=iocc(i,iref,1)
      iquad(3)=iocc(j,iref,1)
      iquad(2)=iocc(k,iref,2)
      iquad(1)=iocc(l,iref,2)
cjp compute the indx ... binary number with digits lkji; p=1, h=0
      indx=0
      do ii=1,4
        indx=2*indx+1-iquad(ii)
      enddo
cjp convert to fermi-vacuum indices
      iquad(4)=iphnum(i,iref,1)
      iquad(3)=iphnum(j,iref,1)
      iquad(2)=iphnum(k,iref,2)
      iquad(1)=iphnum(l,iref,2)
cjp permute the quadruple
      i0=iquad(iperm(4,indx))
      j0=iquad(iperm(3,indx))
      k0=iquad(iperm(2,indx))
      l0=iquad(iperm(1,indx))
      iwhere(ik)=intype(indx)
cjp switch to dirac indexing and store
      ibuf(ik)=ipck(i0,k0,j0,l0)
cjp done
cjp
cjp for bwcc write out integrals with indices in active space just for information
cjp
       if(isactive(i,1).and.isactive(j,1).and.
     &   isactive(k,2).and.isactive(l,2)) then
         if(bwgossip.or.iflags(1).ge.10)
     & write(6,*)'@GMOIAB-I integral of active orbs: (',
     &          i,j,'|',K,L,') = <',i,K,'|',j,L,'> = ',buf(ik)
       endif
cjp 
      else
cjp not bwcc
       I0=I
       J0=J
       K0=K
       L0=L
       IF(I.GT.NOA)I0=I-NOA
       IF(J.GT.NOA)J0=J-NOA
       IF(K.GT.NOB)K0=K-NOB
       IF(L.GT.NOB)L0=L-NOB
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
       IF(J.LE.NOA.AND.L.LE.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=9
       ELSEIF(I.GT.NOA.AND.K.GT.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=1
       ELSEIF(L.LE.NOB.AND.I.GT.NOA)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=5
       ELSEIF(J.LE.NOA.AND.K.GT.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=6
       ELSEIF(J.LE.NOA.AND.L.GT.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=8
       ELSEIF(L.LE.NOB.AND.J.GT.NOA)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=7
       ELSEIF(J.GT.NOA.AND.K.GT.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=2
       ELSEIF(L.GT.NOB.AND.I.GT.NOA)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=3
       ELSEIF(J.GT.NOA.AND.L.GT.NOB)THEN
        IBUF(IK)=IPCK(I0,K0,J0,L0)
        IWHERE(IK)=4
       ENDIF
cjp
      endif
cjp
cjp end bwcc change
C
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
cjp special treatment of hf2 deleting for mrbwcc
cjp also take into account parallelization
      IF(IFLAGS(84).NE.2.and. (.not.isbwcc.or.iref+1.gt.nref)) THEN 
cjp
       CLOSE(UNIT=25,STATUS='DELETE')
      ELSE
       CLOSE(UNIT=25,STATUS='KEEP')
      ENDIF
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
cjp
cjp debug info
cjp
      if(isbwcc.and.bwgossip) then
      write(6,*) ' gmoiab-fock A contrib, reference no. ',iref
      do i=1,norb
      write(6,888)(fockcontr(findex(i,j),1),j=1,i)
888   format(10f12.7)
      enddo
      write(6,*) ' gmoiab-fock A contrib, reference no. ',iref
      do i=1,norb
      write(6,888)(fockcontr(findex(i,j),2),j=1,i)
      enddo
      endif

      RETURN
999   WRITE(LUOUT,1400)
1400  FORMAT(T3,'@GMOIAB-F, Unexpected end-of-file on HF2.')
      CALL ERREX
      END
