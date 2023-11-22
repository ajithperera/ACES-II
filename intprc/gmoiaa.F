      SUBROUTINE GMOIAA(IORDER,BUCK,IBUCK,BUF,IBUF,NO,NBKINT,
     &                   ISPIN,NIRREP,NSTO,IWHERE,ILOOKUP,IUHF)
cjp
cjp process integrals for different fermi vacua of Hilbert space BWCC
cjp besides that, accumulate one triangle of fock matrix difference
cjp and fermi vacuum energy
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
cjp
       parameter (nrecord=1024, mask255=255)
cjp data for integral sorting  wrto general fermi vacuum
       integer iquad,indx,iperm,intype
       dimension iquad(4),iperm(4,0:15),intype(0:15)
c     lkji    type permutation to get aces ordering
c                  (must be element of symmetry group of index permutations
c                   of tei in mullikan conv.(1234,2134,1243,2143,3412,4312,3421,4321)),
c                   still is not unique
c     hhhh    1    1234
c     hhhp    2    4312
c     hhph    2    3412
c     hhpp    3    3412
c     hphh    2    2134
c     hphp    4    2143
c     hpph    4    2134
c     hppp    5    4321
c     phhh    2    1234
c     phhp    4    1243
c     phph    4    1234
c     phpp    5    3412     conservatively with original aces
c     pphh    3    1234
c     pphp    5    1243
c     ppph    5    1234
c     pppp    7    1234
      data intype/1,2,2,3,2,4,4,5,2,4,4,5,3,5,5,7/
      data iperm/
     +1,2,3,4,
     +4,3,1,2,
     +3,4,1,2,
     +3,4,1,2,
     +2,1,3,4,
     +2,1,4,3,
     +2,1,3,4,
     +4,3,2,1,
     +1,2,3,4,
     +1,2,4,3,
     +1,2,3,4,
     +3,4,1,2,
     +1,2,3,4,
     +1,2,4,3,
     +1,2,3,4,
     +1,2,3,4/
cjp for fock contribution generation
      integer fperm,fockind,nfperm
      parameter(nfperm=8)
      dimension fperm(4,nfperm)
      data fperm/
     +1,2,3,4,
     +2,1,3,4,
     +1,2,4,3,
     +2,1,4,3,
     +3,4,1,2,
     +4,3,1,2,
     +3,4,2,1,
     +4,3,2,1/
      dimension fockind(4,nfperm)
C
C IGRNBF IS THE LENGTH OF INTEGRAL RECORDS IN THE HF2 FILE WRITTEN BY
C  GRNFNC.  NBKINT IS THE LENGTH OF RECORDS WRITTEN BY GMOIAA.
C
      PARAMETER (IGRNBF =600)
      DIMENSION BUCK(NBKINT,7),IBUCK(NBKINT,7),ILOOKUP(IGRNBF)
      DIMENSION BUF(IGRNBF),IBUF(IGRNBF),IORDER(NSTO),IWHERE(IGRNBF)
      DIMENSION NFULL(7),IU(7),NINBCK(7)
cjp
cjp for suggestion of orbital transformation
      real*8 hhhh,hhhl,hhll,hlhl,hlll,llll
cjp
      real*8 xx
      integer *4 npacked, ijklx
      dimension xx(nrecord), ijklx(nrecord)
cjp
cjp for fock contributions
      logical b1,b2,b3,b4,b5,b6
      real*8 ruhffock
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      COMMON /ABCD/ ABCDTYPE
cjp
#include "bwcc.com"

c      real*8 enermo
c      dimension enermo(maxorb)
cjp
      DATA SPCASE /'AA  ','BB  ','AB  '/
      DATA SP /'A','B'/
      DATA IU /21,20,19,18,17,17,16/
cjp
      integer findex
      findex(i,j)=(i*(i-1))/2+j

      IPCK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),
     &ISHFT(K,2*IBITWD)),ISHFT(L,3*IBITWD))
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      indp(i)=(i*(i+1))/2
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
cjp
c       iword1(iarg)=iand(ishft(iarg,-24),mask255)
c       iword2(iarg)=iand(ishft(iarg,-16),mask255)
c       iword3(iarg)=iand(ishft(iarg, -8),mask255)
c       iword4(iarg)=iand(ishft(iarg,  0),mask255)

cjp faktor for integral contribution to fock matrix
      if(iuhf.eq.0) then
       ruhffock=2.0
      else
       ruhffock=1.0
      endif

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
cjp
cjp interface to masik's program
       if(masik) then
       open(unit=33,file='ftn03',status='unknown', form='unformatted')
       npacked=0
       endif
cjp zero out fock contributions
       if(isbwcc) then
       do j=1,nsto*(nsto+1)/2
       fockcontr(j,ispin)=0d0
       enddo
       endif
cjp zero out contributions for scf energy
       do j=1,nsto
       fock2elcontr(j,ispin)=0d0
       enddo
C      
cjp reading loop, we also accumulate contributions to fock matrix differencies
cjp
CSSS      Write(6,*) "Reading Ints"
 2    READ(22,END=999)BUF,IBUF,NUT
CSS      Write(6,*) "What is going on"
CSS      Write(6,*) NUT
      IF(NUT.EQ.-1) GO TO 103
      DO 11 IK=1,NUT
       I=IUPKI(IBUF(IK))
       J=IUPKJ(IBUF(IK))
       K=IUPKK(IBUF(IK))
       L=IUPKL(IBUF(IK))
cjp
cjpccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjp fock contributions
cjp
       if(isbwcc) then
cjp due to vacuum energy calc. must be done also for iref=1
cjp check for all possibilities of equal at least two indices-only such integrals can contribute
       b1=i.eq.j
       b2= i.eq.k
       b3= i.eq.l 
       b4=j.eq.k
       b5=j.eq.l
       b6=k.eq.l
       if(b1.or.b2.or.b3.or.b4.or.b5.or.b6) then
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
cjp as coulomb or exchange integral appropriatelly for closed and open shells
cjp these integrals here contribute only to either alpha or beta fock matrix
c
cjp here all such contributions are considered, addfock then itself selects
cjp which of them belong to the difference due to different occupation in
cjp reference of SCF and iref
c coulomb
       if(myi.ge.myj .and. myk.eq.myl)
     &     call addfock(myi,myj,myk,ruhffock*buf(ik),ispin,ispin,iuhf)
c exchange
       if(myi.ge.myl .and. myk.eq.myj) 
     &     call addfock(myi,myl,myk,-buf(ik),ispin,ispin,iuhf)
       enddo
       endif
       endif
cjp end fock contribution
cjpccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       IF(NIRREP.NE.1)THEN
        I=IUPKI(IBUF(IK))
        J=IUPKJ(IBUF(IK))
        K=IUPKK(IBUF(IK))
        L=IUPKL(IBUF(IK))
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
       ELSE
        I=IUPKI(IBUF(IK))
        J=IUPKJ(IBUF(IK))
        K=IUPKK(IBUF(IK))
        L=IUPKL(IBUF(IK))
       ENDIF
C
cjp write records to the masik's file
       if(masik) then
       npacked=npacked+1
       xx(npacked)=buf(ik)
C insure canonical ordering of indices for masik
        if(j.gt.i)then
         it=i
         i=j
         j=it
        endif
        if(l.gt.k)then
         it=k
         k=l
         l=it
        endif
        if(indp(i)+j.lt.indp(k)+l)then
         it=i
         i=k
         k=it
         it=j
         j=l
         l=it
        endif
       ijklx(npacked)=
     &    ior(ior(ior(ishft(i,24),ishft(j,16)),ishft(k,8)),l)
       if(npacked.eq.nrecord) then
         write(33)xx,ijklx,npacked
         npacked=0
       endif
       endif
cjp end masik

       if(.not.isbwcc) then
cjp this algorithm is valid only for the ground state reference

       I0=I
       J0=J
       K0=K
       L0=L
       IF(I.GT.NO)I0=I-NO
       IF(J.GT.NO)J0=J-NO
       IF(K.GT.NO)K0=K-NO
       IF(L.GT.NO)L0=L-NO
C
       IBUF(IK)=IPCK(I0,K0,J0,L0)
C
C  Designations above refer to Dirac indices.
C
       IF(I.GT.NO.AND.K.GT.NO)THEN
        IWHERE(IK)=7 
       ELSEIF(K.GT.NO.AND.J.GT.NO.AND.L.GT.NO)THEN
        IWHERE(IK)=5
       ELSEIF(I.GT.NO.AND.J.GT.NO.AND.L.GT.NO)THEN
        IWHERE(IK)=5
        IBUF(IK)=IPCK(K0,I0,L0,J0)
       ELSEIF(L.LE.NO)THEN
         IWHERE(IK)=1
       ELSEIF(K.LE.NO.AND.J.LE.NO)THEN
         IWHERE(IK)=2
       ELSEIF(K.GT.NO.AND.J.LE.NO)THEN
         IWHERE(IK)=3
       ELSEIF(K.LE.NO.AND.J.GT.NO.AND.I.LE.NO)THEN
         IWHERE(IK)=4
       ENDIF
cjp   
      else
cjp general algorithm for mr-bwcc case, in fact, it does not require canonical ordering at all
c
cjp this algorithm should be applicable to UHF AA, BB cases
      iquad(4)=iocc(i,iref,ispin)
      iquad(3)=iocc(j,iref,ispin)
      iquad(2)=iocc(k,iref,ispin)
      iquad(1)=iocc(l,iref,ispin)
cjp compute the indx ... binary number with digits lkji; p=1, h=0
      indx=0
      do ii=1,4
        indx=2*indx+1-iquad(ii)
      enddo
cjp permute the quadruple
      iquad(4)=iphnum(i,iref,ispin)
      iquad(3)=iphnum(j,iref,ispin)
      iquad(2)=iphnum(k,iref,ispin)
      iquad(1)=iphnum(l,iref,ispin)
      i0=iquad(iperm(4,indx))
      j0=iquad(iperm(3,indx))
      k0=iquad(iperm(2,indx))
      l0=iquad(iperm(1,indx))
      iwhere(ik)=intype(indx)
cjp switch to dirac at the end and store it
      ibuf(ik)=ipck(i0,k0,j0,l0)
      endif
cjp
cjp for bwcc write out integrals with indices in active space just for information
      if(isbwcc.and.isactive(i,ispin).and.isactive(j,ispin).and.
     &   isactive(k,ispin).and.isactive(l,ispin)) then
        if(bwgossip.or.iflags(1).ge.10)
     &  write(6,*)'@GMOIAA-I spin ',ispin,' integral of active orbs: (',
     &          i,j,'|',k,l,') = <',i,k,'|',j,l,'> = ',buf(ik)
       if(iuhf.eq.0 .and. nactive(1).eq.2) then
cjp    store integrals
       if(i+j+k+l.eq.4*numactive(1,1)) hhhh=buf(ik)
       if(i+j+k+l.eq.4*numactive(2,1)) llll=buf(ik)
       if(i+j.eq.2*numactive(1,1).and.k+l.eq.2*numactive(2,1).or.
     &  i+j.eq.2*numactive(2,1).and.k+l.eq.2*numactive(1,1))
     &     hhll=buf(ik)
       if(i+j.eq.k+l.and.i+j.eq.numactive(1,1)+numactive(2,1))
     &      hlhl=buf(ik)
       if(i+j+k+l.eq.3*numactive(1,1)+numactive(2,1)) hhhl=buf(ik)
       if(i+j+k+l.eq.3*numactive(2,1)+numactive(1,1)) hlll=buf(ik)
       endif
      endif
#ifdef _FAIL_COMPILE
cjp debug
      if(bwgossip) then
      if(isbwcc) then
       write(6,*)'<',i,k,'|',j,l,'> = ',buf(ik)
       write(6,*)'gmoiaa: ik, old ijkl, new ijkl, type',ik,':',i,j,k,l
     &       ,':',i0,j0,k0,l0,':',iwhere(ik)
      else
chandle the swap to have strightforward debug output
      if (iwhere(ik).eq.5 .and. I.GT.NO.AND.J.GT.NO.AND.L.GT.NO ) then
         write(6,*)'gmoiaa: ik, old ijkl, new ijkl, type',ik,':',i,j,k,
     &               l,':',k0,l0,i0,j0,':',iwhere(ik)
      else
         write(6,*)'gmoiaa: ik, old ijkl, new ijkl, type',ik,':',i,j,k,
     &               l,':',i0,j0,k0,l0,':',iwhere(ik)
      endif
      endif
      endif
#endif /* _FAIL_COMPILE */
 11   CONTINUE
C
C DEAL WITH STUFF NOW
cjp this should be already independent on the method of sorting
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
C        
      GOTO 2
C
103   CONTINUE
#ifdef _DEBUG_LVLM	
      Write(6,*)
      Write(6,*) "Writing in GMOIAA"
      Write(*,'(4(1X,F12.7)))'), (BUF(I), I=1, IGRNBF)
#endif
      if(masik) then
      if(npacked.gt.0) write(33)xx,ijklx,npacked
      npacked=0
      write(33)xx,ijklx,npacked
cjp
c nsto is total no. of (not dropped) orbitals 
c      call getrec(20,'JOBARC','SCFEVALA',nsto*iintfp,enermo)
cjp this has to be replaced by diagonal fock elements in general
c      write(33) (enermo(i),i=1,nsto)
c
c      call getrec(20,'JOBARC','SCFENEG ',iintfp,escf)
c      enuc=0d0
c      ehf=0d0
c      write(33) enuc,ehf,escf
cjp 
c     close(33)
       write(6,*)' masik data part 1 prepared'
c     stop 'masik data prepared'
      endif
c end reading loop
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
cjp special treatment of hf2 deleting for mrbwcc
cjp also has to take into account parallelization
      IF(IFLAGS(84).NE.2 .and. (.not.isbwcc .or. masik .or.
     &           iref+1 .gt. nref)) THEN
       CLOSE(UNIT=22,STATUS='DELETE')
      ELSE 
       CLOSE(UNIT=22,STATUS='KEEP')
      ENDIF
C
C WRITE SUMMARY.
C
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
cjp
      if(isbwcc.and.bwgossip) then
      write(6,*) ' gmoiaa-fock contrib, reference no. ',
     &      iref,' spin ',ispin
      do i=1,nsto
      write(6,888)(fockcontr(findex(i,j),ispin),j=1,i)
888   format(10f12.7)
      enddo
      endif
c
cjp suggest orbital transformation to minimize exchange of active orbitals
      if(isbwcc .and. iuhf.eq.0 .and. nactive(1).eq.2) then
         call actorbtransf(hhhh,hhhl,hhll,hlhl,hlll,llll)
      endif
      RETURN
999   WRITE(LUOUT,1400)
1400  FORMAT(T3,'@GMOIAA-F, Unexpected end-of-file on HF2.')
      CALL ERREX
      END
