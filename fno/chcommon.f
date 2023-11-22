      subroutine chcommon(nbfirr,droppop,scr,maxcor,uhf)
c
c this subroutine is entered only when there are dropped occupied orbitals 
c it resets all variable (i hope) from their dropmo dimensions to their 
c full dimensions
c
c NOCCORB(spin) total number of occupied orbitals by spin
c SYMPOPOA(irrep) number of occuppied orbitals per irrep
c SYMPOPVA(irrep) number of virtual orbitals per irrep
c OCCUPYA(irrep) OCCUPYB(irrep) occupation by irrep
c NUMBASIR(irrep) total number of functions by irrep
c IRREPALP(nbas) IRREPBET(nbas) irrep for individual orbitals  
c
C COMMON BLOCK /SYM/ is initialized by initpop
C
C  DESCRIPTION OF /SYM/
C
C         POP    ...... NUMBER OF OCCUPIED ORBITALS WITHIN EACH IRREP
C         VRT    ...... NUMBER OF VIRTUAL ORBITALS WITHIN EACH IRREP
C         NTAA   ...... LENGTH OF T1(A,I) ALPHA
C         NTBB   ...... LENGTH OF T1(A,I) BETA
C         NF1AA  ...... LENGTH OF F(M,I) ALPHA ( NOTE ALL Fs ARE NOT SYMMETRIC)
C         NF1BB  ...... LENGTH OF F(M,I) BETA
C         NF2AA  ...... LENGTH OF F(A,E) ALPHA
C         NF2BB  ...... LENGTH OF F(A,E) BETA
c
c REORDERA/REORDERB are changed back to its orginal scf to post scf form  
c============variable declarations and common blocks=====================
c      implicit double precision (a-h,o-z)
      implicit none
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
C     Input Variables
      integer uhf,maxcor
C     Pre-allocated Local Variables
      integer nbfirr(8),droppop(8),scr(maxcor)
C     Local Variables
      integer itwo,ivrt,nmo,spin,nbas,irrep,icnt,iocc,i,j
      character*8 coccupy(2)
      data coccupy /'OCCUPYA ','OCCUPYB '/
      parameter(itwo=2)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c trap cases where virtuals have been dropped
c
      call getrec(20,'JOBARC','NDROPVRT',nirrep,scr)
      do 10 irrep=1,nirrep
       if(scr(irrep).ne.0) then
        write(6,*)' a virtual orbital was specified with the ',
     &   'dropmo keyword '
        write(6,*)' use the keyword fno_keep to ',
     &   'reduce the virtual space'
        call errex
       endif
10    continue
      call getrec(20,'JOBARC','NDROPPOP',nirrep,droppop)
c
c take care of pop(irrep,spin) and nocco(spin)
c
      do 20 irrep=1,nirrep
        pop(irrep,1)=pop(irrep,1)+droppop(irrep)
        pop(irrep,2)=pop(irrep,2)+droppop(irrep)
        nocco(1)=nocco(1)+droppop(irrep)
        nocco(2)=nocco(2)+droppop(irrep)
20    continue
      call putrec(20,'JOBARC','NOCCORB ',itwo,nocco)
      call putrec(20,'JOBARC','SYMPOPOA',nirrep,pop(1,1))
      if(uhf.ne.0)then 
        call putrec(20,'JOBARC','SYMPOPOB',nirrep,pop(1,2))
      endif
c
c take care of occupya and if necessary occupyb 
c
      do 30 spin=1,uhf+1
        call getrec(20,'JOBARC',coccupy(spin),nirrep,scr)
        do 35 irrep=1,nirrep
          scr(irrep)=scr(irrep)+droppop(irrep)
35      continue
        call putrec(20,'JOBARC',coccupy(spin),nirrep,scr) 
30    continue
c
c take care of numbasir
c take care of irrepalp and irrepbet which are initially defined in
c the subroutine dmpjob of program vscf
c
      do 40 irrep=1,nirrep
        nbfirr(irrep)=pop(irrep,1)+vrt(irrep,1)
40    continue
      call putrec(20,'JOBARC','NUMBASIR',nirrep,nbfirr)

      nbas=nocco(1)+nvrto(1)
      ICNT=0
      CALL IZERO(SCR,NBAS)
      DO 50 I=1,NIRREP
        DO 60 J=1,NBFIRR(I)
          ICNT=ICNT+1
          SCR(ICNT)=I
60      CONTINUE
50    CONTINUE
      CALL PUTREC(20,'JOBARC','IRREPALP',NBAS,SCR)
      IF(UHF.EQ.1)CALL PUTREC(20,'JOBARC','IRREPBET',NBAS,SCR)
c
c now take care of all other variable within sym that are dependent on
c pop and vrt 
c
      nt(1)=0
      nt(2)=0
      nd1(1)=0
      nd1(2)=0
      DO 90 IRREP=1,NIRREP         
         nt(1)=nt(1)+pop(irrep,1)*vrt(irrep,1)
         nt(2)=nt(2)+pop(irrep,2)*vrt(irrep,2)
         nd1(1)=nd1(1)+pop(irrep,1)**2
         nd1(2)=nd1(2)+pop(irrep,2)**2
90    CONTINUE
c
c take care of reordera reorderb
c these are written by vtran and seem to be altered in cases of dropped
c core
c
       i=1
       iocc=1
       ivrt=nocco(1)+1
       do 100 irrep=1,nirrep
         do 101 j=iocc,iocc+pop(irrep,1)-1
           scr(i)=j
           i=i+1
101      continue
         do 102 j=ivrt,ivrt+vrt(irrep,1)-1
           scr(i)=j
           i=i+1
102      continue
         iocc=iocc+pop(irrep,1)
         ivrt=ivrt+vrt(irrep,1)
100    continue
       nmo=nocco(1)+nvrto(1)
       call putrec(20,'JOBARC','REORDERA',nmo,scr)  
       if(uhf.ne.0)call putrec(20,'JOBARC','REORDERB',nmo,scr)  
       return
1001   format((15(2x,I3)))
       end
