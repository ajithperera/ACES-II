      subroutine redvrt(occ,list,fnodrop,scr,uhf,dropmo)
c determine the number of virtuals out of each irrep that 
c will be dropped this is based on the lowest occupation numbers
c
c occ		natural orbital occupation numbers
c list		reordering vector for occupation numbers 
c fnodrop	number of orbitals out of each irrep to be dropped
c scr		integer scratch space dimension 16
c
c input can come from several different places
c *ACES2 namelist
c 	FNO_KEEP=	percent of the virtual space to keep
c				iflags2(h_IFLAGS2_fno_keep)
c	NATURAL=	natural virtual orbitals will be generated 
c				non-hf reference
c *fno_gen namelist
c	keep_sym	number of orbitals from each symmetry to keep
c                       if uhf alpha/beta
c	keep_thresh	double precision natural orbital occupation threshold 
c			for virtual space to keep	
c================variable declarations and common blocks==================  
      implicit none
#include "flags.h"
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iflags2(500)
      common/flags2/iflags2
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      logical rohf,qrhf,semi,natorb
      common/fnoparam/rohf,qrhf,semi,natorb
C     Input variables
      integer uhf,dropmo
      double precision occ(*)
C     Output Variables
      integer fnodrop(8,2)
C     Pre-allocated Local Variables
      integer scr(16),list(*)
C     Local Variables
      integer err,correct,flagcor,reducedrop,i,ikeep,ii,iln,ioff,idrop,
     &   spin,iocc,i2,irrep,i1,vrto
      double precision k_thresh,tolerance,zilch,samethresh,threshold
      character*5 csptype(2)
      character*70 ck_sym
      data csptype /'alpha','beta '/
      parameter(tolerance=1.0d-20)
      parameter(zilch=0.d0)
      parameter(samethresh=1.0d-12)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      iocc=1
      call izero(fnodrop,16)
c     open the ZMAT file and search for the fno_gen namelist
      call nl_init('fno_gen',err,.true.)
      call nl_real('keep_thresh',0.d0,k_thresh)
      call nl_str('keep_sym',' ',ck_sym)
      Write(6,*) "The FNO keep string", iflags2(h_IFLAGS2_fno_keep)
      if (iflags2(h_IFLAGS2_fno_keep).lt.10) 
     &   write(6,1) iflags2(h_IFLAGS2_fno_keep)
      if ((iflags2(h_IFLAGS2_fno_keep).ge.10).and.
     &   (iflags2(h_IFLAGS2_fno_keep).lt.100)) 
     &   write(6,2) iflags2(h_IFLAGS2_fno_keep)
      if (iflags2(h_IFLAGS2_fno_keep).eq.100) 
     &   write(6,3) iflags2(h_IFLAGS2_fno_keep)
 1    format(12x,'FNO_KEEP',5x,'integer',13x,i1,1x,'percent',7x,
     &   '*ACES2 namelist')
 2    format(12x,'FNO_KEEP',5x,'integer',12x,i2,1x,'percent',7x,
     &   '*ACES2 namelist')
 3    format(12x,'FNO_KEEP',5x,'integer',11x,i3,1x,'percent',7x,
     &   '*ACES2 namelist')
      if(natorb)then
         write(6,4) 
 4       format(13x,'NATURAL',5x,'integer',19x,' ON',7x,
     &      '*ACES2 namelist')
      else
         write(6,5) 
 5       format(13x,'NATURAL',5x,'integer',19x,'OFF',7x,
     &      '*ACES2 namelist')
      endif
      write(6,6)
 6    format(76('='))
c     loop over alpha and beta spin cases if non rhf 
      do 10 spin=1,uhf+1
         idrop=0
c sort occupation numbers
         do 20 i=1,nvrto(spin)
            list(i)=i
 20      continue
         call piksr2(nvrto(spin),occ(iocc),list)
               
c     check different possible inputs and deal with them 
         if(abs(k_thresh-zilch).gt.tolerance)then
            do 30 i=iocc,iocc+nvrto(spin)-1
               if(k_thresh.gt.occ(i)) idrop=idrop+1
 30         continue
         else if(ck_sym(1:1).ne.' ')then
            if (spin.eq.1) then
               i=1
               do 50 while(i.ne.0)
                  i=index(ck_sym,'-')
                  ck_sym(i:i)=' '
 50            continue
c     alpha spin
               i1=1 
               do 55 irrep=1,nirrep
                  i2=index(ck_sym,' ')
                  iln=i2-i1
                  if(iln.eq.1)then
                     read(ck_sym(i1:i2),'(I1)') scr(irrep)
                  else if(iln.eq.2)then
                     read(ck_sym(i1:i2),'(I2)') scr(irrep)
                  else if(iln.eq.3)then
                     read(ck_sym(i1:i2),'(I3)') scr(irrep)
                  else
                     write(6,*)' problem with the keep_sym string'
                     write(6,*)' iln ',iln
                     call errex
                  endif
                  fnodrop(irrep,1)=vrt(irrep,1)-scr(irrep)
                  ck_sym(i2:i2)='x'
                  i1=i2+1
 55            continue
            else
               do 56 irrep=1,nirrep
                  fnodrop(irrep,2)=fnodrop(irrep,1)
 56            continue
            endif
            goto 900 
         else if(iflags2(h_IFLAGS2_fno_keep).ne.0)then
            vrto=nvrto(spin)
            if (rohf) vrto=nvrto(1)
            idrop=nint(vrto*((100-iflags2(h_IFLAGS2_fno_keep))/100.0d0))
     &         +dropmo
            threshold=occ(iocc+idrop)
         else
            write(6,*)' no input was made to specify ',
     &         'virtual space partition'
            call errex 
            
         endif
c     
c if we've got here idrop should be defined and is used for orbital selection
c
         ikeep=nvrto(spin)-idrop 
c     loop over orbitals to be dropped 
         do 200 i=1,idrop
c     figure out which irrep this orbital belongs to
            ioff=0
            do 210 irrep=1,nirrep
               if((list(i).gt.ioff).and.
     &            (list(i).le.(ioff+vrt(irrep,spin))))then
                  fnodrop(irrep,spin)=fnodrop(irrep,spin)+1
                  goto 200
               endif
               ioff=ioff+vrt(irrep,spin)
 210        continue
 200     continue 
c the rest is just printing
 900     continue

         reducedrop=0
C     Checking for near-degeneracies among occupation numbers
         do ii=idrop,1,-1
            if(abs(occ(iocc+ii-1)-occ(iocc+idrop)).lt.samethresh) then
               reducedrop=reducedrop+1
               ioff=0
               do irrep=1,nirrep
                  if ((list(ii).gt.ioff).and.(list(ii).le.ioff+
     &               vrt(irrep,spin)))
     &               fnodrop(irrep,spin)=fnodrop(irrep,spin)-1
                  ioff=ioff+vrt(irrep,spin)
               end do
            endif
         end do
         idrop=idrop-reducedrop
      
         write(6,*)' Virtual orbitals with the lowest natural '
         write(6,*)' orbital occupation numbers will be dropped.'
         write(6,2001) occ(iocc+idrop)
         write(6,*)' Breakdown for ',csptype(spin),' spin follows: '
         write(6,1001)
 2001    format('  Occupation Number cutoff:',t30,f13.10)
 1001    format(t5,'irrep',t17,'virtuals',t30,
     &      'virtuals to be dropped')
         write(6,1000)
         do 220 irrep=1,nirrep
            write(6,1003) irrep,vrt(irrep,spin),fnodrop(irrep,spin)
 1003       format(t5,i3,t18,i3,t43,i3)
 220     continue
         write(6,1000)
c/////////////////////////////////////////////////////////////////////////
c     check to see is fnodrop is valid
         do 230 irrep=1,nirrep
            if(fnodrop(irrep,spin).lt.0) then 
               write(6,*)' problem with irrep ',irrep
               write(6,*)' number of virtuals present ',vrt(irrep,spin)
               write(6,*)' number requested to be kept',
     &            vrt(irrep,spin)-fnodrop(irrep,spin)
               call errex
            endif
 230     continue
         iocc=iocc+nvrto(1)
 10   continue
      flagcor=0
      if (uhf.ne.0) then
         do 260 irrep=1,nirrep
            correct=fnodrop(irrep,1)-fnodrop(irrep,2)
            if (correct.gt.0) then
               fnodrop(irrep,1)=fnodrop(irrep,1)-correct
               flagcor=1
            elseif (correct.lt.0) then
               fnodrop(irrep,2)=fnodrop(irrep,2)+correct
               flagcor=1
            endif
 260     continue
         if (flagcor.eq.1) then
            write(6,*) 'Imbalance in dropped orbitals found.  In at'
     &         //' least one irrep, differing'
            write(6,*) 'numbers of orbitals were dropped in the alpha'
     &         //' and beta spins.'
            write(6,*) 'Fewer orbitals will be dropped to correct this'
     &         //' imbalance.'
            write(6,*) 'Imbalance corrected results:'
            write(6,*)' Virtual orbitals with the lowest natural '
            write(6,*)' orbital occupation numbers will be dropped.'
            write(6,*)' Breakdown for ',csptype(1),' spin follows: '
            write(6,1001)
            write(6,1000)
            do 270 irrep=1,nirrep
               write(6,1003) irrep,vrt(irrep,1),fnodrop(irrep,1)
 270        continue
            write(6,1000)
c     take care of beta spin
            write(6,*)' Virtual orbitals with the lowest natural '
            write(6,*)' orbital occupation numbers will be dropped.'
            write(6,*)' Breakdown for ',csptype(2),' spin follows: '
            write(6,1001)
            write(6,1000)
            do 280 irrep=1,nirrep
               write(6,1003) irrep,vrt(irrep,2),fnodrop(irrep,2)
 280        continue
            write(6,1000)
         endif
      endif
      
      return
 1000 format(T3,70('-'))
 2000 format((10(2x,I3)))
      end 

