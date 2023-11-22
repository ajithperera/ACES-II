      subroutine shku(u,u2,fnodrop,scr,maxcor,uhf)
c
c this routine simply reduces the size of u such that 
c
c	kept		dropped
c	    no-d 	    d
c       vrt u		vrt u2
c
c===============variable declarations and common blocks=======================
      implicit none
C     common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
C     Input variables
      integer fnodrop(8,2),uhf,maxcor
C     Output variables
      double precision u2(*)
C     Input/Output variables
      double precision u(*)
C     Pre-allocated local variables
      double precision scr(maxcor)
C     Local variables
      integer i000,i010,iu,iu2,spin,irrep,virt,vrtdrop,drop,ioff
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i000=1
      iu=1
      iu2=1
      do 10 spin=1,uhf+1
         i010=i000+nd2(spin)
         if (i010.gt.maxcor) call insmem('shku',i010,maxcor)
         call dcopy(nd2(spin),u(1+(spin-1)*nd2(1)),1,scr(i000),1)
         ioff=1
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            drop=fnodrop(irrep,spin)
            vrtdrop=virt-drop
            call dcopy(virt*vrtdrop,scr(ioff),1,u(iu),1)
            ioff=ioff+virt*vrtdrop
            iu=iu+virt*vrtdrop
            call dcopy(virt*drop,scr(ioff),1,u2(iu2),1)
            ioff=ioff+virt*drop
            iu2=iu2+virt*drop
 20      continue
 10   continue
      return
      end
