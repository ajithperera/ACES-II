
c This function loads the T amplitudes and calculates <0|adj(Tn)*Tn|0>
c after each call to getlst. This routine DOES NOT calculate all of the
c overlap.

c INPUT
c    double precision dcore(*) : core space for doubles

c OUTPUT
c    integer get_t : the exit status of the routine

#include "disttype.h"
#include "listnum.h"

      function get_t(dcore)
      implicit none

c ARGUMENT LIST
      integer get_t
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irrep, i
      integer listnum, sublist, isize, itmp
      double precision dtmp, dtmp2

c EXTERNAL FUNCTIONS
      double precision ddot

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c   o assume the worst
      get_t=-1

c ----------------------------------------------------------------------

C T1

c dtmp  : the cumulative contribution to the overlap from T1
c dtmp2 : the individual dot products from each spin case

      dtmp=0.0D0

c   o ALPHA spin
      listnum=h_T1_IA_AI
      sublist=1
      do irrep=1,nirrep
         ndx_t1(irrep,1)=freecore
         freecore=freecore+(pop(irrep,1)*vrt(irrep,1))
      end do
      call getlst(dcore(ndx_t1(1,1)),1,1,1,sublist,listnum)
      isize=irpdpd(1,h_DistType_AI)
      dtmp2=ddot(isize,dcore(ndx_t1(1,1)),1,
     &                 dcore(ndx_t1(1,1)),1)
      dtmp=dtmp+dtmp2

c   o BETA spin
      if (iuhf.eq.0) then
         dtmp=dtmp+dtmp2
         do irrep=1,nirrep
            ndx_t1(irrep,2)=ndx_t1(irrep,1)
         end do
      else
         listnum=h_T1_ia_ai
         sublist=2
         do irrep=1,nirrep
            ndx_t1(irrep,2)=freecore
            freecore=freecore+(pop(irrep,2)*vrt(irrep,2))
         end do
         call getlst(dcore(ndx_t1(1,2)),1,1,1,sublist,listnum)
         isize=irpdpd(1,h_DistType_ai)
         dtmp2=ddot(isize,dcore(ndx_t1(1,2)),1,
     &                    dcore(ndx_t1(1,2)),1)
         dtmp=dtmp+dtmp2
      end if

      overlap=overlap+dtmp

c END T1

c ----------------------------------------------------------------------

c T2

c dtmp  : the cumulative contribution to the overlap from T2
c dtmp2 : the individual dot products from each spin case

      dtmp=0.0D0

c   o T2(IJ,AB); A<B, I<J
      listnum=h_T2_IJAB_AltB_IltJ
      do irrep=1,nirrep
         ndx_t2(irrep,1)=freecore
         freecore=   freecore
     &             + (   irpdpd(irrep,h_DistType_AltB)
     &                 * irpdpd(irrep,h_DistType_IltJ) )
      end do
      isize=freecore-ndx_t2(1,1)
      call getall(dcore(ndx_t2(1,1)),isize*iintfp,1,listnum)
      dtmp2=ddot(isize,dcore(ndx_t2(1,1)),1,
     &                 dcore(ndx_t2(1,1)),1)
      dtmp=dtmp+dtmp2

c   o T2(ij,ab); a<b, i<j
      if (iuhf.eq.0) then
         dtmp=dtmp+dtmp2
         do irrep=1,nirrep
            ndx_t2(irrep,2)=ndx_t2(irrep,1)
         end do
      else
         listnum=h_T2_ijab_altb_iltj
         do irrep=1,nirrep
            ndx_t2(irrep,2)=freecore
            freecore=   freecore
     &                + (   irpdpd(irrep,h_DistType_altb)
     &                    * irpdpd(irrep,h_DistType_iltj) )
         end do
         isize=freecore-ndx_t2(1,2)
         call getall(dcore(ndx_t2(1,2)),isize*iintfp,1,listnum)
         dtmp2=ddot(isize,dcore(ndx_t2(1,2)),1,
     &                    dcore(ndx_t2(1,2)),1)
         dtmp=dtmp+dtmp2
      end if

c   o T2(Ij,Ab); A,b,I,j
      listnum=h_T2_IjAb_Ab_Ij
      do irrep=1,nirrep
         ndx_t2(irrep,3)=freecore
         freecore=   freecore
     &             + (   irpdpd(irrep,h_DistType_Ab)
     &                 * irpdpd(irrep,h_DistType_Ij) )
      end do
      isize=freecore-ndx_t2(1,3)
      call getall(dcore(ndx_t2(1,3)),isize*iintfp,1,listnum)
      dtmp2=ddot(isize,dcore(ndx_t2(1,3)),1,
     &                 dcore(ndx_t2(1,3)),1)
      dtmp=dtmp+dtmp2

      overlap=overlap+dtmp

c END T2

c ----------------------------------------------------------------------

      get_t=0

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

