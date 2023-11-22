#define _MINE
c#define _DO_TRANSPOSE

c This routine allocates space at dcore(freecore) and forms the tau 2-particle
c intermediate:
c       AB     AB    A B    B A
c    tau   := T   + T T  - T T
c       IJ     IJ    I J    I J

c INPUT
c double dcore(*) : address space for ft2w.com

#include "disttype.h"
#include "parallel_aces.h" /* for _PA_MAX_PROCS */

#define _DO_T2AA
#define _DO_T2BB
#define _DO_T2AB
#define _PRINT_AMPLITUDES

      subroutine make_tau2(dcore)
      implicit none

#include "mpif.h"

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irrep, nRows, nCols, nBatch, iOff, i, iErr
      integer off_t1(8,2)
      integer recvcounts(_PA_MAX_PROCS), displs(_PA_MAX_PROCS)

c COMMON BLOCKS
#include "ft2w.com"
#include "method.com" /* for iuhf */
#include "sym.com" /* for pop(,) and vrt(,) */
#include "syminf.com" /* for nirrep */
#include "sympop.com" /* for irpdpd(,) */
#include "parallel_aces.com" /* for parallel info */

c ----------------------------------------------------------------------

c   o assign space for each tau2 array and copy T2 batches into it

c   o tau2(A<B,I<J)
      do irrep = 1, nirrep
         ndx_tau2(irrep,1) = freecore
         i =   irpdpd(irrep,h_DistType_AltB)
     &       * irpdpd(irrep,h_DistType_IltJ)
         freecore = freecore + i
         if (freecore.gt.dcoresize) then
            print *, '@MAKE_TAU2: Not enough memory.'
            call aces_exit(1)
         end if
         call dcopy(i,dcore(ndx_t2  (irrep,1)),1,
     &                dcore(ndx_tau2(irrep,1)),1)
      end do

c   o tau2(a<b,i<j)
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_tau2(irrep,2) = ndx_tau2(irrep,1)
         end do
      else
      do irrep = 1, nirrep
         ndx_tau2(irrep,2) = freecore
         i =   irpdpd(irrep,h_DistType_altb)
     &       * irpdpd(irrep,h_DistType_iltj)
         freecore = freecore + i
         if (freecore.gt.dcoresize) then
            print *, '@MAKE_TAU2: Not enough memory.'
            call aces_exit(1)
         end if
         call dcopy(i,dcore(ndx_t2  (irrep,2)),1,
     &                dcore(ndx_tau2(irrep,2)),1)
      end do
      end if

c   o tau2(Ab,Ij)
      do irrep = 1, nirrep
         ndx_tau2(irrep,3) = freecore
         i =   irpdpd(irrep,h_DistType_Ab)
     &       * irpdpd(irrep,h_DistType_Ij)
         freecore = freecore + i
         if (freecore.gt.dcoresize) then
            print *, '@MAKE_TAU2: Not enough memory.'
            call aces_exit(1)
         end if
         call dcopy(i,dcore(ndx_t2  (irrep,3)),1,
     &                dcore(ndx_tau2(irrep,3)),1)
      end do

c ----------------------------------------------------------------------

c   o create the antisymmetrized products of T1*T1 added into tau2

#ifdef _DO_TRANSPOSE
      call transp_dist(dcore(ndx_t1(1,1)),vrt(1,1),pop(1,1),1,1,.false.)
      call transp_dist(dcore(ndx_t1(1,2)),vrt(1,2),pop(1,2),1,1,.false.)
      do irrep = 1, nirrep
         call dmat_xpose(dcore(ndx_tau2(irrep,1)),
     &                   irpdpd(irrep,h_DistType_AltB),
     &                   irpdpd(irrep,h_DistType_IltJ))
         if (iuhf.eq.1) then
         call dmat_xpose(dcore(ndx_tau2(irrep,2)),
     &                   irpdpd(irrep,h_DistType_altb),
     &                   irpdpd(irrep,h_DistType_iltj))
         end if
         call dmat_xpose(dcore(ndx_tau2(irrep,3)),
     &                   irpdpd(irrep,h_DistType_Ab),
     &                   irpdpd(irrep,h_DistType_Ij))
      end do
#endif

c   o hack in the T1 offsets
      off_t1(1,1) = 0
      off_t1(1,2) = 0
      if (nirrep.gt.1) then
      do irrep = 1, nirrep-1
         off_t1(irrep+1,1) = off_t1(irrep,1) + vrt(irrep,1)*pop(irrep,1)
         off_t1(irrep+1,2) = off_t1(irrep,2) + vrt(irrep,2)*pop(irrep,2)
      end do
      end if
      if (irank.eq.0) then
         print *, 'off_t1(*,1) = ',(off_t1(i,1),i=1,nirrep)
         print *, 'off_t1(*,2) = ',(off_t1(i,2),i=1,nirrep)
      end if

#ifdef _DO_T2AA
c   o tau2(A<B,I<J) += T1(A,I)*T1(B,J) - T1(B,I)*T1(A,J)
      do irrep = 1, nirrep
#ifdef _MINE
#ifdef _DO_TRANSPOSE
         nRows = irpdpd(irrep,h_DistType_IltJ)
         nCols = irpdpd(irrep,h_DistType_AltB)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call antisym_s2pack_s1pq_sqr(
     &                dcore(ndx_tau2(irrep,1)+nRows*iOff),
     &                     irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,1)),
     &                     1,pop(1,1),vrt(1,1),off_t1(1,1)
     &        )
#else
         nRows = irpdpd(irrep,h_DistType_AltB)
         nCols = irpdpd(irrep,h_DistType_IltJ)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call antisym_s2pack_s1pq_sqr(
     &                dcore(ndx_tau2(irrep,1)+nRows*iOff),
     &                     irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,1)),
     &                     1,vrt(1,1),pop(1,1),off_t1(1,1)
     &        )
#endif
         do i = 1, nprocs
            call paces_batch_stat(i-1,nprocs,0,nCols,
     &                            displs(i),recvcounts(i))
            recvcounts(i) = nRows*recvcounts(i)
            displs(i)     = nRows*displs(i)
         end do
         call MPI_Allgatherv(
     &            dcore(ndx_tau2(irrep,1)+nRows*iOff),nRows*nBatch,
     &                 MPI_DOUBLE_PRECISION,
     &            dcore(ndx_tau2(irrep,1)),recvcounts,displs,
     &                 MPI_DOUBLE_PRECISION,
     &            MPI_COMM_WORLD,iErr)
         if (iErr.ne.0) then
            print *, '@MAKE_TAU2: unable to gather tau2(A<B,I<J)'
            call aces_exit(iErr)
         end if
#else
         call ftau(dcore(ndx_tau2(irrep,1)),
     &             dcore(ndx_t1(1,1)),
     &             dcore(ndx_t1(1,1)),
     &             irpdpd(irrep,h_DistType_AltB),
     &             irpdpd(irrep,h_DistType_IltJ),
     &             pop(1,1),pop(1,1),
     &             vrt(1,1),vrt(1,1),
     &             irrep,1,1.d0)
#endif
      end do
#endif /* _DO_T2AA */

#ifdef _DO_T2BB
c   o tau2(a<b,i<j) += T1(a,i)*T1(b,j) - T1(b,i)*T1(a,j)
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
#ifdef _MINE
#ifdef _DO_TRANSPOSE
         nRows = irpdpd(irrep,h_DistType_iltj)
         nCols = irpdpd(irrep,h_DistType_altb)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call antisym_s2pack_s1pq_sqr(
     &                dcore(ndx_tau2(irrep,2)+nRows*iOff),
     &                     irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,2)),
     &                     1,pop(1,2),vrt(1,2),off_t1(1,2)
     &        )
#else
         nRows = irpdpd(irrep,h_DistType_altb)
         nCols = irpdpd(irrep,h_DistType_iltj)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call antisym_s2pack_s1pq_sqr(
     &                dcore(ndx_tau2(irrep,2)+nRows*iOff),
     &                     irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,2)),
     &                     1,vrt(1,2),pop(1,2),off_t1(1,2)
     &        )
#endif
         do i = 1, nprocs
            call paces_batch_stat(i-1,nprocs,0,nCols,
     &                            displs(i),recvcounts(i))
            recvcounts(i) = nRows*recvcounts(i)
            displs(i)     = nRows*displs(i)
         end do
         call MPI_Allgatherv(
     &            dcore(ndx_tau2(irrep,2)+nRows*iOff),nRows*nBatch,
     &                 MPI_DOUBLE_PRECISION,
     &            dcore(ndx_tau2(irrep,2)),recvcounts,displs,
     &                 MPI_DOUBLE_PRECISION,
     &            MPI_COMM_WORLD,iErr)
         if (iErr.ne.0) then
            print *, '@MAKE_TAU2: unable to gather tau2(a<b,i<j)'
            call aces_exit(iErr)
         end if
#else
         call ftau(dcore(ndx_tau2(irrep,2)),
     &             dcore(ndx_t1(1,2)),
     &             dcore(ndx_t1(1,2)),
     &             irpdpd(irrep,h_DistType_altb),
     &             irpdpd(irrep,h_DistType_iltj),
     &             pop(1,2),pop(1,2),
     &             vrt(1,2),vrt(1,2),
     &             irrep,2,1.d0)
#endif
      end do
      end if
#endif /* _DO_T2BB */

#ifdef _DO_T2AB
c   o tau2(Ab,Ij) += T1(A,I)*T1(b,j)
      do irrep = 1, nirrep
#ifdef _MINE
#ifdef _DO_TRANSPOSE
         nRows = irpdpd(irrep,h_DistType_Ij)
         nCols = irpdpd(irrep,h_DistType_Ab)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call s2pqrs_pe_s1pr_s1qs(
     &                dcore(ndx_tau2(irrep,3)+nRows*iOff),
     &                     irrep,irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,1)),
     &                     1,pop(1,1),vrt(1,1),off_t1(1,1),
     &                dcore(ndx_t1(1,2)),
     &                     1,pop(1,2),vrt(1,2),off_t1(1,2)
     &        )
#else
         nRows = irpdpd(irrep,h_DistType_Ab)
         nCols = irpdpd(irrep,h_DistType_Ij)
         call paces_batch_stat(irank,nprocs,0,nCols,iOff,nBatch)
         call s2pqrs_pe_s1pr_s1qs(
     &                dcore(ndx_tau2(irrep,3)+nRows*iOff),
     &                     irrep,irrep,1+iOff,nBatch,
     &           1.d0,dcore(ndx_t1(1,1)),
     &                     1,vrt(1,1),pop(1,1),off_t1(1,1),
     &                dcore(ndx_t1(1,2)),
     &                     1,vrt(1,2),pop(1,2),off_t1(1,2)
     &        )
#endif
         do i = 1, nprocs
            call paces_batch_stat(i-1,nprocs,0,nCols,
     &                            displs(i),recvcounts(i))
            recvcounts(i) = nRows*recvcounts(i)
            displs(i)     = nRows*displs(i)
         end do
         call MPI_Allgatherv(
     &            dcore(ndx_tau2(irrep,3)+nRows*iOff),nRows*nBatch,
     &                 MPI_DOUBLE_PRECISION,
     &            dcore(ndx_tau2(irrep,3)),recvcounts,displs,
     &                 MPI_DOUBLE_PRECISION,
     &            MPI_COMM_WORLD,iErr)
         if (iErr.ne.0) then
            print *, '@MAKE_TAU2: unable to gather tau2(Ab,Ij)'
            call aces_exit(iErr)
         end if
#else
         call ftau(dcore(ndx_tau2(irrep,3)),
     &             dcore(ndx_t1(1,1)),
     &             dcore(ndx_t1(1,2)),
     &             irpdpd(irrep,h_DistType_Ab),
     &             irpdpd(irrep,h_DistType_Ij),
     &             pop(1,1),pop(1,2),
     &             vrt(1,1),vrt(1,2),
     &             irrep,3,1.d0)
#endif
      end do
#endif /* _DO_T2AB */

#ifdef _DO_TRANSPOSE
      call transp_dist(dcore(ndx_t1(1,1)),pop(1,1),vrt(1,1),1,1,.false.)
      call transp_dist(dcore(ndx_t1(1,2)),pop(1,2),vrt(1,2),1,1,.false.)
      do irrep = 1, nirrep
         call dmat_xpose(dcore(ndx_tau2(irrep,1)),
     &                   irpdpd(irrep,h_DistType_IltJ),
     &                   irpdpd(irrep,h_DistType_AltB))
         if (iuhf.eq.1) then
         call dmat_xpose(dcore(ndx_tau2(irrep,2)),
     &                   irpdpd(irrep,h_DistType_iltj),
     &                   irpdpd(irrep,h_DistType_altb))
         end if
         call dmat_xpose(dcore(ndx_tau2(irrep,3)),
     &                   irpdpd(irrep,h_DistType_Ij),
     &                   irpdpd(irrep,h_DistType_Ab))
      end do
#endif

c ----------------------------------------------------------------------

#ifdef _PRINT_AMPLITUDES
      if (irank.eq.0) then
      do irrep = 1, nirrep
         nRows = irpdpd(irrep,h_DistType_AltB)
         nCols = irpdpd(irrep,h_DistType_IltJ)
         print *, 'There are ',nRows*nCols,
     &            ' tau2AA amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_tau2(irrep,1)),nRows,nCols,nRows)
      end do
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
         nRows = irpdpd(irrep,h_DistType_altb)
         nCols = irpdpd(irrep,h_DistType_iltj)
         print *, 'There are ',nRows*nCols,
     &            ' tau2BB amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_tau2(irrep,2)),nRows,nCols,nRows)
      end do
      end if
      do irrep = 1, nirrep
         nRows = irpdpd(irrep,h_DistType_Ab)
         nCols = irpdpd(irrep,h_DistType_Ij)
         print *, 'There are ',nRows*nCols,
     &            ' tau2AB amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_tau2(irrep,3)),nRows,nCols,nRows)
      end do
c     end if (irank.eq.0)
      end if
#endif /* _PRINT_AMPLITUDES */

c ----------------------------------------------------------------------

      return
c     end subroutine make_tau2
      end

