
c This routine drives the calculation.

c INPUT
c double dcore(n) : address space for ft2w.com

#if defined(_DEBUG)
#   define _INIT_DCORE
#endif

      subroutine ft2w(dcore,n)
      implicit none

c ARGUMENT LIST
      integer n
      double precision dcore(n)

c INTERNAL VARIABLES
      integer i

c COMMON BLOCKS
#include "ft2w.com"
#include "machsp.com" /* for ifltln */
#include "parallel_aces.com" /* for irank */
c#include "sym.com"
c#include "method.com"

c ----------------------------------------------------------------------

#ifdef _ASSERT
      i = 0
c   o assert n > 0
      if (n.lt.1) then
         print *, '@FT2W: Assertion failed.'
         print *, '   n = ',n
         i = 1
      end if
      if (i.ne.0) call aces_exit(i)
#endif

c ----------------------------------------------------------------------

c   o initialize dcore, the pointer to free memory, and the overlap
      dcoresize = n
#ifdef _INIT_DCORE
      call c_memset(dcore,255,dcoresize*ifltln)
#endif
      freecore = 1

cc   o fill in the standard tables of offsets
c      call offset_ind_pltp(pop(1,1),off_pltp(1,h_Offsets_IltJ))
c      call offset_ind_pltp(pop(1,2),off_pltp(1,h_Offsets_iltj))
c      call offset_ind_pltp(vrt(1,1),off_pltp(1,h_Offsets_AltB))
c      call offset_ind_pltp(vrt(1,2),off_pltp(1,h_Offsets_altb))
c      call offset_ind_pq(pop(1,1),pop(1,1),off_pq(1,h_Offsets_IJ))
c      call offset_ind_pq(pop(1,1),pop(1,2),off_pq(1,h_Offsets_Ij))
c      call offset_ind_pq(pop(1,1),vrt(1,1),off_pq(1,h_Offsets_IA))
c      call offset_ind_pq(pop(1,1),vrt(1,2),off_pq(1,h_Offsets_Ia))
c      call offset_ind_pq(pop(1,2),pop(1,1),off_pq(1,h_Offsets_iJ))
c      call offset_ind_pq(pop(1,2),pop(1,2),off_pq(1,h_Offsets_ij))
c      call offset_ind_pq(pop(1,2),vrt(1,1),off_pq(1,h_Offsets_iA))
c      call offset_ind_pq(pop(1,2),vrt(1,2),off_pq(1,h_Offsets_ia))
c      call offset_ind_pq(vrt(1,1),pop(1,1),off_pq(1,h_Offsets_AI))
c      call offset_ind_pq(vrt(1,1),pop(1,2),off_pq(1,h_Offsets_Ai))
c      call offset_ind_pq(vrt(1,1),vrt(1,1),off_pq(1,h_Offsets_AB))
c      call offset_ind_pq(vrt(1,1),vrt(1,2),off_pq(1,h_Offsets_Ab))
c      call offset_ind_pq(vrt(1,2),pop(1,1),off_pq(1,h_Offsets_aI))
c      call offset_ind_pq(vrt(1,2),pop(1,2),off_pq(1,h_Offsets_ai))
c      call offset_ind_pq(vrt(1,2),vrt(1,1),off_pq(1,h_Offsets_aB))
c      call offset_ind_pq(vrt(1,2),vrt(1,2),off_pq(1,h_Offsets_ab))

c NOTE:
c    All of the get_* routines load data into dcore and point to them
c with addresses stored in ft2w_com. Depending on whether the data is
c new, freecore will shift to preserve the information.

c   o print chemical system info
      if (irank.eq.0) call print_header

c   o load the T1 amplitudes into dcore
      call get_t1(dcore)

c   o load the T2 amplitudes into dcore
      call get_t2(dcore)

c   o make the tau 2-particle intermediate
      call make_tau2(dcore)

c   o verify the freecore pointer is -NaN
      print *, 'dcore(freecore-1) = ', dcore(freecore-1)
      print *, 'dcore(freecore)   = ', dcore(freecore)

c ----------------------------------------------------------------------

      return
c     end subroutine ft2w
      end

