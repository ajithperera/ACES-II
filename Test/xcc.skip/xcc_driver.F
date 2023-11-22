
c This routine drives the calculation.

c INPUT
c    double precision dcore(n) : address space for xcc.com

#ifdef _DEBUG
#define _INIT_DCORE
#endif /* _DEBUG */

      subroutine xcc_driver(dcore,n)
      implicit none

c ARGUMENT LIST
      integer n
      double precision dcore(n)

c INTERNAL VARIABLES
c      integer i

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "machsp.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o initialize dcore, the pointer to free memory, and the overlap
      dcoresize = n
#ifdef _INIT_DCORE
      call init_dcore(dcore,dcoresize*ifltln)
#endif
      freecore = 1
      overlap = 1.0d0

c   o fill in the standard tables of offsets
      call offset_ind_pltp(pop(1,1),off_pltp(1,h_Offsets_IltJ))
      call offset_ind_pltp(pop(1,2),off_pltp(1,h_Offsets_iltj))
      call offset_ind_pltp(vrt(1,1),off_pltp(1,h_Offsets_AltB))
      call offset_ind_pltp(vrt(1,2),off_pltp(1,h_Offsets_altb))
      call offset_ind_pq(pop(1,1),pop(1,1),off_pq(1,h_Offsets_IJ))
      call offset_ind_pq(pop(1,1),pop(1,2),off_pq(1,h_Offsets_Ij))
      call offset_ind_pq(pop(1,1),vrt(1,1),off_pq(1,h_Offsets_IA))
      call offset_ind_pq(pop(1,1),vrt(1,2),off_pq(1,h_Offsets_Ia))
      call offset_ind_pq(pop(1,2),pop(1,1),off_pq(1,h_Offsets_iJ))
      call offset_ind_pq(pop(1,2),pop(1,2),off_pq(1,h_Offsets_ij))
      call offset_ind_pq(pop(1,2),vrt(1,1),off_pq(1,h_Offsets_iA))
      call offset_ind_pq(pop(1,2),vrt(1,2),off_pq(1,h_Offsets_ia))
      call offset_ind_pq(vrt(1,1),pop(1,1),off_pq(1,h_Offsets_AI))
      call offset_ind_pq(vrt(1,1),pop(1,2),off_pq(1,h_Offsets_Ai))
      call offset_ind_pq(vrt(1,1),vrt(1,1),off_pq(1,h_Offsets_AB))
      call offset_ind_pq(vrt(1,1),vrt(1,2),off_pq(1,h_Offsets_Ab))
      call offset_ind_pq(vrt(1,2),pop(1,1),off_pq(1,h_Offsets_aI))
      call offset_ind_pq(vrt(1,2),pop(1,2),off_pq(1,h_Offsets_ai))
      call offset_ind_pq(vrt(1,2),vrt(1,1),off_pq(1,h_Offsets_aB))
      call offset_ind_pq(vrt(1,2),vrt(1,2),off_pq(1,h_Offsets_ab))

c NOTE:
c    All of the get_* routines load data into dcore and point to them
c with addresses stored in xcc_com. Depending on whether the data is
c new, freecore will shift to preserve the information.
c    It is important to track the values held in operator intermediates
c like s1 and s2. The information will change depending on which diagram
c is to be processed next.

c   o load the T1 amplitudes into dcore and compute their overlap
      call get_t1(dcore)

c   o calculate the 4th order overlap term <0|adj(T1)^2*T1^2|0>
      call ovrlp_t1t1_t1t1(dcore,dcore(freecore))

c   o load the T2 amplitudes into dcore and compute their overlap
      call get_t2(dcore)

c   o calculate the 3rd order overlap term(s) <0|adj(T2)*T1^2|0>
      call ovrlp_t2_t1t1(dcore)

c   o load the SCF eigenvalues into dcore
      call get_eval(dcore)

c   o load the T2(1) amplitudes into dcore for CCSD(TQf)
c      call get_t2p1(dcore)

c   o load the intermediates for Qf from triples
c      call get_qf3_s1vo(dcore)

c   o reset the I/O cache since everything is in dcore
      call aces_cache_reset

c   o load the T3 amplitudes into dcore and compute their overlap
      call get_t3(dcore)

c     MORE TO COME...

      print '(A,F15.12,A)', 'The XCC overlap is: ',overlap,' a.u.'
c      call dump_xcc

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

