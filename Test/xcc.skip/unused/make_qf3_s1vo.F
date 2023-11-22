
c This routine builds the 1-particle intermediate needed to compute the
c Qf energy from triples.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(*)  : the scratch array for doubles
c    integer dim_dscr : the size of the dscr array
c    integer ispin1   : the spin case to build

#include "disttype.h"
#include "listnum.h"
#include "xcc.h"

c#define _DEBUG_MAKE_QF3_S1VO

      subroutine make_qf3_s1vo(dcore,dscr,dim_dscr,ispin1)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*)
      integer dim_dscr, ispin1

c INTERNAL VARIABLES
      integer ispin2, irrep, itmp, start_scr
      integer irp_ket
      integer start_t2p1, start_w, start_s1
      integer a, max_a
      integer i, max_i
      integer k
      integer len_bra, len_ket1, len_ket2
      integer off_a1i1
      integer off_i1ltj1, off_i1j2, off_i1a1, off_i2a1
      integer off_a1ltb1, off_a1b2, off_a2b1
      integer dist_i1ltj1, dist_i1j2
      integer dist_i1a1, dist_i1a2
      integer dist_i2a1
      integer dist_a1ltb1, dist_a1b2
      integer list_i1j1k1a1, list_i1j2k1a2, list_a1b1c1i1, list_a1b2c1i2

c PARAMETERS
      double precision zero, one, minus
      parameter (zero=0.0d0, one=1.0d0, minus=-1.0d0)

c EXTERNAL FUNCTIONS
      double precision ddot

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert ispin1 is 1 or 2
      if ((ispin1.lt.1).or.(2.lt.ispin1)) then
         write(*,*) '@MAKE_QF3_S1VO: Assertion failed.'
         write(*,*) '   ispin1 = ',ispin1
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

c   o set handles for these dependents
      if (ispin1.eq.1) then
         ispin2        = 2
         off_i1ltj1    = h_Offsets_IltJ
         off_a1ltb1    = h_Offsets_AltB
         off_i1a1      = h_Offsets_IA
         off_i1j2      = h_Offsets_Ij
         off_i2a1      = h_Offsets_iA
         off_a1i1      = h_Offsets_AI
         off_a1b2      = h_Offsets_Ab
         off_a2b1      = h_Offsets_aB
         list_i1j1k1a1 = h_HH_IJKA_IltJ_KA
         list_i1j2k1a2 = h_GS_IjKa_Ij_Ka
         list_a1b1c1i1 = h_HH_ABCI_AltB_CI
         list_a1b2c1i2 = h_GS_AbCi_Ab_Ci
         dist_i1ltj1   = h_DistType_IltJ
         dist_i1a1     = h_DistType_IA
         dist_i1a2     = h_DistType_Ia
         dist_i2a1     = h_DistType_Ai
         dist_a1ltb1   = h_DistType_AltB
      else
         ispin2        = 1
         off_i1ltj1    = h_Offsets_iltj
         off_a1ltb1    = h_Offsets_altb
         off_i1a1      = h_Offsets_ia
         off_i1j2      = h_Offsets_iJ
         off_i2a1      = h_Offsets_Ia
         off_a1i1      = h_Offsets_ai
         off_a1b2      = h_Offsets_aB
         off_a2b1      = h_Offsets_Ab
         list_i1j1k1a1 = h_HH_ijka_iltj_ka
         list_i1j2k1a2 = h_GS_IjAk_Ij_Ak
         list_a1b1c1i1 = h_HH_abci_altb_ci
         list_a1b2c1i2 = h_GS_AbIc_Ab_Ic
         dist_i1ltj1   = h_DistType_iltj
         dist_i1a1     = h_DistType_ia
         dist_i1a2     = h_DistType_Ai
         dist_i2a1     = h_DistType_Ia
         dist_a1ltb1   = h_DistType_altb
      end if
      dist_i1j2 = h_DistType_Ij
      dist_a1b2 = h_DistType_Ab

c ----------------------------------------------------------------------

c   o initialize the S1vo intermediates
      do i = 0, irpdpd(1,dist_i1a1)-1
         dcore(ndx_s1vo(1,ispin1)+i) = zero
      end do

      do irp_ket = 1, nirrep

c ----------------------------------------------------------------------

c ABCI integrals

c   o alpha:   <Ef|Am> * T2(1)(Im,Ef)
c   o beta:    <Ef|Ma> * T2(1)(Mi,Ef)
         len_bra  = irpdpd(irp_ket,dist_a1b2)
         len_ket1 = irpdpd(irp_ket,dist_i2a1)
         len_ket2 = irpdpd(irp_ket,dist_i1j2)
      if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0)) then

c      o <Ef|Am>; Ef, Am at dscr(1)
         itmp = len_bra * len_ket1
         if (itmp.gt.dim_dscr) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr,1,len_ket1,2,irp_ket,list_a1b2c1i2)

c      o S1(I,A); A, I at dcore(ndx_s1vo(1,ispin1))
         if (ispin1.eq.1) then
            call kmm_pq_rq(
     &              len_bra,vrt(1,ispin1),pop(1,ispin1),pop(1,ispin2),
     &              one,dscr,                      irp_ket,
     &                  dcore(ndx_t2p1(irp_ket,3)),irp_ket,
     &                  dcore(ndx_s1vo(1,ispin1)), off_pq(1,off_a1i1)
     &           )
         else
            call kmm_qp_qr(
     &              len_bra,vrt(1,ispin1),pop(1,ispin1),pop(1,ispin2),
     &              one,dscr,irp_ket,off_pq(1,off_i2a1),
     &                  dcore(ndx_t2p1(irp_ket,3)),irp_ket,
     &                  dcore(ndx_s1vo(1,ispin1))
     &           )
         end if

c     end if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0))
      end if

c   o - <EF||AM> * T2(1)(MI,EF)
         len_bra  = irpdpd(irp_ket,dist_a1ltb1)
         len_ket1 = irpdpd(irp_ket,dist_i1a1)
         len_ket2 = irpdpd(irp_ket,dist_i1ltj1)
      if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0)) then

c      o <EF||AM>; E<F, AM at dscr(1)
         if (iuhf.eq.0) then
c         o antisymmetrize the integrals in memory at dscr(1)
            call antisym_pp_pltp_in_situ(
     &              dscr,irpdpd(irp_ket,dist_a1b2),
     &              vrt(1,1),len_ket1,irp_ket,off_pq(1,off_a1b2)
     &           )
         else
            itmp = len_bra * len_ket1
            if (itmp.gt.dim_dscr) then
               write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dscr,1,len_ket1,2,irp_ket,list_a1b1c1i1)
         end if

c      o S1(I,A); A, I at dcore(ndx_s1vo(1,ispin1))
         call kmm_pq_qltq(
     &           len_bra,vrt(1,ispin1),pop(1,ispin1),
     &           minus,dscr,                           irp_ket,
     &                    off_pq(1,off_a1i1),
     &                 dcore(ndx_t2p1(irp_ket,ispin1)),irp_ket,
     &                    off_pltp(1,off_i1ltj1),
     &                 dcore(ndx_s1vo(1,ispin1))
     &        )

c     end if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0))
      end if

c ----------------------------------------------------------------------

c IJKA integrals

c   o alpha: - T2(1)(Mn,Ae) * <Mn|Ie>
c   o beta:  - T2(1)(Mn,Ea) * <Mn|Ei>
         len_bra  = irpdpd(irp_ket,dist_i1j2)
         len_ket1 = irpdpd(irp_ket,dist_a1b2)
         len_ket2 = irpdpd(irp_ket,dist_i1a2)
      if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0)) then

c      o <Mn|Ie>; Mn, Ie at dscr(1)
         itmp = len_bra * len_ket2
         if (itmp.gt.dim_dscr) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr,1,len_ket2,2,irp_ket,list_i1j2k1a2)
         start_scr = 1 + itmp

c      o T2(1)(Mn,Ae); Mn, Ae at dscr(start_t2p1)
         start_t2p1 = start_scr
         itmp = len_bra * len_ket1
         if (itmp.gt.1+dim_dscr-start_t2p1) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         do k = 0, itmp-1
            dscr(start_t2p1+k) = dcore(ndx_t2p1(irp_ket,3)+k)
         end do
         if ((len_bra.gt.1).and.(len_ket1.gt.1)) then
            call dmat_xpose(dscr(start_t2p1),len_ket1,len_bra)
         end if
         start_scr = start_scr + itmp

c      o S1(I,A); A, I at dcore(ndx_s1vo(1,ispin1))
         if (ispin1.eq.1) then
            call kmm_pq_rq(
     &              len_bra,vrt(1,ispin1),pop(1,ispin1),vrt(1,ispin2),
     &              minus,dscr(start_t2p1),         irp_ket,
     &                    dscr,                     irp_ket,
     &                    dcore(ndx_s1vo(1,ispin1)),off_pq(1,off_a1i1)
     &           )
         else
            call kmm_qp_qr(
     &              len_bra,vrt(1,ispin1),pop(1,ispin1),vrt(1,ispin2),
     &              minus,dscr(start_t2p1),irp_ket,off_pq(1,off_a2b1),
     &                    dscr,            irp_ket,
     &                    dcore(ndx_s1vo(1,ispin1))
     &           )
         end if

c     end if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0))
      end if

c   o   <MN||IE> * T2(1)(MN,EA)
         len_bra  = irpdpd(irp_ket,dist_i1ltj1)
         len_ket1 = irpdpd(irp_ket,dist_i1a1)
         len_ket2 = irpdpd(irp_ket,dist_a1ltb1)
      if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0)) then

         start_scr = 1

c      o <MN||IE>; M<N, IE at dscr(start_w)
         start_w = start_scr
         itmp = len_bra * len_ket1
         if (itmp.gt.1+dim_dscr-start_w) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(start_w),1,len_ket1,2,irp_ket,list_i1j1k1a1)
         start_scr = start_scr + itmp

c      o T2(1)(MN,EA); M<N, E<A at dscr(start_t2p1)
         start_t2p1 = start_scr
         itmp = len_bra * len_ket2
         if (itmp.gt.1+dim_dscr-start_t2p1) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         do k = 0, itmp-1
            dscr(start_t2p1+k) = dcore(ndx_t2p1(irp_ket,ispin1)+k)
         end do
         if ((len_bra.gt.1).and.(len_ket2.gt.1)) then
            call dmat_xpose(dscr(start_t2p1),len_ket2,len_bra)
         end if
         start_scr = start_scr + itmp

c      o S1(I,A); I, A at dscr(start_s1)
         start_s1 = start_scr
         itmp = irpdpd(1,dist_i1a1)
         if (itmp.gt.1+dim_dscr-start_s1) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         do k = 0, itmp-1
            dscr(start_s1+k) = zero
         end do
         call kmm_pq_qltq(
     &           len_bra,pop(1,ispin1),vrt(1,ispin1),
     &           one,dscr(start_w),   irp_ket,off_pq(1,off_i1a1),
     &               dscr(start_t2p1),irp_ket,off_pltp(1,off_a1ltb1),
     &               dscr(start_s1)
     &        )

c      o S1(I,A); A, I += S1(I,A); I, A (use start_s1 as a pointer)
         do irrep = 1, nirrep
            max_i = pop(irrep,ispin1)
            max_a = vrt(irrep,ispin1)
            do i = 0, max_i-1
               do a = 0, max_a-1
                  dcore(ndx_s1vo(irrep,ispin1)+i*max_a+a) =
     &            dcore(ndx_s1vo(irrep,ispin1)+i*max_a+a)
     &          + dscr(start_s1+i+a*max_i)
               end do
            end do
            start_s1 = start_s1 + ( max_i * max_a )
         end do

c     end if ((len_bra.ne.0).and.(len_ket1.ne.0).and.(len_ket2.ne.0))
      end if

c ----------------------------------------------------------------------

c     end do irp_ket = 1, nirrep
      end do

cc   o take the dot product as a check
c      itmp = irpdpd(1,dist_i1a1)
c      dscr(1) = ddot(itmp,dcore(ndx_s1vo(1,ispin1)),1,
c     &                    dcore(ndx_s1vo(1,ispin1)),1)
c      write(*,*) 'S1*S1 in spin ',ispin1,' for Qf is ',dscr(1)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

