
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
      integer ispin2, irrep, itmp, start_scr, irp_ket
      integer start_t2p1, add_t2p1, start_w, ndx_w
      integer irp_a, a, max_a, num_ae
      integer irp_i, i, max_i, num_ie
      integer m, k
      integer max_mn, max_ie, max_ae
      integer irp_e, e, max_e
      logical do_irp_e
      integer dist_i1ltj1, dist_i1a1, dist_a1ltb1
      integer list_i1j1k1a1

c PARAMETERS
      double precision zero, one
      parameter (zero=0.0d0, one=1.0d0)

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _DEBUG
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
#endif /* _DEBUG */

c   o set handles for these dependents
      if (ispin1.eq.1) then
         ispin2        = 2
c         off_i1ltj1    = h_Offsets_IltJ
c         off_a1ltb1    = h_Offsets_AltB
c         off_i1j2      = h_Offsets_Ij
c         off_i2j1      = h_Offsets_iJ
c         off_a1b2      = h_Offsets_Ab
         list_i1j1k1a1 = h_HH_IJKA_IltJ_KA
c         list_i1j2k1a2 = h_GS_IjKa_Ij_Ka
c         list_i2j1k2a1 = h_GS_IjAk_Ij_Ak
c         list_a1b1c1i1 = h_HH_ABCI_AltB_CI
c         list_a1b2c1i2 = h_GS_AbCi_Ab_Ci
c         list_a2b1c2i1 = h_GS_AbIc_Ab_Ic
         dist_i1ltj1   = h_DistType_IltJ
         dist_i1a1     = h_DistType_IA
c         dist_i1a2     = h_DistType_Ia
c         dist_i2a1     = h_DistType_Ai
         dist_a1ltb1   = h_DistType_AltB
      else
         ispin2        = 1
c         off_i1ltj1    = h_Offsets_iltj
c         off_a1ltb1    = h_Offsets_altb
c         off_i1j2      = h_Offsets_iJ
c         off_i2j1      = h_Offsets_Ij
c         off_a1b2      = h_Offsets_aB
         list_i1j1k1a1 = h_HH_ijka_iltj_ka
c         list_i1j2k1a2 = h_GS_IjAk_Ij_Ak
c         list_i2j1k2a1 = h_GS_IjKa_Ij_Ka
c         list_a1b1c1i1 = h_HH_abci_altb_ci
c         list_a1b2c1i2 = h_GS_AbIc_Ab_Ic
c         list_a2b1c2i1 = h_GS_AbCi_Ab_Ci
         dist_i1ltj1   = h_DistType_iltj
         dist_i1a1     = h_DistType_ia
c         dist_i1a2     = h_DistType_Ai
c         dist_i2a1     = h_DistType_Ia
         dist_a1ltb1   = h_DistType_altb
      end if
c      dist_i1j2 = h_DistType_Ij
c      dist_a1b2 = h_DistType_Ab

c   o initialize the S1 intermediates
      itmp = irpdpd(1,h_DistType_AI)
      do i = 0, itmp-1
         dcore(ndx_s1vo(1,ispin1)+i) = zero
      end do

c ----------------------------------------------------------------------

c   o <MN||IE>*T2(1)(MN,AE)
      do irp_ket = 1, nirrep
         max_mn = irpdpd(irp_ket,dist_i1ltj1)
         max_ie = irpdpd(irp_ket,dist_i1a1)
         max_ae = irpdpd(irp_ket,dist_a1ltb1)
      if ((max_ae.ne.0).and.(max_ie.ne.0).and.(max_mn.ne.0)) then

         start_scr = 1
c      o T2(1)(MN,AE); M<N, E<A at dscr(start_t2p1)
         start_t2p1 = start_scr
         itmp = max_ae * max_mn
         if (itmp.gt.1+dim_dscr-start_t2p1) then
            write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         do k = 0, itmp-1
            dscr(start_t2p1+k) = dcore(ndx_t2p1(irp_ket,ispin1)+k)
         end do
         if ((max_ae.gt.1).and.(max_mn.gt.1)) then
            call dmat_xpose(dscr(start_t2p1),max_ae,max_mn)
         end if
         start_scr = start_scr + itmp
c      o loop over irp_e and make S1(I,A); I, A but write out A, I
         if (irp_ket.gt.1) then
            add_t2p1 = start_t2p1
            start_w  = start_scr
            ndx_w    = 1
            do irp_e = 1, nirrep
               irp_i = dirprd(irp_e,irp_ket)
               if (irp_i.lt.irp_e) then
                  irp_a = irp_i
                  max_i = pop(irp_i,ispin1)
                  max_a = vrt(irp_a,ispin1)
                  max_e = vrt(irp_e,ispin1)
                  num_ie = max_i * max_e
                  num_ae = max_a * max_e
                  do_irp_e = ((num_ie.ne.0).and.(num_ae.ne.0))
               else
                  num_ie = 0
                  num_ae = 0
                  do_irp_e = .false.
               end if
            if (do_irp_e) then
c            o <MN||IE>; M<N, IE at dscr(start_w)
               itmp = num_ie * max_mn
               if (itmp.gt.1+dim_dscr-start_w) then
                  write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
                  call xcc_errex
               end if
               start_scr = start_w + itmp
               if (num_ie.gt.1+dim_dscr-start_scr) then
                  write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dscr(start_w),ndx_w,num_ie,2,irp_ket,
     &                     list_i1j1k1a1)
               do i = 0, num_ie-1
                  dscr(start_scr+i) = zero
               end do
               call omm_pq_rq_ddot(
     &                  max_mn,max_i,max_a,max_e,
     &                  one,dscr(start_w),  max_i,
     &                      dscr(add_t2p1), max_a,
     &                      dscr(start_scr),max_i
     &              )
               do i = 0, max_i-1
                  do a = 0, max_a-1
      write(*,*) dcore(ndx_s1vo(irp_i,ispin1)+i*max_a+a),' += ',
     &           dscr(start_scr+i+a*max_i)
                     dcore(ndx_s1vo(irp_i,ispin1)+i*max_a+a) =
     &               dcore(ndx_s1vo(irp_i,ispin1)+i*max_a+a)
     &             + dscr(start_scr+i+a*max_i)
                  end do
               end do
c           end if (do_irp_e)
            end if
               add_t2p1 = add_t2p1 + num_ae
               ndx_w    = ndx_w    + num_ie
c           end do irp_e = 1, nirrep
            end do
c        else if (irp_ket.eq.0) then
         else
            add_t2p1 = start_t2p1
            start_w  = start_scr
            ndx_w    = 1
            do irp_e = 1, nirrep
               max_i = pop(irp_e,ispin1)
               max_e = vrt(irp_e,ispin1)
               num_ie = max_i * max_e
            if ((max_i.ne.0).and.(max_e.gt.1)) then
c            o <MN||IE>; M<N, IE at dscr(start_w)
               itmp = num_ie * max_mn
               if (itmp.gt.1+dim_dscr-start_w) then
                  write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
                  call xcc_errex
               end if
               start_scr = start_w + itmp
               if (num_ie.gt.1+dim_dscr-start_scr) then
                  write(*,*) '@MAKE_QF3_S1VO: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dscr(start_w),ndx_w,num_ie,2,irp_ket,
     &                     list_i1j1k1a1)
               do i = 0, num_ie-1
                  dscr(start_scr+i) = zero
               end do
               call omm_pq_qltq_ddot(
     &                  max_mn,max_i,max_e,
     &                  one,dscr(start_w),  max_i,0,max_i-1,
     &                      dscr(add_t2p1),
     &                      dscr(start_scr),max_i
     &              )
               do i = 0, max_i-1
                  do a = 0, max_e-1
      write(*,*) dcore(ndx_s1vo(irp_e,ispin1)+i*max_e+a),' += ',
     &           dscr(start_scr+i+a*max_i)
                     dcore(ndx_s1vo(irp_e,ispin1)+i*max_e+a) =
     &               dcore(ndx_s1vo(irp_e,ispin1)+i*max_e+a)
     &             + dscr(start_scr+i+a*max_i)
                  end do
               end do
c           end if ((max_i.ne.0).and.(max_e.gt.1))
            end if
               add_t2p1 = add_t2p1 + rshift(max_e*(max_e-1),1)
               ndx_w    = ndx_w    + num_ie
c           end do irp_e = 1, nirrep
            end do
c        end if (irp_ket.gt.1)
         end if

c     end if ((max_ae.ne.0).and.(max_ie.ne.0).and.(max_mn.ne.0))
      end if
c     end do irp_ket = 1, nirrep
      end do

c ----------------------------------------------------------------------

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

