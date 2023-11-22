
c This routine builds <RM{}JK> intermediates used in <adj(T2(1)^2)|W|T3>.
c The actual contractions are W*T2(1) -> S2(+1).

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(dim_dscr) : the scratch array for doubles

c    integer ispin1 : the spin case to handle

c    integer j, irp_j, max_j : j statistics
c    integer k, irp_k, max_k : k statistics
c    integer ioff_jk : the common offset for the JK pair

c    integer irp_rm, irp_jk : the left and right symmetries
c                             NOTE: this is not fully tested for
c                                   cases when irp_rm != irp_jk

c OUTPUT
c    double precision s2rmjk(*) : the destination array
c                                 NOTE: this is initialized to 0.0d0

c ----------------------------------------------------------------------

#include "disttype.h"
#include "listnum.h"

c#define _DEBUG_MAKE_QF3_HH_AIJK
#define _DO_HH_EFRM
#define _DO_GS_eKMn
#define _DO_HH_EKMN
#define _DO_GS_eJMn
#define _DO_HH_EJMN

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine make_qf3_hh_aijk(
     &              dcore, dscr, dim_dscr,
     &              ispin1, s2rmjk,
     &              j, irp_j, max_j, irp_jk, ioff_jk,
     &              k, irp_k, max_k, irp_rm
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr), s2rmjk(*)
      integer ispin1, ispin2
      integer j, irp_j, max_j
      integer k, irp_k, max_k
      integer irp_rm, irp_jk, ioff_jk

c INTERNAL VARIABLES
      integer list_abci, list_ijka, list_i1j2k1a2
      integer dist_altb, dist_iltj, dist_ia, dist_i1j2
      integer off_iltj, off_ia, off_ai, off_altb
      integer e, irp_e, max_e
      integer n, irp_n, max_n, from_n, num_n
      integer    irp_m, max_m
      integer    irp_r, max_r
      integer irp_mn, max_mn
      logical do_irp_n
      integer start_scr
      integer start_t2p1, add_t2p1, inc_t2p1
      integer             add_w,    inc_w,    ldw
      integer             add_s2
      integer itmp, ioff, max_dists
      double precision s2_factor

c PARAMETERS
      double precision one, minus
      parameter (one=1.0d0, minus=-1.0d0)

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
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Assertion failed.'
         write(*,*) '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin1 is 1 or 2
      if ((ispin1.lt.1).or.(2.lt.ispin1)) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Assertion failed.'
         write(*,*) '   ispin1 = ',ispin1
         itmp = -1
      end if
c   o assert j and k are well-behaved
      if ((max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Assertion failed.'
         write(*,*) '   j = ',j,', max_j = ',max_j
         write(*,*) '   k = ',k,', max_k = ',max_k
         itmp = -1
      end if
c   o assert irp_j and irp_k are in [1,nirrep]
      if ((irp_j .lt.1).or.(nirrep.lt.irp_j ).or.
     &    (irp_k .lt.1).or.(nirrep.lt.irp_k ).or.
     &    (irp_jk.lt.1).or.(nirrep.lt.irp_jk).or.
     &    (irp_rm.lt.1).or.(nirrep.lt.irp_rm)    ) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Assertion failed.'
         write(*,*) '   irp_j  = ',irp_j
         write(*,*) '   irp_k  = ',irp_k
         write(*,*) '   irp_jk = ',irp_jk
         write(*,*) '   irp_rm = ',irp_rm
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

c ----------------------------------------------------------------------

c The plan is:
c   - initialize S2 with <EF||RM>*T2(1)(JK,EF)
c   o loop over irreps of Mn/MN
c      ( irp_k -> irp_e )
c      - load all <Mn|Ke> and convert to <eK|Mn>
c      o loop over irreps of n
c         o loop over orbital n
c            - load/transpose T2(1)(nJ,eR)
c            - add in  <eK|Mn>*T2(1)(nJ,eR) for all RM
c      - if RHF, antisymmetrize <eK|Mn> to <EK||MN>
c      - if UHF, load all <MN||KE> and convert to <EK||MN>
c      o loop over irreps of N
c         o loop over orbital N
c            - load/transpose T2(1)(NJ,ER)
c            - add in -<EK||MN>*T2(1)(NJ,ER) for all RM
c      ( irp_j -> irp_e )
c      - load all <Mn|Je> and convert to <eJ|Mn>
c      o loop over irreps of n
c         o loop over orbital n
c            - load/transpose T2(1)(nK,eR)
c            - add in -<eJ|Mn>*T2(1)(nK,eR) for all RM
c      - if RHF, antisymmetrize <eJ|Mn> to <EJ||MN>
c      - if UHF, load all <MN||JE> and convert to <EJ||MN>
c      o loop over irreps of N
c         o loop over orbital N
c            - load/transpose T2(1)(NK,ER)
c            - add in  <EJ||MN>*T2(1)(NK,ER) for all RM

c   o set handles
      if (ispin1.eq.1) then
         ispin2 = 2
         list_i1j2k1a2 = h_GS_IjKa_Ij_Ka
         list_abci = h_HH_ABCI_AltB_CI
         list_ijka = h_HH_IJKA_IltJ_KA
         dist_altb = h_DistType_AltB
         dist_iltj = h_DistType_IltJ
         dist_ia   = h_DistType_IA
         off_ia    = h_Offsets_IA
         off_ai    = h_Offsets_AI
         off_iltj  = h_Offsets_IltJ
         off_altb  = h_Offsets_AltB
      else
         ispin2 = 1
         list_i1j2k1a2 = h_GS_IjAk_Ij_Ak
         list_abci = h_HH_abci_altb_ci
         list_ijka = h_HH_ijka_iltj_ka
         dist_altb = h_DistType_altb
         dist_iltj = h_DistType_iltj
         dist_ia   = h_DistType_ia
         off_ia    = h_Offsets_ia
         off_ai    = h_Offsets_ai
         off_iltj  = h_Offsets_iltj
         off_altb  = h_Offsets_altb
      end if
      dist_i1j2 = h_DistType_Ij

c   o initialize S2 with <EF||RM>*T2(1)(JK,EF)
      n = irpdpd(irp_rm,dist_ia)
      if (n.eq.0) return
      e = irpdpd(irp_rm,dist_altb)
      if (e.ne.0) then
c      o we need at least 2 <EF| distributions
         if (dim_dscr.lt.(e*2)) then
            write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
            call xcc_errex
         end if
c      o load T2(1)(JK,EF)
         ioff = ndx_t2p1(irp_jk,ispin1) + ( e * ioff_jk )
         do itmp = 0, e-1
            dscr(1+itmp) = dcore(ioff+itmp)
         end do
         start_scr = 1 + e
         max_dists = ( dim_dscr / e ) - 1
c      o Load as many <EF||RM> lists as possible,
c        dot with T2(1), and write into S2.
c        NOTE: these *_n variables are NOT OFFSETS
         max_n  = n + 1
         from_n = 1
         do while (from_n.lt.max_n)
            num_n = max_n-from_n
            if (iuhf.eq.0) then
c            o RHF must antisymmetrize the integrals on the fly
               call get_hh_pltp_rhf(
     &                 dscr(start_scr),dim_dscr-e,from_n,num_n,
     &                 irp_rm,h_GS_AbCi_Ab_Ci,
     &                 irpdpd(irp_rm,h_DistType_Ab),
     &                 irpdpd(irp_rm,h_DistType_AltB),
     &                 vrt(1,1),off_pq(1,h_Offsets_Ab)
     &              )
               if (num_n.eq.0) then
                  write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
                  write(*,*) '                   At least ',
     &                       lshift(( irpdpd(irp_rm,h_DistType_Ab)
     &                               -irpdpd(irp_rm,h_DistType_AltB)),
     &                              3),' Bytes are needed.'
                  call xcc_errex
               end if
            else
               num_n = min(num_n,max_dists)
               call getlst(dscr(start_scr),from_n,num_n,
     &                     2,irp_rm,list_abci)
            end if
            itmp = start_scr
            do n = from_n, (from_n-1+num_n)
#ifdef _DO_HH_EFRM
               s2rmjk(n) = ddot(e,dscr(1),1,dscr(itmp),1)
               itmp = itmp + e
#else /* _DO_HH_EFRM */
               s2rmjk(n) = 0.0d0
#endif /* _DO_HH_EFRM */
            end do
            from_n = from_n + num_n
c        end do while (from_n.lt.max_n)
         end do
c     else if (e.eq.0) then
      else
         do itmp = 1, n
            s2rmjk(itmp) = 0.0d0
         end do
      end if

c ----------------------------------------------------------------------

c   o loop over irreps of Mn/MN
      do irp_mn = 1, nirrep

c ----------------------------------------------------------------------

c   ( irp_k -> irp_e )
      irp_e = dirprd(irp_mn,irp_k)

      max_mn = irpdpd(irp_mn,dist_i1j2)
      max_e  = vrt(irp_e,ispin2)
      if ((max_e.ne.0).and.(max_mn.ne.0)) then

      start_scr = 1 + ( max_mn * max_e )
      if (start_scr.gt.dim_dscr) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
         call xcc_errex
      end if

c   o load all <Mn|Ke> and convert to <eK|Mn>
      if (ispin1.eq.1) then
c      o <Mn|Ke>;Mn,Ke
         itmp = nirrep*(irp_e-1) + irp_k
         itmp = 1 + off_pq(itmp,h_Offsets_Ia) + k
         start_scr = 1
         do e = 1, max_e
            call getlst(dscr(start_scr),itmp,1,2,irp_mn,list_i1j2k1a2)
            start_scr = start_scr + max_mn
            itmp = itmp + max_k
         end do
      else
c      o <Nm|Ek>;Nm,Ek
         itmp = nirrep*(irp_k-1) + irp_e
         itmp = 1 + off_pq(itmp,h_Offsets_Ai) + ( max_e * k )
         call getlst(dscr,itmp,max_e,2,irp_mn,list_i1j2k1a2)
         call transp_dist(dscr,pop(1,1),pop(1,2),irp_mn,max_e,.false.)
      end if
      call dmat_xpose(dscr,max_mn,max_e)

#ifdef _DO_GS_eKMn
      start_t2p1 = start_scr
      add_w = 1
c   o loop over irreps of n
      do irp_n = 1, nirrep
         irp_m = dirprd(irp_n,irp_mn)
         irp_r = dirprd(irp_m,irp_rm)
         max_m = pop(irp_m,ispin1)
         max_n = pop(irp_n,ispin2)
         max_r = vrt(irp_r,ispin1)
         do_irp_n = ((max_n.ne.0).and.(max_m.ne.0).and.(max_r.ne.0))
      if (do_irp_n) then

c      o point to the RM block of S2
         itmp = nirrep*(irp_m-1) + irp_r
         add_s2 = 1 + off_pq(itmp,off_ai)

c      o point to the first T2(1)(nJ,*) and measure the increment
         start_scr = start_t2p1 + ( max_e * max_r )
         if (start_scr.gt.(dim_dscr+1)) then
            write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
            call xcc_errex
         end if
         itmp = dirprd(irp_n,irp_j)
         inc_t2p1 = irpdpd(itmp,h_DistType_Ab)
         add_t2p1 = ndx_t2p1(itmp,3)
         if (ispin1.eq.1) then
            itmp = nirrep*(irp_n-1) + irp_j
            itmp = off_pq(itmp,h_Offsets_Ij) + j
            add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
            inc_t2p1 = inc_t2p1 * max_j
            s2_factor = minus
         else
            itmp = nirrep*(irp_j-1) + irp_n
            itmp = off_pq(itmp,h_Offsets_Ij) + ( max_n * j )
            add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
            s2_factor = one
         end if

c      o loop over orbital n (not an offset, here)
         inc_w = max_e * max_m
         do n = 1, max_n

c         o load/transpose/increment T2(1)(nJ,eR)
            call load_t2_virts(
     &              dcore(add_t2p1), dscr(start_t2p1),
     &              irp_e, max_e,             ispin2,
     &              irp_r, max_r, 0, max_r-1, ispin1, h_Offsets_Ab
     &           )
            add_t2p1 = add_t2p1 + inc_t2p1

c         o add in <eK|Mn>*T2(1)(nJ,eR) for all RM
            call B_GEMM('T','N',max_r,max_m,max_e,
     &                  s2_factor,dscr(start_t2p1),max_e,
     &                            dscr(add_w),     max_e,
     &                  one,      s2rmjk(add_s2),  max_r)

c      o increment and close the loop over orbital n
            add_w = add_w + inc_w
c        end do n = 1, max_n
         end do

c   o close the loop over irreps of n
c     else if (.not.do_irp_n) then
      else
         add_w = add_w + ( max_e * max_m * max_n )
c     end if (do_irp_n)
      end if
c     end do irp_n = 1, nirrep
      end do
#endif /* _DO_GS_eKMn */

c     end if ((max_e.ne.0).and.(max_mn.ne.0))
      end if

      max_mn = irpdpd(irp_mn,dist_iltj)
      max_e  = vrt(irp_e,ispin1)
      if ((max_e.ne.0).and.(max_mn.ne.0)) then

      start_scr = 1 + ( max_mn * max_e )
      if (start_scr.gt.dim_dscr) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
         call xcc_errex
      end if

      if (iuhf.eq.0) then
c      o antisymmetrize <eK|Mn> to <EK||MN>
         call antisym_ket_pp_pltp_in_situ(
     &           dscr(1),max_e,
     &           pop(1,1),irp_mn,off_pq(1,h_Offsets_Ij)
     &        )
      else
c      o load all <MN||KE> and convert to <EK||MN>
         itmp = nirrep*(irp_e-1) + irp_k
         itmp = 1 + off_pq(itmp,off_ia) + k
         start_scr = 1
         do e = 1, max_e
            call getlst(dscr(start_scr),itmp,1,2,irp_mn,list_ijka)
            start_scr = start_scr + max_mn
            itmp = itmp + max_k
         end do
         call dmat_xpose(dscr,max_mn,max_e)
c     end if (iuhf.eq.0)
      end if

#ifdef _DO_HH_EKMN
      start_t2p1 = start_scr
c   o loop over irreps of N
      do irp_n = 1, nirrep
         irp_m = dirprd(irp_n,irp_mn)
         irp_r = dirprd(irp_m,irp_rm)
         max_m = pop(irp_m,ispin1)
         max_n = pop(irp_n,ispin1)
         max_r = vrt(irp_r,ispin1)
         do_irp_n = ((max_n.ne.0).and.(max_m.ne.0).and.(max_r.ne.0))
      if (do_irp_n) then

c      o point to the RM block of S2
         itmp = nirrep*(irp_m-1) + irp_r
         add_s2 = 1 + off_pq(itmp,off_ai)

c      o allow for unpacked T2(1)(NJ,*)
         start_scr = start_t2p1 + ( max_e * max_r )
         if (start_scr.gt.(dim_dscr+1)) then
            write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
            call xcc_errex
         end if

c      o four cases: irp_m (==,!=) irp_n, and irp_j (==,!=) irp_n
c        If irp_m and irp_n are the same, then the effects of looping
c        over orbital n are much worse; therefore, that is the highest
c        criterion.
         if (irp_m.ne.irp_n) then

c         o find the first <*||MN> and calculate the metrics
c           (Remember, this term is -<EK||MN>*T2(1)...)
            if (irp_m.lt.irp_n) then
               s2_factor = minus
               inc_w = max_e * max_m
               ldw   = max_e
               itmp  = rshift(irp_n*(irp_n-1),1) + irp_m
               itmp  = off_pltp(itmp,off_iltj)
               add_w = 1 + ( max_e * itmp )
            else
               s2_factor = one
               inc_w = max_e
               ldw   = max_e * max_n
               itmp  = rshift(irp_m*(irp_m-1),1) + irp_n
               itmp  = off_pltp(itmp,off_iltj)
               add_w = 1 + ( max_e * itmp )
            end if

c         o find the first T2(1)(NJ,*), calculate the metrics, and
c           loop over orbital N
            itmp = dirprd(irp_n,irp_j)
            inc_t2p1 = irpdpd(itmp,dist_altb)
            add_t2p1 = ndx_t2p1(itmp,ispin1)
            if (irp_n.ne.irp_j) then
               if (irp_n.lt.irp_j) then
                  itmp = rshift(irp_j*(irp_j-1),1) + irp_n
                  itmp = off_pltp(itmp,off_iltj) + ( max_n * j )
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               else
                  itmp = rshift(irp_n*(irp_n-1),1) + irp_j
                  itmp = off_pltp(itmp,off_iltj) + j
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
                  inc_t2p1 = inc_t2p1 * max_j
                  s2_factor = -s2_factor
               end if
               do n = 1, max_n
                  call load_t2_virts(
     &                    dcore(add_t2p1), dscr(start_t2p1),
     &                    irp_e, max_e,             ispin1,
     &                    irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                 )
                  add_t2p1 = add_t2p1 + inc_t2p1
                  call B_GEMM('T','N',max_r,max_m,max_e,
     &                        s2_factor,dscr(start_t2p1),max_e,
     &                                  dscr(add_w),     ldw,
     &                        one,      s2rmjk(add_s2),  max_r)
                  add_w = add_w + inc_w
               end do
c           else if (irp_n.eq.irp_j) then
            else
c               itmp = rshift(irp_j*(irp_j-1),1) + irp_j
               itmp = rshift(irp_j*(irp_j+1),1)
               itmp = off_pltp(itmp,off_iltj) + rshift(j*(j-1),1)
               add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               if (j.ne.0) then
                  do n = 1, j
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + inc_t2p1
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           s2_factor,dscr(start_t2p1),max_e,
     &                                     dscr(add_w),     ldw,
     &                           one,      s2rmjk(add_s2),  max_r)
                     add_w = add_w + inc_w
                  end do
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * j )
               end if
               add_w = add_w + inc_w
               if ((j+1).ne.max_j) then
                  do n = 1, (max_n-1)-j
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + ( inc_t2p1 * (n+j) )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           -s2_factor,dscr(start_t2p1),max_e,
     &                                      dscr(add_w),     ldw,
     &                           one,       s2rmjk(add_s2),  max_r)
                     add_w = add_w + inc_w
                  end do
               end if
c           end if (irp_n.ne.irp_j)
            end if

c        else if (irp_m.eq.irp_n) then
         else

c         o allow for unpacked <*||MN>
            if ((start_scr+(max_e*max_m)).gt.(dim_dscr+1)) then
               write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
               call xcc_errex
            end if

c         o find the first <*||MN>
            s2_factor = minus
c            itmp  = rshift(irp_n*(irp_n-1),1) + irp_n
            itmp  = rshift(irp_n*(irp_n+1),1)
            itmp  = off_pltp(itmp,off_iltj)
            add_w = 1 + ( max_e * itmp )

c         o find the first T2(1)(NJ,*), calculate the metrics, and
c           loop over orbital N
            itmp = dirprd(irp_n,irp_j)
            inc_t2p1 = irpdpd(itmp,dist_altb)
            add_t2p1 = ndx_t2p1(itmp,ispin1)
            if (irp_n.ne.irp_j) then
               if (irp_n.lt.irp_j) then
                  itmp = rshift(irp_j*(irp_j-1),1) + irp_n
                  itmp = off_pltp(itmp,off_iltj) + ( max_n * j )
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               else
                  itmp = rshift(irp_n*(irp_n-1),1) + irp_j
                  itmp = off_pltp(itmp,off_iltj) + j
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
                  inc_t2p1 = inc_t2p1 * max_j
                  s2_factor = -s2_factor
               end if
               do n = 0, max_n-1
                  call load_t2_virts(
     &                    dcore(add_t2p1), dscr(start_t2p1),
     &                    irp_e, max_e,             ispin1,
     &                    irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                 )
                  add_t2p1 = add_t2p1 + inc_t2p1
                  call extract_s2_dists(
     &                    dscr(add_w), dscr(start_scr), max_e,
     &                    irp_m, max_m, 0, max_m-1, max_m,
     &                    irp_n, max_n, n,
     &                    .true., .false.
     &                 )
                  call B_GEMM('T','N',max_r,max_m,max_e,
     &                        s2_factor,dscr(start_t2p1),max_e,
     &                                  dscr(start_scr), max_e,
     &                        one,      s2rmjk(add_s2),  max_r)
               end do
c           else if (irp_n.eq.irp_j) then
            else
c               itmp = rshift(irp_j*(irp_j-1),1) + irp_j
               itmp = rshift(irp_j*(irp_j+1),1)
               itmp = off_pltp(itmp,off_iltj) + rshift(j*(j-1),1)
               add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               if (j.ne.0) then
                  do n = 0, j-1
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + inc_t2p1
                     call extract_s2_dists(
     &                       dscr(add_w), dscr(start_scr), max_e,
     &                       irp_m, max_m, 0, max_m-1, max_m,
     &                       irp_n, max_n, n,
     &                       .true., .false.
     &                    )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           s2_factor,dscr(start_t2p1),max_e,
     &                                     dscr(start_scr), max_e,
     &                           one,      s2rmjk(add_s2),  max_r)
                  end do
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * j )
               end if
               if ((j+1).ne.max_j) then
                  do n = 1, (max_n-1)-j
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + ( inc_t2p1 * (n+j) )
                     call extract_s2_dists(
     &                       dscr(add_w), dscr(start_scr), max_e,
     &                       irp_m, max_m, 0, max_m-1, max_m,
     &                       irp_n, max_n, j+n,
     &                       .true., .false.
     &                    )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           -s2_factor,dscr(start_t2p1),max_e,
     &                                      dscr(start_scr), max_e,
     &                           one,       s2rmjk(add_s2),  max_r)
                  end do
               end if
c           end if (irp_n.ne.irp_j)
            end if

c        end if (irp_m.ne.irp_n)
         end if

c   o close the loop over irreps of N
c     end if (do_irp_n)
      end if
c     end do irp_n = 1, nirrep
      end do
#endif /* _DO_HH_EKMN */

c     end if ((max_e.ne.0).and.(max_mn.ne.0))
      end if

c ----------------------------------------------------------------------

c   ( irp_j -> irp_e )
      irp_e = dirprd(irp_mn,irp_j)

      max_mn = irpdpd(irp_mn,dist_i1j2)
      max_e  = vrt(irp_e,ispin2)
      if ((max_e.ne.0).and.(max_mn.ne.0)) then

      start_scr = 1 + ( max_mn * max_e )
      if (start_scr.gt.dim_dscr) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
         call xcc_errex
      end if

c   o load all <Mn|Je> and convert to <eJ|Mn>
      if (ispin1.eq.1) then
c      o <Mn|Je>;Mn,Je
         itmp = nirrep*(irp_e-1) + irp_j
         itmp = 1 + off_pq(itmp,h_Offsets_Ia) + j
         start_scr = 1
         do e = 1, max_e
            call getlst(dscr(start_scr),itmp,1,2,irp_mn,list_i1j2k1a2)
            start_scr = start_scr + max_mn
            itmp = itmp + max_j
         end do
      else
c      o <Nm|Ej>;Nm,Ej
         itmp = nirrep*(irp_j-1) + irp_e
         itmp = 1 + off_pq(itmp,h_Offsets_Ai) + ( max_e * j )
         call getlst(dscr,itmp,max_e,2,irp_mn,list_i1j2k1a2)
         call transp_dist(dscr,pop(1,1),pop(1,2),irp_mn,max_e,.false.)
      end if
      call dmat_xpose(dscr,max_mn,max_e)

#ifdef _DO_GS_eJMn
      start_t2p1 = start_scr
      add_w = 1
c   o loop over irreps of n
      do irp_n = 1, nirrep
         irp_m = dirprd(irp_n,irp_mn)
         irp_r = dirprd(irp_m,irp_rm)
         max_m = pop(irp_m,ispin1)
         max_n = pop(irp_n,ispin2)
         max_r = vrt(irp_r,ispin1)
         do_irp_n = ((max_n.ne.0).and.(max_m.ne.0).and.(max_r.ne.0))
      if (do_irp_n) then

c      o point to the RM block of S2
         itmp = nirrep*(irp_m-1) + irp_r
         add_s2 = 1 + off_pq(itmp,off_ai)

c      o point to the first T2(1)(nK,*) and measure the increment
         start_scr = start_t2p1 + ( max_e * max_r )
         if (start_scr.gt.(dim_dscr+1)) then
            write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
            call xcc_errex
         end if
         itmp = dirprd(irp_n,irp_k)
         inc_t2p1 = irpdpd(itmp,h_DistType_Ab)
         add_t2p1 = ndx_t2p1(itmp,3)
         if (ispin1.eq.1) then
            itmp = nirrep*(irp_n-1) + irp_k
            itmp = off_pq(itmp,h_Offsets_Ij) + k
            add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
            inc_t2p1 = inc_t2p1 * max_k
            s2_factor = one
         else
            itmp = nirrep*(irp_k-1) + irp_n
            itmp = off_pq(itmp,h_Offsets_Ij) + ( max_n * k )
            add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
            s2_factor = minus
         end if

c      o loop over orbital n (not an offset, here)
         inc_w = max_e * max_m
         do n = 1, max_n

c         o load/transpose/increment T2(1)(nK,eR)
            call load_t2_virts(
     &              dcore(add_t2p1), dscr(start_t2p1),
     &              irp_e, max_e,             ispin2,
     &              irp_r, max_r, 0, max_r-1, ispin1, h_Offsets_Ab
     &           )
            add_t2p1 = add_t2p1 + inc_t2p1

c         o add in -<eJ|Mn>*T2(1)(nK,eR) for all RM
            call B_GEMM('T','N',max_r,max_m,max_e,
     &                  s2_factor,dscr(start_t2p1),max_e,
     &                            dscr(add_w),     max_e,
     &                  one,      s2rmjk(add_s2),  max_r)

c      o increment and close the loop over orbital n
            add_w = add_w + inc_w
c        end do n = 1, max_n
         end do

c   o close the loop over irreps of n
c     else if (.not.do_irp_n) then
      else
         add_w = add_w + ( max_e * max_m * max_n )
c     end if (do_irp_n)
      end if
c     end do irp_n = 1, nirrep
      end do
#endif /* _DO_GS_eJMn */

c     end if ((max_e.ne.0).and.(max_mn.ne.0))
      end if

      max_mn = irpdpd(irp_mn,dist_iltj)
      max_e  = vrt(irp_e,ispin1)
      if ((max_e.ne.0).and.(max_mn.ne.0)) then

      start_scr = 1 + ( max_mn * max_e )
      if (start_scr.gt.dim_dscr) then
         write(*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
         call xcc_errex
      end if

      if (iuhf.eq.0) then
c      o antisymmetrize <eJ|Mn> to <EJ||MN>
         call antisym_ket_pp_pltp_in_situ(
     &           dscr(1),max_e,
     &           pop(1,1),irp_mn,off_pq(1,h_Offsets_Ij)
     &        )
      else
c      o load all <MN||JE> and convert to <EJ||MN>
         itmp = nirrep*(irp_e-1) + irp_j
         itmp = 1 + off_pq(itmp,off_ia) + j
         start_scr = 1
         do e = 1, max_e
            call getlst(dscr(start_scr),itmp,1,2,irp_mn,list_ijka)
            start_scr = start_scr + max_mn
            itmp = itmp + max_j
         end do
         call dmat_xpose(dscr,max_mn,max_e)
c     end if (iuhf.eq.0)
      end if

#ifdef _DO_HH_EJMN
      start_t2p1 = start_scr
c   o loop over irreps of N
      do irp_n = 1, nirrep
         irp_m = dirprd(irp_n,irp_mn)
         irp_r = dirprd(irp_m,irp_rm)
         max_m = pop(irp_m,ispin1)
         max_n = pop(irp_n,ispin1)
         max_r = vrt(irp_r,ispin1)
         do_irp_n = ((max_n.ne.0).and.(max_m.ne.0).and.(max_r.ne.0))
      if (do_irp_n) then

c      o point to the RM block of S2
         itmp = nirrep*(irp_m-1) + irp_r
         add_s2 = 1 + off_pq(itmp,off_ai)

c      o allow for unpacked T2(1)(NK,*)
         start_scr = start_t2p1 + ( max_e * max_r )
         if (start_scr.gt.(dim_dscr+1)) then
            write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
            call xcc_errex
         end if

c      o four cases: irp_m (==,!=) irp_n, and irp_k (==,!=) irp_n
c        If irp_m and irp_n are the same, then the effects of looping
c        over orbital n are much worse; therefore, that is the highest
c        criterion.
         if (irp_m.ne.irp_n) then

c         o find the first <*||MN> and calculate the metrics
            if (irp_m.lt.irp_n) then
               s2_factor = one
               inc_w = max_e * max_m
               ldw   = max_e
               itmp  = rshift(irp_n*(irp_n-1),1) + irp_m
               itmp  = off_pltp(itmp,off_iltj)
               add_w = 1 + ( max_e * itmp )
            else
               s2_factor = minus
               inc_w = max_e
               ldw   = max_e * max_n
               itmp  = rshift(irp_m*(irp_m-1),1) + irp_n
               itmp  = off_pltp(itmp,off_iltj)
               add_w = 1 + ( max_e * itmp )
            end if

c         o find the first T2(1)(NK,*), calculate the metrics, and
c           loop over orbital N
            itmp = dirprd(irp_n,irp_k)
            inc_t2p1 = irpdpd(itmp,dist_altb)
            add_t2p1 = ndx_t2p1(itmp,ispin1)
            if (irp_n.ne.irp_k) then
               if (irp_n.lt.irp_k) then
                  itmp = rshift(irp_k*(irp_k-1),1) + irp_n
                  itmp = off_pltp(itmp,off_iltj) + ( max_n * k )
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               else
                  itmp = rshift(irp_n*(irp_n-1),1) + irp_k
                  itmp = off_pltp(itmp,off_iltj) + k
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
                  inc_t2p1 = inc_t2p1 * max_k
                  s2_factor = -s2_factor
               end if
               do n = 1, max_n
                  call load_t2_virts(
     &                    dcore(add_t2p1), dscr(start_t2p1),
     &                    irp_e, max_e,             ispin1,
     &                    irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                 )
                  add_t2p1 = add_t2p1 + inc_t2p1
                  call B_GEMM('T','N',max_r,max_m,max_e,
     &                        s2_factor,dscr(start_t2p1),max_e,
     &                                  dscr(add_w),     ldw,
     &                        one,      s2rmjk(add_s2),  max_r)
                  add_w = add_w + inc_w
               end do
c           else if (irp_n.eq.irp_k) then
            else
c               itmp = rshift(irp_k*(irp_k-1),1) + irp_k
               itmp = rshift(irp_k*(irp_k+1),1)
               itmp = off_pltp(itmp,off_iltj) + rshift(k*(k-1),1)
               add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               if (k.ne.0) then
                  do n = 1, k
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + inc_t2p1
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           s2_factor,dscr(start_t2p1),max_e,
     &                                     dscr(add_w),     ldw,
     &                           one,      s2rmjk(add_s2),  max_r)
                     add_w = add_w + inc_w
                  end do
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * k )
               end if
               add_w = add_w + inc_w
               if ((k+1).ne.max_k) then
                  do n = 1, (max_n-1)-k
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + ( inc_t2p1 * (n+k) )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           -s2_factor,dscr(start_t2p1),max_e,
     &                                      dscr(add_w),     ldw,
     &                           one,       s2rmjk(add_s2),  max_r)
                     add_w = add_w + inc_w
                  end do
               end if
c           end if (irp_n.ne.irp_k)
            end if

c        else if (irp_m.eq.irp_n) then
         else

c         o allow for unpacked <*||MN>
            if ((start_scr+(max_e*max_m)).gt.(dim_dscr+1)) then
               write (*,*) '@MAKE_QF3_HH_AIJK: Not enough memory.'
               call xcc_errex
            end if

c         o find the first <*||MN>
            s2_factor = one
c            itmp  = rshift(irp_n*(irp_n-1),1) + irp_n
            itmp  = rshift(irp_n*(irp_n+1),1)
            itmp  = off_pltp(itmp,off_iltj)
            add_w = 1 + ( max_e * itmp )

c         o find the first T2(1)(NK,*), calculate the metrics, and
c           loop over orbital N
            itmp = dirprd(irp_n,irp_k)
            inc_t2p1 = irpdpd(itmp,dist_altb)
            add_t2p1 = ndx_t2p1(itmp,ispin1)
            if (irp_n.ne.irp_k) then
               if (irp_n.lt.irp_k) then
                  itmp = rshift(irp_k*(irp_k-1),1) + irp_n
                  itmp = off_pltp(itmp,off_iltj) + ( max_n * k )
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               else
                  itmp = rshift(irp_n*(irp_n-1),1) + irp_k
                  itmp = off_pltp(itmp,off_iltj) + k
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
                  inc_t2p1 = inc_t2p1 * max_k
                  s2_factor = -s2_factor
               end if
               do n = 0, max_n-1
                  call load_t2_virts(
     &                    dcore(add_t2p1), dscr(start_t2p1),
     &                    irp_e, max_e,             ispin1,
     &                    irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                 )
                  add_t2p1 = add_t2p1 + inc_t2p1
                  call extract_s2_dists(
     &                    dscr(add_w), dscr(start_scr), max_e,
     &                    irp_m, max_m, 0, max_m-1, max_m,
     &                    irp_n, max_n, n,
     &                    .true., .false.
     &                 )
                  call B_GEMM('T','N',max_r,max_m,max_e,
     &                        s2_factor,dscr(start_t2p1),max_e,
     &                                  dscr(start_scr), max_e,
     &                        one,      s2rmjk(add_s2),  max_r)
               end do
c           else if (irp_n.eq.irp_k) then
            else
c               itmp = rshift(irp_k*(irp_k-1),1) + irp_k
               itmp = rshift(irp_k*(irp_k+1),1)
               itmp = off_pltp(itmp,off_iltj) + rshift(k*(k-1),1)
               add_t2p1 = add_t2p1 + ( inc_t2p1 * itmp )
               if (k.ne.0) then
                  do n = 0, k-1
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + inc_t2p1
                     call extract_s2_dists(
     &                       dscr(add_w), dscr(start_scr), max_e,
     &                       irp_m, max_m, 0, max_m-1, max_m,
     &                       irp_n, max_n, n,
     &                       .true., .false.
     &                    )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           s2_factor,dscr(start_t2p1),max_e,
     &                                     dscr(start_scr), max_e,
     &                           one,      s2rmjk(add_s2),  max_r)
                  end do
                  add_t2p1 = add_t2p1 + ( inc_t2p1 * k )
               end if
               if ((k+1).ne.max_k) then
                  do n = 1, (max_n-1)-k
                     call load_t2_virts(
     &                       dcore(add_t2p1), dscr(start_t2p1),
     &                       irp_e, max_e,             ispin1,
     &                       irp_r, max_r, 0, max_r-1, ispin1, off_altb
     &                    )
                     add_t2p1 = add_t2p1 + ( inc_t2p1 * (n+k) )
                     call extract_s2_dists(
     &                       dscr(add_w), dscr(start_scr), max_e,
     &                       irp_m, max_m, 0, max_m-1, max_m,
     &                       irp_n, max_n, k+n,
     &                       .true., .false.
     &                    )
                     call B_GEMM('T','N',max_r,max_m,max_e,
     &                           -s2_factor,dscr(start_t2p1),max_e,
     &                                      dscr(start_scr), max_e,
     &                           one,       s2rmjk(add_s2),  max_r)
                  end do
               end if
c           end if (irp_n.ne.irp_k)
            end if

c        end if (irp_m.ne.irp_n)
         end if

c   o close the loop over irreps of N
c     end if (do_irp_n)
      end if
c     end do irp_n = 1, nirrep
      end do
#endif /* _DO_HH_EJMN */

c     end if ((max_e.ne.0).and.(max_mn.ne.0))
      end if

c ----------------------------------------------------------------------

c   o close the loop over irreps of MN/Mn
c     end do irp_mn = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

