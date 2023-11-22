
c This routine builds <PQ{}EI> intermediates used in <adj(T2(1)^2)|W|T3>.
c The actual contractions are W*T2(1) -> S2(+1).

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(dim_dscr) : the scratch array for doubles

c    integer ispin1 : the spin case to handle

c    integer from_e, to_e, irp_e, max_e  : e statistics
c    integer i,            irp_i, max_i  : i statistics
c    integer                      max_pq : the PQ distribution length

c    integer irp_pq, irp_ei : the left and right symmetries
c                             NOTE: this is not fully tested for
c                                   cases when irp_pq != irp_ei

c OUTPUT
c    double precision s2pqei(*) : the destination array
c                                 NOTE: this is initialized to 0.0d0

c ----------------------------------------------------------------------

#include "disttype.h"
#include "listnum.h"

c#define _DEBUG_MAKE_QF3_HH_ABCI
#define _DO_HH_MNEI
#define _DO_GS_EfQm
#define _DO_HH_EFQM

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine make_qf3_hh_abci(
     &              dcore, dscr, dim_dscr,
     &              ispin1, s2pqei,
     &              from_e, to_e, irp_e,  max_e,  irp_ei,
     &              i,            irp_i,  max_i,
     &                            irp_pq, max_pq
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr), s2pqei(*)
      integer ispin1, ispin2
      integer    irp_e,  max_e, from_e, to_e
      integer i, irp_i,  max_i
      integer    irp_pq, max_pq, irp_ei

c INTERNAL VARIABLES
      integer list_ijka, list_a1b2c1i2, list_abci
      integer dist_iltj, dist_a1i2, dist_a1b2, dist_ai, dist_altb
      integer off_ia, off_a2b1, off_altb, off_ai, off_iltj
      integer e, num_e
      integer m, irp_m, max_m
      integer    irp_f, max_f
      integer p, irp_p, max_p
      integer q, irp_q, max_q
      integer    irp_qm, max_qm
      integer    irp_pf, max_pf
      integer            max_ef, ndx_ef, rhf_ndx_fe
      logical do_irp_q
      integer start_scr, start_product
      integer start_t2p1, add_t2p1, inc_t2p1, unp_t2p1
      integer             add_w,              unp_w
      integer add_s2
      integer itmp, max_dists
      double precision s2_factor

c PARAMETERS
      double precision zero, one, minus
      parameter (zero=0.0d0, one=1.0d0, minus=-1.0d0)

c EXTERNAL FUNCTIONS
      double precision ddot

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com" /* for pop/vrt */
#include "method.com" /* for iuhf */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin1 is 1 or 2
      if ((ispin1.lt.1).or.(2.lt.ispin1)) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   ispin1 = ',ispin1
         itmp = -1
      end if
c   o assert max_pq is well-behaved
      if (max_pq.lt.0) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   max_pq = ',max_pq
         itmp = -1
      end if
c   o assert e is well-behaved
      if ((max_e.le.to_e).or.(to_e.lt.from_e).or.(from_e.lt.0)) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   from_e = ',from_e
         write(*,*) '   to_e   = ',to_e
         write(*,*) '   max_e  = ',max_e
         itmp = -1
      end if
c   o assert i is well-behaved
      if ((max_i.le.i).or.(i.lt.0)) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   i = ',i,', max_i = ',max_i
         itmp = -1
      end if
c   o assert irreps are in [1,nirrep]
      if ((irp_e .lt.1).or.(nirrep.lt.irp_e ).or.
     &    (irp_i .lt.1).or.(nirrep.lt.irp_i ).or.
     &    (irp_pq.lt.1).or.(nirrep.lt.irp_pq).or.
     &    (irp_ei.lt.1).or.(nirrep.lt.irp_ei)    ) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   irp_e  = ',irp_e
         write(*,*) '   irp_i  = ',irp_i
         write(*,*) '   irp_pq = ',irp_pq
         write(*,*) '   irp_ei = ',irp_ei
         itmp = -1
      end if
c   o assert only alpha intermediates are built if RHF
      if ((ispin1.eq.2).and.(iuhf.eq.0)) then
         write(*,*) '@MAKE_QF3_HH_ABCI: Assertion failed.'
         write(*,*) '   ispin1 = ',ispin1
         write(*,*) '   iuhf   = ',iuhf
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      if (max_pq.lt.1) return

c ----------------------------------------------------------------------

c The plan is:
c   - initialize S2 with T2(1)(MN,PQ)*<MN||EI>
c   o loop over irreps of Qm/QM ( -> irp_f )
c      o if alpha
c         o loop over irreps of m
c            o loop over orbital m
c               - if RHF, load T2(1)(MI,PF)
c               o loop over orbital Q
c                  - load <Ef|Qm>
c                  - add in -T2(1)(Im,Pf)*<Ef|Qm>
c                  - if RHF, add in T2(1)(MI,PF)*<EF||QM>
c      o if beta
c         o loop over irreps of q
c            o loop over orbital q
c               o loop over orbital M
c                  - load <Fe|Mq>
c                  - add in -T2(1)(Mi,Fp)*<Fe|Mq>
c      o if UHF
c         o loop over irreps of M
c            o loop over orbital M
c               o loop over orbital Q
c                  - load <FE||QM>
c                  - add in T2(1)(MI,FP)*<FE||QM>

c   o set handles
      if (ispin1.eq.1) then
         ispin2 = 2
         list_a1b2c1i2 = h_GS_AbCi_Ab_Ci
         list_abci = h_HH_ABCI_AltB_CI
         list_ijka = h_HH_IJKA_IltJ_KA
         dist_altb = h_DistType_AltB
         dist_iltj = h_DistType_IltJ
         dist_a1i2 = h_DistType_Ai
         dist_ai   = h_DistType_AI
         off_ia    = h_Offsets_IA
         off_ai    = h_Offsets_AI
         off_iltj  = h_Offsets_IltJ
         off_altb  = h_Offsets_AltB
         off_a2b1  = h_Offsets_aB
      else
         ispin2 = 1
         list_a1b2c1i2 = h_GS_AbIc_Ab_Ic
         list_abci = h_HH_abci_altb_ci
         list_ijka = h_HH_ijka_iltj_ka
         dist_altb = h_DistType_altb
         dist_iltj = h_DistType_iltj
         dist_a1i2 = h_DistType_Ia
         dist_ai   = h_DistType_ai
         off_ia    = h_Offsets_ia
         off_ai    = h_Offsets_ai
         off_iltj  = h_Offsets_iltj
         off_altb  = h_Offsets_altb
         off_a2b1  = h_Offsets_Ab
      end if
      dist_a1b2 = h_DistType_Ab

#ifdef _DEBUG_MAKE_QF3_HH_ABCI
      print '(a,i1,a,i1,a,i3, a,i1,a,i3,a,i3)',
     &         '@MAKE_QF3_HH_ABCI: ',ispin1,': [',irp_i,'] ',i,
     &         ' with [',irp_e,'] ',from_e,'-',to_e
#endif /* _DEBUG_MAKE_QF3_HH_ABCI */

c   o   T2(1)(MN,PQ) * <MN||EI>
      m = irpdpd(irp_ei,dist_iltj)
      if (m.ne.0) then
#ifdef _DO_HH_MNEI
      if (iuhf.eq.0) then
c      o we need at least 1 <Mn| distribution
         if (dim_dscr.lt.irpdpd(irp_ei,h_DistType_Ij)) then
            write(*,*) '@MAKE_QF3_HH_ABCI: Not enough memory.'
            call xcc_errex
         end if
c      o point to T2(1)(MN,PQ)
         add_t2p1 = ndx_t2p1(irp_pq,ispin1)
c      o load and antisymmetrize <Mn|Ie>,
c        gemm with T2(1), and write into S2
         add_s2 = 1
         itmp = nirrep*(irp_e-1) + irp_i
         add_w = 1 + off_pq(itmp,h_Offsets_Ia) + (max_i*from_e) + i
         do e = from_e, to_e
            call getlst(dscr,add_w,1,2,irp_ei,h_GS_IjKa_Ij_Ka)
            call antisym_pp_pltp_in_situ(
     &              dscr,irpdpd(irp_ei,h_DistType_Ij),
     &              pop(1,1),1,irp_ei,off_pq(1,h_Offsets_Ij)
     &           )
            do itmp = 0, max_pq-1
               s2pqei(add_s2+itmp) =
     &          - ddot(m,dcore(add_t2p1+itmp),max_pq,dscr,1)
            end do
            add_s2 = add_s2 + max_pq
            add_w  = add_w  + max_i
c        end do e = from_e, to_e
         end do
c     else if (iuhf.eq.1) then
      else
c      o we need at least 1 <MN| distribution
         if (dim_dscr.lt.m) then
            write(*,*) '@MAKE_QF3_HH_ABCI: Not enough memory.'
            call xcc_errex
         end if
c      o point to T2(1)(MN,PQ)
         add_t2p1 = ndx_t2p1(irp_pq,ispin1)
c      o load many <MN||EI> lists, gemm with T2(1), and write into S2
         max_dists = dim_dscr / m
         itmp = nirrep*(irp_e-1) + irp_i
         add_w = 1 + off_pq(itmp,off_ia) + (max_i*from_e) + i
         add_s2 = 1
         e = from_e
         do while (e.lt.(to_e+1))
            start_scr = 1
            num_e = to_e+1 - e
            num_e = min(num_e,max_dists)
            do itmp = 1, num_e
               call getlst(dscr(start_scr),add_w,1,
     &                     2,irp_ei,list_ijka)
               start_scr = start_scr + m
               add_w     = add_w     + max_i
            end do
            call B_GEMM('N','N',max_pq,num_e,m,
     &                  minus,dcore(add_t2p1),max_pq,
     &                        dscr,           m,
     &                  zero, s2pqei(add_s2), max_pq)
            add_s2 = add_s2 + (max_pq*num_e)
            e      = e      + num_e
c        end do while (e.lt.(to_e+1))
         end do
c     end if (iuhf.eq.0)
      end if
#else /* _DO_HH_MNEI */
         do itmp = 1, max_pq*(to_e+1-from_e)
            s2pqei(itmp) = zero
         end do
#endif /* _DO_HH_MNEI */
c     else if (m.eq.0) then
      else
         do itmp = 1, max_pq*(to_e+1-from_e)
            s2pqei(itmp) = zero
         end do
c     end if (m.ne.0)
      end if

c ----------------------------------------------------------------------

      num_e = to_e+1 - from_e

c   o loop over irreps of Qm/QM ( -> irp_f )
      do irp_qm = 1, nirrep
         irp_f  = dirprd(irp_qm,irp_e)

c ----------------------------------------------------------------------

      if (ispin1.eq.1) then

c   o - T2(1)(Im,Pf) * <Ef|Qm>
      max_qm = irpdpd(irp_qm,h_DistType_Ai)
      max_f  = vrt(irp_f,2)
      if ((max_f.ne.0).and.(max_qm.ne.0)) then

c   o memory allocation will be checked for each irp_q
      max_ef = irpdpd(irp_qm,h_DistType_Ab)
      start_scr = 1 + max_ef

c   o <Ef|Qm>;Ef,Qm
c     The actual <Ef| offset will be constant for all |Qm>
c     provided each distribution is loaded one at a time and into
c     the same place...
      itmp = nirrep*(irp_f-1) + irp_e
      ndx_ef = 1 + off_pq(itmp,h_Offsets_Ab) + from_e

c   o point to <FE| for RHF
      if (iuhf.eq.0) then
         itmp = nirrep*(irp_e-1) + irp_f
         rhf_ndx_fe = 1 + off_pq(itmp,h_Offsets_Ab) + (max_f*from_e)
      end if

c   o loop over irreps of m
      do irp_m = 1, nirrep
         irp_q = dirprd(irp_m,irp_qm)
         irp_p = dirprd(irp_q,irp_pq)
         max_m = pop(irp_m,2)
         max_p = vrt(irp_p,1)
         max_q = vrt(irp_q,1)
         do_irp_q = ((max_m.ne.0).and.
     &               (max_p.ne.0).and.
     &               (max_q.ne.0).and.
     &               ((irp_p.ne.irp_q).or.(max_p.gt.1))
     &              )
      if (do_irp_q) then

         irp_pf = dirprd(irp_p,irp_f)

c      o If this is an RHF calculation, then unp_t2p1 will replace
c        start_scr for the unpacked T2(1) amplitudes. At the end of this
c        Q irrep, start_scr must recover that space.
         if (iuhf.eq.0) then
            unp_t2p1  = start_scr
            start_scr = start_scr + (max_p*max_f)
         end if

c      o we need all the PQ*E contributions and start_product will track Q
         itmp = max_p * max_q * num_e
         if ((start_scr+itmp).gt.(dim_dscr+1)) then
            write(*,*) '@MAKE_QF3_HH_ABCI: Not enough memory.'
            call xcc_errex
         end if
         do q = 0, itmp-1
            dscr(start_scr+q) = zero
         end do

c      o pointers/addresses
c        <Ef|Qm>;Ef,Qm
         itmp = nirrep*(irp_m-1) + irp_q
         add_w = 1 + off_pq(itmp,h_Offsets_Ai)
c        (Pf,Im);Pf,Im
         itmp = nirrep*(irp_m-1) + irp_i
         itmp = off_pq(itmp,h_Offsets_Ij) + i
         add_t2p1 =   ndx_t2p1(irp_pf,3)
     &              + (irpdpd(irp_pf,h_DistType_Ab)*itmp)
         itmp = nirrep*(irp_f-1) + irp_p
         add_t2p1 = add_t2p1 + off_pq(itmp,h_Offsets_Ab)
c        (MI,PF);P<F,M<I
         if (iuhf.eq.0) then
            if (irp_m.lt.irp_i) then
               itmp = rshift(irp_i*(irp_i-1),1) + irp_m
               itmp = off_pltp(itmp,h_Offsets_IltJ) + (max_m*i)
               inc_t2p1 = irpdpd(irp_pf,h_DistType_AltB)
               s2_factor = one
            else
            if (irp_i.lt.irp_m) then
               itmp = rshift(irp_m*(irp_m-1),1) + irp_i
               itmp = off_pltp(itmp,h_Offsets_IltJ) + i
               inc_t2p1 = irpdpd(irp_pf,h_DistType_AltB)*max_i
               s2_factor = minus
            else
               itmp = rshift(irp_i*(irp_i+1),1)
               itmp = off_pltp(itmp,h_Offsets_IltJ)
               itmp = itmp + rshift(i*(i-1),1)
c           end if (irp_i.lt.irp_m)
            end if
c           end if (irp_m.lt.irp_i)
            end if
            start_t2p1 =   ndx_t2p1(irp_pf,1)
     &                   + (irpdpd(irp_pf,h_DistType_AltB)*itmp)
c        end if (iuhf.eq.0)
         end if

c      o loop over orbitals m and Q
         do m = 0, max_m-1
c         o unpack T2(1)(MI,PF)
            if (iuhf.eq.0) then
               if (irp_m.ne.irp_i) then
                  itmp = start_t2p1
                  start_t2p1 = start_t2p1 + inc_t2p1
               else
                  if (m.lt.i) then
                     itmp = m
                     s2_factor = one
                  else
                     itmp = i + rshift((m*(m-1)+i*(-i+1)),1)
                     s2_factor = minus
                  end if
                  itmp =   start_t2p1
     &                   + (irpdpd(irp_pf,h_DistType_AltB)*itmp)
               end if
               call load_t2_virts(
     &                 dcore(itmp), dscr(unp_t2p1),
     &                 irp_p, max_p,             1,
     &                 irp_f, max_f, 0, max_f-1, 1, h_Offsets_AltB
     &              )
            end if
c         o reset the PQ*E product pointer
            start_product = start_scr
         do q = 0, max_q-1

c         o load <Ef|Qm>
            call getlst(dscr,add_w,1,2,irp_qm,list_a1b2c1i2)

#ifdef _DO_GS_EfQm
c         o - T2(1)(Im,Pf) * <Ef|Qm>
            call B_GEMM('N','T',max_p,num_e,max_f,
     &                  minus,dcore(add_t2p1),    max_p,
     &                        dscr(ndx_ef),       max_e,
     &                  one,  dscr(start_product),max_p*max_q)
#endif /* _DO_GS_EfQm */

#ifdef _DO_HH_EFQM
c         o RHF: T2(1)(MI,PF)*(<EF|QM>-<FE|QM>)
            if (iuhf.eq.0) then
            if ((m.ne.i).or.(irp_m.ne.irp_i)) then
c            o   T2(1)(MI,PF) * <EF|QM>
               call B_GEMM('N','T',max_p,num_e,max_f,
     &                     s2_factor,dscr(unp_t2p1),     max_p,
     &                               dscr(ndx_ef),       max_e,
     &                     one,      dscr(start_product),max_p*max_q)
c            o - T2(1)(MI,PF) * <FE|QM>
               call B_GEMM('N','N',max_p,num_e,max_f,
     &                     -s2_factor,dscr(unp_t2p1),     max_p,
     &                                dscr(rhf_ndx_fe),   max_f,
     &                     one,       dscr(start_product),max_p*max_q)
c           end if ((m.ne.i).or.(irp_m.ne.irp_i))
            end if
c           end if (iuhf.eq.0)
            end if
#endif /* _DO_HH_EFQM */

c      o increment and close the loop over orbitals m and Q
            start_product = start_product + max_p
            add_w         = add_w         + 1
c        end do q = 0, max_q-1
         end do
            add_t2p1 = add_t2p1 + (irpdpd(irp_pf,h_DistType_Ab)*max_i)
c        end do m = 0, max_m-1
         end do

#include "make_qf3_hh_abci.F.src"

c   o recover scr and close the loop over irreps of m
         if (iuhf.eq.0) start_scr = unp_t2p1
c     end if (do_irp_q)
      end if
c     end do irp_m = 1, nirrep
      end do

c     end if ((max_f.ne.0).and.(max_qm.ne.0))
      end if

c ----------------------------------------------------------------------

c     else if (ispin1.eq.2) then
      else

c   o - T2(1)(Mi,Fp) * <Fe|Mq>
      max_qm = irpdpd(irp_qm,h_DistType_Ia)
      max_f  = vrt(irp_f,1)
      if ((max_f.ne.0).and.(max_qm.ne.0)) then

c   o memory allocation will be checked for each irp_q
      max_ef = irpdpd(irp_qm,h_DistType_Ab)
      start_scr = 1 + max_ef

c   o <Fe|Mq>;Fe,Mq
c     The actual <Fe| offset will be constant for all |Mq>
c     provided each distribution is loaded one at a time and into
c     the same place...
      itmp = nirrep*(irp_e-1) + irp_f
      ndx_ef = 1 + off_pq(itmp,h_Offsets_Ab) + (max_f*from_e)

c   o loop over irreps of q
      do irp_q = 1, nirrep
         irp_m = dirprd(irp_q,irp_qm)
         irp_p = dirprd(irp_q,irp_pq)
         max_m = pop(irp_m,1)
         max_p = vrt(irp_p,2)
         max_q = vrt(irp_q,2)
         do_irp_q = ((max_m.ne.0).and.
     &               (max_p.ne.0).and.
     &               (max_q.ne.0).and.
     &               ((irp_p.ne.irp_q).or.(max_p.gt.1))
     &              )
      if (do_irp_q) then

         irp_pf = dirprd(irp_p,irp_f)

c      o we need all the PQ*E contributions and start_product will track Q
         itmp = max_p * max_q * num_e
         if ((start_scr+itmp).gt.(dim_dscr+1)) then
            write(*,*) '@MAKE_QF3_HH_ABCI: Not enough memory.'
            call xcc_errex
         end if
         do q = 0, itmp-1
            dscr(start_scr+q) = zero
         end do

c      o pointers/addresses
c        <Fe|Mq>;Fe,Mq
         itmp = nirrep*(irp_q-1) + irp_m
         add_w = 1 + off_pq(itmp,h_Offsets_Ia)
c        (Fp,Mi);Fp,Mi
         itmp = nirrep*(irp_i-1) + irp_m
         itmp = off_pq(itmp,h_Offsets_Ij) + (max_m*i)
         start_t2p1 =   ndx_t2p1(irp_pf,3)
     &                + (irpdpd(irp_pf,h_DistType_Ab)*itmp)
         itmp = nirrep*(irp_p-1) + irp_f
         start_t2p1 = start_t2p1 + off_pq(itmp,h_Offsets_Ab)

c      o loop over orbitals q and M
            start_product = start_scr
         do q = 0, max_q-1
            add_t2p1 = start_t2p1
         do m = 0, max_m-1

c         o load <Fe|Mq>
            call getlst(dscr,add_w,1,2,irp_qm,list_a1b2c1i2)

#ifdef _DO_GS_EfQm
c         o - T2(1)(Mi,Fp) * <Fe|Mq>
            call B_GEMM('T','N',max_p,num_e,max_f,
     &                  minus,dcore(add_t2p1),    max_f,
     &                        dscr(ndx_ef),       max_f,
     &                  one,  dscr(start_product),max_p*max_q)
#endif /* _DO_GS_EfQm */

c      o increment and close the loop over orbitals q and M
            add_t2p1 = add_t2p1 + irpdpd(irp_pf,h_DistType_Ab)
            add_w    = add_w    + 1
c        end do m = 0, max_m-1
         end do
            start_product = start_product + max_p
c        end do q = 0, max_q-1
         end do

#include "make_qf3_hh_abci.F.src"

c   o close the loop over irreps of q
c     end if (do_irp_q)
      end if
c     end do irp_q = 1, nirrep
      end do

c     end if ((max_f.ne.0).and.(max_qm.ne.0))
      end if

c     end if (ispin1.eq.1)
      end if

c ----------------------------------------------------------------------

#ifdef _DO_HH_EFQM
      if (iuhf.eq.1) then

c   o   T2(1)(MI,FP) * <FE||QM>
      max_qm = irpdpd(irp_qm,dist_ai)
      max_f  = vrt(irp_f,ispin1)
      if ((max_f.ne.0).and.(max_qm.ne.0).and.
     &    ((irp_f.ne.irp_e).or.(max_f.gt.1))
     &   ) then

c   o memory allocation will be checked for each irp_q
      max_ef = irpdpd(irp_qm,dist_altb)
      start_scr = 1 + max_ef

c   o loop over irreps of M
      do irp_m = 1, nirrep
         irp_q = dirprd(irp_m,irp_qm)
         irp_p = dirprd(irp_q,irp_pq)
         max_m = pop(irp_m,ispin1)
         max_p = vrt(irp_p,ispin1)
         max_q = vrt(irp_q,ispin1)
         do_irp_q = ((max_m.ne.0).and.
     &               (max_p.ne.0).and.
     &               (max_q.ne.0).and.
     &               ((irp_p.ne.irp_q).or.(max_p.gt.1)).and.
     &               ((irp_m.ne.irp_i).or.(max_m.gt.1))
     &              )
      if (do_irp_q) then

         irp_pf = dirprd(irp_p,irp_f)
         max_pf = irpdpd(irp_pf,dist_altb)

c      o we need <FE||QM> and T2(1)(FP,MI) unpacked as well as
c        all the PQ*E contributions and start_product will track Q
         unp_w = start_scr
         unp_t2p1 = unp_w + (max_f*num_e)
         start_scr = unp_t2p1 + (max_f*max_p)
         itmp = max_p * max_q * num_e
         if ((start_scr+itmp).gt.(dim_dscr+1)) then
            write(*,*) '@MAKE_QF3_HH_ABCI: Not enough memory.'
            call xcc_errex
         end if
         do q = 0, itmp-1
            dscr(start_scr+q) = zero
         end do

c      o pointers/addresses
c        <FE||QM>;F<E,QM
         itmp = nirrep*(irp_m-1) + irp_q
         add_w = 1 + off_pq(itmp,off_ai)
c        (FP,MI);F<P,M<I
         if (irp_m.lt.irp_i) then
            itmp = rshift(irp_i*(irp_i-1),1) + irp_m
            itmp = off_pltp(itmp,off_iltj) + (max_m*i)
            inc_t2p1 = max_pf
            s2_factor = one
         else
         if (irp_i.lt.irp_m) then
            itmp = rshift(irp_m*(irp_m-1),1) + irp_i
            itmp = off_pltp(itmp,off_iltj) + i
            inc_t2p1 = max_pf*max_i
            s2_factor = minus
         else
            itmp = rshift(irp_i*(irp_i+1),1)
            itmp = off_pltp(itmp,off_iltj)
            itmp = itmp + rshift(i*(i-1),1)
c        end if (irp_i.lt.irp_m)
         end if
c        end if (irp_m.lt.irp_i)
         end if
         start_t2p1 = ndx_t2p1(irp_pf,ispin1) + (max_pf*itmp)

c      o loop over orbitals M and Q
         do m = 0, max_m-1
         if ((m.ne.i).or.(irp_m.ne.irp_i)) then
c         o unpack T2(1)(MI,FP)
            if (irp_m.ne.irp_i) then
               itmp = start_t2p1
               start_t2p1 = start_t2p1 + inc_t2p1
            else
               if (m.lt.i) then
                  itmp = m
                  s2_factor = one
               else
                  itmp = i + rshift((m*(m-1)+i*(-i+1)),1)
                  s2_factor = minus
               end if
               itmp = start_t2p1 + (max_pf*itmp)
            end if
            call load_t2_virts(
     &              dcore(itmp), dscr(unp_t2p1),
     &              irp_f, max_f,             ispin1,
     &              irp_p, max_p, 0, max_p-1, ispin1, off_altb
     &           )
c         o reset the PQ*E pointer
            start_product = start_scr
         do q = 0, max_q-1

c         o load <FE||QM> and unpack
            call getlst(dscr,add_w,1,2,irp_qm,list_abci)
            call load_t2_virts(
     &              dscr, dscr(unp_w),
     &              irp_f, max_f,               ispin1,
     &              irp_e, max_e, from_e, to_e, ispin1, off_altb
     &           )

c         o   T2(1)(MI,FP) * <FE||QM>
            call B_GEMM('T','N',max_p,num_e,max_f,
     &                  s2_factor,dscr(unp_t2p1),     max_f,
     &                            dscr(unp_w),        max_f,
     &                  one,      dscr(start_product),max_p*max_q)

c      o increment and close the loop over orbitals M and Q
            start_product = start_product + max_p
            add_w         = add_w         + 1
c        end do q = 0, max_q-1
         end do
         else
            add_w = add_w + max_q
c        end if ((m.ne.i).or.(irp_m.ne.irp_i))
         end if
c        end do m = 0, max_m-1
         end do

#include "make_qf3_hh_abci.F.src"

c   o recover scr and close the loop over irreps of M
         start_scr = unp_w
c     end if (do_irp_q)
      end if
c     end do irp_m = 1, nirrep
      end do

c     end if ((max_f.ne.0).and.(max_qm.ne.0))
      end if

c     end if (iuhf.eq.1)
      end if
#endif /* _DO_HH_EFQM */

c ----------------------------------------------------------------------

c   o close the loop over irreps of Qm/QM
c     end do irp_qm = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

