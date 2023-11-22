
c This routine processes single distributions of T3AAA amplitudes.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(dim_dscr) : scratch array for doubles
c    integer          ispin : spin case to handle
c    double precision t3aaa(max_abc) : T3(IJK,*) distribution
c    integer off_altbltc(1+nirrep*(nirrep+1)/2) : irrep offsets for T3AAA/BBB

c    [ orbital info ]

c    double precision e5staaa : the 5-th order energy correction from T3AAA
c                               <0|adj(T1)*W(-2)*T3AAA|0>

#include "disttype.h"
#include "listnum.h"

c#define _DEBUG_PROCESS_T3AAA
#define _DO_E5ST_T1I
#define _DO_E5ST_T1J
#define _DO_E5ST_T1K
c#define _DO_E5QfT_S1I
c#define _DO_E5QfT_S1J
c#define _DO_E5QfT_S1K
c#define _DO_E5QfT_S2I
c#define _DO_E5QfT_S2J
c#define _DO_E5QfT_S2K
c#define _DO_E5QfT_S2IJ
c#define _DO_E5QfT_S2IK
c#define _DO_E5QfT_S2JK

      subroutine process_t3aaa(
     &              dcore,dscr,dim_dscr,
     &              ispin,
     &              t3aaa,max_abc,off_altbltc,
     &              i,irp_i,max_i,irp_ij,add_t2ij,
     &              j,irp_j,max_j,irp_jk,add_t2jk,
     &              k,irp_k,max_k,irp_ik,add_t2ik,irp_ijk,
     &              e5staaa
     &           )
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*), t3aaa(*)
      integer ispin, max_abc, off_altbltc(*), dim_dscr
      integer irp_i, i, max_i, irp_ij, add_t2ij
      integer irp_j, j, max_j, irp_jk, add_t2jk
      integer irp_k, k, max_k, irp_ik, add_t2ik, irp_ijk
      double precision e5staaa

c INTERNAL VARIABLES
      integer max_a, max_bc
      integer ioff_jk, ioff_ik, ioff_ij
      integer add_s1, add_s2
      integer start_scr, size_scr
      integer itmp
      integer list_abij, list_abci
      integer off_iltj, off_altb
      integer dist_altb, dist_ia
      double precision dtmp

c PARAMETERS
      double precision zero, one, minus
      parameter (zero=0.0d0, one=1.0d0, minus=-1.0d0)

c EXTERNAL FUNCTIONS
      double precision ddot

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         print *, '@PROCESS_T3AAA: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = 1
      end if
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(2.lt.ispin)) then
         print *, '@PROCESS_T3AAA: Assertion failed.'
         print *, '   ispin = ',ispin
         itmp = 1
      end if
c   o assert max_abc > 0
      if (max_abc.lt.1) then
         print *, '@PROCESS_T3AAA: Assertion failed.'
         print *, '   max_abc = ',max_abc
         itmp = 1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@PROCESS_T3AAA: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         print *, '   j = ',j,', max_j = ',max_j
         print *, '   k = ',k,', max_k = ',max_k
         itmp = 1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_ij .lt.1).or.(nirrep.lt.irp_ij ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ik .lt.1).or.(nirrep.lt.irp_ik ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         print *, '@PROCESS_T3AAA: Assertion failed.'
         print *, '   irp_i   = ',irp_i
         print *, '   irp_j   = ',irp_j
         print *, '   irp_k   = ',irp_k
         print *, '   irp_ij  = ',irp_ij
         print *, '   irp_jk  = ',irp_jk
         print *, '   irp_ik  = ',irp_ik
         print *, '   irp_ijk = ',irp_ijk
         itmp = 1
      end if
      if (itmp.ne.0) call xcc_errex
#endif /* _ASSERT */

c ----------------------------------------------------------------------

#ifdef _DEBUG_PROCESS_T3AAA
      print *, '@',__FILE__,': T3( ',
     &         1+i,'[',irp_i,'] ',
     &         1+j,'[',irp_j,'] ',
     &         1+k,'[',irp_k,'] , * )'
#endif

c   o set handles for these dependents
      if (ispin.eq.1) then
         list_abij = h_HH_ABIJ_AltB_IltJ
         list_abci = h_HH_ABCI_AltB_CI
         dist_altb = h_DistType_AltB
         dist_ia   = h_DistType_IA
         off_iltj  = h_Offsets_IltJ
         off_altb  = h_Offsets_AltB
      else
         list_abij = h_HH_abij_altb_iltj
         list_abci = h_HH_abci_altb_ci
         dist_altb = h_DistType_altb
         dist_ia   = h_DistType_ia
         off_iltj  = h_Offsets_iltj
         off_altb  = h_Offsets_altb
      end if

c ----------------------------------------------------------------------

c   o initialize the JK, IK, and IJ offsets
c     These represent the number of occ-occ pairs to skip in order to
c     address the current occupied orbitals. This should probably be
c     passed in.

      itmp    = rshift(irp_k*(irp_k-1),1) + irp_j
      ioff_jk = off_pltp(itmp,off_iltj)
      if (irp_j.eq.irp_k) then
         ioff_jk = ioff_jk + rshift(k*(k-1),1) + j
      else
         ioff_jk = ioff_jk + ( k * max_j ) + j
      end if

      itmp    = rshift(irp_k*(irp_k-1),1) + irp_i
      ioff_ik = off_pltp(itmp,off_iltj)
      if (irp_i.eq.irp_k) then
         ioff_ik = ioff_ik + rshift(k*(k-1),1) + i
      else
         ioff_ik = ioff_ik + ( k * max_i ) + i
      end if

      itmp    = rshift(irp_j*(irp_j-1),1) + irp_i
      ioff_ij = off_pltp(itmp,off_iltj)
      if (irp_i.eq.irp_j) then
         ioff_ij = ioff_ij + rshift(j*(j-1),1) + i
      else
         ioff_ij = ioff_ij + ( j * max_i ) + i
      end if

c ----------------------------------------------------------------------

c E[5]ST ( = adj(T1)*W(-2)*T3 )

c   o make the disconnected operator T1*W(+2) at dscr(1)
c     The W(+2) integrals will start at dscr(add_s2).
      if (max_abc.gt.dim_dscr) then
         print *, '@PROCESS_T3AAA: Not enough memory.'
         call xcc_errex
      end if
      do itmp = 1, max_abc
         dscr(itmp) = zero
      end do
      add_s2 = 1 + max_abc

#ifdef _DO_E5ST_T1I
c   o   <PQ||JK> * T1(I,R)
      max_a  = vrt(irp_i,ispin)
      max_bc = irpdpd(irp_jk,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_i,ispin) + ( i * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAA: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_jk,1,2,irp_jk,list_abij)
         call cont_s3pltpltp_s2pltp_s1p(
     &           dscr(add_s2),  irp_jk, max_bc,
     &           dcore(add_s1), irp_i,  max_a,
     &           one, dscr(1),  off_altbltc,
     &           vrt(1,ispin),  dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1I */

#ifdef _DO_E5ST_T1J
c   o - <PQ||IK> * T1(J,R)
      max_a  = vrt(irp_j,ispin)
      max_bc = irpdpd(irp_ik,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_j,ispin) + ( j * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAA: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_ik,1,2,irp_ik,list_abij)
         call cont_s3pltpltp_s2pltp_s1p(
     &           dscr(add_s2),   irp_ik, max_bc,
     &           dcore(add_s1),  irp_j,  max_a,
     &           minus, dscr(1), off_altbltc,
     &           vrt(1,ispin),   dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1J */

#ifdef _DO_E5ST_T1K
c   o   <PQ||IJ> * T1(K,R)
      max_a  = vrt(irp_k,ispin)
      max_bc = irpdpd(irp_ij,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_k,ispin) + ( k * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAA: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_ij,1,2,irp_ij,list_abij)
         call cont_s3pltpltp_s2pltp_s1p(
     &           dscr(add_s2),  irp_ij, max_bc,
     &           dcore(add_s1), irp_k,  max_a,
     &           one, dscr(1),  off_altbltc,
     &           vrt(1,ispin),  dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1K */

c   o add in the E[5]ST contribution
      e5staaa = e5staaa + ddot(max_abc,t3aaa,1,dscr,1)

c ----------------------------------------------------------------------

c E[5]QfT

c   o initialize the 3-particle intermediate operator S3 at dscr(1)
      start_scr = 1 + max_abc
      size_scr  = dim_dscr - max_abc
      do itmp = 1, max_abc
         dscr(itmp) = zero
      end do

c   o The first three blocks create the disconnected operator:
c     S1vo*T2(1)

#ifdef _DO_E5QfT_S1I
c   o   <PQ{}JK> * S1(I,R)
      max_a  = vrt(irp_i,ispin)
      max_bc = irpdpd(irp_jk,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_i,ispin) + ( i * max_a )
         add_s2 = ndx_t2p1(irp_jk,ispin) + ( ioff_jk * max_bc )
         call cont_s3pltpltp_s2pltp_s1p(
     &           dcore(add_s2), irp_jk, max_bc,
     &           dcore(add_s1), irp_i,  max_a,
     &           one, dscr(1),  off_altbltc,
     &           vrt(1,ispin),  dscr(start_scr), size_scr
     &        )
      end if
#endif /* _DO_E5QfT_S1I */

#ifdef _DO_E5QfT_S1J
c   o - <PQ{}IK> * S1(J,R)
      max_a  = vrt(irp_j,ispin)
      max_bc = irpdpd(irp_ik,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_j,ispin) + ( j * max_a )
         add_s2 = ndx_t2p1(irp_ik,ispin) + ( ioff_ik * max_bc )
         call cont_s3pltpltp_s2pltp_s1p(
     &           dcore(add_s2),  irp_ik, max_bc,
     &           dcore(add_s1),  irp_j,  max_a,
     &           minus, dscr(1), off_altbltc,
     &           vrt(1,ispin),   dscr(start_scr), size_scr
     &        )
      end if
#endif /* _DO_E5QfT_S1J */

#ifdef _DO_E5QfT_S1K
c   o   <PQ{}IJ> * S1(K,R)
      max_a  = vrt(irp_k,ispin)
      max_bc = irpdpd(irp_ij,dist_altb)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_k,ispin) + ( k * max_a )
         add_s2 = ndx_t2p1(irp_ij,ispin) + ( ioff_ij * max_bc )
         call cont_s3pltpltp_s2pltp_s1p(
     &           dcore(add_s2), irp_ij, max_bc,
     &           dcore(add_s1), irp_k,  max_a,
     &           one, dscr(1),  off_altbltc,
     &           vrt(1,ispin),  dscr(start_scr), size_scr
     &        )
      end if
#endif /* _DO_E5QfT_S1K */

c   o The next three blocks create the connected operator:
c     S2vvvo*T2(1)

#ifdef _DO_E5QfT_S2I
c   o   <PQ{}EI> * T2(1)(JK,ER)
      add_s2 =   ndx_t2p1(irp_jk,ispin)
     &         + ( ioff_jk * irpdpd(irp_jk,dist_altb) )
      call cont_qf3s2_abci_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        list_abci, dist_altb,
     &        dcore(add_s2),    off_altb,
     &        i, irp_i,  max_i, irp_ij,
     &        j, irp_j,  max_j, irp_jk,
     &        k, irp_k,  max_k, irp_ijk,
     &        one, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2I */

#ifdef _DO_E5QfT_S2J
c   o - <PQ{}EJ> * T2(1)(IK,ER)
      add_s2 =   ndx_t2p1(irp_ik,ispin)
     &         + ( ioff_ik * irpdpd(irp_ik,dist_altb) )
      call cont_qf3s2_abci_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        list_abci, dist_altb,
     &        dcore(add_s2),    off_altb,
     &        j, irp_j,  max_j, irp_ij,
     &        i, irp_i,  max_i, irp_ik,
     &        k, irp_k,  max_k, irp_ijk,
     &        minus, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2J */

#ifdef _DO_E5QfT_S2K
c   o   <PQ{}EK> * T2(1)(IJ,ER)
      add_s2 =   ndx_t2p1(irp_ij,ispin)
     &         + ( ioff_ij * irpdpd(irp_ij,dist_altb) )
      call cont_qf3s2_abci_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        list_abci, dist_altb,
     &        dcore(add_s2),    off_altb,
     &        k, irp_k,  max_k, irp_ik,
     &        i, irp_i,  max_i, irp_ij,
     &        j, irp_j,  max_j, irp_ijk,
     &        one, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2K */

c   o The next three blocks create the connected operator:
c     S2ooov*T2(1)

#ifdef _DO_E5QfT_S2IJ
c   o   T2(1)(MK,PQ) * <MR{}IJ>
      call cont_qf3s2_iajk_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        off_iltj, dist_altb, dist_ia,
     &        k, irp_k, max_k, irp_ik,
     &        i, irp_i, max_i, irp_ij, ioff_ij,
     &        j, irp_j, max_j, irp_ijk,
     &        one, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2IJ */

#ifdef _DO_E5QfT_S2IK
c   o - T2(1)(MJ,PQ) * <MR{}IK>
      call cont_qf3s2_iajk_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        off_iltj, dist_altb, dist_ia,
     &        j, irp_j, max_j, irp_ij,
     &        i, irp_i, max_i, irp_ik, ioff_ik,
     &        k, irp_k, max_k, irp_ijk,
     &        minus, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2IK */

#ifdef _DO_E5QfT_S2JK
c   o   T2(1)(MI,PQ) * <MR{}JK>
      call cont_qf3s2_iajk_t2p1aa(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin,
     &        off_iltj, dist_altb, dist_ia,
     &        i, irp_i, max_i, irp_ij,
     &        j, irp_j, max_j, irp_jk, ioff_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, dscr(1), off_altbltc
     &     )
#endif /* _DO_E5QfT_S2JK */

c   o add in the E[5]QfT contribution
#ifdef _DEBUG_PROCESS_T3AAA
      if (ispin.eq.1) then
         dtmp = ddot(max_abc,t3aaa,1,dscr,1)
         print *, 'E[5]QfT, ',ispin,': ',e5qftaaa,' += ',dtmp
      end if
#endif /* _DEBUG_PROCESS_T3AAA */
      e5qftaaa = e5qftaaa + ddot(max_abc,t3aaa,1,dscr,1)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

