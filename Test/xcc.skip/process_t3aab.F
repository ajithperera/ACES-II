
c This routine processes single distributions of T3AAB amplitudes.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(dim_dscr) : scratch array for doubles
c    integer          ispin1 : primary spin case to handle
c    double precision t3aab(max_abc) : T3(IJk,*) distribution
c    integer off_a1ltb1c2(1+nirrep*nirrep) : irrep offsets for T3AAB/BBA

c    [ orbital info ]

c    double precision e5staab : the 5-th order energy correction from T3AAB
c                               <0|adj(T1)*W(-2)*T3AAB|0>

#include "disttype.h"
#include "listnum.h"

c#define _DEBUG_PROCESS_T3AAB
#define _DO_E5ST_T1I
#define _DO_E5ST_T1J
#define _DO_E5ST_T1k
c#define _DO_E5QfT_S1I
c#define _DO_E5QfT_S1J
c#define _DO_E5QfT_S1k

      subroutine process_t3aab(
     &              dcore,dscr,dim_dscr,
     &              ispin1,
     &              t3aab,max_abc,off_a1ltb1c2,
     &              i,irp_i,max_i,irp_ij,
     &              j,irp_j,max_j,irp_jk,
     &              k,irp_k,max_k,irp_ik,irp_ijk,
     &              e5staab
     &           )
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*), t3aab(*)
      integer ispin1, max_abc, off_a1ltb1c2(*), dim_dscr
      integer irp_i, i, max_i, irp_ij
      integer irp_j, j, max_j, irp_jk
      integer irp_k, k, max_k, irp_ik, irp_ijk
      double precision e5staab

c INTERNAL VARIABLES
      integer ispin2
      integer max_a, bc, max_bc
      integer ioff_jk, ioff_ik, ioff_ij
      integer add_s1, add_s2, start_scr
      integer itmp
      integer list_gs_abij, list_hh_abij
      integer off_i1ltj1, off_i1j2, off_i2j1
      integer dist_a1ltb1, dist_a1b2

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
         print *, '@PROCESS_T3AAB: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = 1
      end if
c   o assert ispin1 is 1 or 2
      if ((ispin1.lt.1).or.(2.lt.ispin1)) then
         print *, '@PROCESS_T3AAB: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         itmp = 1
      end if
c   o assert max_abc > 0
      if (max_abc.lt.1) then
         print *, '@PROCESS_T3AAB: Assertion failed.'
         print *, '   max_abc = ',max_abc
         itmp = 1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@PROCESS_T3AAB: Assertion failed.'
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
         print *, '@PROCESS_T3AAB: Assertion failed.'
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

#ifdef _DEBUG_PROCESS_T3AAB
      print *, '@',__FILE__,': T3( ',
     &         1+i,'[',irp_i,'] ',
     &         1+j,'[',irp_j,'] ',
     &         1+k,'[',irp_k,'] , * )'
#endif

c   o set handles for these dependents
      if (ispin1.eq.1) then
         ispin2 = 2
         list_hh_abij = h_HH_ABIJ_AltB_IltJ
         dist_a1ltb1 = h_DistType_AltB
         off_i1ltj1  = h_Offsets_IltJ
         off_i1j2    = h_Offsets_Ij
         off_i2j1    = h_Offsets_iJ
c         off_altb  = h_Offsets_AltB
      else
         ispin2 = 1
         list_hh_abij = h_HH_abij_altb_iltj
         dist_a1ltb1 = h_DistType_altb
         off_i1ltj1  = h_Offsets_iltj
         off_i1j2    = h_Offsets_iJ
         off_i2j1    = h_Offsets_Ij
      end if
      list_gs_abij = h_GS_AbIj_Ab_Ij
      dist_a1b2 = h_DistType_Ab

c ----------------------------------------------------------------------

c   o initialize the Jk, Ik, and IJ offsets
c     These represent the number of occ-occ pairs to skip in order to
c     address the current occupied orbitals. This should probably be
c     passed in. Note also that these are for addressing integral lists
c     stored <alpha-beta|alpha-beta>.

      if (ispin1.eq.1) then
         itmp    = nirrep*(irp_k-1) + irp_j
         ioff_jk = off_pq(itmp,h_Offsets_Ij) + ( k * max_j ) + j
         itmp    = nirrep*(irp_k-1) + irp_i
         ioff_ik = off_pq(itmp,h_Offsets_Ij) + ( k * max_i ) + i
      else
         itmp    = nirrep*(irp_j-1) + irp_k
         ioff_jk = off_pq(itmp,h_Offsets_Ij) + ( j * max_k ) + k
         itmp    = nirrep*(irp_i-1) + irp_k
         ioff_ik = off_pq(itmp,h_Offsets_Ij) + ( i * max_k ) + k
      end if

      itmp    = rshift(irp_j*(irp_j-1),1) + irp_i
      ioff_ij = off_pltp(itmp,off_i1ltj1)
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
         print *, '@PROCESS_T3AAB: Not enough memory.'
         call xcc_errex
      end if
      do itmp = 1, max_abc
         dscr(itmp) = zero
      end do
      add_s2 = 1 + max_abc

#ifdef _DO_E5ST_T1I
c   o   <Pq|Jk> * T1(I,R)
      max_a  = vrt(irp_i,ispin1)
      max_bc = irpdpd(irp_jk,dist_a1b2)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_i,ispin1) + ( i * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAB: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_jk,1,2,irp_jk,list_gs_abij)
         if (ispin1.eq.2) then
            call transp_dist(dscr(add_s2),vrt(1,1),vrt(1,2),
     &                       irp_jk,1,.false.)
         end if
         call cont_s3pltpq_s2pq_s1p(
     &           dscr(add_s2),  irp_jk, max_bc,
     &           dcore(add_s1), irp_i,  max_a,
     &           one, dscr(1),  off_a1ltb1c2,
     &           vrt(1,ispin1), vrt(1,ispin2),
     &           dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1I */

#ifdef _DO_E5ST_T1J
c   o - <Pq|Ik> * T1(J,R)
      max_a  = vrt(irp_j,ispin1)
      max_bc = irpdpd(irp_ik,dist_a1b2)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_j,ispin1) + ( j * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAB: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_ik,1,2,irp_ik,list_gs_abij)
         if (ispin1.eq.2) then
            call transp_dist(dscr(add_s2),vrt(1,1),vrt(1,2),
     &                       irp_ik,1,.false.)
         end if
         call cont_s3pltpq_s2pq_s1p(
     &           dscr(add_s2),   irp_ik, max_bc,
     &           dcore(add_s1),  irp_j,  max_a,
     &           minus, dscr(1), off_a1ltb1c2,
     &           vrt(1,ispin1),  vrt(1,ispin2),
     &           dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1J */

#ifdef _DO_E5ST_T1k
c   o   <PQ||IJ> * T1(k,r)
      max_a  = vrt(irp_k,ispin2)
      max_bc = irpdpd(irp_ij,dist_a1ltb1)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_t1(irp_k,ispin2) + ( k * max_a )
         start_scr = add_s2 + max_bc
         if (start_scr.gt.dim_dscr+1) then
            print *, '@PROCESS_T3AAB: Not enough memory.'
            call xcc_errex
         end if
         call getlst(dscr(add_s2),1+ioff_ij,1,2,irp_ij,list_hh_abij)
         call cont_s3pltpq_s2pltp_s1q(
     &           dscr(add_s2),  irp_ij, max_bc,
     &           dcore(add_s1), irp_k,  max_a,
     &           one, dscr(1),  off_a1ltb1c2,
     &           vrt(1,ispin1), dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5ST_T1k */

c   o add in the E[5]ST contribution
      e5staab = e5staab + ddot(max_abc,t3aab,1,dscr,1)

c ----------------------------------------------------------------------

c E[5]QfT

c   o initialize the 3-particle intermediate operator S3 at dscr(1)
      start_scr = 1 + max_abc
      do itmp = 1, max_abc
         dscr(itmp) = zero
      end do

c   o The first three blocks create the disconnected operator:
c     S1vo*T2(1)

#ifdef _DO_E5QfT_S1k
c   o   <PQ{}IJ> * S1(k,r) ( <- the others muck with start_scr )
      max_a  = vrt(irp_k,ispin2)
      max_bc = irpdpd(irp_ij,dist_a1ltb1)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_k,ispin2) + ( k * max_a )
         add_s2 = ndx_t2p1(irp_ij,ispin1) + ( ioff_ij * max_bc )
         call cont_s3pltpq_s2pltp_s1q(
     &           dcore(add_s2), irp_ij, max_bc,
     &           dcore(add_s1), irp_k,  max_a,
     &           one, dscr(1),  off_a1ltb1c2,
     &           vrt(1,ispin1), dscr(start_scr), 1+dim_dscr-start_scr
     &        )
      end if
#endif /* _DO_E5QfT_S1k */

#ifdef _DO_E5QfT_S1J
c   o - <Pq{}Ik> * S1(J,R)
      max_a  = vrt(irp_j,ispin1)
      max_bc = irpdpd(irp_ik,dist_a1b2)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_j,ispin1) + ( j * max_a )
         if (ispin1.eq.1) then
            add_s2 = ndx_t2p1(irp_ik,3) + ( ioff_ik * max_bc )
            call cont_s3pltpq_s2pq_s1p(
     &              dcore(add_s2),  irp_ik, max_bc,
     &              dcore(add_s1),  irp_j,  max_a,
     &              minus, dscr(1), off_a1ltb1c2,
     &              vrt(1,ispin1),  vrt(1,ispin2),
     &              dscr(start_scr), 1+dim_dscr-start_scr
     &           )
         else
            add_s2 = max_abc + 1
            start_scr = add_s2 + max_bc
            if (start_scr.gt.dim_dscr+1) then
               print *, '@PROCESS_T3AAB: Not enough memory.'
               call xcc_errex
            end if
            itmp = ndx_t2p1(irp_ik,3) + ( ioff_ik * max_bc )
            do bc = 0, max_bc-1
               dscr(add_s2+bc) = dcore(itmp+bc)
            end do
            call transp_dist(dscr(add_s2),vrt(1,1),vrt(1,2),
     &                       irp_ik,1,.false.)
            call cont_s3pltpq_s2pq_s1p(
     &              dscr(add_s2),   irp_ik, max_bc,
     &              dcore(add_s1),  irp_j,  max_a,
     &              minus, dscr(1), off_a1ltb1c2,
     &              vrt(1,ispin1),  vrt(1,ispin2),
     &              dscr(start_scr), 1+dim_dscr-start_scr
     &           )
         end if
      end if
#endif /* _DO_E5QfT_S1J */

#ifdef _DO_E5QfT_S1I
c   o   <Pq{}Jk> * S1(I,R)
      max_a  = vrt(irp_i,ispin1)
      max_bc = irpdpd(irp_jk,dist_a1b2)
      if ((max_a.ne.0).and.(max_bc.ne.0)) then
         add_s1 = ndx_s1vo(irp_i,ispin1) + ( i * max_a )
         if (ispin1.eq.1) then
            add_s2 = ndx_t2p1(irp_jk,3) + ( ioff_jk * max_bc )
            call cont_s3pltpq_s2pq_s1p(
     &              dcore(add_s2), irp_jk, max_bc,
     &              dcore(add_s1), irp_i,  max_a,
     &              one, dscr(1),  off_a1ltb1c2,
     &              vrt(1,ispin1), vrt(1,ispin2),
     &              dscr(start_scr), 1+dim_dscr-start_scr
     &           )
         else
            add_s2 = max_abc + 1
            start_scr = add_s2 + max_bc
            if (start_scr.gt.dim_dscr+1) then
               print *, '@PROCESS_T3AAB: Not enough memory.'
               call xcc_errex
            end if
            itmp = ndx_t2p1(irp_jk,3) + ( ioff_jk * max_bc )
            do bc = 0, max_bc-1
               dscr(add_s2+bc) = dcore(itmp+bc)
            end do
            call transp_dist(dscr(add_s2),vrt(1,1),vrt(1,2),
     &                       irp_jk,1,.false.)
            call cont_s3pltpq_s2pq_s1p(
     &              dscr(add_s2),  irp_jk, max_bc,
     &              dcore(add_s1), irp_i,  max_a,
     &              one, dscr(1),  off_a1ltb1c2,
     &              vrt(1,ispin1), vrt(1,ispin2),
     &              dscr(start_scr), 1+dim_dscr-start_scr
     &           )
         end if
      end if
#endif /* _DO_E5QfT_S1I */

c   o The next three blocks create the connected operator:
c     S2vvvo*T2(1)

cc   o <Pq{}Ek> * T2(1)(IJ,ER)
cc   o <pQ{}eI> * T2(1)(kJ,eR)
cc   o <pQ{}eJ> * T2(1)(kI,eR)
cc   o <PQ{}EI> * T2(1)(Jk,Er)
cc   o <PQ{}EJ> * T2(1)(Ik,Er)
cc   o T2(Mk,Pq)*<MR||IJ>
cc   o T2(mI,pQ)*<mR|kJ>
cc   o T2(mJ,pQ)*<mR|kI>
cc   o T2(MI,PQ)*<Mr|Jk>
cc   o T2(MJ,PQ)*<Mr|Ik>

C#ifdef _DO_E5QfT_S2I
Cc   o <PQ{}EI> * T2(1)(JK,ER)
C      add_s2 =   ndx_t2p1(irp_jk,ispin)
C     &         + ( ioff_jk * irpdpd(irp_jk,dist_altb) )
C      call cont_qf3s2_abci_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        list_abci, dist_altb,
C     &        dcore(add_s2),    off_altb,
C     &        i, irp_i,  max_i, irp_ij,
C     &        j, irp_j,  max_j, irp_jk,
C     &        k, irp_k,  max_k, irp_ijk,
C     &        one, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2I */

C#ifdef _DO_E5QfT_S2J
Cc   o - <PQ{}EJ> * T2(1)(IK,ER)
C      add_s2 =   ndx_t2p1(irp_ik,ispin)
C     &         + ( ioff_ik * irpdpd(irp_ik,dist_altb) )
C      call cont_qf3s2_abci_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        list_abci, dist_altb,
C     &        dcore(add_s2),    off_altb,
C     &        j, irp_j,  max_j, irp_ij,
C     &        i, irp_i,  max_i, irp_ik,
C     &        k, irp_k,  max_k, irp_ijk,
C     &        minus, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2J */

C#ifdef _DO_E5QfT_S2K
Cc   o <PQ{}EK> * T2(1)(IJ,ER)
C      add_s2 =   ndx_t2p1(irp_ij,ispin)
C     &         + ( ioff_ij * irpdpd(irp_ij,dist_altb) )
C      call cont_qf3s2_abci_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        list_abci, dist_altb,
C     &        dcore(add_s2),    off_altb,
C     &        k, irp_k,  max_k, irp_ik,
C     &        i, irp_i,  max_i, irp_ij,
C     &        j, irp_j,  max_j, irp_ijk,
C     &        one, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2K */

c   o The next three blocks create the connected operator:
c     S2ooov*T2(1)

C#ifdef _DO_E5QfT_S2IJ
Cc   o <MR{}IJ> * T2(1)(MK,PQ)
C      call cont_qf3s2_iajk_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        off_iltj, dist_altb, dist_ia,
C     &        k, irp_k, max_k, irp_ik,
C     &        i, irp_i, max_i, irp_ij, ioff_ij,
C     &        j, irp_j, max_j, irp_ijk,
C     &        one, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2IJ */

C#ifdef _DO_E5QfT_S2IK
Cc   o - <MR{}IK> * T2(1)(MJ,PQ)
C      call cont_qf3s2_iajk_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        off_iltj, dist_altb, dist_ia,
C     &        j, irp_j, max_j, irp_ij,
C     &        i, irp_i, max_i, irp_ik, ioff_ik,
C     &        k, irp_k, max_k, irp_ijk,
C     &        minus, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2IK */

C#ifdef _DO_E5QfT_S2JK
Cc   o <MR{}JK> * T2(1)(MI,PQ)
C      call cont_qf3s2_iajk_t2p1aa(
C     &        dcore, dscr(start_scr), size_scr,
C     &        ispin,
C     &        off_iltj, dist_altb, dist_ia,
C     &        i, irp_i, max_i, irp_ij,
C     &        j, irp_j, max_j, irp_jk, ioff_jk,
C     &        k, irp_k, max_k, irp_ijk,
C     &        one, dscr(1), off_altbltc
C     &     )
C#endif /* _DO_E5QfT_S2JK */

c   o add in the E[5]QfT contribution
#ifdef _DEBUG_PROCESS_T3AAB
      dtmp = ddot(max_abc,t3aab,1,dscr,1)
      print *, 'E[5]QfT, ',ispin,': ',e5qftaab,' += ',dtmp
#endif /* _DEBUG_PROCESS_T3AAB */
      e5qftaab = e5qftaab + ddot(max_abc,t3aab,1,dscr,1)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

