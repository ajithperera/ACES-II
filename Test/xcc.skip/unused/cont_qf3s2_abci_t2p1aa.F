
c This routine contracts <AB{}EI> from QfT with T2(1)(JK,EC) for all
c A, B, and C under a given I, J, and K.

c NOTE:
c    Orbitals I, J, and K may be arbitrarily ordered. It is up to the
c calling routine to scale the contributions by -1 (via s3_factor)
c if J<I<K.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_CONT_QF3S2_ABCI_T2P1AA

      subroutine cont_qf3s2_abci_t2p1aa(
     &              dcore, dscr, dim_dscr,
     &              ispin,
     &              list_abci, dist_altb,
     &              t2p1jk,    off_altb,
     &              i, irp_i,  max_i, irp_ij,
     &              j, irp_j,  max_j, irp_jk,
     &              k, irp_k,  max_k, irp_ijk,
     &              s3_factor, s3,    off_altbltc
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr), t2p1jk(*), s3(*)
      integer ispin, list_abci, dist_altb, off_altb, off_altbltc(*)
      integer i, irp_i, max_i, irp_ij
      integer j, irp_j, max_j, irp_jk
      integer k, irp_k, max_k, irp_ijk
      double precision s3_factor

c INTERNAL VARIABLES
      integer     irp_pq, max_pq
      integer     irp_e,  max_e, from_e,       num_e
      integer     irp_r,  max_r, from_r, to_r, num_r
      logical do_irp_r
      integer itmp, size_scr, ioff
      integer max_dists, start_t2p1, ndx
      double precision dtmp_s3_factor

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Assertion failed.'
         write(*,*) '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(2.lt.ispin)) then
         write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Assertion failed.'
         write(*,*) '   ispin = ',ispin
         itmp = -1
      end if
c   o assert something about dist_altb, list_abci, off_altb?
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Assertion failed.'
         write(*,*) '   i = ',i,', max_i = ',max_i
         write(*,*) '   j = ',j,', max_j = ',max_j
         write(*,*) '   k = ',k,', max_k = ',max_k
         itmp = -1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_ij .lt.1).or.(nirrep.lt.irp_ij ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Assertion failed.'
         write(*,*) '   irp_i   = ',irp_i
         write(*,*) '   irp_j   = ',irp_j
         write(*,*) '   irp_k   = ',irp_k
         write(*,*) '   irp_ij  = ',irp_ij
         write(*,*) '   irp_jk  = ',irp_jk
         write(*,*) '   irp_ijk = ',irp_ijk
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
c   o warn if s3_factor = 0
      if (s3_factor.eq.0.0d0) then
         write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Warning.'
         write(*,*) '   (s3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (s3_factor.eq.0.0d0) return

c ----------------------------------------------------------------------

      dtmp_s3_factor = -s3_factor

c The plan is:
c   o loop over irreps of R in T2(1)(JK,ER)
c      o loop over orbital E
c         - load as many <PQ{}EI> lists as possible (all P<Q for each E)
c         o loop over orbital R in T2(1)(JK,ER)
c            - load/unpack T2(1)(JK,ER) (all E for each R)
c            - contract <PQ{}EI> with T2(1)(JK,ER) and add into S3

c   o loop over irreps of R in T2(1)(JK,ER)
      do irp_r = 1, nirrep
c      o check for T2 amplitudes
         irp_e = dirprd(irp_r,irp_jk)
         max_r = vrt(irp_r,ispin)
         if (irp_e.eq.irp_r) then
            max_e = max_r
            do_irp_r = (max_e.gt.1)
         else
            max_e = vrt(irp_e,ispin)
            do_irp_r = ((max_e.ne.0).and.(max_r.ne.0))
         end if
         if (do_irp_r) then
c         o check for intermediates
            irp_pq = dirprd(irp_e,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            do_irp_r = (max_pq.ne.0)
         end if
      if (do_irp_r) then

c      o initialize values for these irreps
c        size_scr : scratch space minus at least one S3(IJK,PQR) product
c                   matrix and T2(1)(JK,ER) vector
         size_scr = dim_dscr - ( max_pq + max_e )
         itmp = mod(size_scr,max_pq)
         max_dists = ( size_scr - itmp ) / max_pq
         if (max_dists.eq.0) then
            write(*,*) '@CONT_QF3S2_ABCI_T2P1AA: Not enough memory.'
            write(*,*) '                         At least ',
     &                 lshift((max_pq-size_scr),3),
     &                 ' more bytes are needed.'
            call xcc_errex
         end if

c      o loop over orbital E in <PQ{}EI>
         from_e = 0
         do while (from_e.lt.max_e)

c         o reset the T2(1)(JK,ER) pointer
c           This is the scratch address of the temporary factor.
            num_e = min(max_dists,(max_e-from_e))
            itmp  = max_pq * num_e
            start_t2p1 = 1        + itmp
            size_scr   = dim_dscr - itmp

c         o get as many <PQ{}EI> lists as possible
            call make_qf3_hh_abci(
     &              dcore, dscr(start_t2p1), size_scr,
     &              ispin, dscr(1),
     &              from_e, from_e-1+num_e, irp_e,  max_e,  irp_pq,
     &              i,                      irp_i,  max_i,
     &                                      irp_pq, max_pq
     &           )

c         o loop over orbital R in T2(1)(JK,ER)
            from_r = 0
            do while (from_r.lt.max_r)

c            o We know there is enough memory to process at least one R orbital
c              (from size_scr), so we do not need to check for memory.
c              ndx : (temp) a T2(1), S3 column pair (factor, product)
               ndx   = max_e + max_pq
               itmp  = 1 + dim_dscr - start_t2p1
               ioff  = mod(itmp,ndx)
c              to_r : (temp) the number of column pairs that can fit in dscr
               to_r  = ( itmp - ioff ) / ndx
               num_r = min(to_r,(max_r-from_r))
               to_r  = from_r + num_r - 1

c            o unpack/load T2(1)(JK,ER) (all E for each R)
               call load_t2_virts(
     &                 t2p1jk, dscr(start_t2p1),
     &                 irp_e, max_e,               ispin,
     &                 irp_r, max_r, from_r, to_r, ispin, off_altb
     &              )

c ----------------------------------------------------------------------

      itmp = start_t2p1 + from_e
      ndx  = start_t2p1 + ( max_e * num_r )
      call cont_s3pltpltp_s2pltp_s2xp(
     &        dscr(1),    irp_pq, max_pq, num_e,
     &        dscr(itmp), irp_r,  max_e,  from_r, to_r, max_r,
     &        s3_factor,  s3,     off_altbltc,
     &        vrt(1,ispin), dscr(ndx), 1+dim_dscr-ndx
     &     )

c ----------------------------------------------------------------------

c         o increment and close the loop over orbital R in T2(1)(JK,ER)
               from_r = from_r + num_r
c           end do while (from_r.lt.max_r)
            end do

c      o increment and close the loop over orbital E in <PQ{}EI>
            from_e = from_e + num_e
c        end do while (from_e.lt.max_e)
         end do

c   o close the loop over irreps of R in T2(1)(JK,ER)
c     end if (do_irp_r)
      end if
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

