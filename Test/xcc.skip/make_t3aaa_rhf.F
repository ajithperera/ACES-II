
c This routine builds the T3AAA amplitudes for RHF calculations from the
c T3AAB pre-amplitudes. For a given distribution of T3(I<J<K;A<B,C) pre-
c amplitudes, the T3AAA pre-amplitudes are given by:
c    (A<B,C)-(A<C,B)+(B<C,A)

c INPUT
c    double precision dcore(*)        : address space for xcc.com
c    integer ndx_t3aab(nirrep*nirrep) : addresses of T3AAB irreps
c    double precision t3aaa(*)        : the destination array for T3(IJK,*)
c    integer off_altbltc(1+nirrep*(nirrep+1)/2) : offsets for T3AAA irreps

c    double precision dscr(*) : scratch array for doubles
c    integer dim_dscr         : the scratch dimension

c    [occupied orbital info]

c    double precision eps_ijk  : the first half of the epsilon denominator
c    double precision e4t      : the E[4]T correlation energy (not zeroed)
c    double precision t3_ovrlp : <0|adj(T3AAA)*T3AAA|0> (not zeroed)

c#define _DEBUG_MAKE_T3AAA_RHF
#define _DO_EPS

      subroutine make_t3aaa_rhf(
     &              dcore,ndx_t3aab,
     &              t3aaa,off_altbltc,
     &              dscr, dim_dscr,
     &              i,irp_i,max_i,irp_ij,
     &              j,irp_j,max_j,irp_jk,
     &              k,irp_k,max_k,irp_ik,irp_ijk,
     &              eps_ijk,e4t,t3_ovrlp
     &           )
      implicit none

c ARGUMENT LIST
      double precision dcore(*), t3aaa(*), dscr(*)
      double precision eps_ijk, e4t, t3_ovrlp
      integer ndx_t3aab(*), off_altbltc(*), dim_dscr
      integer irp_i, i, max_i, irp_ij
      integer irp_j, j, max_j, irp_jk
      integer irp_k, k, max_k, irp_ik, irp_ijk

c INTERNAL VARIABLES
      integer a, irp_a, max_a,         to_a, start_bca, add_bca, num_bc
      integer b, irp_b, max_b, from_b, to_b, start_acb, add_acb, num_ac
      integer c, irp_c, max_c, from_c,       start_abc, add_abc, num_ab
      integer irp_bc
      logical do_irp_b
      integer ndx
      integer itmp

c PARAMETERS
      integer ispin
      parameter (ispin=1)

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
         print *, '@MAKE_T3AAA_RHF: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = 1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@MAKE_T3AAA_RHF: Assertion failed.'
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
         print *, '@MAKE_T3AAA_RHF: Assertion failed.'
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

#ifdef _DEBUG_MAKE_T3AAA_RHF
      print *, '@',__FILE__,': T3( ',
     &         1+i,'[',irp_i,'] ',
     &         1+j,'[',irp_j,'] ',
     &         1+k,'[',irp_k,'] , * )'
#endif

      ndx = 1

c   o loop over irreps of virtuals in T3(IJK,ABC)
      do irp_c = 1, nirrep
         max_c = vrt(irp_c,ispin)
      if (max_c.ne.0) then
      do irp_b = 1, irp_c
         irp_bc = dirprd(irp_b,irp_c)
         irp_a  = dirprd(irp_bc,irp_ijk)
         if (irp_a.lt.irp_b) then
            max_a = vrt(irp_a,ispin)
            max_b = vrt(irp_b,ispin)
            if (irp_b.lt.irp_c) then
               do_irp_b = ((max_b.ne.0).and.(max_a.ne.0))
            else
               do_irp_b = ((max_b.ne.1).and.(max_a.ne.0))
            end if
         else
            if (irp_a.gt.irp_b) then
               do_irp_b = .false.
            else
               max_a = vrt(irp_a,ispin)
               max_b = max_a
               do_irp_b = ((max_b.gt.2).or.
     &                     ((irp_b.lt.irp_c).and.(max_b.eq.2)))
            end if
         end if
      if (do_irp_b) then

c      o set orbital boundary conditions (as offsets)
         to_a = max_a-1
         to_b = max_b-1
         if (irp_a.lt.irp_b) then
            from_b = 0
         else
            from_b = 1
         end if
         if (irp_b.lt.irp_c) then
               from_c = 0
         else
            if (irp_a.lt.irp_b) then
               from_c = 1
            else
               from_c = 2
            end if
         end if

c      o find starting dcore addresses for T3AAB irreps and skip values
c      o ABC
         if (irp_a.eq.irp_b) then
            num_ab = rshift(max_b*(max_b-1),1)
         else
            num_ab = max_a * max_b
         end if
         itmp      = nirrep*(irp_c-1) + irp_b
         start_abc = ndx_t3aab(itmp) + ( from_c * num_ab )
c      o ACB
         itmp      = nirrep*(irp_b-1) + irp_c
         start_acb = ndx_t3aab(itmp)
         if (irp_a.eq.irp_c) then
            num_ac = rshift(max_c*(max_c-1),1)
            itmp   = 1
         else
            num_ac = max_a * max_c
            itmp   = max_a * from_c
         end if
         start_acb = start_acb + ( from_b * num_ac ) + itmp
c      o BCA
         itmp      = nirrep*(irp_a-1) + irp_c
         start_bca = ndx_t3aab(itmp)
         if (irp_b.eq.irp_c) then
            num_bc = rshift(max_c*(max_c-1),1)
            itmp   = rshift(from_c*(from_c-1),1) + from_b
         else
            num_bc = max_b * max_c
            itmp   = ( max_b * from_c ) + from_b
         end if
         start_bca = start_bca + itmp

c      o loop over orbitals A<B<C
            add_bca = start_bca
         do c = from_c, max_c-1
            add_abc = start_abc
            add_acb = start_acb
            if (irp_b.eq.irp_c) to_b = c-1
         do b = from_b, to_b
            if (irp_a.eq.irp_b) to_a = b-1
         do a = 0,      to_a

c ----------------------------------------------------------------------

      t3aaa(ndx) =   dcore(add_abc+a)
     &             - dcore(add_acb+a)
     &             + dcore(add_bca+(a*num_bc))

c ----------------------------------------------------------------------

c         o increment the T3 pointer
            ndx = ndx + 1

c      o close the loop over orbitals A<B<C
c        end do a = 0, to_a
         end do
            add_abc = add_abc + to_a + 1
            add_acb = add_acb + num_ac
            add_bca = add_bca + 1
c        end do b = from_b, to_b
         end do
            start_abc = start_abc + num_ab
            if (irp_a.eq.irp_c) then
               start_acb = start_acb + c
            else
               start_acb = start_acb + max_a
            end if
            if (irp_a.eq.irp_b) then
               add_bca = add_bca + 1
            end if
c        end do c = from_c, max_c-1
         end do

c   o close the loop over irreps of virtuals in T3(IJK,ABC)
c     end if (do_irp_b)
      end if
c     end do irp_b = 1, irp_c
      end do
c     end if (max_c.ne.0)
      end if
c     end do irp_c = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _DO_EPS
c   o divide the pre-amplitudes with epsilon(IJK,ABC) and return
c     E[4]T and the T3 overlap contribution
      call cont_eps_t3aaa(dcore,t3aaa,
     &                    eps_ijk,irp_ijk,
     &                    ispin,
     &                    e4t,t3_ovrlp)
#endif /* _DO_EPS */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

