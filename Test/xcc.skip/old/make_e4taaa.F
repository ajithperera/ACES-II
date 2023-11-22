
c NOTE:
c    cont_eps_t3aaa and make_e4taaa have EXACTLY the same loop structures.
c Any debugging here must be carried over to the other file.

c This routine calculates the E[4]T contribution from
c T3(IJK,ABC)*T3(IJK,ABC)*epsilon(IJK,ABC)

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision t3aaa(*) : T3(IJK,ABC) amplitudes
c    integer          ispin    : the spin case to handle
c    double precision e4taaa   : the correlation energy
c                                [the calling routine is responsible
c                                 for passing in 0.0d0 or the actual
c                                 correlation energy]

c OUTPUT
c    double precision e4taaa : the correlation energy plus
c                              the E[4]T contribution from T3(IJK,ABC)

c#define _DEBUG_MAKE_E4TAAA

      subroutine make_e4taaa(dcore,t3aaa,ispin,e4taaa)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), t3aaa(*)
      integer ispin
      double precision e4taaa

c INTERNAL VARIABLES
      integer i, irp_i, max_i,         to_i, start_i
      integer j, irp_j, max_j, from_j, to_j, start_j
      integer k, irp_k, max_k, from_k,       start_k
      integer a, irp_a, max_a,         to_a, start_a
      integer b, irp_b, max_b, from_b, to_b, start_b
      integer c, irp_c, max_c, from_c,       start_c
      integer irp_ijk, irp_jk,   irp_bc
      logical          do_irp_j, do_irp_b
      integer ndx
      double precision eps_k
      double precision eps_jk
      double precision eps_ijk
      double precision eps_ijk_c
      double precision eps_ijk_bc
      double precision eps_ijk_abc

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _DEBUG
      ndx = 0
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(ispin.gt.2)) then
         write(*,*) '@MAKE_E4TAAA: Assertion failed.'
         write(*,*) '   ispin = ',ispin
         ndx = -1
      end if
      if (ndx.ne.0) then
         call xcc_errex
      end if
#endif /* _DEBUG */

      ndx = 1
      do irp_ijk = 1, nirrep

c   o loop over irreps of occupieds in T3(IJK,ABC)
      do irp_k = 1, nirrep
         max_k = pop(irp_k,ispin)
      if (max_k.ne.0) then
      do irp_j = 1, irp_k
         irp_jk = dirprd(irp_j,irp_k)
         irp_i  = dirprd(irp_jk,irp_ijk)
         if (irp_i.eq.irp_j) then
            max_i = pop(irp_i,ispin)
            max_j = max_i
            if (irp_j.eq.irp_k) then
               do_irp_j = (max_j.gt.2)
            else
               do_irp_j = (max_j.gt.1)
            end if
         else
            if (irp_i.lt.irp_j) then
               max_i = pop(irp_i,ispin)
               max_j = pop(irp_j,ispin)
               if (irp_j.eq.irp_k) then
                  do_irp_j = ((max_j.ne.1).and.(max_i.ne.0))
               else
                  do_irp_j = ((max_j.ne.0).and.(max_i.ne.0))
               end if
            else
               do_irp_j = .false.
            end if
         end if
      if (do_irp_j) then

c      o retrieve starting addresses once
         start_i = ndx_eval(irp_i,1,ispin)
         start_j = ndx_eval(irp_j,1,ispin)
         start_k = ndx_eval(irp_k,1,ispin)

c      o set orbital boundary conditions (as offsets)
         to_i = max_i-1
         to_j = max_j-1
         if (irp_i.eq.irp_j) then
            from_j = 1
         else
            from_j = 0
         end if
         if (irp_j.eq.irp_k) then
            if (irp_i.eq.irp_j) then
               from_k = 2
            else
               from_k = 1
            end if
         else
               from_k = 0
         end if

c      o loop over orbitals I<J<K
         do k = from_k, max_k-1
            eps_k   =          dcore(start_k + k)
            if (irp_j.eq.irp_k) to_j = k-1
         do j = from_j, to_j
            eps_jk  = eps_k  + dcore(start_j + j)
            if (irp_i.eq.irp_j) to_i = j-1
         do i = 0,      to_i
            eps_ijk = eps_jk + dcore(start_i + i)

c ----------------------------------------------------------------------

c   o loop over irreps of virtuals in T3(IJK,ABC)
      do irp_c = 1, nirrep
         max_c = vrt(irp_c,ispin)
      if (max_c.ne.0) then
      do irp_b = 1, irp_c
         irp_bc = dirprd(irp_b,irp_c)
         irp_a  = dirprd(irp_bc,irp_ijk)
         if (irp_a.eq.irp_b) then
            max_a = vrt(irp_a,ispin)
            max_b = max_a
            if (irp_b.eq.irp_c) then
               do_irp_b = (max_b.gt.2)
            else
               do_irp_b = (max_b.gt.1)
            end if
         else
            if (irp_a.lt.irp_b) then
               max_a = vrt(irp_a,ispin)
               max_b = vrt(irp_b,ispin)
               if (irp_b.eq.irp_c) then
                  do_irp_b = ((max_b.ne.1).and.(max_a.ne.0))
               else
                  do_irp_b = ((max_b.ne.0).and.(max_a.ne.0))
               end if
            else
               do_irp_b = .false.
            end if
         end if
      if (do_irp_b) then

c      o retrieve starting addresses once
         start_a = ndx_eval(irp_a,2,ispin)
         start_b = ndx_eval(irp_b,2,ispin)
         start_c = ndx_eval(irp_c,2,ispin)

c      o set orbital boundary conditions (as offsets)
         to_a = max_a-1
         to_b = max_b-1
         if (irp_a.eq.irp_b) then
            from_b = 1
         else
            from_b = 0
         end if
         if (irp_b.eq.irp_c) then
            if (irp_a.eq.irp_b) then
               from_c = 2
            else
               from_c = 1
            end if
         else
               from_c = 0
         end if

c      o loop over orbitals A<B<C
         do c = from_c, max_c-1
            eps_ijk_c   = eps_ijk    - dcore(start_c + c)
            if (irp_b.eq.irp_c) to_b = c-1
         do b = from_b, to_b
            eps_ijk_bc  = eps_ijk_c  - dcore(start_b + b)
            if (irp_a.eq.irp_b) to_a = b-1
         do a = 0,      to_a
            eps_ijk_abc = eps_ijk_bc - dcore(start_a + a)

c ----------------------------------------------------------------------

      e4taaa = e4taaa + ( t3aaa(ndx) * t3aaa(ndx) * eps_ijk_abc )
#ifdef _DEBUG_MAKE_E4TAAA
c      write(*,*) 'e4taaa += ',eps_ijk_abc,' * (',t3aaa(ndx),')^2'
      write(*,*) 'trps14.F : e4trip += ',eps_ijk_abc,
     &           ' * (',t3aaa(ndx),')^2'
#endif

c ----------------------------------------------------------------------

c         o increment the T3 pointer
            ndx = ndx + 1

c      o close the loop over orbitals A<B<C
c        end do a = 0, to_a
         end do
c        end do b = from_b, to_b
         end do
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

c      o close the loop over orbitals I<J<K
c        end do i = 0, to_i
         end do
c        end do j = from_j, to_j
         end do
c        end do k = from_k, max_k-1
         end do

c   o close the loop over irreps of occupieds in T3(IJK,ABC)
c     end if (do_irp_j)
      end if
c     end do irp_j = 1, irp_k
      end do
c     end if (max_k.ne.0)
      end if
c     end do irp_k = 1, nirrep
      end do

c     end do irp_ijk = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

