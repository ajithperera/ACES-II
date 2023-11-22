
      subroutine tmp_cont_eps_t3aaa_dist(dcore,start_t3,
     &                                   eps_ijk,irp_ijk,ispin)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), eps_ijk
      integer start_t3, irp_ijk, ispin

c INTERNAL VARIABLES
      integer a, irp_a, max_a,         to_a, start_a
      integer b, irp_b, max_b, from_b, to_b, start_b
      integer c, irp_c, max_c, from_c,       start_c
      integer irp_bc
      logical do_irp_b
      integer ndx
      double precision dtmp
      double precision eps_ijk_c
      double precision eps_ijk_bc
      double precision eps_ijk_abc

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

      ndx = 0

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

      dtmp = dcore(start_t3+ndx) * ( 1.0d0 / eps_ijk_abc )
         write(*,'(a,f20.10,a,f20.10,a)')
     &      'trps14.F : e4trip += ',eps_ijk_abc,' * (',dtmp,')^2'

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

      return
      end

