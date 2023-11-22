
c This routine returns the individual offsets from element (1,1) of
c dirprd(irp_1,irp_2) needed to address any combination of symmetry-UNpacked
c orbitals of a 2-particle operator.
c EXAMPLE: If we want to jump to (irp_I,irp_j) then offsets(irp_I,irp_j)
c          would equal the number of occupied pairs from irreps with a
c          total symmetry of dirprd(irp_I,irp_j) you would need to skip from
c          element 1 of the block (dirprd(irp_I,irp_j),1) provided
c          pop1()=pop(*,1) and pop2()=pop(*,2).

c INPUT
c    integer pop1(nirrep) : the symmetry population vector of index 1
c    integer pop2(nirrep) : the symmetry population vector of index 2

c OUTPUT
c    integer offsets(nirrep,nirrep) : the offsets for irp_1 and irp_2
c                                     Internally, this is treated as
c                                     a 1-dimensional array. Since
c                                     the syminf common block holds nirrep.

      subroutine offset_ind_pq(pop1,pop2,offsets)
      implicit none

c ARGUMENT LIST
      integer pop1(*), pop2(*), offsets(*)

c INTERNAL VARIABLES
      integer irp_p, irp_q, irp_pq
      integer irp_1, irp_2, ioff,   index

c COMMON BLOCKS
#include "syminf.com" /* sets nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o the first blocks (irp_q=1) never have an offset
      do irp_p = 1, nirrep
         offsets(irp_p) = 0
      end do

      if (nirrep.ne.1) then
         index = 1 + nirrep
         do irp_q = 2, nirrep
            do irp_p = 1, nirrep
               irp_pq = dirprd(irp_p,irp_q)
               ioff = 0
               do irp_2 = 1, irp_q-1
                  irp_1 = dirprd(irp_2,irp_pq)
                  ioff = ioff + ( pop1(irp_1) * pop2(irp_2) )
               end do
               offsets(index) = ioff
               index = index + 1
c           end do irp_p = 1, nirrep
            end do
c        end do irp_q = 2, nirrep
         end do
c     end if (nirrep.ne.1)
      end if

c#define _DEBUG_OFFSET_IND_PQ
#ifdef _DEBUG_OFFSET_IND_PQ
      print *, '@OFFSET_IND_PQ: For the population vectors'
      print *, '                   ',(pop1(irp_p),irp_p=1,nirrep)
      print *, '                   ',(pop2(irp_p),irp_p=1,nirrep)
      print *, '                the pq individual offset vector is'
      do irp_p = 1, nirrep
         print *, '                   ',
     &              (offsets(irp_p+nirrep*irp_q),irp_q=0,nirrep-1)
      end do
#endif /* _DEBUG_OFFSET_IND_PQ */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

