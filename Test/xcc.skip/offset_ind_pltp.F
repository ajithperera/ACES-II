
c This routine returns the individual offsets from element (1,1) of
c dirprd(irp_p,irp_q) needed to address any combination of symmetry-packed
c orbitals of a 2-particle operator.
c EXAMPLE: If we want to jump to (irp_i=3,irp_j=3) then offsets(6)
c          would equal the number of occupied pairs from irreps 1 and 2
c          you would need skip from element 1 of total irrep 1
c          provided pop()=pop(1,ispin).
c NOTE: A concrete example of this follows. The first element in the
c       block of orbital pairs with total symmetry of 1 would need no
c       offset, so offsets(1) is 0. The first element in the block of
c       orbital pairs with total symmetry of 2 would also need no offset
c       provided you had the address of LIST(irp_p=1,irp_q=2); therefore,
c       offsets(2) is also 0. offsets(7), may not be 0 since
c       irp_p=2,irp_q=3 comes before irp_p=1,irp_q=4.

c INPUT
c    integer pop(nirrep) : the symmetry population vector

c OUTPUT
c    integer offsets(nirrep*(nirrep+1)/2) :
c                         the offsets for irp_p and irp_q
c                         in upper-triangular packed form
c                         addressable by irp_p+(irp_q*(irp_q-1)/2)

      subroutine offset_ind_pltp(pop,offsets)
      implicit none

c ARGUMENT LIST
      integer pop(*), offsets(*)

c INTERNAL VARIABLES
      integer irp_p, irp_q, irp_pq
      integer p,     q,     ioff,   index

c COMMON BLOCKS
#include "syminf.com" /* sets nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      index = 1
      do irp_q = 1, nirrep
         do irp_p = 1, irp_q
            irp_pq = dirprd(irp_p,irp_q)
            if (irp_p.eq.irp_q) then
               if (irp_q.ne.1) then
                  ioff = 0
                  do q = 1, irp_q-1
                     ioff = ioff + ( pop(q) * ( pop(q) - 1 ) / 2 )
                  end do
                  offsets(index) = ioff
               else
                  offsets(index) = 0
               end if
            else
               ioff = 0
               do q = 1, irp_q-1
                  p = dirprd(q,irp_pq)
                  if (p.lt.q) then
                     ioff = ioff + ( pop(p) * pop(q) )
                  end if
               end do
               offsets(index) = ioff
            end if
            index = index + 1
c        end do irp_p = 1, irp_q
         end do
c     end do irp_q = 1, nirrep
      end do

c#define _DEBUG_OFFSET_IND_PLTP
#ifdef _DEBUG_OFFSET_IND_PLTP
      print *, '@OFFSET_IND_PLTP: For the population vector'
      print *, '                     ',(pop(p),p=1,nirrep)
      print *, '                  the p<q individual offset vector ',
     &           '(upper triangle) is'
      print *, '                     ',
     &           (offsets(p),p=1,(nirrep*(nirrep+1)/2))
#endif /* _DEBUG_OFFSET_IND_PLTP */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

