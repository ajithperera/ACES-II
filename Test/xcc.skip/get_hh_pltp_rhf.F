
c This routine loads lists and antisymmetrizes them in-place. It DOES NOT
c stop if there is insufficient room to load one unpacked distribution.
c This means the caller must handle the error.
c EXAMPLE:
c    <Ab|Ci> -> <AB||CI> (A<B) for RHF.

c INPUT
c    double precision dest(*) : the destination array

c    integer dim_dest    : the size of the destination array
c    integer list1       : the index of the first list to retrieve
c    integer nlists      : the number of distributions to retrieve (request)
c    integer irp_dist    : the symmetry of the distributions
c    integer listnum     : the ACES list number
c    integer pq_size     : the size of the full   <pq| distribution
c    integer pltp_size   : the size of the packed <pq| distribution
c    integer orb(nirrep) : the orbital occupation vector
c    integer offsets(nirrep*nirrep) : the matrix/array of offsets
c                                     to address blocks of unpacked irreps
c                                     (for any irp_dist)

c OUTPUT
c    integer nlists : the actual number of lists retrieved

c ----------------------------------------------------------------------

      subroutine get_hh_pltp_rhf(dest,dim_dest,list1,nlists,
     &                           irp_dist,listnum,pq_size,pltp_size,
     &                           orb,offsets)
      implicit none

c ARGUMENT LIST
      double precision dest(*)
      integer dim_dest, list1, nlists, irp_dist, listnum
      integer pq_size, pltp_size, orb(*), offsets(*)

c INTERNAL VARIABLES
      integer add_dest
      integer list, ioff
      integer dim_left, loaded

c COMMON BLOCKS
#include "syminf.com" /* for dirprd and nirrep */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      ioff = 0
c   o assert dim_dest is a natural number
      if (dim_dest.lt.1) then
         print *, '@GET_HH_PLTP_RHF: Assertion failed.'
         print *, '   dim_dest = ',dim_dest
         ioff = 1
      end if
c   o assert the getlst arguments are valid
      if ((list1.lt.1).or.(nlists.lt.0).or.(listnum.lt.1)) then
         print *, '@GET_HH_PLTP_RHF: Assertion failed.'
         print *, '   list1   = ',list1
         print *, '   nlists  = ',nlists
         print *, '   listnum = ',listnum
         ioff = 1
      end if
c   o assert irp_dist is in [1,nirrep]
      if ((irp_dist.lt.1).or.(nirrep.lt.irp_dist)) then
         print *, '@GET_HH_PLTP_RHF: Assertion failed.'
         print *, '   irp_dist = ',irp_dist
         ioff = 1
      end if
c   o assert distribution sizes are non-negative
      if ((pq_size.lt.0).or.(pltp_size.lt.0)) then
         print *, '@GET_HH_PLTP_RHF: Assertion failed.'
         print *, '   pq_size   = ',pq_size
         print *, '   pltp_size = ',pltp_size
         ioff = 1
      end if
      if (ioff.ne.0) call xcc_errex
#endif /* _ASSERT */

      loaded = 0

      if (nlists.gt.0) then

         add_dest = 1
         dim_left = dim_dest

         do list = list1, (list1+nlists)-1
         if (dim_left.ge.pq_size) then

            call getlst(dest(add_dest),list,1,2,irp_dist,listnum)
            call antisym_pp_pltp_in_situ(dest(add_dest),pq_size,
     &                                   orb,1,irp_dist,offsets)
            add_dest = add_dest + pltp_size
            dim_left = dim_left - pltp_size
            loaded   = loaded + 1

c        end if (dim_left.ge.pq_size)
         end if
c        end do list = list1, (list1+nlists)-1
         end do

c     end if (nlists.gt.0)
      end if

      nlists = loaded

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

