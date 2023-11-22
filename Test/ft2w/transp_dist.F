
c This routine transposes the two indices for EACH DISTRIBUTION
c (that's all irreps and all orbitals) across lists of distributions
c provided each distribution has the same total symmetry.

c WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c    This routine does not change the sign of the matrix elements.
c    The caller is responsible for doing that independently.

c EXAMPLE:
c    <pq|rs>   -> <qp|rs>
c    T2(IJ,PQ) -> T2(IJ,QP)

c ARGUMENTS:
c
c double dist(*) : the list(s) of distributions of irp_dist
c
c integer pop_p(*) : the symmetry population vector for index 1
c integer pop_q(*) : the symmetry population vector for index 2
c integer irp_dist : the total symmetry of the distribution(s)
c integer nlists   : the number of distributions to transpose
c
c logical pltq : a flag specifying whether the elements are packed p<q
c                [Usually .true. if pop_p(*) and pop_q(*) point to
c                 the same population vectors. This must be .false.
c                 if irp_q can be less than irp_p, e.g., <ka|ij>.]

c ----------------------------------------------------------------------

      subroutine transp_dist(
     &              dist,
     &              pop_p,    pop_q,
     &              irp_dist, nlists,
     &              pltq
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*)
      integer pop_p(*), pop_q(*), irp_dist, nlists
      logical pltq

c INTERNAL VARIABLES
      integer iadd, list, ierr, idim
      integer irp_q, pops(2,8), shifts(3,7)

c COMMON BLOCKS
#include "syminf.com"

c ----------------------------------------------------------------------

c ASSERTIONS
#ifdef _ASSERT
      ierr = 0
c   o assert nirrep is at most 8
      if (nirrep.gt.8) then
         print *, '@TRANSP_BRA: Assertion failed.'
         print *, '   nirrep = ',nirrep
         ierr = 1
      end if
c   o assert irp_dist is in [1,nirrep]
      if ((irp_dist.lt.1).or.(nirrep.lt.irp_dist)) then
         print *, '@TRANSP_BRA: Assertion failed.'
         print *, '   irp_dist = ',irp_dist
         ierr = 1
      end if
c   o assert pop_p(1:nirrep) and pop_q(1:nirrep) are non-negative
      do iadd = 1, nirrep
         if (pop_p(iadd).lt.0) then
            print *, '@TRANSP_BRA: Assertion failed.'
            print *, '   pop_p(',iadd,') = ',pop_p(iadd)
            ierr = 1
         end if
      end do
      do iadd = 1, nirrep
         if (pop_q(iadd).lt.0) then
            print *, '@TRANSP_BRA: Assertion failed.'
            print *, '   pop_q(',iadd,') = ',pop_q(iadd)
            ierr = 1
         end if
      end do
      if (ierr.ne.0) call aces_exit(ierr)
#endif /* _ASSERT */

      if ((irp_dist.eq.1).and.(pltq)) return
      if (nlists.lt.1) return

c ----------------------------------------------------------------------

c   o initialize the pops(,) array
      do irp_q = 1, 8
         pops(1,irp_q) = 0
         pops(2,irp_q) = 0
      end do

c   o fill in the symmetry info for allowed combinations
      if (pltq) then
         do irp_q = 1, nirrep
            ierr = dirprd(irp_q,irp_dist)
            if (ierr.lt.irp_q) then
               pops(1,irp_q) = pop_p(ierr)
               pops(2,irp_q) = pop_q(irp_q)
            end if
         end do
      else
         do irp_q = 1, nirrep
            ierr = dirprd(irp_q,irp_dist)
            pops(1,irp_q) = pop_p(ierr)
            pops(2,irp_q) = pop_q(irp_q)
         end do
      end if

c   o transpose the individual irreps for each list (sorting irreps later)
      iadd = 1
      do list = 1, nlists
         do irp_q = 1, nirrep
            if ((pops(1,irp_q).gt.1).and.(pops(2,irp_q).gt.1)) then
               call dmat_xpose(dist(iadd),pops(1,irp_q),pops(2,irp_q))
            end if
            iadd = iadd + ( pops(1,irp_q) * pops(2,irp_q) )
         end do
      end do

c   o sort irreps of p by circular shifting
c     (only necessary if irp_dist is not 1)
      if (irp_dist.ne.1) then

c      o pops(2,irp) is the number of q*p orbitals in irp
         do irp_q = 1, 8
            pops(2,irp_q) = pops(1,irp_q) * pops(2,irp_q)
         end do

c      o pops(1,irp) is irp_p
         do irp_q = 1, nirrep
            pops(1,irp_q) = dirprd(irp_q,irp_dist)
         end do

c      o initialize the shifts(,) array
         do irp_q = 1, 7
            shifts(1,irp_q) = 0
            shifts(2,irp_q) = 0
            shifts(3,irp_q) = 0
         end do

c      o hard code as many shifts as possible
         if (nirrep.eq.8) then
            if (8.eq.pops(1,1)) then
c            o 8-7-6-5-4-3-2-1
               shifts(1,1) = 0
               shifts(2,1) = pops(2,1) + pops(2,2)
               shifts(3,1) = pops(2,1)
               shifts(1,2) = shifts(2,1)
               shifts(2,2) = pops(2,3) + pops(2,4)
               shifts(3,2) = pops(2,3)
               shifts(1,3) = 0
               shifts(2,3) = shifts(2,1) + shifts(2,2)
               shifts(3,3) = shifts(2,1)
               shifts(1,4) = shifts(2,3)
               shifts(2,4) = pops(2,5) + pops(2,6)
               shifts(3,4) = pops(2,5)
               shifts(1,5) = shifts(2,3) + shifts(2,4)
               shifts(2,5) = pops(2,7) + pops(2,8)
               shifts(3,5) = pops(2,7)
               shifts(1,6) = shifts(2,3)
               shifts(2,6) = shifts(2,4) + shifts(2,5)
               shifts(3,6) = shifts(2,4)
               shifts(1,7) = 0
               shifts(2,7) = shifts(1,5) + shifts(2,5)
               shifts(3,7) = shifts(2,1) + shifts(2,2)
               idim = shifts(2,7)
            else
               if (4.eq.pops(1,1)) then
c               o 4-3-2-1-8-7-6-5
                  shifts(1,1) = 0
                  shifts(2,1) = pops(2,1) + pops(2,2)
                  shifts(3,1) = pops(2,1)
                  shifts(1,2) = shifts(2,1)
                  shifts(2,2) = pops(2,3) + pops(2,4)
                  shifts(3,2) = pops(2,3)
                  shifts(1,3) = 0
                  shifts(2,3) = shifts(2,1) + shifts(2,2)
                  shifts(3,3) = shifts(2,1)
                  shifts(1,4) = shifts(2,3)
                  shifts(2,4) = pops(2,5) + pops(2,6)
                  shifts(3,4) = pops(2,5)
                  shifts(1,5) = shifts(2,3) + shifts(2,4)
                  shifts(2,5) = pops(2,7) + pops(2,8)
                  shifts(3,5) = pops(2,7)
                  shifts(1,6) = shifts(2,3)
                  shifts(2,6) = shifts(2,4) + shifts(2,5)
                  shifts(3,6) = shifts(2,4)
                  idim = shifts(1,5) + shifts(2,5)
               else
                  if (6.eq.pops(1,1)) then
c                  o 6-5-8-7-2-1-4-3
                     shifts(1,1) = 0
                     shifts(2,1) = pops(2,1) + pops(2,2)
                     shifts(3,1) = pops(2,1)
                     shifts(1,2) = shifts(2,1)
                     shifts(2,2) = pops(2,3) + pops(2,4)
                     shifts(3,2) = pops(2,3)
                     shifts(1,3) = shifts(2,1) + shifts(2,2)
                     shifts(2,3) = pops(2,5) + pops(2,6)
                     shifts(3,3) = pops(2,5)
                     shifts(1,4) = shifts(1,3) + shifts(2,3)
                     shifts(2,4) = pops(2,7) + pops(2,8)
                     shifts(3,4) = pops(2,7)
                     shifts(1,5) = 0
                     shifts(2,5) = shifts(1,4) + shifts(2,4)
                     shifts(3,5) = shifts(1,3)
                     idim = shifts(2,5)
                  else
                     if (2.eq.pops(1,1)) then
c                     o 2-1-4-3-6-5-8-7
                        shifts(1,1) = 0
                        shifts(2,1) = pops(2,1) + pops(2,2)
                        shifts(3,1) = pops(2,1)
                        shifts(1,2) = shifts(2,1)
                        shifts(2,2) = pops(2,3) + pops(2,4)
                        shifts(3,2) = pops(2,3)
                        shifts(1,3) = shifts(2,1) + shifts(2,2)
                        shifts(2,3) = pops(2,5) + pops(2,6)
                        shifts(3,3) = pops(2,5)
                        shifts(1,4) = shifts(1,3) + shifts(2,3)
                        shifts(2,4) = pops(2,7) + pops(2,8)
                        shifts(3,4) = pops(2,7)
                        idim = shifts(1,4) + shifts(2,4)
                     else
                        if (7.eq.pops(1,1)) then
c                        o 7-8-5-6-3-4-1-2
                           shifts(1,1) = 0
                           shifts(3,1) = pops(2,1) + pops(2,2)
                           shifts(2,1) =   shifts(3,1)
     &                                   + pops(2,3) + pops(2,4)
                           shifts(1,2) = shifts(2,1)
                           shifts(3,2) = pops(2,5) + pops(2,6)
                           shifts(2,2) =   shifts(3,2)
     &                                   + pops(2,7) + pops(2,8)
                           shifts(1,3) = 0
                           shifts(2,3) = shifts(2,1) + shifts(2,2)
                           shifts(3,3) = shifts(2,1)
                           idim = shifts(2,3)
                        else
                           if (3.eq.pops(1,1)) then
c                           o 3-4-1-2-7-8-5-6
                              shifts(1,1) = 0
                              shifts(3,1) = pops(2,1) + pops(2,2)
                              shifts(2,1) =   shifts(3,1)
     &                                      + pops(2,3) + pops(2,4)
                              shifts(1,2) = shifts(2,1)
                              shifts(3,2) = pops(2,5) + pops(2,6)
                              shifts(2,2) =   shifts(3,2)
     &                                      + pops(2,7) + pops(2,8)
                              idim = shifts(2,1) + shifts(2,2)
                           else
                              if (5.eq.pops(1,1)) then
c                              o 5-6-7-8-1-2-3-4
                                 shifts(1,1) = 0
                                 shifts(3,1) =   pops(2,1) + pops(2,2)
     &                                         + pops(2,3) + pops(2,4)
                                 shifts(2,1) =   shifts(3,1)
     &                                         + pops(2,5) + pops(2,6)
     &                                         + pops(2,7) + pops(2,8)
                                 idim = shifts(2,1)
c                             end if (5.eq.pops(1,1))
                              end if
c                          end if (3.eq.pops(1,1))
                           end if
c                       end if (7.eq.pops(1,1))
                        end if
c                    end if (2.eq.pops(1,1))
                     end if
c                 end if (6.eq.pops(1,1))
                  end if
c              end if (4.eq.pops(1,1))
               end if
c           end if (8.eq.pops(1,1))
            end if
         else
            if (nirrep.eq.4) then
               if (4.eq.pops(1,1)) then
c               o 4-3-2-1
                  shifts(1,1) = 0
                  shifts(2,1) = pops(2,1) + pops(2,2)
                  shifts(3,1) = pops(2,1)
                  shifts(1,2) = shifts(2,1)
                  shifts(2,2) = pops(2,3) + pops(2,4)
                  shifts(3,2) = pops(2,3)
                  shifts(1,3) = 0
                  shifts(2,3) = shifts(2,1) + shifts(2,2)
                  shifts(3,3) = shifts(2,1)
                  idim = shifts(2,3)
               else
                  if (2.eq.pops(1,1)) then
c                  o 2-1-4-3
                     shifts(1,1) = 0
                     shifts(2,1) = pops(2,1) + pops(2,2)
                     shifts(3,1) = pops(2,1)
                     shifts(1,2) = shifts(2,1)
                     shifts(2,2) = pops(2,3) + pops(2,4)
                     shifts(3,2) = pops(2,3)
                     idim = shifts(2,1) + shifts(2,2)
                  else
c                  o 3-4-1-2
                     shifts(1,1) = 0
                     shifts(3,1) = pops(2,1) + pops(2,2)
                     shifts(2,1) =    shifts(3,1)
     &                              + pops(2,3) + pops(2,4)
                     idim = shifts(2,1)
c                 end if (2.eq.pops(1,1))
                  end if
c              end if (4.eq.pops(1,1))
               end if
            else
               if (nirrep.eq.2) then
c               o 2-1
                  shifts(1,1) = 0
                  shifts(2,1) = pops(2,1) + pops(2,2)
                  shifts(3,1) = pops(2,1)
                  idim = shifts(2,1)
               else
                  print *, '@TRANSP_BRA: Code this point group.'
                  call aces_exit(1)
c              end if (nirrep.eq.2)
               end if
c           end if (nirrep.eq.4)
            end if
c        end if (nirrep.eq.8)
         end if

c      o compress the shifts to remove "empty" shifts
         iadd = 0
         do ierr = 1, 7
            if (iadd.eq.0) then
               if (shifts(3,ierr).eq.0) iadd = ierr
            else
               if (shifts(3,ierr).ne.0) then
                  shifts(1,iadd) = shifts(1,ierr)
                  shifts(2,iadd) = shifts(2,ierr)
                  shifts(3,iadd) = shifts(3,ierr)
                  iadd = iadd + 1
               end if
            end if
         end do
         if (iadd.ne.0) then
            do ierr = iadd, 7
               shifts(1,ierr) = 0
               shifts(2,ierr) = 0
               shifts(3,ierr) = 0
            end do
         end if

c      o shift all the lists
         iadd = 1
         do list = 1, nlists
            if (shifts(3,1).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,1)),
     &                          shifts(2,1),shifts(3,1))
            if (shifts(3,2).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,2)),
     &                          shifts(2,2),shifts(3,2))
            if (shifts(3,3).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,3)),
     &                          shifts(2,3),shifts(3,3))
            if (shifts(3,4).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,4)),
     &                          shifts(2,4),shifts(3,4))
            if (shifts(3,5).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,5)),
     &                          shifts(2,5),shifts(3,5))
            if (shifts(3,6).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,6)),
     &                          shifts(2,6),shifts(3,6))
            if (shifts(3,7).ne.0) then
               call dvec_cshift(dist(iadd+shifts(1,7)),
     &                          shifts(2,7),shifts(3,7))
c           end if (shifts(3,7).ne.0)
            end if
c           end if (shifts(3,6).ne.0)
            end if
c           end if (shifts(3,5).ne.0)
            end if
c           end if (shifts(3,4).ne.0)
            end if
c           end if (shifts(3,3).ne.0)
            end if
c           end if (shifts(3,2).ne.0)
            end if
c           end if (shifts(3,1).ne.0)
            end if
            iadd = iadd + idim
         end do

c     end if (irp_dist.ne.1)
      end if

c ----------------------------------------------------------------------

      return
      end

