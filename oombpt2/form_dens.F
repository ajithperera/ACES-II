
      subroutine form_dens(dens,totlen,spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
C Input variables
      integer totlen,spin
C Output variables
      double precision dens(totlen)
C Local variables
      integer idens,irrep,mos,ii,occ
      double precision one
      data one /1.0d0/

      idens=1
      do irrep = 1,nirrep
        occ = pop_full(irrep,spin)
        mos=occ+vrt_full(irrep,spin)
        if (mos .gt. 0) then
          call zero(dens(idens),mos*mos)
          do ii = 1,occ
            dens(idens+(ii-1)*mos + ii-1) = one
          end do
        endif
        idens = idens + mos*mos
      end do

      return
      end
