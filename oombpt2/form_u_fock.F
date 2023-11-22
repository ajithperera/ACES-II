
      subroutine form_u_fock(fock,umat,totlen,spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
C Input variables
      integer totlen,spin
      double precision fock(totlen)
C Output variables
      double precision umat(totlen)
C Local variables
      integer ioff,irrep,occ,virt,mos,ii
      double precision one
      data one /1.0d0/

      ioff = 1
      do irrep = 1,nirrep
        occ=pop_full(irrep,spin)
        virt=vrt_full(irrep,spin)
        mos=occ+virt
        if (min(occ,virt) .gt. 0) then
          call eig(fock(ioff),umat(ioff),1,mos,1)
        else
          if (mos .gt. 0) then
            call zero(umat(ioff),mos*mos)
            do ii = 1,mos
              umat(ioff+(ii-1)*mos+ii-1) = one
            end do
          endif
        endif
        ioff = ioff + mos*mos
      end do

      return
      end
