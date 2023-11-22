
      subroutine form_err(fock,dens,err,totlen,spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
C Input variables
      integer totlen,spin
      double precision fock(totlen),dens(totlen)
C Output variables
      double precision err(totlen)
C Local variables
      integer ioff,irrep,mos
      double precision one,zilch
      data one,zilch /1.0d0,0.0d0/

C Transform density to original basis
      ioff = 1
      do irrep = 1,nirrep
        mos=pop_full(irrep,spin)+vrt_full(irrep,spin)
        if (mos .gt. 0) then
C Calculate commutator with the fock matrix
C F*D
          call xgemm('n','n',mos,mos,mos,one,fock(ioff),mos,dens(ioff),
     &               mos,zilch,err(ioff),mos)
C -D*F
          call xgemm('n','n',mos,mos,mos,-one,dens(ioff),mos,fock(ioff),
     &               mos,one,err(ioff),mos)
        endif
        ioff = ioff + mos*mos
      end do

      return
      end
