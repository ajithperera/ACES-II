
      subroutine transform_obj(obj_in,umat,obj_out,trans,side,totlen,
     &                         spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
C Input variables
      character*1 trans,side
      integer totlen,spin
      double precision obj_in(totlen),umat(totlen),obj_out(totlen)
C Local variables
      integer ioff,irrep,mos
      double precision one,zilch
      data one,zilch /1.0d0,0.0d0/

      ioff = 1
      if ((side .eq. 'r').or.(side .eq. 'R')) then
        do irrep = 1,nirrep
          mos = pop_full(irrep,spin) + vrt_full(irrep,spin)
          if (mos .gt. 0)
     &      call xgemm('n',trans,mos,mos,mos,one,obj_in(ioff),mos,
     &                 umat(ioff),mos,zilch,obj_out(ioff),mos)
          ioff = ioff + mos*mos
        end do
      else if ((side .eq. 'l').or.(side .eq. 'L')) then
        do irrep = 1,nirrep
          mos = pop_full(irrep,spin) + vrt_full(irrep,spin)
          if (mos .gt. 0)
     &      call xgemm(trans,'n',mos,mos,mos,one,umat(ioff),mos,
     &                 obj_in(ioff),mos,zilch,obj_out(ioff),mos)
          ioff = ioff + mos*mos
        end do
      else
        write(6,*) '@TRANSFORM_OBJ: side must be in [r,R,l,L]'
        call errex
      endif
      return
      end
