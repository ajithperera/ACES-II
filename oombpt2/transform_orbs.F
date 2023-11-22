
      subroutine transform_orbs(orbs,umat,scr,nao,nmo,totlen,mxirr2,
     &                          spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
C Input variables
      integer totlen,spin,nao,nmo,mxirr2
      double precision orbs(nao*nmo),umat(totlen),scr(mxirr2)
C Local variables
      integer iorbs_c,iorbs_r,ioff,irrep,mos,aos
      double precision one,zilch
      data one,zilch /1.0d0,0.0d0/

      iorbs_c = 1
      iorbs_r = 1
      ioff = 1
      do irrep = 1,nirrep
        mos = pop_full(irrep,spin) + vrt_full(irrep,spin)
        aos = mos
        if (mos .gt. 0) then
          call xgemm('n','n',aos,mos,mos,one,
     &               orbs((iorbs_c-1)*nao+iorbs_r),nao,umat(ioff),mos,
     &               zilch,scr,aos)
          call blkcpy(scr,aos,mos,orbs,nao,nmo,iorbs_r,iorbs_c)
        endif
        iorbs_c = iorbs_c + mos
        iorbs_r = iorbs_r + aos
        ioff = ioff + mos*mos
      end do

      return
      end
