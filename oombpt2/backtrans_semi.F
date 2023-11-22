
      subroutine backtrans_semi(wmat,umat,totlen,spin)
      implicit none
C Common blocks
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2),nt_full(2),nd1_full(2),
     &        nd2_full(2)
      common/sym_full/pop_full,vrt_full,nt_full,nd1_full,nd2_full
C Input variables
      integer totlen,spin
C Output variables
      double precision umat(totlen)
C Pre-allocated Local variables
      double precision wmat(totlen)
C Local variables
      integer iwmat,iumat,irrep,occ,virt,mos
      character*8 wlabel(2,2)
      data wlabel /'WMATRXOA','WMATRXVA','WMATRXOB','WMATRXVB'/

      call getrec(20,'JOBARC',wlabel(1,spin),nd1_full(spin)*iintfp,wmat)
      iwmat = 1
      iumat = 1
      call zero(umat,totlen)
      do irrep = 1,nirrep
        occ = pop_full(irrep,spin)
        mos = occ + vrt_full(irrep,spin)
        if (occ .gt. 0)
     &    call blkcpy(wmat(iwmat),occ,occ,umat(iumat),mos,mos,1,1)
        iwmat = iwmat + occ*occ
        iumat = iumat + mos*mos
      end do
      call getrec(20,'JOBARC',wlabel(2,spin),nd2_full(spin)*iintfp,wmat)
      iwmat = 1
      iumat = 1
      do irrep = 1,nirrep
        virt = vrt_full(irrep,spin)
        occ = pop_full(irrep,spin)
        mos = occ + virt
        if (virt .gt. 0)
     &    call blkcpy(wmat(iwmat),virt,virt,umat(iumat),mos,mos,
     &                1+occ,1+occ)
        iwmat = iwmat + virt*virt
        iumat = iumat + mos*mos
      end do

      return
      end

