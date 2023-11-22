
      subroutine form_dropmo_length(nmo,totlen,mxirr2)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer pop_full(8,2),vrt_full(8,2),nt_full(2),nd1_full(2),
     &        nd2_full(2)
      common/sym_full/pop_full,vrt_full,nt_full,nd1_full,nd2_full
C Output variables
      integer nmo,totlen,mxirr2
C Local variables
      integer spin,irrep,irr2,ndrop(8)

      call getrec(-1,'JOBARC','NDROPPOP',nirrep,ndrop)
      do spin = 1,2
        nd1_full(spin) = 0
        do irrep = 1,nirrep
          pop_full(irrep,spin) = pop(irrep,spin)+ndrop(irrep)
          nd1_full(spin) = nd1_full(spin) + pop_full(irrep,spin)**2
        end do
      end do

      call getrec(-1,'JOBARC','NDROPVRT',nirrep,ndrop)
      do spin = 1,2
        nt_full(spin) = 0
        nd2_full(spin) = 0
        do irrep = 1,nirrep
          vrt_full(irrep,spin) = vrt(irrep,spin)+ndrop(irrep)
          nt_full(spin) = nt_full(spin) +
     &                         pop_full(irrep,spin)*vrt_full(irrep,spin)
          nd2_full(spin) = nd2_full(spin) + vrt_full(irrep,spin)**2
        end do
      end do

      totlen = 0
      mxirr2 = 0
      nmo = 0
      do irrep = 1,nirrep
        nmo = nmo + pop_full(irrep,1)+vrt_full(irrep,1)
        irr2 = (pop_full(irrep,1)+vrt_full(irrep,1))**2
        totlen = totlen+irr2
        mxirr2 = max(mxirr2,irr2)
      end do

      return
      end
