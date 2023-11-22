
      subroutine expand_vai(vai_drop,vai_full,spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer pop_full(8,2),vrt_full(8,2),nt_full(2),nd1_full(2),
     &        nd2_full(2)
      common/sym_full/pop_full,vrt_full,nt_full,nd1_full,nd2_full
C Input variables
      integer spin,len_full,len_drop, nocc, nvrt
      double precision vai_drop(nt(spin))
C Output variables
      double precision vai_full(nt_full(spin))
C Local variables
      integer ivai_f,ivai_d,irrep,occ_f,virt_f,virt_d,occ_d,odiff,ii,aa

      ivai_f = 1
      ivai_d = 1
      call zero(vai_full,nt_full(spin))
      do irrep = 1,nirrep
        occ_f = pop_full(irrep,spin)
        virt_f = vrt_full(irrep,spin)
        occ_d = pop(irrep,spin)
        virt_d = vrt(irrep,spin)
        odiff = occ_f - occ_d
        do ii = 1+odiff,occ_f
          call dcopy(virt_d,vai_drop(ivai_d+(ii-1-odiff)*virt_d),1,
     &               vai_full(ivai_f+(ii-1)*virt_f),1)
        end do
        ivai_f = ivai_f + occ_f*virt_f
        ivai_d = ivai_d + occ_d*virt_d
      end do

      write(6,*) "VAI"
      nocc= pop_full(1,1)
      nvrt= vrt_full(1,1)
      call output(vai_full, 1, nvrt, 1, nocc, nvrt, nocc,1)
      return
      end
