
      subroutine formfock_d(vai,fock,evals,totlen,nmo,spin)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer nmo,spin,totlen
C Pre-allocated local variables
      double precision vai(totlen),fock(totlen),evals(nmo)
C Local variables
      integer ifock,ivai,ievalo,ievalv,irrep,occ,virt,mos,aa,ii,dropocc
      character*8 cscfevl0(2)
      data cscfevl0 /'SCFEVLA0','SCFEVLB0'/

      ifock = 1
      ivai = 1
      ievalo = 1
      call getrec(20,'JOBARC','NDROTPOP',1,dropocc)
      ievalv=nocco(spin)+1+dropocc
      call getrec(20,'JOBARC',cscfevl0(spin),nmo*iintfp,evals)
      do irrep=1,nirrep
        occ=pop_full(irrep,spin)
        virt=vrt_full(irrep,spin)
        mos=occ+virt
        call zero(fock(ifock),mos*mos)
C Diagonal OO/VV part
        if (occ .gt. 0)
     &    call dcopy(occ,evals(ievalo),1,fock(ifock),mos+1)
        if (virt .gt. 0)
     &    call dcopy(virt,evals(ievalv),1,fock(occ*mos+occ+ifock),mos+1)

        if (min(occ,virt) .gt. 0) then
C Add OV part Fock + potential part
          call blkcpy(vai(ivai),virt,occ,fock(ifock),mos,mos,occ+1,1)
C symmetrize ov block
          do ii=1,occ
            do aa=occ+1,mos
              fock((aa-1)*mos+ii-1+ifock)=fock((ii-1)*mos+aa-1+ifock)
            end do
          end do
        endif
        ifock = ifock + mos*mos
        ivai = ivai + occ*virt
        ievalo = ievalo + occ
        ievalv = ievalv + virt
      end do


      return
      end
