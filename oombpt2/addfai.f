
      subroutine addfai2(vai,fai,doo,dvv,lenvo,lenoo,lenvv,spin,listv)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C Input variables
      integer lenvo,lenoo,lenvv,spin,listv
C Preallocated local variables
      double precision vai(lenvo),fai(lenvo),doo(lenoo),dvv(lenvv)
C Local variables
      integer iai,ioo,ivv,occ,virt,irrep,ii,aa,offai
      double precision ddoo

      call getlst(vai,1,1,1,spin,listv)
      call getlst(fai,1,1,1,2+spin,93)
      call getlst(doo,1,1,1,spin,160)
      call getlst(dvv,1,1,1,2+spin,160)

      iai=0
      ioo=0
      ivv=0
      do irrep = 1,nirrep
        occ=pop(irrep,spin)
        virt=vrt(irrep,spin)
        if (min(occ,virt) .gt. 0) then
          do ii=1,occ
            ddoo = doo(ioo+(ii-1)*occ+ii)
            do aa=1,virt
              offai=iai+(ii-1)*virt+aa
              vai(offai) = vai(offai)
     &                    + fai(offai)*(ddoo-dvv(ivv+(aa-1)*virt+aa))
            end do
          end do
          iai=iai+occ*virt
          ioo=ioo+occ*occ
          ivv=ivv+virt*virt
        endif
      end do
      call putlst(vai,1,1,1,spin,listv)

      return
      end


      subroutine addfai(vai,fai,doo,dvv,lenvo,lenoo,lenvv,spin,listv)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C Input variables
      integer lenvo,lenoo,lenvv,spin,listv
C Preallocated local variables
      double precision vai(lenvo),fai(lenvo),doo(lenoo),dvv(lenvv)
C Local variables
      integer iai,ioo,ivv,nocc,nvrt,irrep
      double precision one
      data one /1.0d0/

      call getlst(vai,1,1,1,spin,listv)
      call getlst(fai,1,1,1,2+spin,93)
      call getlst(doo,1,1,1,spin,160)
      call getlst(dvv,1,1,1,2+spin,160)
       
      iai=1
      ioo=1
      ivv=1
      do irrep = 1,nirrep
        nocc = pop(irrep,spin)
        nvrt = vrt(irrep,spin)
csss        Write(6,*) "Print fai in addfai"
csss        call output(fai, 1, nvrt, 1, nocc, nvrt, nocc, 1)

        if (min(nocc,nvrt) .gt. 0) then
          call xgemm('n','n',nvrt,nocc,nocc,-one,fai(iai),nvrt,doo(ioo),
     &               nocc,one,vai(iai),nvrt)
          call xgemm('n','n',nvrt,nocc,nvrt,+one,dvv(ivv),nvrt,fai(iai),
     &               nvrt,one,vai(iai),nvrt)
        endif
        iai = iai + nvrt*nocc
        ivv = ivv + nvrt*nvrt
        ioo = ioo + nocc*nocc
      end do

      call putlst(vai,1,1,1,spin,listv)

      return
      end
