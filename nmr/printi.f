   

      subroutine printi(a,ipert)
C
      implicit double precision (a-h,o-z)
      integer pop,vrt,dirprd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),ISYTYP(2,500),ntot(18)
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
      common/dsym/irrepx,i1pert,ndt(2),ndf1(2),ndf2(2),ioff(8,6)
      dimension a(1)
C
      write(*,*) 'pop',(pop(i,1),i=1,8)
      write(*,*) 'vrt',(vrt(i,1),i=1,8)
C
C dump d f(ij)/dchi
C
      call getlst(a,ipert,1,1,irrepx,176)
      write(*,*) 'd f ij '
      write(*,1111) (a(i),i=1,ndf1(1))
1111  format(8F10.7)
C
C dump d f(ab)/dchi
C
      call getlst(a,ipert,1,1,irrepx,178)
      write(*,*) 'd f ab '
      write(*,1111) (a(i),i=1,ndf2(1))
C
C dump d S(ab)/dchi
C
      call getlst(a,ipert,1,1,irrepx,172)
      write(*,*) 'S ab '
      write(*,1111) (a(i),i=1,ndf2(1))
C
C dump d S(ij)/dchi
C
      call getlst(a,ipert,1,1,irrepx,170)
      write(*,*) 'S ij '
      write(*,1111) (a(i),i=1,ndf1(1))
C
C dump d S(ai)/dchi
C
      call getlst(a,ipert,1,1,irrepx,174)
      write(*,*) 'S ai '
      write(*,1111) (a(i),i=1,ndt(1))
C
C dump d U(ai)/dchi
C
      call getlst(a,ipert,1,1,irrepx,182)
      write(*,*) 'U ai '
      write(*,1111) (a(i),i=1,ndt(1))
C
C dump D(ij) 
C
      call getlst(a,1,1,1,1,160)
      write(*,*) 'D ij '
      write(*,1111) (a(i),i=1,nf1(1))
C
C dump D(ab) 
C
      call getlst(a,1,1,1,3,160)
      write(*,*) 'D ab '
      write(*,1111) (a(i),i=1,nf2(1))
C
C dump D(ai) 
C
      call getlst(a,1,1,1,5,160)
      write(*,*) 'D ai '
      write(*,1111) (a(i),i=1,nt(1))
C
C
C dump integrals ijab
C
      sum=0.0
      sum1=0.0
      sum2=0.0
      do 1233 irrep1=1,nirrep
      irrep2=dirprd(irrep1,irrepx)
      n1=irpdpd(irrep1,isytyp(2,16))
      n2=irpdpd(irrep2,isytyp(1,16))
      write(*,*) n1,n2
      call getlst(a,1,n1,1,irrep1,316)
C
      call checksm2('dabij ',a,n1,n2,sum,sum1,sum2)
      write(*,*) 'ijab '
      do 1232 j=1,n1
      write(*,*) 'irrep2',irrep2,'irrep1',irrep1,'dis',j
      write(*,1111) (a(i),i=(j-1)*n2+1,j*n2)
1232  continue
1233  continue
C
C
C dump integrals d<ij||kl>/dchi
C
      write(*,*) 'ijkl'
      sum=0.0
      sum1=0.0
      sum2=0.0
      do 1234 irrep1=1,nirrep
      irrep2=dirprd(irrep1,irrepx)
      n1=irpdpd(irrep1,isytyp(2,13))
      n2=irpdpd(irrep2,isytyp(1,13))
      call getlst(a,1,n1,1,irrep1,313)
C
      call checksm2('dijkl ',a,n1,n2,sum,sum1,sum2)
      write(*,*) 'd ijkl ',irrepx,irrep1,irrep2
      write(*,1111) (a(i),i=1,n1*n2)
1234  continue
C
C
      sum=0.0
      sum1=0.0
      sum2=0.0
      do 1235 irrep1=1,nirrep
      irrep2=dirprd(irrep1,irrepx)
      n1=irpdpd(irrep1,isytyp(2,10))
      n2=irpdpd(irrep2,isytyp(1,10))
      call getlst(a,1,n1,1,irrep1,310)
C
      call checksm2('dijka ',a,n1,n2,sum,sum1,sum2)
c      write(*,*) 'd ijka ',irrepx,irrep1,irrep2
c      write(*,1111) (a(i),i=1,n1*n2)
C
1235  continue
C
      do 1238 irrep1=1,nirrep
      irrep2=dirprd(irrep1,irrepx)
      n1=irpdpd(irrep1,isytyp(2,25))
      n2=irpdpd(irrep2,isytyp(1,25))
      call getlst(a,1,n1,1,irrep1,325)
      call checksm2('diajb ',a,n1,n2,sum,sum1,sum2)
1238  continue
C
c      write(*,*) 'iabc '
c      write(*,1111) (a(i),i=1,n1*n2)
C
      do 1239 irrep1=1,nirrep
      irrep2=dirprd(irrep1,irrepx)
      n1=irpdpd(irrep1,isytyp(2,30))
      n2=irpdpd(irrep2,isytyp(1,30))
      call getlst(a,1,n1,1,irrep1,330)
      call checksm2('diabc ',a,n1,n2,sum,sum1,sum2)
1239  continue
C
      return
      write(*,*) 'abci '
      write(*,1111) (a(i),i=1,n1*n2)
C
c      n1=irpdpd(1,isytyp(2,33))
c      n2=irpdpd(1,isytyp(1,33))
c      call getlst(a,1,n1,1,1,233)
C
c      write(*,*) 'abcd '
c      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,16))
      n2=irpdpd(1,isytyp(1,16))
      call getlst(a,1,n1,1,1,316)
C
      write(*,*) 'd abij '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,10))
      n2=irpdpd(1,isytyp(1,10))
      call getlst(a,1,n1,1,1,310)
C
      write(*,*) 'd ijka '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,30))
      n2=irpdpd(1,isytyp(1,30))
      call getlst(a,1,n1,1,1,330)
C
      write(*,*) 'd abci '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,25))
      n2=irpdpd(1,isytyp(1,25))
      call getlst(a,1,n1,1,1,325)
      write(*,*) 'd aibj '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,21))
      n2=irpdpd(1,isytyp(1,21))
      call getlst(a,1,n1,1,1,321)
      write(*,*) 'd iabj '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,16))
      n2=irpdpd(1,isytyp(1,16))
      call getlst(a,1,n1,1,1,446)
C
      write(*,*) 'd t(abij) '
      write(*,1111) (a(i),i=1,n1*n2)
C
      n1=irpdpd(1,isytyp(2,16))
      n2=irpdpd(1,isytyp(1,16))
      call getlst(a,1,n1,1,1,50)
C
      write(*,*) 'den(abij) '
      write(*,1111) (a(i),i=1,n1*n2)
      return
      end
