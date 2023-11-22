      subroutine setbasHFKS
c      nbas ..... number of basis function per irrep
c      nbasis ... total number of basis function
c      nbassq ... size of a symmetry packed square matrix
c      nbastt ... size of a symmetry packed triangular matrix
      
      implicit integer(a-z)
      common/bassym/nbas(8),nbasis,nbassq,nbastt
      common/lsym/nlenq(8),nlent(8)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
      
      do 100 irrep=1,nirrep
           nbas(irrep)=pop(irrep,1)+vrt(irrep,1)
100   continue
           nbasis=0
           nbassq=0
           nbastt=0
      do 200 irrep=1,nirrep
            nbasis=nbasis+nbas(irrep)
            nbassq=nbassq+nbas(irrep)*nbas(irrep)
            nbastt=nbastt+(nbas(irrep)*(nbas(irrep)+1))/2
200   continue
             nlenq(1)=nbassq
             nlent(1)=nbastt
      do 300 irrep=2,nirrep
             nlenq(irrep)=0
             nlent(irrep)=0
         do 400 irrep1=1,nirrep
                irrep2=dirprd(irrep,irrep1)
                if(irrep2.gt.irrep1) then
      nlent(irrep)=nlent(irrep)+nbas(irrep1)*nbas(irrep2)
                endif
      nlenq(irrep)=nlenq(irrep)+nbas(irrep1)*nbas(irrep2)
400      continue
300   continue
      return
      end
