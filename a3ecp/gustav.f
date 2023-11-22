      subroutine gustav(ltdrv,jx,jy,jz,iadr,imasc,maxnft,mxlt,imin,
     &                  imax,nftdrv)


      implicit double precision (a-h,o-z)
      dimension jx(maxnft),jy(maxnft),jz(maxnft),
     &          iadr(0:mxlt,0:mxlt,0:mxlt),imin(ltdrv),imax(ltdrv)
      logical imasc(maxnft),odd
c
c     set data needed for calculation of derivative integrals over
c     gaussian basis functions up to l quantum number ltmax
c     ltdrv = ltmax+1    cartesian first derivatives
c                  +2    cartesian second derivatives or
c                        derivatives with respect to basis set exponents
c
c --- generate jx,y,z() and iadr() arrays ---
c     jx,y,z : monomial degrees of gaussians
c     iadr   : integral addresses labelled with respect to
c                 monomial degrees
c     nftdrv : number of reducible gaussians with l quantum number ltdrv
c
c     call monom(ltdrv,nftdrv,jx,jy,jz,iadr,maxnft,mxlt,1)
c
c --- initialize logical masc for derivatives of cartesian gaussians ---
c     imasc = .true.  for  s,d,g, ... gaussians (i.e. functions with
c                                                even l-quntum numbers)
c             .false. for  p,f,h, ... gaussians
c
      dimension nx(20),ny(20),nz(20),iorder(20)
      data nx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,1,1,1,2,2,3/
      data ny/0,0,1,0,0,1,2,0,1,0,0,1,2,3,0,1,2,0,1,0/
      data nz/0,1,0,0,2,1,0,1,0,0,3,2,1,0,2,1,0,1,0,0/
cch      data iorder/1,4,3,2,10,7,5,9,8,6,20,14,11,19,18,17,13,15,12,16/
CCH   iorder changed for ACES 2
      data iorder/1,4,3,2,10,9,8,7,6,5,20,19,18,17,16,15,14,13,12,11/

      nft=0
      do 100 i=1,20
        nft=nft+1
        jx(i)=nx(iorder(i))
        jy(i)=ny(iorder(i))
        jz(i)=nz(iorder(i))
        iadr(jx(i),jy(i),jz(i))=nft-1
  100 continue
 
      k=0
      do 200 i=1,ltdrv
        odd=mod(i,2).ne.0
        do 300 j=imin(i),imax(i)
          k=k+1
          imasc(k)=odd
  300   continue
  200 continue
 
      return
      end
