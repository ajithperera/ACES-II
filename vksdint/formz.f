      subroutine formz(xia,eval,pop,vrt,nocc)

      implicit double precision(a-h,o-z)
      integer pop,vrt,dirprd
      dimension xia(1),eval(1),pop(8),vrt(8)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)

      ind=0
      ioffp=0
      ioffv=nocc
      do 10 irrep=1,nirrep
       nocci=pop(irrep)
       nvrti=vrt(irrep)
       do 5 i=1,nocci
        indi=i+ioffp
        do 5 ia=1,nvrti
        inda=ia+ioffv
        ind=ind+1
        xia(ind)=xia(ind)/(eval(indi)-eval(inda))
5      continue 
       ioffp=ioffp+nocci
       ioffv=ioffv+nvrti
10    continue
      return
      end
