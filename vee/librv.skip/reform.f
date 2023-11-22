
      SUBROUTINE REFORM(R,SCR,NDIM,NDIMVEC,IRREPX,MAXCOR,ISIDE,IUHF)
C
C THIS ROUTINE CALCULATES THE COMPLETE R MATRIX FROM THE
C EXPANSION VECTORS AND MATRIX VECTOR PRODUCTS ON DISK.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(*),R(*),SDOT,X
      PARAMETER (MAXORD=100)
      COMMON/EXTINF/ITER,IOLDEST
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
      INDXF(I,J,N)=I+(J-1)*N
C
C SEE IF WE CAN GET ALL OLD AND NEW VECTORS INTO CORE
C SIMULTANEOUSLY.  THIS IS A GREAT SITUATION!
C
      IOFF1=1
      IOFF2=IOFF1+NDIMVEC
      do 10 i=1,ndim
       call getlst(scr(ioff1),iget(i),1,1,iside,470)
       if(iuhf.eq.0)then
        call spntsing(irrepx,scr(ioff1),scr(ioff2),maxcor-ioff2+1)
       endif
       do 11 j=1,ndim
        call getlst(scr(ioff2),iget(j),1,1,iside,471)
        x=sdot(ndimvec,scr(ioff1),1,scr(ioff2),1)
        ipos=indxf(i,j,maxord)
        r(ipos)=x
11     continue
10    continue
c
      return
      end
