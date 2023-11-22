      subroutine lineqz2(amat,xia,xianew,xupdate,asmall,asquare,
     &                   EVALA,EVALB,
     &                   conv,na,nb,kmax,nocca,noccB
     &                   )
      implicit double precision(a-h,o-z)
      INTEGER DIRPRD,POP,VRT
      Logical semi
      DIMENSION DET(2)
      dimension amat(1),xia(na+nb), xianew(na+nb),xupdate(na+nb,kmax),
     &          asmall(kmax,kmax),asquare(kmax+1),evala(1),
     &          evalb(1)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common /sym/ pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common /flags/ iflags(100)

      data azero,one,onem,two/0.d+0,1.0d+0,-1.0d+0,2.d+0/

      if(iflags(39).eq.1) semi=.true.
      cutoff=conv*conv
      tol=azero
      n=na+nb

      ioff=1
      call zero(asmall,kmax*kmax)
      anorm=sdot(n,xia,1,xia,1)
      ascale=sqrt(anorm)

      if(ascale/n.lt.cutoff) then
       write(6,*) ' Warning form lineqz2 : the initial vector is zero.'
       ascale=one
      endif

      scale=one/ascale
      call sscal(n,scale,xia,1)

      asquare(1)=one
      call dcopy(n,xia,1,xupdate(1,1),1)

      do 1000 k=1,kmax
       do 90 ispin=1,2
       if(ispin.eq.1) then
        n1=na 
       else
        n1=nb
       endif
       listw=18+ispin  
       call getlst(amat,1,n1,2,1,listw)
       listw=22+ispin  
       ioffa=1+n1*n1
       call getlst(amat(ioffa),1,n1,2,1,listw)
       call saxpy(n1*n1,onem,amat(ioffa),1,amat,1)
       call xgemm('N','N',1,n1,n1,one,xia(1+(ispin-1)*na),1,amat,n1,
     &            azero,xianew(1+(ispin-1)*na),1)
90     continue
       LISTW=18 
       call getlst(amat,1,nb,2,1,listw)
       call xgemm('N','T',1,na,nb,two,xia(1+na),1,amat,
     &            na,one,xianew,1)
       call xgemm('N','N',1,nb,na,two,xia,1,amat,na,
     &            one,xianew(1+na),1)
        CALL FORMZ(XIANEW,EVALA,POP(1,1),VRT(1,1),NOCCA)
        CALL FORMZ(XIANEW(1+NA),EVALB,POP(1,2),VRT(1,2),NOCCB)
      DO 100 L=1,K
      ASMALL(L,K)=SDOT(N,XUPDATE(1,L),1,XIANEW,1)
100   CONTINUE
      IF(K.GT.1) ASMALL(K,K-1)=ASQUARE(K)
      DO 200 L=1,K
      SCALE=-ASMALL(L,K)/ASQUARE(L)
      CALL SAXPY(N,SCALE,XUPDATE(1,L),1,XIANEW,1)
200   CONTINUE
      ASQUARE(K+1)=SDOT(N,XIANEW,1,XIANEW,1)
      CALL DCOPY(N,XIANEW,1,XIA,1)
      CALL DCOPY(N,XIANEW,1,XUPDATE(1,K+1),1)
      TEST=ASQUARE(K+1)/N
      IF(TEST.LE.CUTOFF) THEN
       ASQUARE(K+1)=ONE
       GO TO 300
      ENDIF
1000  CONTINUE
      WRITE(6,3000)
3000  FORMAT('  The iterative expansion of D(ai) did not ',
     &           'converge, abort !')
      WRITE(6,3001)TEST
3001  FORMAT('  The convergence is : ',E20.10) 
      CALL ERREX
  300 CONTINUE
      WRITE(6,3002) K 
3002  FORMAT('  The iterative expansion of D(ai) ',
     &           'converged after ',I3,' iterations.')
      DO 310 I=1,K
       ASMALL(I,I)=ASMALL(I,I)-ASQUARE(I)
310   CONTINUE
      CALL MINV(ASMALL,K,KMAX,XIANEW,DET,TOL,0,1)
      DO 330 I=1,K
       ASQUARE(I)=-ASCALE*ASMALL(I,1)
330   CONTINUE
      CALL ZERO(XIA,N)
      DO 340 I=1,K
       CALL SAXPY(N,ASQUARE(I),XUPDATE(1,I),1,XIA,1)
340   CONTINUE
      RETURN
      END
