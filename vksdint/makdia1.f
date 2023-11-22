      subroutine makdia1(xia,icore,maxcor,iuhf,bredundant)
      
      implicit double precision(a-h,o-z)

      LOGICAL INCORE,DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,
     &        TRIP1,TRIP2,GABCD,RELAXED,TRULY_NONHF,ONLYM,ONLYP

      logical bRedundant
      INTEGER DIRPRD,POP,VRT,POPRHF,VRTRHF,POPDOC,VRTDOC
      DIMENSION XIA(1),ICORE(MAXCOR)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common/flags/ iflags(100)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
       common /sym/ pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)    
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2,GABCD,RELAXED,TRULY_NONHF
      common /dropgeo/ ndrgeo
      
      data tol/0.0d+0/
      kmax=iflags(31)
      conv=10.0**(-iflags(30))
      mxcor=maxcor

       norba=nocco(1)+nvrto(1)
       ievaa=mxcor+1-iintfp*norba
       mxcor=mxcor-iintfp*norba
       call getrec(20,'JOBARC','SCFEVLA0',iintfp*norba,icore(ievaa))
       if(iuhf.eq.1) then
        norbb=nocco(2)+nvrto(2)
        ievbb=ievaa-iintfp*norbb
        mxcor=mxcor-iintfp*norbb
        call getrec(20,'JOBARC','SCFEVALB',iintfp*norbb,icore(ievbb))
       endif

       if(iuhf.eq.0) then

       numsyw=irpdpd(1,isytyp(1,19)) 
       i010=1
       i020=i010+iintfp*numsyw*numsyw
       i030=i020+iintfp*numsyw*numsyw
       if(i030.lt.mxcor) then
        call mkarhf1(icore(i010),icore(i020),icore(i030),mxcor-i030+1,
     &               bredundant)
       else
        call insmem('MKARHF',i030,mxcor)
       endif 
       call formz(xia,icore(ievaa),pop(1,1),vrt(1,1),nocco(1))
c       call kkk(nt(1),xia)
      i030=i020+iintfp*max(numsyw,2*kmax)
      i040=i030+iintfp*(kmax+1)*numsyw
      i050=i040+iintfp*kmax*kmax
      i060=i050+iintfp*(kmax+1)
      if(i060.lt.mxcor) then
       call lineqz1(icore(i010),xia,icore(i020),icore(i030),icore(i040),
     &             icore(i050),icore(ievaa),conv,numsyw,kmax,nocco(1))
      else
       call insmem('lineqz1',i060,mxcor)
      endif
      else
      numsywa=irpdpd(1,isytyp(2,19)) 
      numsywb=irpdpd(1,isytyp(2,20))
      i010=1
      i020=i010+iintfp*(numsywa*numsywa+numsywb*numsywb+numsywa*numsywb)
      iend1=i020+iintfp*max(numsywa*numsywa,numsywb*numsywb)           
      i030=i020+iintfp*max(numsywa+numsywb,2*kmax)
      i040=i030+iintfp*(kmax+1)*(numsywa+numsywb)
      i050=i040+iintfp*kmax*kmax
      iend2=i050+iintfp*(kmax+1)
      if(max(iend1,iend2).lt.mxcor) then
       incore=.true.
       call mkauhf(icore(i010),icore(i020),ICORE(i030),mxcor-i030+1,
     &             bRedundant)
      else
       incore=.false.
       i020=i010+iintfp*2*max(numsywa*numsywa,numsywb*numsywb)
       i030=i020+iintfp*max(numsywa+numsywb,2*kmax)
       i040=i030+iintfp*(kmax+1)*(numsywa+numsywb)
       i050=i040+iintfp*kmax*kmax
       iend=i050+iintfp*(kmax+1)
       if(iend.ge.mxcor) call insmem('lineqz2',iend,mxcor)
      endif
        call formz(xia,icore(ievaa),pop(1,1),vrt(1,1),nocco(1))
        call formz(xia(1+numsywa),icore(ievbb),pop(1,2),vrt(1,2),
     &             nocco(2))
        call lineqz2(icore(i010),xia,icore(i020),icore(i030),
     &               icore(i040),icore(i050),icore(ievaa),
     &               icore(ievbb),
     &               conv,numsywa,numsywb,kmax,
     &               nocco(1),nocco(2))
      endif
      return
      ENd
