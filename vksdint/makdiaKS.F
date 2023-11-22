
      subroutine makdiaKS(xia,icore,maxcor,iuhf,bredundant)
      
      implicit double precision(a-h,o-z)
      LOGICAL INCORE,DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,
     &        TRIP1,TRIP2,GABCD,RELAXED,TRULY_NONHF,ONLYM,ONLYP
      logical bRedundant
      INTEGER DIRPRD,POP,VRT,POPRHF,VRTRHF,POPDOC,VRTDOC
      DIMENSION XIA(1),ICORE(MAXCOR)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common/flags/ iflags(100)
c below common block sympop is strange

      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)

c      common/sympop2/irpdpd(8,22)
c      common/sympop/irp_dm(8,22),isytyp(2,500),ntot(18)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
c little problem with below common block
c  sym2 doesn't make any sence
c nf1,nf2 weird nonsence
c used in linear2 so need to be changed
c i am replaing the sym2 common block with sym common block
       common /sym/ pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)    
c      common/sym2/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
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
       call getrec(20,'JOBARC','SCFEVALA',iintfp*norba,icore(ievaa))
       if(iuhf.eq.1) then
        norbb=nocco(2)+nvrto(2)
        ievbb=ievaa-iintfp*norbb
        mxcor=mxcor-iintfp*norbb
        call getrec(20,'JOBARC','SCFEVALB',iintfp*norbb,icore(ievbb))
       endif

       call formz(xia,icore(ievaa),pop(1,1),vrt(1,1),nocco(1))
      return
      end
