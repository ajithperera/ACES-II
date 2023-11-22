      SUBROUTINE RECALAT(IUHF,NIRREP,NDROP0,EVAL,EVEC,IINTFP,ifl77)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SYM2/NPOP2(8,2),NVRT2(8,2),NT1,NT2,NI1,NI2,NA1,NA2
      COMMON /SYM/NPOP(8,2),NVRT(8,2),NT10,NT20,NI10,NI20,NA10,NA20
      COMMON /SYMDROP/NDRPOP(8),NDRVRT(8)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      common /dropgeo/ ndrgeo
      dimension eval(2),evec(2)
c---------------------------------------------
      ndrgeo = 0
      ione=1
      call getrec(20,'JOBARC','NUMDROPA',ione,NDROP0)
      if (ndrop0.eq.0) go to 333
c
      write (6,111) ndrop0
 111  format (/2x,70('=')
     x    /2x,'[  Drop-MO Analytic Gradient with Total No of',
     x     1x,'Drop-MO = ',i3,' ,',3x,'KKB  ]'/2x,70('=')/)
c
      ndrgeo = 1
      call getrec(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
      call getrec(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT)
      call putrec(20,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
      call putrec(20,'JOBARC','NDROPVRT',NIRREP,NDRVRT)
c-----------
      nbasa = 0
      nbasd = 0
      DO 5 I=1,NIRREP
       NPOP2(I,1) = NPOP(I,1) + NDRPOP(I)
       NVRT2(I,1) = NVRT(I,1) + NDRVRT(I)
       nbasa = nbasa + npop2(I,1) + nvrt2(I,1)
       nbasd = nbasd + npop(I,1) + nvrt(I,1)
       IF (IUHF.EQ.0) THEN
        NPOP2(I,2) = NPOP2(I,1)
        NVRT2(I,2) = NVRT2(I,1)
       ELSE
        NPOP2(I,2) = NPOP(I,2) + NDRPOP(I)
        NVRT2(I,2) = NVRT(I,2) + NDRVRT(I)
       ENDIF
   5  CONTINUE
c---------------------------------------------------------
      call getrec(20,'JOBARC','NBASTOT ',ione,nbas)
      call putrec(20,'JOBARC','NBASTOT ',ione,nbasa)
      call putrec(20,'JOBARC','NBASTOTD',ione,nbasd)
c
      call getrec(20,'JOBARC','SCFEVALA',NBASd*IINTFP,EVAL)
      call putrec(20,'JOBARC','DCFEVALA',NBASd*IINTFP,EVAL)
      call getrec(20,'JOBARC','SCFEVLA0',NBASa*IINTFP,EVAL)
      call putrec(20,'JOBARC','SCFEVALA',NBASa*IINTFP,EVAL)
      call getrec(20,'JOBARC','SCFEVECA',NBASd*NBAS*IINTFP,EVEC)
      call putrec(20,'JOBARC','DCFEVECA',NBASd*NBAS*IINTFP,EVEC)
      call getrec(20,'JOBARC','SCFEVCA0',NBASa*NBAS*IINTFP,EVEC)
      call putrec(20,'JOBARC','SCFEVECA',NBASa*NBAS*IINTFP,EVEC)
      if (iuhf.eq.1) then
       call getrec(20,'JOBARC','SCFEVALB',NBASd*IINTFP,EVAL)
       call putrec(20,'JOBARC','DCFEVALB',NBASd*IINTFP,EVAL)
       call getrec(20,'JOBARC','SCFEVLB0',NBASa*IINTFP,EVAL)
       call putrec(20,'JOBARC','SCFEVALB',NBASa*IINTFP,EVAL)
       call getrec(20,'JOBARC','SCFEVECB',NBASd*NBAS*IINTFP,EVEC)
       call putrec(20,'JOBARC','DCFEVECB',NBASd*NBAS*IINTFP,EVEC)
       call getrec(20,'JOBARC','SCFEVCB0',NBASa*NBAS*IINTFP,EVEC)
       call putrec(20,'JOBARC','SCFEVECB',NBASa*NBAS*IINTFP,EVEC)
      endif
c---------------------------------------------------------
      if (ifl77.ne.0) then
       call aces_ja_fin
       istate = ishell('mv JOBARC JOBARC_DM')
       istate = ishell('mv JAINDX JAINDX_DM')
       istate = ishell('mv JOBARC_AM JOBARC')
       istate = ishell('mv JAINDX_AM JAINDX')
       call aces_ja_init
c
       call getrec(20,'JOBARC','SCFEVALA',NBASa*IINTFP,EVAL)
       call getrec(20,'JOBARC','SCFEVECA',NBASa*NBAS*IINTFP,EVEC)
       if (iuhf.eq.1) then
        call getrec(20,'JOBARC','SCFEVALB',NBASa*IINTFP,
     &              EVAL(NBASa+1))
        call getrec(20,'JOBARC','SCFEVECB',NBASa*NBAS*IINTFP,
     &              EVEC(NBASa*NBAS+1))
       endif
c
       call aces_ja_fin
       istate = ishell('mv JOBARC JOBARC_AM')
       istate = ishell('mv JAINDX JAINDX_AM')
       istate = ishell('mv JOBARC_DM JOBARC')
       istate = ishell('mv JAINDX_DM JAINDX')
       call aces_ja_init
c
       call putrec(20,'JOBARC','SCFEVALA',NBASa*IINTFP,EVAL)
       call putrec(20,'JOBARC','SCFEVECA',NBASa*NBAS*IINTFP,EVEC)
       if (iuhf.eq.1) then
        call putrec(20,'JOBARC','SCFEVALB',NBASa*IINTFP,
     &              EVAL(NBASa+1))
        call putrec(20,'JOBARC','SCFEVECB',NBASa*NBAS*IINTFP,
     &              EVEC(NBASa*NBAS+1))
       endif
c
      endif
c---------------------------------------------------------
      NT1 = 0
      NI1 = 0
      NA1 = 0
      DO 10 I=1,NIRREP
       NI = NPOP(I,1) + NDRPOP(I)
       NA = NVRT(I,1) + NDRVRT(I)
       NT1 = NT1 + NI * NA
       NI1 = NI1 + NI * NI
       NA1 = NA1 + NA * NA
  10  CONTINUE
c
      IF (IUHF.EQ.0) THEN
       NT2 = NT1
       NI2 = NI1
       NA2 = NA1
      ELSE
c
       NT2 = 0
       NI2 = 0
       NA2 = 0
       DO 20 I=1,NIRREP
        NI = NPOP(I,2) + NDRPOP(I)
        NA = NVRT(I,2) + NDRVRT(I)
        NT2 = NT2 + NI * NA
        NI2 = NI2 + NI * NI
        NA2 = NA2 + NA * NA
  20   CONTINUE
      ENDIF
c--------------
      return
c----------------------------------------------
c      for the case  ndrop = 0
c----------------------------------------------
 333  continue
      do 400 i=1,nirrep
       NPOP2(I,1) = NPOP(I,1)
       NPOP2(I,2) = NPOP(I,2)
       NVRT2(I,1) = NVRT(I,1)
       NVRT2(I,2) = NVRT(I,2)
 400  continue
      nt1 = nt10
      nt2 = nt20
      ni1 = ni10
      ni2 = ni20
      na1 = na10
      na2 = na20
c
      do 500 i=1,22
       do 490 j=1,nirrep
        irpdpd(j,i) = irp_dm(j,i)
 490   continue
 500  continue
      return
      end
