      SUBROUTINE RESETSCF(IUHF,NIRREP,NDROP0,EVEC,IINTFP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
c     &               MOIODS(10,500),MOIOFL(10,500)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      DIMENSION EVEC(2)
      ione = 1
      call GETrec(20,'JOBARC','NBASTOT ',ione,nbasa)
      call GETrec(20,'JOBARC','NBASTOTD',ione,nbasd)
      call putrec(20,'JOBARC','NBASTOT ',ione,nbasd)
      call getrec(20,'JOBARC','DCFEVALA',NBASd*IINTFP,EVEC)
      call putrec(20,'JOBARC','SCFEVALA',NBASd*IINTFP,EVEC)
      call getrec(20,'JOBARC','DCFEVECA',NBASd*NBASa*IINTFP,EVEC)
      call putrec(20,'JOBARC','SCFEVECA',NBASd*NBASa*IINTFP,EVEC)
c-------------
      if (iuhf.eq.1) then
       call getrec(20,'JOBARC','DCFEVALB',NBASd*IINTFP,EVEC)
       call putrec(20,'JOBARC','SCFEVALB',NBASd*IINTFP,EVEC)
       call getrec(20,'JOBARC','DCFEVECB',NBASd*NBASa*IINTFP,EVEC)
       call putrec(20,'JOBARC','SCFEVECB',NBASd*NBASa*IINTFP,EVEC)
      endif
c-------------
c YAU : old
c     do 510 i=300,400
c      isytyp(1,i) = 0
c      isytyp(2,i) = 0
c      do 505 k=1,10
c       moio(k,i)   = 0
c       moiowd(k,i) = 0
c       moiosz(k,i) = 0
c       moiods(k,i) = 0
c       moiofl(k,i) = 0
c505   continue
c510  continue
c YAU : new
      call aces_io_remove(53,'DERINT')
c YAU : end
c---------------------------------------------
c     call aces_ja_fin
      istate = ishell('rm JOBARC_AM')
      istate = ishell('rm JAINDX_AM')
c     call aces_ja_init
c---------------------------------------------
      return
      end
