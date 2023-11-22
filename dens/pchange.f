      SUBROUTINE PCHANGE (IUHF,NIRREP,ienter,ifl77)
      IMPLICIT INTEGER (P,V)
      CHARACTER*6 SYMLST(22)
      COMMON/SYM/NPVTF(38)
      COMMON/SYM2/MPVTF(38)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON /SYMPOP/ IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
      DATA SYMLST /'SVAVA0','SVBVB0','SOAOA0','SOBOB0',
     &             'SVAVA1','SVBVB1','SOAOA1','SOBOB1',
     &             'SVAOA2','SVBOB2','SOBVA2','SVBOA2',
     &             'SVAVB2','SOAOB2','SVAVB2','SOAVA2',
     &             'SOBVB2','SOAVB2','SVAVA2','SVBVB2',
     &             'SOAOA2','SOBOB2'/
      if (ienter.eq.2) go to 222
      if (ienter.ne.1) then
       write (6,10) ienter
 10    format (4x,'==== ienter should be 1 or 2, ienter=',i4)
       call errex
      endif
c---------------------------------------------------------
C---    RE-INITIALIZES COMMON BLOCK /INFO/, /SYMPOP2/
c---------------------------------------------------------
      CALL GETREC(20,'JOBARC','NDROTPOP',1,ndrpopt)
      CALL GETREC(20,'JOBARC','NDROTVRT',1,ndrvrtt)
      NOCCO(1) = NOCCO(1) + ndrpopt
      NOCCO(2) = NOCCO(2) + ndrpopt
      NVRTO(1) = NVRTO(1) + ndrvrtt
      NVRTO(2) = NVRTO(2) + ndrvrtt
c---------------------------------------------------------
      CALL IZERO(IRPDPD,176)
c---
      call aces_ja_fin
      istate = ishell('mv JOBARC JOBARC_DM')
      istate = ishell('mv JAINDX JAINDX_DM')
      istate = ishell('mv JOBARC_AM JOBARC')
      istate = ishell('mv JAINDX_AM JAINDX')
      call aces_ja_init
c---
      DO 30 ITYPE=1,22
       CALL GETREC(20,'JOBARC',SYMLST(ITYPE)//'X ',NIRREP,
     &             IRPDPD(1,ITYPE))
  30  CONTINUE
c---
      call aces_ja_fin
      istate = ishell('mv JOBARC JOBARC_AM')
      istate = ishell('mv JAINDX JAINDX_AM')
      istate = ishell('mv JOBARC_DM JOBARC')
      istate = ishell('mv JAINDX_DM JAINDX')
      call aces_ja_init
c---
      if (ifl77.ne.0) call QRHFSET2
      return
c---------------------------------------------------------
c---        for the case ientry = 2
C---    RE-INITIALIZES COMMON BLOCK /SYM/
c---------------------------------------------------------
 222  continue
      do 225 i=1,38
       npvtf(i) = mpvtf(i)
 225  continue
      return
      end
