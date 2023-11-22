      SUBROUTINE GT2RING(CORE,MAXCOR,IUHF,LZOFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD,DSZZ,DSZTA,DSZTB,DSZWA,DSZWB
      DIMENSION CORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      DATA ONE   / 1.0D+00/
      DATA ZILCH / 0.0D+00/
      DATA ONEM  /-1.0D+00/
C
C     Compute contribution
C
C     W(ncdk) = - \sum_{og} t(cg,ok) * <no || dg>
C
C     Write to lists LZOFF+1,...,LZOFF+6. Note that we negate for lists
C     LZOFF+1, LZOFF+2, LZOFF+5, LZOFF+6 (a well known quirk of the MBEJ
C     lists). t(cg,ok) and <no || dg> are
C     assumed to be symmetric quantities. Target lists are assumed NOT to
C     have been spin-adapted. Target is stored dn;ck.
C
C                               Integral         T
C
C NCDK    - <DG||NO> * T(GC,KO)   DN;GO List 19; GO;CK List 34
C ncdk    - <dg||no> * T(gc,ko)   dn;go List 20; go;ck List 35
C
C NCDK    + <Dg||No> * T(Cg,Ko)   DN;go List 18; go;CK List 36 (UHF)
C                      T(gC,oK)                  GO;ck List 37 (RHF)
C ncdk    + <Gd||On> * T(Gc,Ok)   dn;GO List 17; GO;ck List 37
C
C-----------------------------------------------------------------------
C
C NcDk    - <Dg||No> * T(gc,ko)   DN;go List 18; go;ck List 35
C nCdK    - <Gd||On> * T(GC,KO)   dn;GO List 17; GO;CK List 34
C
C NcDk      <DG||NO> * T(Gc,Ok)   DN;GO List 19; GO;ck List 37
C nCdK      <dg||no> * T(Cg,Ko)   dn;go List 20; GO;ck List 37 (RHF)
C                                                go;CK List 36 (UHF)
C
C-----------------------------------------------------------------------
C
C nCDk      <Dg||On> * T(Cg,Ok)   Dn;gO List 21; gO;Ck List 38 (UHF)
C                      T(gC,kO)                  Go;Ck List 39 (RHF)
C NcdK      <Gd||No> * T(Gc,Ko)   dN;Go List 22; Go;cK List 39
C
C-----------------------------------------------------------------------
C
      DO 20 ISPIN=1,IUHF+1
       LISTZ  = LZOFF + ISPIN
       LISTTA =    33 + ISPIN
       LISTWA =    18 + ISPIN
       LISTTB =    35 + ISPIN
       IF(IUHF.EQ.0) LISTTB = 37
       LISTWB =    19 - ISPIN
C
      DO 10 IRREP=1,NIRREP
C
       DSZZ  = IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NDSZ  = IRPDPD(IRREP,ISYTYP(2,LISTZ))
       DSZTA = IRPDPD(IRREP,ISYTYP(1,LISTTA))
       NDSTA = IRPDPD(IRREP,ISYTYP(2,LISTTA))
       DSZTB = IRPDPD(IRREP,ISYTYP(1,LISTTB))
       NDSTB = IRPDPD(IRREP,ISYTYP(2,LISTTB))
       DSZWA = IRPDPD(IRREP,ISYTYP(1,LISTWA))
       NDSWA = IRPDPD(IRREP,ISYTYP(2,LISTWA))
       DSZWB = IRPDPD(IRREP,ISYTYP(1,LISTWB))
       NDSWB = IRPDPD(IRREP,ISYTYP(2,LISTWB))
C
       I000 = 1
       I010 = I000 + DSZZ  * NDSZ
       I020 = I010 + DSZWA * NDSWA
       I030 = I020 + DSZTA * NDSTA
       NEED = I030 * IINTFP
       IF(NEED.GT.MAXCOR)THEN
        CALL INSMEM('GT2RING',NEED,MAXCOR)
       ENDIF
C
       CALL GETLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
       CALL GETLST(CORE(I010),1,NDSWA,2,IRREP,LISTWA)       
       CALL GETLST(CORE(I020),1,NDSTA,2,IRREP,LISTTA)       
       CALL XGEMM('N','N',DSZZ,NDSZ,NDSWA,ONE ,
     &            CORE(I010),DSZWA,CORE(I020),DSZTA,ONE,
     &            CORE(I000),DSZZ)
C
       I000 = 1
       I010 = I000 + DSZZ  * NDSZ
       I020 = I010 + DSZWB * NDSWB
       I030 = I020 + DSZTB * NDSTB
       NEED = I030 * IINTFP
       IF(NEED.GT.MAXCOR)THEN
        CALL INSMEM('GT2RING',NEED,MAXCOR)
       ENDIF
C
       CALL GETLST(CORE(I010),1,NDSWB,2,IRREP,LISTWB)       
       CALL GETLST(CORE(I020),1,NDSTB,2,IRREP,LISTTB)       
       CALL XGEMM('N','N',DSZZ,NDSZ,NDSWB,ONEM,
     &            CORE(I010),DSZWB,CORE(I020),DSZTB,ONE,
     &            CORE(I000),DSZZ)
       CALL PUTLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
   10 CONTINUE
   20 CONTINUE
C
      DO 40 ISPIN=1,IUHF+1
       LISTZ  = LZOFF + ISPIN + 2
       LISTTA =    36 - ISPIN
       IF(IUHF.EQ.0) LISTTA = 34
       LISTWA =    19 - ISPIN
       LISTTB =    38 - ISPIN
       LISTWB =    18 + ISPIN
C
      DO 30 IRREP=1,NIRREP
C
       DSZZ  = IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NDSZ  = IRPDPD(IRREP,ISYTYP(2,LISTZ))
       DSZTA = IRPDPD(IRREP,ISYTYP(1,LISTTA))
       NDSTA = IRPDPD(IRREP,ISYTYP(2,LISTTA))
       DSZTB = IRPDPD(IRREP,ISYTYP(1,LISTTB))
       NDSTB = IRPDPD(IRREP,ISYTYP(2,LISTTB))
       DSZWA = IRPDPD(IRREP,ISYTYP(1,LISTWA))
       NDSWA = IRPDPD(IRREP,ISYTYP(2,LISTWA))
       DSZWB = IRPDPD(IRREP,ISYTYP(1,LISTWB))
       NDSWB = IRPDPD(IRREP,ISYTYP(2,LISTWB))
C
       I000 = 1
       I010 = I000 + DSZZ  * NDSZ
       I020 = I010 + DSZWA * NDSWA
       I030 = I020 + DSZTA * NDSTA
       NEED = I030 * IINTFP
       IF(NEED.GT.MAXCOR)THEN
        CALL INSMEM('GT2RING',NEED,MAXCOR)
       ENDIF
C
       CALL GETLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
       CALL GETLST(CORE(I010),1,NDSWA,2,IRREP,LISTWA)       
       CALL GETLST(CORE(I020),1,NDSTA,2,IRREP,LISTTA)       
       CALL XGEMM('N','N',DSZZ,NDSZ,NDSWA,ONEM,
     &            CORE(I010),DSZWA,CORE(I020),DSZTA,ONE,
     &            CORE(I000),DSZZ)
C
       I000 = 1
       I010 = I000 + DSZZ  * NDSZ
       I020 = I010 + DSZWB * NDSWB
       I030 = I020 + DSZTB * NDSTB
       NEED = I030 * IINTFP
       IF(NEED.GT.MAXCOR)THEN
        CALL INSMEM('GT2RING',NEED,MAXCOR)
       ENDIF
C
       CALL GETLST(CORE(I010),1,NDSWB,2,IRREP,LISTWB)       
       CALL GETLST(CORE(I020),1,NDSTB,2,IRREP,LISTTB)       
       CALL XGEMM('N','N',DSZZ,NDSZ,NDSWB,ONE ,
     &            CORE(I010),DSZWB,CORE(I020),DSZTB,ONE,
     &            CORE(I000),DSZZ)
       CALL PUTLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
   30 CONTINUE
   40 CONTINUE
C
      DO 60 ISPIN=1,IUHF+1
       LISTZ  = LZOFF + ISPIN + 4
       LISTTA =    37 + ISPIN
       IF(IUHF.EQ.0) LISTTA = 39
       LISTWA =    20 + ISPIN
C
      DO 50 IRREP=1,NIRREP
C
       DSZZ  = IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NDSZ  = IRPDPD(IRREP,ISYTYP(2,LISTZ))
       DSZTA = IRPDPD(IRREP,ISYTYP(1,LISTTA))
       NDSTA = IRPDPD(IRREP,ISYTYP(2,LISTTA))
       DSZWA = IRPDPD(IRREP,ISYTYP(1,LISTWA))
       NDSWA = IRPDPD(IRREP,ISYTYP(2,LISTWA))
C
       I000 = 1
       I010 = I000 + DSZZ  * NDSZ
       I020 = I010 + DSZWA * NDSWA
       I030 = I020 + DSZTA * NDSTA
       NEED = I030 * IINTFP
       IF(NEED.GT.MAXCOR)THEN
        CALL INSMEM('GT2RING',NEED,MAXCOR)
       ENDIF
C
       CALL GETLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
       CALL GETLST(CORE(I010),1,NDSWA,2,IRREP,LISTWA)       
       CALL GETLST(CORE(I020),1,NDSTA,2,IRREP,LISTTA)       
       CALL XGEMM('N','N',DSZZ,NDSZ,NDSWA,ONEM,
     &            CORE(I010),DSZWA,CORE(I020),DSZTA,ONE,
     &            CORE(I000),DSZZ)
       CALL PUTLST(CORE(I000),1,NDSZ ,2,IRREP,LISTZ )       
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
