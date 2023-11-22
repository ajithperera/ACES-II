      SUBROUTINE T1BINT(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
C
      WRITE(6,1000)
 1000 FORMAT(' @T1BINT-I, Creating lists for CCSDT-1b Calculations ')
C
C     Subroutine which forms the abci and ijka intermediates which
C     are needed for CCSDT-1B and higher calculations. Part of the
C     work is already handled in CCSD calculations (the FOV-like
C     intermediate). This routine makes use of the MBPT(4) gamma
C     code.
C
C     Assume lists for intermediates have been created. Copy integrals
C     to these lists.
C
      iNeed = -1
      do list = 10-3*iUHF, 10
      do irrep = 1, nirrep
         nRows = aces_list_rows(irrep,list)
         nCols = aces_list_cols(irrep,list)
         iNeed = max(iNeed,nRows*nCols)
      end do
      end do
      do list = 30-3*iUHF, 30
      do irrep = 1, nirrep
         nRows = aces_list_rows(irrep,list)
         nCols = aces_list_cols(irrep,list)
         iNeed = max(iNeed,nRows*nCols)
      end do
      end do
      iNeed = IINTFP*iNeed
      if (iNeed.gt.MAXCOR) call insmem('T1BINT',iNeed,MAXCOR)

      DO   10 IRPKA=1,NIRREP
      IF(IUHF.GT.0)THEN
      IF(IRPDPD(IRPKA,16).GT.0)THEN
      CALL GETLST(ICORE(1),1,IRPDPD(IRPKA,16),2,IRPKA,     7)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPKA,16),2,IRPKA,LWIC21)
      ENDIF
      IF(IRPDPD(IRPKA,17).GT.0)THEN
      CALL GETLST(ICORE(1),1,IRPDPD(IRPKA,17),2,IRPKA,     8)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPKA,17),2,IRPKA,LWIC22)
      ENDIF
      IF(IRPDPD(IRPKA,11).GT.0)THEN
      CALL GETLST(ICORE(1),1,IRPDPD(IRPKA,11),2,IRPKA,     9)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPKA,11),2,IRPKA,LWIC23)
      ENDIF
      ENDIF
      IF(IRPDPD(IRPKA,18).GT.0)THEN
      CALL GETLST(ICORE(1),1,IRPDPD(IRPKA,18),2,IRPKA,    10)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPKA,18),2,IRPKA,LWIC24)
      ENDIF
   10 CONTINUE
C
      DO   20 IRPCI=1,NIRREP
C
      IF(IUHF.GT.0)THEN
      CALL GETLST(ICORE(1),1,IRPDPD(IRPCI, 9),2,IRPCI,    27)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPCI, 9),2,IRPCI,LWIC25)
      CALL GETLST(ICORE(1),1,IRPDPD(IRPCI,10),2,IRPCI,    28)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPCI,10),2,IRPCI,LWIC26)
      CALL GETLST(ICORE(1),1,IRPDPD(IRPCI,18),2,IRPCI,    29)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPCI,18),2,IRPCI,LWIC27)
      ENDIF
      CALL GETLST(ICORE(1),1,IRPDPD(IRPCI,11),2,IRPCI,    30)
      CALL PUTLST(ICORE(1),1,IRPDPD(IRPCI,11),2,IRPCI,LWIC28)
   20 CONTINUE
C
C     Now add W*T1 pieces to form intermediate.
C
      CALL GAMMA5(ICORE,MAXCOR,IUHF)
      CALL GAMMA6(ICORE,MAXCOR,IUHF)
C
      RETURN
      END
