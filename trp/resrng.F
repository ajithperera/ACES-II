      SUBROUTINE RESRNG(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
      LOGICAL CC3FLG,HACKUP
      DIMENSION ICORE(MAXCOR)
      DIMENSION LSTRNG(6),LENRNG(7)
      DIMENSION LNIJKA(4),LNABCI(4),LNIJKL(3),LNABCD(3),LNMBEJ(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /FILES/  LUOUT,MOINTS
C
cc      IF(  IFLAGS(2).EQ.33 .OR.
cc     &    (IFLAGS(2).EQ.22 .AND. IFLAGS2(124).EQ.5) )THEN
cc        HACKUP = .TRUE.
cc      ELSE
c         HACKUP = .FALSE.
cc      ENDIF
c         HACKUP = IFLAGS(2) .EQ. 34
cC
cC     HACKUP = .TRUE. means linearized CC3 approach. If we set HACKUP
cC     to .FALSE. we should get normal stuff (we better). See also
cC     modifications in TRPINT and possibly EOMNT3.
cC     Lists 107-110, 127-130 get W*(1+T1), copied to 307-310, 327-330.
cC     Lists 331-333 get plain <ab||cd> integrals.
cC     Lists 351-353 get plain <ij||kl> integrals.
cC     Lists 354-359 get plain <mb||ej> integrals.
cC     Lists 364-369 get plain <mb||ej> integrals.
C
      IF(IFLAGS(2).NE.10.AND.IFLAGS(2).NE.22.AND.IFLAGS(2).NE.13.AND.
     &   IFLAGS(2).NE.14.AND.IFLAGS(2).NE.16.AND.IFLAGS(2).NE.33.AND.
     &   IFLAGS(2).NE.5 .AND.IFLAGS(2).NE.6 .AND.IFLAGS(2).NE.16)THEN
      WRITE(LUOUT,1010)
 1010 FORMAT(' @RESRNG-F, Code not implemented for method ',I3)
      CALL ERREX
      ENDIF
C
      CALL IZERO(LNIJKA,4)
      CALL IZERO(LNABCI,4)
      CALL IZERO(LNIJKL,3)
      CALL IZERO(LNABCD,3)
      CALL IZERO(LNMBEJ,6)
      DO   50 IRREP=1,NIRREP
      LNIJKA(1)=LNIJKA(1) + IRPDPD(IRREP, 3) * IRPDPD(IRREP, 9)
      LNIJKA(2)=LNIJKA(2) + IRPDPD(IRREP, 4) * IRPDPD(IRREP,10)
      LNIJKA(3)=LNIJKA(3) + IRPDPD(IRREP,14) * IRPDPD(IRREP,11)
      LNIJKA(4)=LNIJKA(4) + IRPDPD(IRREP,14) * IRPDPD(IRREP,12)
C
      LNABCI(1)=LNABCI(1) + IRPDPD(IRREP, 1) * IRPDPD(IRREP, 9)
      LNABCI(2)=LNABCI(2) + IRPDPD(IRREP, 2) * IRPDPD(IRREP,10)
      LNABCI(3)=LNABCI(3) + IRPDPD(IRREP,13) * IRPDPD(IRREP,11)
      LNABCI(4)=LNABCI(4) + IRPDPD(IRREP,13) * IRPDPD(IRREP,12)
C
      LNIJKL(1)=LNIJKL(1) + IRPDPD(IRREP, 3) * IRPDPD(IRREP, 3)
      LNIJKL(2)=LNIJKL(2) + IRPDPD(IRREP, 4) * IRPDPD(IRREP, 4)
      LNIJKL(3)=LNIJKL(3) + IRPDPD(IRREP,14) * IRPDPD(IRREP,14)
C
      LNABCD(1)=LNABCD(1) + IRPDPD(IRREP, 1) * IRPDPD(IRREP, 1)
      LNABCD(2)=LNABCD(2) + IRPDPD(IRREP, 2) * IRPDPD(IRREP, 2)
      LNABCD(3)=LNABCD(3) + IRPDPD(IRREP,13) * IRPDPD(IRREP,13)
C
      LNMBEJ(1)=LNMBEJ(1) + IRPDPD(IRREP, 9) * IRPDPD(IRREP, 9)
      LNMBEJ(2)=LNMBEJ(2) + IRPDPD(IRREP,10) * IRPDPD(IRREP,10)
      LNMBEJ(3)=LNMBEJ(3) + IRPDPD(IRREP, 9) * IRPDPD(IRREP,10)
      LNMBEJ(4)=LNMBEJ(4) + IRPDPD(IRREP,10) * IRPDPD(IRREP, 9)
      LNMBEJ(5)=LNMBEJ(5) + IRPDPD(IRREP,11) * IRPDPD(IRREP,11)
      LNMBEJ(6)=LNMBEJ(6) + IRPDPD(IRREP,12) * IRPDPD(IRREP,12)
   50 CONTINUE
C
      CALL IZERO(LENRNG,7)
      DO   10 IRREP=1,NIRREP
      LENRNG(1) = LENRNG(1) + IRPDPD(IRREP, 9) * IRPDPD(IRREP, 9)
      LENRNG(2) = LENRNG(2) + IRPDPD(IRREP,10) * IRPDPD(IRREP,10)
      LENRNG(3) = LENRNG(3) + IRPDPD(IRREP, 9) * IRPDPD(IRREP,10)
      LENRNG(4) = LENRNG(4) + IRPDPD(IRREP,10) * IRPDPD(IRREP, 9)
      LENRNG(5) = LENRNG(5) + IRPDPD(IRREP,11) * IRPDPD(IRREP,11)
      LENRNG(6) = LENRNG(6) + IRPDPD(IRREP,12) * IRPDPD(IRREP,12)
      LENRNG(7) = LENRNG(7) + IRPDPD(IRREP,19) * IRPDPD(IRREP,9)
   10 CONTINUE
C
      LSTRNG(1) = 54
      LSTRNG(2) = 55
      LSTRNG(3) = 56
      LSTRNG(4) = 57
      LSTRNG(5) = 58
      LSTRNG(6) = 59
C
C-----------------------------------------------------------------------
C     Logic for CCSDT-1a, CCSDT-1b, CCSD, NCCSDT-1, F-NCCSDT-1 (includes
C     old EOM-CCSD(T)), NCCSDT-1b, F-NCCSDT-1b.
C
C     Put raw ijka and abci integrals on lists 307-310 and 327-330.
C-----------------------------------------------------------------------
C
#ifdef _DEBUG_LVL0

      Print*, "Entering iflags(2)", iflags(2)
#endif
      IF( IFLAGS(2).EQ.13. OR.IFLAGS(2).EQ.14.OR.IFLAGS(2).EQ.10.OR.
     &    IFLAGS(2).EQ.5.OR.IFLAGS(2).EQ.6.OR.
     &   (IFLAGS(2).EQ.22.AND.
     &    IFLAGS2(124).GE.1.AND.IFLAGS2(124).LE.4))THEN
C
       IF(IFLAGS(2).NE.14)THEN
C#ifdef _EOMTRPS
C         IMODE = 1
c#else
c        IMODE=1
        call aces_io_remove(53,'DERINT')
        IMODE=0
C#endif
C
#ifdef _DEBUG_LVL0
        Print*, "Entering iflags(2)", iflags(2)
#endif 
        if(iuhf.ne.0)then
         call inipck(1,3,16,307,imode,0,1)
         call inipck(1,4,17,308,imode,0,1)
         call inipck(1,14,11,309,imode,0,1)
        endif
        call inipck(1,14,18,310,imode,0,1)
c
        if(iuhf.ne.0)then
         call inipck(1,1,9,327,imode,0,1)
         call inipck(1,2,10,328,imode,0,1)
         call inipck(1,13,18,329,imode,0,1)
        endif
        call inipck(1,13,11,330,imode,0,1)
       ENDIF
C
       CALL GETALL(ICORE,LNIJKA(4),1, 10)
       CALL PUTALL(ICORE,LNIJKA(4),1, 10 + 300)
       IF(IUHF.GT.0)THEN
        CALL GETALL(ICORE,LNIJKA(1),1,  7)
        CALL PUTALL(ICORE,LNIJKA(1),1,  7 + 300)
        CALL GETALL(ICORE,LNIJKA(2),1,  8)
        CALL PUTALL(ICORE,LNIJKA(2),1,  8 + 300)
        CALL GETALL(ICORE,LNIJKA(3),1,  9)
        CALL PUTALL(ICORE,LNIJKA(3),1,  9 + 300)
       ENDIF
C
       CALL GETALL(ICORE,LNABCI(4),1, 30)
       CALL PUTALL(ICORE,LNABCI(4),1, 30 + 300)
       IF(IUHF.GT.0)THEN
        CALL GETALL(ICORE,LNABCI(1),1, 27)
        CALL PUTALL(ICORE,LNABCI(1),1, 27 + 300)
        CALL GETALL(ICORE,LNABCI(2),1, 28)
        CALL PUTALL(ICORE,LNABCI(2),1, 28 + 300)
        CALL GETALL(ICORE,LNABCI(3),1, 29)
        CALL PUTALL(ICORE,LNABCI(3),1, 29 + 300)
       ENDIF
C
C     Note : Need memory checks in above !!!
C
      Write(6,*) "Returned from here for non-it methods"
      RETURN
      ENDIF
C
C-----------------------------------------------------------------------
C     Logic for iterative methods CCSDT-T1T2, CC3, CCSDT-3, CCSDT and
C     non-iterative methods CCSD(T'), CCSDR(T), CCSDR(3), NCCSDT-3.
C-----------------------------------------------------------------------
C
C     For CCSDT-T1T2, CC3, CCSDT-3, and CCSDT W(abci) and W(mcjk)
C     intermediates are on lists 107-100 and 127-130. For CCSDT-T1T2
C     and CC3 these are exactly the quantities which will be contracted
C     with R2 in vee. For CCSDT-3 and CCSDT we
C     must add F*T2 to W(abci). The contents of 107-110 and 127-130
C     are then copied to 307-310 and 327-330.
C
C     For the CCSD(T'), CCSDR(T), CCSDR(3), and NCCSDT-3 methods, there is
C     nothing on lists 107-110, 127-130 at this time. They are created
C     by TRPINT.
C
C     CCSD(T') : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) =  5 (CCSD[TPR])
C     CCSD(T') : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) =  6 (F-CCSD[TPR])
C     CCSDR(T) : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) =  7 (CCSDR[T])
C     CCSDR(T) : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) =  8 (F-CCSDR[T])
C     CCSDR(3) : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) =  9 (CCSDR[3])
C     CCSDR(3) : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) = 10 (F-CCSDR[3])
C     NCCSDT-3 : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) = 11 (NCCSDT-3)
C     NCCSDT-3 : IFLAGS(2) = 22 (CCSD[T]); IFLAGS2(124) = 12 (F-NCCSDT-3)
C
      IF(IFLAGS(2).EQ.22.AND.
     &   (IFLAGS2(124).EQ. 5 .OR.IFLAGS2(124).EQ. 6))THEN
          IFLAGS(2) = 34
          CALL TRPINT(ICORE,MAXCOR,IUHF)
          IFLAGS(2) = 22
      ENDIF
C
      IF(IFLAGS(2).EQ.22.AND.
     &   (IFLAGS2(124).GE. 7.AND.IFLAGS2(124).LE.10))THEN
          IFLAGS(2) = 33
          CALL TRPINT(ICORE,MAXCOR,IUHF)
          IFLAGS(2) = 22
      ENDIF
C
      IF(IFLAGS(2).EQ.22.AND.
     &   (IFLAGS2(124).EQ.11.OR.IFLAGS2(124).EQ.12))THEN
          IFLAGS(2) = 16
          CALL TRPINT(ICORE,MAXCOR,IUHF)
          IFLAGS(2) = 22
      ENDIF
C
C-----------------------------------------------------------------------
C     Ugly digression :
C
C     Logic primarily for handling ring lists in CCSDT-3,
C     CCSDT-4, CCSDT, CC3. The lists are put back so that lambda gets
C     what it expects on 54-59.
C-----------------------------------------------------------------------
C
cjdw 3/31/96. reset is done by trpint ?????? how did this work before ?
      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.17.OR.IFLAGS(2).EQ.18.OR.
     &   IFLAGS(2).EQ.33.OR.IFLAGS(2).EQ.34)THEN
c      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.17.OR.IFLAGS(2).EQ.18.OR.
c     &   IFLAGS(2).EQ.33.OR.(IFLAGS(2).EQ.22.AND.IFLAGS2(124).GE.5))THEN
C
C     In RHF move list 56 to list 354.
C
        IF(IUHF.EQ.0)THEN
          CALL GETALL(ICORE,LENRNG(1),1,LSTRNG(3))
          Print*, "First getall 1" 
          CALL PUTALL(ICORE,LENRNG(1),1,LSTRNG(1) + 300)
        ENDIF
C
        DO   20 ILIST=1,6,2-IUHF
        CALL GETALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST) + 300)
        Print*, "First getall 2" 
        CALL PUTALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST))
   20   CONTINUE
C
C     Now copy list 354 to 356.
C
        IF(IUHF.EQ.0)THEN
          CALL GETALL(ICORE,LENRNG(1),1,LSTRNG(1)+300)
          CALL PUTALL(ICORE,LENRNG(1),1,LSTRNG(3)+300)
C 130-->330 WHY ??????????????????????
          CALL GETALL(ICORE,LENRNG(7),1,130)
          CALL PUTALL(ICORE,LENRNG(7),1,330)
        ENDIF
C
C     In RHF we must return spin-adapted rings to UHF form and
C     generate list 54.
C
        IF(IUHF.EQ.0) CALL RESET(ICORE,MAXCOR,IUHF)
C
      ENDIF
C
      IF(IFLAGS(2).EQ.22.AND.IFLAGS2(124).GE.5)THEN
        DO   40 ILIST=1,6,2-IUHF
        CALL GETALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST) + 300)
        CALL PUTALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST))
   40   CONTINUE
        IF(IUHF.EQ.0) CALL RESET(ICORE,MAXCOR,IUHF)
      ENDIF
C-----------------------------------------------------------------------
C     End of Ugly Digression (but certainly not the end of the ugly).
C-----------------------------------------------------------------------
C
      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.33.OR.IFLAGS(2).EQ.34.OR.
     &   IFLAGS(2).EQ.22)THEN
C
C     For the CCSD(T) variants, set IFLAGS(2) to corresponding iterative
C     method. It will be reset at the end.
C
      IF(IFLAGS(2).EQ.22)THEN
       IFLOLD = 22
       IF(    IFLAGS2(124).GE. 5.AND.IFLAGS2(124).LE. 6)THEN
        IFLAGS(2)=34
       ELSEIF(IFLAGS2(124).GE. 7.AND.IFLAGS2(124).LE.10)THEN
        IFLAGS(2)=33
       ELSEIF(IFLAGS2(124).GE.11.AND.IFLAGS2(124).LE.12)THEN
        IFLAGS(2)=16
       ENDIF
      ENDIF
C
      CC3FLG = IFLAGS(2).EQ.33
C
      WRITE(6,*) ' @RESRNG-I, Creating lists for CC3/CCSDT-3. '
C
C     W(mcjk) (CC3) already exist on 107-110. Put on 307-310.
C     W(abci) (CC3) already exist on 127-130. Put on 327-330.
C
C     Add additional term to W(abci) for CCSDT-3 or CCSDT.
C
      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.18)THEN
       CALL GWR1T2(ICORE,MAXCOR,IUHF,1,90,126,106,.TRUE.,.FALSE.,
     &             .TRUE.,1.0D+00,0.0D+00)
      ENDIF
C
        IF(IUHF.EQ.0)THEN
          IBOT = 4
        ELSE
          IBOT= 1
        ENDIF
C
      DO 200 ILIST=IBOT,4
C
      LENGTH = LNABCI(ILIST)
      NEED = IINTFP * LENGTH
C
        IF(NEED.GT.MAXCOR)THEN
          CALL INSMEM('RESRNG',NEED,MAXCOR)
        ELSE
          CALL GETALL(ICORE,LENGTH,1,126 + ILIST)
          CALL PUTALL(ICORE,LENGTH,1,326 + ILIST)
        ENDIF
C
      LENGTH = LNIJKA(ILIST)
      NEED = IINTFP * LENGTH
C
        IF(NEED.GT.MAXCOR)THEN
          CALL INSMEM('RESRNG',NEED,MAXCOR)
        ELSE
          CALL GETALL(ICORE,LENGTH,1,106 + ILIST)
          CALL PUTALL(ICORE,LENGTH,1,306 + ILIST)
        ENDIF
C
  200 CONTINUE
      WRITE(6,*) ' @RESRNG-I, mcjk lists written on 307-310. '
      WRITE(6,*) ' @RESRNG-I, abci lists written on 327-330. '
C
C     W(ijkl) already exist on 51-53. Remove T_2 piece and then put
C     on 351-353. At the end put T_2 piece back in lists 51-53. That
C     is, lists 51-53 end up with normal W(ijkl) Hbar element, while
C     lists 351-353 end up with CC3 W(ijkl).
C
      IF(CC3FLG)THEN
        CALL MODIJKL(ICORE,MAXCOR,IUHF,-1.0D+00)
      ENDIF
C
      IMODE = 0
      IF(IUHF.NE.0)THEN
        CALL INIPCK(1, 3, 3,351,imode,0,1)
        CALL INIPCK(1, 4, 4,352,imode,0,1)
      ENDIF
      CALL INIPCK(1,14,14,353,imode,0,1)
C
      IF(IUHF.EQ.0)THEN
        IBOT=3
      ELSE
        IBOT=1
      ENDIF
C
      DO 210 ILIST=IBOT,3
C
      LENGTH = LNIJKL(ILIST)
      NEED = IINTFP * LENGTH
C
        IF(NEED.GT.MAXCOR)THEN
          CALL INSMEM('RESRNG',NEED,MAXCOR)
        ELSE
          IF(IFLAGS(2).EQ.34)THEN
            CALL GETALL(ICORE,LENGTH,1, 10 + ILIST)
            CALL PUTALL(ICORE,LENGTH,1,350 + ILIST)
          ELSE
            CALL GETALL(ICORE,LENGTH,1, 50 + ILIST)
            CALL PUTALL(ICORE,LENGTH,1,350 + ILIST)
          ENDIF
        ENDIF
C
  210 CONTINUE
C
      IF(CC3FLG)THEN
        CALL MODIJKL(ICORE,MAXCOR,IUHF, 1.0D+00)
      ENDIF
      WRITE(6,*) ' @RESRNG-I, ijkl lists written on 351-353. '
C
C-----------------------------------------------------------------------
C     Create CC3 W(abcd) lists.
C-----------------------------------------------------------------------
C
      IF(CC3FLG.OR.IFLAGS(2).EQ.34 .or. iflags2(124).ge.11)THEN
C
       IMODE = 0
       IF(IUHF.NE.0)THEN
         CALL INIPCK(1, 1, 1,331,imode,0,1)
         CALL INIPCK(1, 2, 2,332,imode,0,1)
       ENDIF
       CALL INIPCK(1,13,13,333,imode,0,1)
C
       IF(IUHF.EQ.0)THEN
         IBOT=3
       ELSE
         IBOT=1
       ENDIF
C
       DO 220 ILIST=IBOT,3
C
       LENGTH = LNABCD(ILIST)
       NEED = IINTFP * LENGTH
C
         IF(NEED.GT.MAXCOR)THEN
C
           DO 215 IRREP=1,NIRREP
           NDS = IRPDPD(IRREP,ISYTYP(2,230 + ILIST))
           DO 214 IDIS=1,NDS
            CALL GETLST(ICORE,IDIS,1,2,IRREP,230+ILIST)
            CALL PUTLST(ICORE,IDIS,1,2,IRREP,330+ILIST)
  214      CONTINUE
  215      CONTINUE
C
         ELSE
           CALL GETALL(ICORE,LENGTH,1,230 + ILIST)
           CALL PUTALL(ICORE,LENGTH,1,330 + ILIST)
         ENDIF
C
  220  CONTINUE
C
       IF(CC3FLG .or. iflags2(124).ge.11)THEN
         CALL   GFORMW1(ICORE,MAXCOR,IUHF,.TRUE.,330)
         CALL GHBARABCD(ICORE,MAXCOR,IUHF,       330,CC3FLG)
       ENDIF
       WRITE(6,*) ' @RESRNG-I, abcd lists written on 331-333. '
C
      ENDIF
C
C-----------------------------------------------------------------------
C     Create two types of W(mbej) lists.
C-----------------------------------------------------------------------
C
      IMODE = 0
      IF(IUHF.NE.0)THEN
        CALL INIPCK(1,10,10,355,imode,0,1)
        CALL INIPCK(1,10, 9,357,imode,0,1)
        CALL INIPCK(1,12,12,359,imode,0,1)
        CALL INIPCK(1,10,10,365,imode,0,1)
        CALL INIPCK(1,10, 9,367,imode,0,1)
        CALL INIPCK(1,12,12,369,imode,0,1)
      ENDIF
        CALL INIPCK(1, 9, 9,354,imode,0,1)
        CALL INIPCK(1, 9,10,356,imode,0,1)
        CALL INIPCK(1,11,11,358,imode,0,1)
        CALL INIPCK(1, 9, 9,364,imode,0,1)
        CALL INIPCK(1, 9,10,366,imode,0,1)
        CALL INIPCK(1,11,11,368,imode,0,1)
C
        CALL ZEROLIST(ICORE,MAXCOR,354)
        CALL ZEROLIST(ICORE,MAXCOR,356)
        CALL ZEROLIST(ICORE,MAXCOR,358)
        CALL ZEROLIST(ICORE,MAXCOR,364)
        CALL ZEROLIST(ICORE,MAXCOR,366)
        CALL ZEROLIST(ICORE,MAXCOR,368)
C
      IF(IFLAGS(2).NE.34)THEN
C
        CALL GT1RING(ICORE,MAXCOR,IUHF,.FALSE.,300,.TRUE.,.FALSE.,
     &               .FALSE.,.TRUE.)
        WRITE(6,*) ' @RESRNG-I, T1 contributions written on 354-359 '
        CALL GT1RING(ICORE,MAXCOR,IUHF,.FALSE.,310,.FALSE.,.TRUE.,
     &               .FALSE.,.TRUE.)
        WRITE(6,*) ' @RESRNG-I, T1 contributions written on 364-369 '
C
      ENDIF
C
      DO 230 ISPIN=1,1+IUHF
C
      CALL SAXLST(ICORE,MAXCOR, 22+ISPIN,353+ISPIN,353+ISPIN,
     &            1.0D+00,1.0D+00)
      CALL SAXLST(ICORE,MAXCOR, 22+ISPIN,363+ISPIN,363+ISPIN,
     &            1.0D+00,1.0D+00)
C
      CALL SAXLST(ICORE,MAXCOR, 19-ISPIN,355+ISPIN,355+ISPIN,
     &            1.0D+00,1.0D+00)
      CALL SAXLST(ICORE,MAXCOR, 19-ISPIN,365+ISPIN,365+ISPIN,
     &            1.0D+00,1.0D+00)
C
      CALL SAXLST(ICORE,MAXCOR, 24+ISPIN,357+ISPIN,357+ISPIN,
     &            1.0D+00,1.0D+00)
      CALL SAXLST(ICORE,MAXCOR, 24+ISPIN,367+ISPIN,367+ISPIN,
     &            1.0D+00,1.0D+00)
C
  230 CONTINUE
      WRITE(6,*) ' @RESRNG-I, integrals added to lists 354-359. '
      WRITE(6,*) ' @RESRNG-I, integrals added to lists 364-369. '
C
C     ----- Add T2 contributions to 354-359, 364-369. Not needed -----
C     ----- for CC3. Just testing now prior to CCSDT-3.          -----
C
      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.18)THEN
       CALL GT2RING(ICORE,MAXCOR,IUHF,353)
       CALL GT2RING(ICORE,MAXCOR,IUHF,363)
       WRITE(6,*) ' @RESRNG-I, T2 term added to lists 354-359. '
       WRITE(6,*) ' @RESRNG-I, T2 term added to lists 364-369. '
      ENDIF
C
C     ----- Zero out a series of lists which are not used in CCSD -----
C     ----- energies.                                             -----
C
      IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,117)
      CALL ZEROLIST(ICORE,MAXCOR,118)
      CALL ZEROLIST(ICORE,MAXCOR,123)
      IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,124)
      CALL ZEROLIST(ICORE,MAXCOR,125)
      IF(IUHF.NE.0) CALL ZEROLIST(ICORE,MAXCOR,126)
C
      IF(IUHF.NE.0)THEN
       CALL ZEROLIST(ICORE,MAXCOR,107)
       CALL ZEROLIST(ICORE,MAXCOR,108)
       CALL ZEROLIST(ICORE,MAXCOR,109)
      ENDIF
       CALL ZEROLIST(ICORE,MAXCOR,110)
C
      IF(IUHF.NE.0)THEN
       CALL ZEROLIST(ICORE,MAXCOR,127)
       CALL ZEROLIST(ICORE,MAXCOR,128)
       CALL ZEROLIST(ICORE,MAXCOR,129)
      ENDIF
       CALL ZEROLIST(ICORE,MAXCOR,130)
C
      IF(IFLOLD.EQ.22)THEN
       IFLAGS(2)=22
      ENDIF
C
      RETURN
      ENDIF
C
      RETURN
      END
