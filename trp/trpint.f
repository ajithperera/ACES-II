      SUBROUTINE TRPINT(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE DRIVES THE CALCULATION OF APPROPRIATE
C  INTERMEDIATES FOR CCSDT-n AND CCSDT CALCULATIONS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TERM2,YESNO
      CHARACTER*80 FNAME
      DIMENSION ICORE(MAXCOR)
      DIMENSION LSTRNG(6),LENRNG(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/  IFLAGS(100)
C
C-----------------------------------------------------------------------
C     Return for CCSDT-1a and UCC(4).
C-----------------------------------------------------------------------
      IF(IFLAGS(2).EQ.13.OR.IFLAGS(2).EQ.9) RETURN
C
      IF(IFLAGS(2).GE.14)THEN
C
c       CALL GFNAME('DERINT  ',FNAME,ILENGTH)
c       INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
c       IF(YESNO)THEN
        IMODE = 0
c       ELSE
c        IMODE = 1
c       ENDIF
C
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
C-----------------------------------------------------------------------
C     Return for CCSDT-1b and UCC(4).
C-----------------------------------------------------------------------
      IF(IFLAGS(2).EQ.14) RETURN
C
C
C     Create ring lists on 354-359 in EOM-CCSDT calculations and
C     probably some kinds of gradient calculations too. Copy current
C     contents of 54-59 to these. These will be CCSD MBEJ intermediates
C     (spin adapted in RHF case).
C
      CALL IZERO(LENRNG,6)
      DO   10 IRREP=1,NIRREP
      LENRNG(1) = LENRNG(1) + IRPDPD(IRREP, 9) * IRPDPD(IRREP, 9)
      LENRNG(2) = LENRNG(2) + IRPDPD(IRREP,10) * IRPDPD(IRREP,10)
      LENRNG(3) = LENRNG(3) + IRPDPD(IRREP, 9) * IRPDPD(IRREP,10)
      LENRNG(4) = LENRNG(4) + IRPDPD(IRREP,10) * IRPDPD(IRREP, 9)
      LENRNG(5) = LENRNG(5) + IRPDPD(IRREP,11) * IRPDPD(IRREP,11)
      LENRNG(6) = LENRNG(6) + IRPDPD(IRREP,12) * IRPDPD(IRREP,12)
   10 CONTINUE
C
      LSTRNG(1) = 54
      LSTRNG(2) = 55
      LSTRNG(3) = 56
      LSTRNG(4) = 57
      LSTRNG(5) = 58
      LSTRNG(6) = 59
C
cgeneralize
      IF(IFLAGS(2).EQ.18.or.iflags(2).eq.16.or.iflags(2).eq.33.or.
     &                                         iflags(2).eq.34)THEN
        CALL INIPCK(1, 9, 9,354,IMODE,0,1)
        CALL INIPCK(1, 9,10,356,IMODE,0,1)
        CALL INIPCK(1,11,11,358,IMODE,0,1)
       IF(IUHF.NE.0)THEN
        CALL INIPCK(1,10,10,355,IMODE,0,1)
        CALL INIPCK(1,10, 9,357,IMODE,0,1)
        CALL INIPCK(1,12,12,359,IMODE,0,1)
       ENDIF
        DO   20 ILIST=1,6,2-IUHF
        CALL GETALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST))
        CALL PUTALL(ICORE,LENRNG(ILIST),1,LSTRNG(ILIST) + 300)
   20   CONTINUE
      ENDIF
C
C CREATE INTERMEDIATE LISTS FOR CCSDT-2 AND HIGHER METHODS
C
      IF(IFLAGS(2).GE.15)THEN
C
c       CALL GFNAME('GAMLAM  ',FNAME,ILENGTH)
c       INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
c       IF(YESNO)THEN
        IMODE = 0
c       ELSE
c        IMODE = 1
c       ENDIF
C
       call inipck(1,14,18,110,imode,0,1)
       if(iuhf.ne.0)then
        call inipck(1,3,16,107,imode,0,1)
        call inipck(1,4,17,108,imode,0,1)
        call inipck(1,14,11,109,imode,0,1)
       endif
       call inipck(1,13,11,130,imode,0,1)
       if(iuhf.ne.0)then
        call inipck(1,1,9,127,imode,0,1)
        call inipck(1,2,10,128,imode,0,1)
        call inipck(1,13,18,129,imode,0,1)
       endif
C
       IF(IFLAGS(2).EQ.15)THEN
C-----------------------------------------------------------------------
C      CCSDT-2
C-----------------------------------------------------------------------
        CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.FALSE.,.FALSE.,
     &              .FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.)
        TERM2 = .FALSE.
        CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &              .FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.)
       ELSEIF(IFLAGS(2).GE.16)THEN
C-----------------------------------------------------------------------
C      CCSDT-3, CCSDT-4,CCSDT, CC3, CCSDT-T1T2
C-----------------------------------------------------------------------
        IF(IUHF.NE.0)CALL INIPCK(1,10,9,117,IMODE,0,1)
        CALL INIPCK(1,9,10,118,IMODE,0,1) 
        CALL INIPCK(1,9,9,123,IMODE,0,1)
        IF(IUHF.NE.0)CALL INIPCK(1,10,10,124,IMODE,0,1)
        CALL INIPCK(1,11,11,125,IMODE,0,1)
        IF(IUHF.NE.0)CALL INIPCK(1,12,12,126,IMODE,0,1)
C
C     In RHF we must return spin-adapted rings to UHF form and
C     generate list 54.
C
        IF(IUHF.EQ.0) CALL RESET(ICORE,MAXCOR,IUHF)
        CALL FORMWL(ICORE,MAXCOR,IUHF)
        CALL WTWTW (ICORE,MAXCOR,IUHF)
cjdw 3/17/96        CALL FORMW1(ICORE,MAXCOR,IUHF,.TRUE.)
C
         IF(IFLAGS(2).NE.33 .AND. IFLAGS(2).NE.34)THEN
C-----------------------------------------------------------------------
C      CCSDT-3, CCSDT-4,CCSDT
C-----------------------------------------------------------------------
           CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.TRUE.,.TRUE.,
     &                 .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
           TERM2 = .FALSE.
           CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &              .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
           CALL W5T1ABCD(ICORE,MAXCOR,IUHF)
         ELSEIF(IFLAGS(2).EQ.33)THEN
C-----------------------------------------------------------------------
C      CC3
C-----------------------------------------------------------------------
           CALL MODIJKL(ICORE,MAXCOR,IUHF,-1.0D+00)
           CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.TRUE.,.FALSE.,
     &                 .FALSE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.)
           CALL MODIJKL(ICORE,MAXCOR,IUHF, 1.0D+00)
           TERM2 = .FALSE.
           CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &                 .FALSE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.)
           CALL FORMW1(ICORE,MAXCOR,IUHF,.TRUE.)
           CALL W5T1ABCD(ICORE,MAXCOR,IUHF)
           CALL FORMW1(ICORE,MAXCOR,IUHF,.FALSE.)
           CALL CC3W5RING(ICORE,MAXCOR,IUHF)
           CALL CC3W4RING(ICORE,MAXCOR,IUHF)
         ELSEIF(IFLAGS(2).EQ.34)THEN
C-----------------------------------------------------------------------
C      CCSDT-T1T2
C-----------------------------------------------------------------------
           CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.FALSE.,.FALSE.,
     &                 .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
           CALL CC3W4RING(ICORE,MAXCOR,IUHF)
           CALL DHBIAJK4(ICORE,MAXCOR,IUHF,1,1,90,10,106)
C
           TERM2 = .FALSE.
           CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &                 .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
           CALL W5T1ABCD(ICORE,MAXCOR,IUHF)
           CALL CC3W5RING(ICORE,MAXCOR,IUHF)
C
cjdw 3/17/96       IF(IFLAGS(2).NE.18) CALL FORMW1(ICORE,MAXCOR,IUHF,.FALSE.)
C      
         ENDIF
       ENDIF
      ENDIF
C
      IF(IFLAGS(2).EQ.18) CALL FIXFBAR(ICORE,MAXCOR,IUHF)
      IF(IFLAGS(2).EQ.18) CALL FORMW1(ICORE,MAXCOR,IUHF,.TRUE.)
      IF(IFLAGS(2).EQ.18) CALL HBARCRAP(ICORE,MAXCOR,IUHF,.TRUE.)
C
C     N.B. In CCSDT calculations lists 231-233 are reset by TRPS --- the
C          above three lines assume that plain integrals reside on 231-3.
C
      RETURN
      END
