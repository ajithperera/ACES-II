      SUBROUTINE TRPINT(ICORE,MAXCOR,IUHF,METACT)
C
C THIS ROUTINE DRIVES THE CALCULATION OF APPROPRIATE
C  INTERMEDIATES FOR CCSDT-n AND CCSDT CALCULATIONS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TERM2
      DIMENSION ICORE(MAXCOR)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/ IFLAGS(100)
C
C SIMPLY RETURN FOR CCSDT-1 AND UCC(4) METHODS
C
      IF(IFLAGS(2).EQ.13.OR.IFLAGS(2).EQ.9)THEN
       RETURN
      endif
      IF(IFLAGS(2).GE.14)THEN
c       imode=1
       call aces_io_remove(53,'DERINT')
       imode=0
       call inipck(1,3,16,307,imode,0,1)
       if(iuhf.ne.0)then
        call inipck(1,4,17,308,imode,0,1)
        call inipck(1,14,11,309,imode,0,1)
       endif
       call inipck(1,14,18,310,imode,0,1)
c
       call inipck(1,1,9,327,imode,0,1)
       if(iuhf.ne.0)then
        call inipck(1,2,10,328,imode,0,1)
        call inipck(1,13,18,329,imode,0,1)
       endif
       call inipck(1,13,11,330,imode,0,1)
      endif
c
       if(iflags(2).eq.14)return
C
C CREATE INTERMEDIATE LISTS FOR CCSDT-2 AND HIGHER METHODS
C
      IF(IFLAGS(2).GE.15)THEN
c       imode=1
       call aces_io_remove(51,'GAMLAM')
       imode=0
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
       IF(IFLAGS(2).EQ.15)THEN
        CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.FALSE.,.FALSE.,
     &              .FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.)
CJDW 2/7/95
CAP  5/6/96
        TERM2 = .FALSE.
        CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &              .FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.)
CJDW 3/21/95
C      Specific to this version of TRPINT.
C
        IF(METACT.EQ.29.OR.METACT.EQ.30)THEN
         IF(IUHF.GT.0)THEN
           CALL SAXLST(ICORE,MAXCOR,107, 7,107,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,108, 8,108,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,109, 9,109,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,127,27,127,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,128,28,128,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,129,29,129,1.0D+00,-1.0D+00)
         ENDIF
           CALL SAXLST(ICORE,MAXCOR,110,10,110,1.0D+00,-1.0D+00)
           CALL SAXLST(ICORE,MAXCOR,130,30,130,1.0D+00,-1.0D+00)
        ENDIF

       ELSEIF(IFLAGS(2).GE.16)THEN
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
        CALL FORMW1(ICORE,MAXCOR,IUHF,.TRUE.)
        CALL FORMW4(ICORE,MAXCOR,IUHF,.TRUE.,.TRUE.,.TRUE.,
     &              .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
        TERM2 = .TRUE.
CJDW 2/7/95
CAP  5/6/96
        CALL FORMW5(ICORE,MAXCOR,IUHF,.TRUE.,TERM2,.FALSE.,
     &              .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
       IF(IFLAGS(2).NE.18) CALL FORMW1(ICORE,MAXCOR,IUHF,.FALSE.)
       ENDIF
C      
      ENDIF
C
      IF(IFLAGS(2).EQ.18) CALL FIXFBAR(ICORE,MAXCOR,IUHF)
      IF(IFLAGS(2).EQ.18) CALL HBARCRAP(ICORE,MAXCOR,IUHF,.TRUE.)
C
C     In CCSDT calculations lists 231-233 are reset by TRPS.
C
      RETURN
      END
