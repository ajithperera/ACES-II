      SUBROUTINE HBARDIAG(IRREPX,SCR,MAXCOR,IUHF,FREQ)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,SCR(MAXCOR),FACT,FREQ
      LOGICAL DOUBLE,NONSTD,NONSTD2,SS, SD, DS, DD
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/GUESS/DOUBLE,NONSTD,NONSTD2
      COMMON/DRVHBAR/SS, SD, DS, DD
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      IF (DD) THEN
        IF(IUHF.EQ.0.AND.(IFLAGS(93).NE.2.AND.ISYTYP(1,233).NE.5))THEN
          CALL MKFDENOM(SCR,MAXCOR/IINTFP,IUHF,IRREPX)
        ELSE
          CALL DDENOM(SCR,MAXCOR/IINTFP,IUHF,IRREPX,447)
        ENDIF
      ELSE
        CALL DDENOM(SCR,MAXCOR/IINTFP,IUHF,IRREPX,447)
      ENDIF
C
C SINGLES PART
C
      IOFFC=1
      DO 10 ISPIN=1,1+IUHF
       NDIM=IRPDPD(IRREPX,8+ISPIN)
       ISCR=IOFFC+NDIM
       CALL GETLST(SCR(IOFFC),1,1,1,9,447+ISPIN)
       CALL VMINUS(SCR(IOFFC),NDIM)
       IOFFC=ISCR
10    CONTINUE
C
C DOUBLES PART
C
      DO 20 ISPIN=3,3-2*IUHF,-1
       LISTD=447+ISPIN
       DO 30 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LISTD))
        DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LISTD))
        CALL GETLST(SCR(IOFFC),1,NUMDIS,1,IRREPR,LISTD)
        CALL VMINUS(SCR(IOFFC),NUMDIS*DISSIZ)
        IOFFC=IOFFC+NUMDIS*DISSIZ
30     CONTINUE
20    CONTINUE
C
      IF(FREQ.NE.ZILCH)THEN
       DO 101 I=1,NSIZEC
        SCR(I)=SCR(I)+FREQ
101    CONTINUE
      ENDIF 
      
      CALL PUTLST(SCR,1,1,1,1,472)
C
      RETURN
      END
