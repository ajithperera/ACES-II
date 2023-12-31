      SUBROUTINE T1INT1(ICORE,MAXCOR,IUHF,ISPIN)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTION OF LINEAR T1 IN THE T1
C  EQUATION.
C
C         T1INC(A,I) = - T1(E,M) * <MA||IE> + T1(e,m) * <mI|eA> (AA)
C
C         T1INC(a,i) = - T1(e,m) * <ma||ie> + T1(E,M) * <Mi|Ea> (BB)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL ROHF4,ITRFLG,bRedundant
      LOGICAL IJAB_DEBUG
      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /FLAGS2/ IFLAGS2(500)
      DATA ONE /1.0/ 
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
#ifdef _DEBUG_LVL0
CSSS      IJAB_DEBUG=.FALSE.
CSSS      Write(6,"(a,a)") " @-T1INT1: List 121 and 122 instead of 17 and",
CSSS     +                 " 18 is read for <MB||EJ> Ints"
#endif 
      NTOT1=0
      NTOT2=0
      MAXSIZ=0
      LISTT=90
      LISTW1=22+ISPIN
      LISTW2=16+ISPIN
      IF (IJAB_DEBUG) LISTW2=120+ISPIN 
      bRedundant = IFLAGS2(155).EQ.0
      IF(IUHF.EQ.0)LISTW2=18
CSSS      IF (IJAB_DEBUG .AND. IUHF.EQ.0) LISTW2=122
C
C WE ONLY NEED TO CONSIDER ONE DPD IRREP -- THE TOTALLY SYMMETRIC ONE.
C
      NTOTTAR=IRPDPD(1,8+ISPIN)
      NTOTTA =IRPDPD(1,8+ISPIN)
      NTOTTB =IRPDPD(1,11-ISPIN)
C
C I000 HOLDS T1INC TARGET IN EQUATION ABOVE.
C I010 HOLDS FIRST T IN EQUATION ABOVE.
C I020 HOLDS SECOND T IN EQUATION ABOVE.
C I030 HOLDS INDIVIDUAL A1 IRREP OF BLOCKED INTEGRAL LIST.
C
      I000=1
      I010=I000+NTOTTAR*IINTFP
      IF(IUHF.EQ.0)THEN
       I020=I010
      ELSE
       I020=I010+NTOTTA*IINTFP
      ENDIF
      I030=I020+NTOTTB*IINTFP
      IF(ROHF4.AND..NOT.ITRFLG) THEN
        CALL GETLST(ICORE(I010),1,1,1,8+ISPIN,LISTT)
        IF(IUHF.NE.0)CALL GETLST(ICORE(I020),1,1,1,11-ISPIN,LISTT)
      ELSE
        CALL GETLST(ICORE(I010),1,1,1,ISPIN,LISTT)
        IF(IUHF.NE.0)CALL GETLST(ICORE(I020),1,1,1,3-ISPIN,LISTT)
      ENDIF
      NWDSZ1=IRPDPD(1,ISYTYP(1,LISTW1))
      NWDSZ2=IRPDPD(1,ISYTYP(1,LISTW2))
      NWDIS1=IRPDPD(1,ISYTYP(2,LISTW1))
      NWDIS2=IRPDPD(1,ISYTYP(2,LISTW2))
      I040=I030+NWDSZ1*NWDIS1*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('T1INT1',I040,MAXCOR)
      CALL GETLST(ICORE(I030),1,NWDIS1,2,1,LISTW1)
      CALL XGEMM('N','N',1,NTOTTAR,NTOTTA,ONEM,ICORE(I010),1,
     &           ICORE(I030),NTOTTA,ZILCH,ICORE(I000),1)
      I040=I030+NWDSZ2*NWDIS2*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('T1INT1',I040,MAXCOR)
      IF(bRedundant) THEN
        CALL GETLST(ICORE(I030),1,NWDIS2,2,1,LISTW2)
CSSS        call Izero(ICore(i030),NWDIS2*NWDSZ2)
CSSS        call checksum("t1t1",icore(i030),NWDIS2*NWDSZ2)
CSSS        Write(6,*) "intentionally making t1int1 zero"
      ELSE 
        CALL GETLST_NR(ICORE(I030),ICORE(I040),MAXCOR-I040,
     &                 LISTW2, 1)
      ENDIF    
      CALL XGEMM('N','N',1,NTOTTAR,NTOTTB,ONE,ICORE(I020),1,
     &           ICORE(I030),NTOTTB,ONE,ICORE(I000),1)
      CALL UPDAT1(ICORE(I000),ICORE(I010),NTOTTAR,2+ISPIN,90)
      RETURN
      END
