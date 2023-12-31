C
C THIS ROUTINE COMPUTES THE CONTRIBUTION OF LINEAR L1 IN THE L1
C  EQUATION.
C
C         L1INC(A,I) = - L1(E,M) * <MA||IE> + L1(e,m) * <mI|eA> (AA)
C
C         L1INC(a,i) = - L1(e,m) * <ma||ie> + L1(E,M) * <Mi|Ea> (BB)
C
      SUBROUTINE L1INL1(ICORE,MAXCOR,IUHF,ISPIN)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      DATA ONE /1.0/ 
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
      NTOT1=0
      NTOT2=0
      MAXSIZ=0
      LISTT=190
      LISTW1=53+ISPIN
      LISTW2=55+ISPIN
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
      CALL GETLST(ICORE(I010),1,1,1,ISPIN,LISTT)
      IF(IUHF.NE.0)CALL GETLST(ICORE(I020),1,1,1,3-ISPIN,LISTT)
      NWDSZ1=IRPDPD(1,ISYTYP(1,LISTW1))
      NWDSZ2=IRPDPD(1,ISYTYP(1,LISTW2))
      NWDIS1=IRPDPD(1,ISYTYP(2,LISTW1))
      NWDIS2=IRPDPD(1,ISYTYP(2,LISTW2))
      I040=I030+NWDSZ1*NWDIS1*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('L1INL1',I040,MAXCOR)
      CALL GETLST(ICORE(I030),1,NWDIS1,2,1,LISTW1)
      CALL XGEMM('N','T',1,NTOTTAR,NTOTTA,ONEM,ICORE(I010),1,
     &           ICORE(I030),NTOTTAR,ZILCH,ICORE(I000),1)
      I040=I030+NWDSZ2*NWDIS2*IINTFP
      IF(I040.GT.MAXCOR)CALL INSMEM('L1INL1',I040,MAXCOR)

      CALL GETLST(ICORE(I030),1,NWDIS2,2,1,LISTW2)
      CALL XGEMM('N','T',1,NTOTTAR,NTOTTB,ONE,ICORE(I020),1,
     &           ICORE(I030),NTOTTAR,ONE,ICORE(I000),1)
      CALL UPDAT1(ICORE(I000),ICORE(I010),NTOTTAR,2+ISPIN,90)
      RETURN
      END
