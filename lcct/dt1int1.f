      SUBROUTINE DT1INT1(ICORE,MAXCOR,IUHF,IRREPX,WSPIN,
     &                   LISTW0,LISTT1,ISIDE)
C
C   Z(A,I) = -T1(E,M) * <MA||IE> + T1(e,m) * <mI|eA> (AA)
C
C   Z(a,i) = -T1(e,m) * <ma||ie> + T1(E,M) * <Mi|Ea> (BB)
C
C   Z(A,I) =  T1(E,M) * (2 <Mi|Ea> - <Ma|Ie>)        (RHF)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,TWO
      CHARACTER*1 MATTYP(2)
      LOGICAL WSPIN
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE,ONEM,ZILCH,TWO /1.0,-1.0,0.0,2.0/
      DATA MATTYP/'N','T'/
C
      DO 10 ISPIN=1,1+IUHF
C
C DO SECOND CONTRACTION FIRST
C
       IF(IUHF.EQ.0)THEN
        LISTW=LISTW0-1+3
       ELSE
        IF(ISIDE.EQ.1)THEN
         LISTW=LISTW0-1+5-ISPIN
        ELSE
         LISTW=LISTW0-1+2+ISPIN
        ENDIF
       ENDIF
       DISSYW=IRPDPD(IRREPX,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPX,ISYTYP(2,LISTW))
       LENT1 =IRPDPD(IRREPX,11-ISPIN)
       LENTAR=IRPDPD(IRREPX,8+ISPIN)
       I000=1
       I010=I000+MAX(LENTAR,LENT1)*IINTFP
       I020=I010+MAX(LENTAR,LENT1)*IINTFP
       I030=I020+DISSYW*NUMDSW*IINTFP
C
C INITIALIZE INCREMENTS
C
       CALL ZERO  (ICORE(I000),LENTAR)
C
       IF(IUHF.NE.0)THEN
        CALL GETLST(ICORE(I010),1,1,1,3-ISPIN,LISTT1)
       ELSE
        CALL GETLST(ICORE(I010),1,1,1,1,LISTT1)
       ENDIF
       IF(WSPIN.OR.IUHF.NE.0)THEN
        CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPX,LISTW)
       ELSE
        LISTWA=LISTW0-1+1
        LISTWB=LISTW0-1+5
        I030=I020+DISSYW*NUMDSW*IINTFP
        CALL GETLST(ICORE(I030),1,NUMDSW,1,IRREPX,LISTWA)
        CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPX,LISTWB)
        CALL SAXPY (NUMDSW*DISSYW,-TWO,ICORE(I030),1,ICORE(I020),1)
       ENDIF
       CALL XGEMM ('N',MATTYP(ISIDE),1,LENTAR,LENT1,ONE,ICORE(I010),1,
     &             ICORE(I020),DISSYW,ONE,ICORE(I000),1)
C
C NOW DO SECOND CONTRACTION FOR UHF ONLY
C
       IF(IUHF.NE.0)THEN
        LISTW=LISTW0-1+ISPIN
        DISSYW=IRPDPD(IRREPX,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPX,ISYTYP(2,LISTW))
        LENT1 =DISSYW
        LENTAR=NUMDSW
        I000=1
        I010=I000+LENTAR*IINTFP
        I020=I010+LENT1*IINTFP
        I030=I020+DISSYW*NUMDSW*IINTFP
        CALL GETLST(ICORE(I010),1,1,1,ISPIN,LISTT1)
        CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPX,LISTW)
        CALL XGEMM ('N',MATTYP(ISIDE),1,LENTAR,LENT1,ONEM,ICORE(I010),
     &              1,ICORE(I020),DISSYW,ONE,ICORE(I000),1)
       ENDIF
C
C THIS CONTRACTION INITIALIZES THE INCREMENT LIST
C
       CALL PUTLST(ICORE(I000),1,1,1,ISPIN+2,LISTT1)
C
10    CONTINUE
      RETURN
      END
