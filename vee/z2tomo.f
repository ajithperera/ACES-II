      SUBROUTINE Z2TOMO(ICORE,MAXCOR,IUHF,TAU,IRREPX,LSTMO,
     &                  LSTMOINC,LSTAOINC,NOINC)
C
C TRANSFORMS T2 INCREMENT VECTOR TO AB,IJ REPRESENTATION, ADDS
C TO EXISTING INCREMENT LIST AND CALCULATES LADDER CONTRIBUTION.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SDOT,ABLAD,ZILCH,ONE
      LOGICAL NOINC,TAU
      CHARACTER*2 SPLAB(3)
      DIMENSION ICORE(MAXCOR),I0T(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/EIGPROB/ISIDE
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
C
      IONE=1
      NMO=NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      CALL GETAOINF(IUHF,IRREPX)
      IF(ISIDE.EQ.2)THEN
       I0T(1)=1
       CALL GETLST(ICORE(I0T(1)),1,1,1,3,490)
       IF(IUHF.NE.0)THEN
        I0T(2)=I0T(1)+IRPDPD(IRREPX,9)*IINTFP
        CALL GETLST(ICORE(I0T(2)),1,1,1,4,490)
       ELSE
        I0T(2)=I0T(1)
       ENDIF
       I000=I0T(2)+IRPDPD(IRREPX,10)*IINTFP
      ELSE
       I000=1
      ENDIF
      
C
C TRANSFORM TARGET VECTOR BACK TO MO REPRESENTATION
C
      DO 4 ISPIN=3,3-2*IUHF,-1
       ABLAD=ZILCH
       DO 5 IRREPIJ=1,NIRREP
        IRREPAB=DIRPRD(IRREPX,IRREPIJ)
        NUMAB=IRPDPD(IRREPAB,ISYTYP(1,43+ISPIN))
        NUMIJ=IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
        NUMAO=IRPDPDAO(IRREPAB)
        IF(ISPIN.LE.2)THEN
         NUMABX=IRPDPD(IRREPAB,18+ISPIN)
         NUMIJX=IRPDPD(IRREPIJ,20+ISPIN)
        ELSE
         NUMABX=NUMAB
         NUMIJX=NUMIJ
        ENDIF
        I010=I000+MAX(NUMABX,NUMAO)*NUMIJX*IINTFP
        I020=I010+NUMABX*NUMIJX*IINTFP
        I030=I020+MAX(NT(1),NT(2),NAO*NMO,NUMABX,NUMIJX)*IINTFP
        I040=I030+MAX(NT(1),NT(2),NAO*NMO,NUMABX,NUMIJX)*IINTFP
        I050=I040+MAX(NT(1),NT(2),NAO*NMO)*IINTFP
        CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREPIJ,LSTAOINC+ISPIN)
CSSS        call checksum("R2      :",icore(i000),NUMAO*NUMIJ,s)
        CALL T2TRAN(2,ICORE(I010),ICORE(I000),ICORE(I020),
     &              ICORE(I030),ICORE(I040),NAO,NMO,ISPIN,IUHF,
     &              IRREPIJ,IRREPAB)
        CALL SCOPY(NUMABX*NUMIJ,ICORE(I010),1,ICORE(I000),1)
CSSS        call checksum("R2*ABCD :",icore(i000),NUMABX*NUMIJ,s)
C
C FORM CONTRIBUTION TO L1 INCREMENT
C
        IF(ISPIN.EQ.3.AND.IUHF.EQ.0.AND.ISIDE.EQ.2)THEN
         CALL SPINAD1(IRREPIJ,POP(1,1),NUMAB,ICORE(I000),
     &                ICORE(I020),ICORE(I030))
         CALL GETLST(ICORE(I020),1,1,1,1,90)
         CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &              ICORE(I0T(1)),ICORE(I020),ICORE(I000),
     &              ICORE(I050),NUMAB,VRT(1,1),POP(1,1),
     &              VRT(1,1),VRT(1,1),POP(1,1),POP(1,1),
     &              'STST')
        ELSEIF(ISPIN.LE.2.AND.ISIDE.EQ.2)THEN
         CALL SYMEXP(IRREPIJ,POP(1,ISPIN),NUMABX,ICORE(I000))
         CALL GETLST(ICORE(I020),1,1,1,ISPIN,90)
         CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &              ICORE(I0T(ISPIN)),ICORE(I020),ICORE(I000),
     &              ICORE(I050),NUMABX,VRT(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),
     &              'TSTS')
        ELSEIF(ISPIN.EQ.3.AND.IUHF.NE.0.AND.ISIDE.EQ.2)THEN
         CALL GETLST(ICORE(I020),1,1,1,2,90)
         CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &              ICORE(I0T(1)),ICORE(I020),ICORE(I000),
     &              ICORE(I050),NUMAB,VRT(1,1),POP(1,1),
     &              VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &              'TSTS')
         CALL GETLST(ICORE(I020),1,1,1,1,90)
         CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &              ICORE(I0T(2)),ICORE(I020),ICORE(I000),
     &              ICORE(I050),NUMAB,VRT(1,2),POP(1,2),
     &              VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &              'STST')
        ENDIF
C
C ANTISYMMETRIZE LEFT-HAND INDICES FOR AA AND BB SPIN CASES
C
        IF(ISPIN.LE.2)THEN
         CALL SQSYM(IRREPAB,VRT(1,ISPIN),NUMAB,NUMABX,NUMIJ,
     &              ICORE(I000),ICORE(I010))
         CALL SCOPY(NUMAB*NUMIJ,ICORE(I000),1,ICORE(I010),1)
        ENDIF
        IF(NOINC) THEN
         CALL PUTLST(ICORE(I010),1,NUMIJ,1,IRREPIJ,LSTMOINC+ISPIN)
        ELSE
         CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREPIJ,LSTMOINC+ISPIN)
         CALL SAXPY (NUMIJ*NUMAB,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMIJ,1,IRREPIJ,LSTMOINC+ISPIN)
C
C CALCULATE LADDER CONTRIBUTION
C
         CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREPIJ,LSTMO+ISPIN)
         IF(TAU)THEN
          IF(ISPIN.LE.2)THEN
           CALL GETLST(ICORE(I020),1,1,1,ISPIN,90)
           CALL FTAU  (ICORE(I000),ICORE(I020),ICORE(I020),NUMAB,
     &                 NUMIJ,POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
C SG 1/24/97     &                 VRT(1,ISPIN),IRREP,ISPIN,ONE)
     &                 VRT(1,ISPIN),IRREPIJ,ISPIN,ONE)
          ELSE 
           CALL GETLST(ICORE(I020),1,1,1,1,90)
           CALL GETLST(ICORE(I030),1,1,1,1+IUHF,90)
           CALL FTAU  (ICORE(I000),ICORE(I020),ICORE(I030),NUMAB,
     &                 NUMIJ,POP(1,1),POP(1,2),VRT(1,1),
C SG 1/24/97     &                 VRT(1,2),IRREP,3,ONE)
     &                 VRT(1,2),IRREPIJ,3,ONE)
          ENDIF
         ENDIF
         ABLAD=ABLAD+SDOT(NUMIJ*NUMAB,ICORE(I000),1,ICORE(I010),1)
        ENDIF
5      CONTINUE
CSSS        WRITE(6,5005)SPLAB(ISPIN),ABLAD
CSSS5005    FORMAT(T3,'W(abef) ',A2,' contribution =',F14.10,' a.u.')
C
4     CONTINUE
C
      IF(ISIDE.EQ.2)THEN
       DO 9 ISPIN=1,1+IUHF
        CALL PUTLST(ICORE(I0T(ISPIN)),1,1,1,2+ISPIN,490)
9      CONTINUE
      ENDIF
C
      RETURN
      END
