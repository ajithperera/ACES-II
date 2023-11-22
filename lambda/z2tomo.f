C
      SUBROUTINE Z2TOMO(ICORE,MAXCOR,IUHF,TAU,IRREPX,
     &                  LSTMO,LSTMOINC,LSTAOINC)
C
C TRANSFORMS T2 INCREMENT VECTOR TO AB,IJ REPRESENTATION, ADDS
C TO EXISTING INCREMENT LIST AND CALCULATES LADDER CONTRIBUTION.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION X,SNRM2,SDOT,ABLAD,ZILCH,ONE
      LOGICAL TAU
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      CHARACTER*2 SPLAB(3)
      DIMENSION ICORE(MAXCOR),I0T(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
C
      NNM1O2(I)=(I*(I-1))/2
C
      IONE=1
      LUINT=10
      NMO=NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      CALL GETAOINF(IUHF,IRREPX)
C
C TRANSFORM TARGET VECTOR BACK TO MO REPRESENTATION
C
      IF (CCSD) THEN
      I0T(1)=1
      I0T(2)=I0T(1)+IINTFP*NT(1)
      I000=I0T(2)+IINTFP*NT(2)
      DO 1 ISPIN=1,1+IUHF
       CALL GETLST(ICORE(I0T(ISPIN)),1,1,1,2+ISPIN,90)
1     CONTINUE
      ENDIF 
      DO 4 ISPIN=3,3-2*IUHF,-1
       ABLAD=ZILCH
       DO 5 IRREP=1,NIRREP
        NUMAB=IRPDPD(IRREP,ISYTYP(1,43+ISPIN))
        NUMIJ=IRPDPD(IRREP,ISYTYP(2,43+ISPIN))
        NUMAO=IRPDPDAO(IRREP)
        IF(ISPIN.LE.2)THEN
         NUMABX=IRPDPD(IRREP,18+ISPIN)
         NUMIJX=IRPDPD(IRREP,20+ISPIN)
        ELSE
         NUMABX=NUMAB
         NUMIJX=NUMIJ
        ENDIF
        I010=I000+MAX(NUMABX,NUMAO)*NUMIJX*IINTFP
        I020=I010+NUMABX*NUMIJX*IINTFP
        I030=I020+MAX(NT(1),NT(2),NAO*NMO,NUMABX,NUMIJX)*IINTFP
        I040=I030+MAX(NT(1),NT(2),NAO*NMO,NUMABX,NUMIJX)*IINTFP
        I050=I040+MAX(NT(1),NT(2),NAO*NMO)*IINTFP
        CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREP,LSTAOINC+ISPIN)
        CALL T2TRAN(2,ICORE(I010),ICORE(I000),ICORE(I020),
     &              ICORE(I030),ICORE(I040),NAO,NMO,ISPIN,IUHF,
     &              IRREP,IRREP)
        CALL SCOPY(NUMABX*NUMIJ,ICORE(I010),1,ICORE(I000),1)
C
C FORM CONTRIBUTION TO L1 INCREMENT
C
        IF(CCSD) THEN
         IF(ISPIN.EQ.3.AND.IUHF.EQ.0)THEN
          CALL SPINAD1(IRREP,POP(1,1),NUMAB,ICORE(I000),
     &                 ICORE(I020),ICORE(I030))
          CALL GETLST(ICORE(I020),1,1,1,1,90)   
          CALL DOT24(IRREP,ICORE(I0T(1)),ICORE(I020),ICORE(I000),
     &               ICORE(I050),NUMAB,VRT(1,1),POP(1,1),
     &               VRT(1,1),VRT(1,1),POP(1,1),POP(1,1),
     &               'STST')
         ELSEIF(ISPIN.LE.2)THEN
          CALL SYMEXP(IRREP,POP(1,ISPIN),NUMABX,ICORE(I000))
          CALL GETLST(ICORE(I020),1,1,1,ISPIN,90)   
          CALL DOT24(IRREP,ICORE(I0T(ISPIN)),ICORE(I020),ICORE(I000),
     &               ICORE(I050),NUMABX,VRT(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),
     &               'TSTS')
         ELSEIF(ISPIN.EQ.3.AND.IUHF.NE.0)THEN
          CALL GETLST(ICORE(I020),1,1,1,2,90)   
          CALL DOT24(IRREP,ICORE(I0T(1)),ICORE(I020),ICORE(I000),
     &               ICORE(I050),NUMAB,VRT(1,1),POP(1,1),
     &               VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &               'TSTS')
          CALL GETLST(ICORE(I020),1,1,1,1,90)   
          CALL DOT24(IRREP,ICORE(I0T(2)),ICORE(I020),ICORE(I000),
     &               ICORE(I050),NUMAB,VRT(1,2),POP(1,2),
     &               VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &               'STST')
         ENDIF
        ENDIF
C
C ANTISYMMETRIZE LEFT-HAND INDICES FOR AA AND BB SPIN CASES
C
        IF(ISPIN.LE.2)THEN
         CALL SQSYM(IRREP,VRT(1,ISPIN),NUMAB,NUMABX,NUMIJ,
     &              ICORE(I000),ICORE(I010))
         CALL SCOPY(NUMAB*NUMIJ,ICORE(I000),1,ICORE(I010),1)
        ENDIF
C      
        CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREP,LSTMOINC+ISPIN)
        CALL SAXPY (NUMIJ*NUMAB,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMIJ,1,IRREP,LSTMOINC+ISPIN)
C
C CALCULATE LADDER CONTRIBUTION
C
        CALL GETLST(ICORE(I000),1,NUMIJ,1,IRREP,LSTMO+ISPIN)
        IF(TAU)THEN
         IF(ISPIN.LE.2)THEN
          CALL GETLST(ICORE(I020),1,1,1,ISPIN,90)
          CALL FTAU  (ICORE(I000),ICORE(I020),ICORE(I020),NUMAB,
     &                NUMIJ,POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &                VRT(1,ISPIN),IRREP,ISPIN,ONE)
         ELSE 
          Write(6,*) "I am here", iuhf
          CALL GETLST(ICORE(I020),1,1,1,1,90)
          CALL GETLST(ICORE(I030),1,1,1,1+IUHF,90)
          CALL FTAU  (ICORE(I000),ICORE(I020),ICORE(I030),NUMAB,
     &                NUMIJ,POP(1,1),POP(1,2),VRT(1,1),
     &                VRT(1,2),IRREP,3,ONE)
         ENDIF
        ENDIF
        ABLAD=ABLAD+SDOT(NUMIJ*NUMAB,ICORE(I000),1,ICORE(I010),1)
5      CONTINUE
       WRITE(6,5005)SPLAB(ISPIN),ABLAD
5005   FORMAT(T3,'W(abef) ',A2,' contribution =',F14.10,' a.u.')
C
4     CONTINUE
c
      IF(CCSD) THEN
       DO 9 ISPIN=1,1+IUHF
        CALL PUTLST(ICORE(I0T(ISPIN)),1,1,1,2+ISPIN,90)
9      CONTINUE
      ENDIF
      RETURN
      END