      SUBROUTINE ZIAAO(ICORE,MAXCOR,IUHF,IRREPX,LSTZIA)
C
C THIS ROUTINE CALCULATES THE TERM
C
C            -1/2 * SUM(mef) [t(im,ef) * <ma||ef>] 
C
C WHICH IS A CONTRIBUTION TO
C
C              t(i,a)*D(i,a)
C
CEND
C
C C. HUBER, KARLSRUHE AUG/93

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD,LSTZIA
      CHARACTER*2 SPLAB(3)
      CHARACTER*8 LABELX(2),LABELY(2)
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/ZIA/ IOFFZIA(8,2)
      CHARACTER*7 BLA(4)
C
      DATA BLA /'ispin=1','ispin=2','ispin=3','ispin=4'/
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
      DATA LABELX /'SCFEVECA','SCFEVECB'/
      DATA LABELY /'SCFEVECA','SCFEVECB'/
C
      NMO=NOCCO(1)+NVRTO(1)
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
C
C FILL COMMON /AOSYM/
C
      IF((IFLAGS(2).EQ.3.OR.IFLAGS(2).EQ.3).AND.IFLAGS(11).EQ.2) 
     &    CALL ERREX
C
      CALL GETAOINF(IUHF,IRREPX)
C
C CALCULATE MEMORY OF THE Z(i,a) ARRAY (Z(I,A) AND Z(i,a))
C
      I000=1
      I010=I000+NT(1)*IINTFP
      I020=I010+NT(2)*IUHF*IINTFP
      IF((IFLAGS(2).NE.3.AND.IFLAGS(2).NE.4).OR.IFLAGS(11).EQ.2) THEN
        CALL GETLST(ICORE(I000),1,1,1,3,LSTZIA)
        IF(IUHF.NE.0) CALL GETLST(ICORE(I010),1,1,1,4,LSTZIA)
      ELSE
        CALL UPDMOI(1,NT(1),1,90,0,0)
        CALL ZERO(ICORE(I000),NT(1))
        IF(IUHF.NE.0) THEN
          CALL UPDMOI(1,NT(2),2,90,0,0)
         CALL ZERO(ICORE(I010),NT(2))
        ENDIF
      ENDIF
C
C CALCULATE OFFSETS FOR THE Z(i;a) ARRAY
C
      IOFFZIA(1,1)=1
      IOFFZIA(1,2)=1
      DO 10 ISPIN=1,1+IUHF
        DO 20 IRREP=2,NIRREP
          IOFFZIA(IRREP,ISPIN)=IOFFZIA(IRREP-1,ISPIN)
     &                        +POP(IRREP-1,ISPIN)*VRT(IRREP-1,ISPIN)
 20     CONTINUE
 10   CONTINUE
C
C LOOP OVER SPINCASES AND IRREPS
C
C SPIN CASES AA AND BB (UHF ONLY)
C
      IF(IUHF.EQ.1) THEN
        DO 1100 ISPIN=1,2
          DO 1200 IRREPIJ=1,NIRREP
            IRREPAB=DIRPRD(IRREPIJ,IRREPX)
            NUMIJS=IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
            NUMIJ=IRPDPD(IRREPIJ,20+ISPIN)
            NUMAO=IRPDPDAO(IRREPAB)
C
C DIVIDE UP CORE
C
            I030=I020+NUMAO*NUMIJ*IINTFP
            I040=I030+MAX(NAO*NMO,NT(1))*IINTFP
            I050=I040+MAX(NAO*NMO,NT(2))*IINTFP
            IFIN=I050+NAO*NMO*IINTFP
C
C I020: Z(i,j;nu,mu)
C I030: MO coefficients
C I040: MO coefficients
C I050: scratch array
C
            IF(IFIN.GE.MAXCOR) CALL INSMEM('ZIAAO',IFIN,MAXCOR)
C
C GET Z(i,j ; mu,nu) AND CALL TRANSFORMATION ROUTINE
C
            CALL GETLST(ICORE(I020),1,NUMIJS,1,IRREPIJ,260+ISPIN)
            CALL SYMEXP(IRREPIJ,POP(1,ISPIN),NUMAO,ICORE(I020))
            IF (ISPIN.EQ.1) THEN
              CALL ZIATRAN(ICORE(I000),ICORE(I020),ICORE(I030),
     &                     ICORE(I040),ICORE(I050),NAO,NMO,ISPIN,
     &                     IUHF,IRREPIJ)
            ELSE
              CALL ZIATRAN(ICORE(I010),ICORE(I020),ICORE(I030),
     &                     ICORE(I040),ICORE(I050),NAO,NMO,ISPIN,
     &                     IUHF,IRREPIJ)
            ENDIF
C
1200      CONTINUE
c          if (ispin.eq.1) then
c            call printt1(bla(ispin),icore(i000),nt(1))
c          else
c            call printt1(bla(ispin),icore(i010),nt(2))
c          endif
1100    CONTINUE
      ENDIF
C
C SPIN CASES AB AND BA
C
      DO 100 ISPIN=3,3+IUHF
        DO 200 IRREPIJ=1,NIRREP
          IRREPAB=DIRPRD(IRREPIJ,IRREPX)
          NUMIJ=IRPDPD(IRREPIJ,ISYTYP(2,46))
          NUMAO=IRPDPDAO(IRREPAB)
C
C DIVIDE UP CORE
C
          I030=I020+NUMAO*NUMIJ*IINTFP
          IF ((IUHF.EQ.0).OR.(ISPIN.EQ.4)) THEN
            I040=I030+MAX(NUMAO,NUMIJ)*IINTFP
            I050=I040+MAX(NUMAO,NUMIJ)*IINTFP
            IFIN=I050+MAX(NUMAO,NUMIJ)
          ELSE
            I040=I030+MAX(NAO*NMO,NT(1))*IINTFP
            I050=I040+MAX(NAO*NMO,NT(2))*IINTFP
            IFIN=I050+NAO*NMO*IINTFP
          ENDIF            
C
C I020: Z(i,j;nu,mu)
C I030: MO coefficients
C I040: MO coefficients
C I050: scratch array
C
          IF(IFIN.GE.MAXCOR) CALL INSMEM('ZIAAO',IFIN,MAXCOR)
C
C GET Z(i,j ; mu,nu) AND CALL TRANSFORMATION ROUTINE
C
          IF (ISPIN.EQ.3) THEN
            CALL GETLST(ICORE(I020),1,NUMIJ,1,IRREPIJ,260+ISPIN)
            IF (IUHF.EQ.0) THEN
              CALL SPINAD3(IRREPIJ,IAOPOP,NUMAO,NUMIJ,ICORE(I020),
     &                     ICORE(I030),ICORE(I040))
              I040=I030+MAX(NAO*NMO,NT(1))*IINTFP
              I050=I040+MAX(NAO*NMO,NT(2))*IINTFP
              IFIN=I050+NAO*NMO*IINTFP
              IF(IFIN.GE.MAXCOR) CALL INSMEM('ZIAAO',IFIN,MAXCOR)
            ENDIF
            CALL ZIATRAN(ICORE(I000),ICORE(I020),ICORE(I030),
     &                   ICORE(I040),ICORE(I050),NAO,NMO,ISPIN,IUHF,
     &                   IRREPIJ)
          ELSE
            CALL GETLST(ICORE(I020),1,NUMIJ,1,IRREPIJ,260+3)
            CALL SYMTR3(IRREPIJ,IAOPOP,IAOPOP,NUMAO,NUMIJ,ICORE(I020),
     &                  ICORE(I030),ICORE(I040),ICORE(I050))
            CALL SYMTR1(IRREPIJ,POP(1,1),POP(1,2),NUMAO,ICORE(I020),
     &                  ICORE(I030),ICORE(I040),ICORE(I050))
            I040=I030+MAX(NAO*NMO,NT(1))*IINTFP
            I050=I040+MAX(NAO*NMO,NT(2))*IINTFP
            IFIN=I050+NAO*NMO*IINTFP
            IF(IFIN.GE.MAXCOR) CALL INSMEM('ZIAAO',IFIN,MAXCOR)
            CALL ZIATRAN(ICORE(I010),ICORE(I020),ICORE(I030),
     &                   ICORE(I040),ICORE(I050),NAO,NMO,ISPIN,
     &                   IUHF,IRREPIJ)
          ENDIF
C
 200    CONTINUE
c        if (ispin.eq.1.or.ispin.eq.3) then
c          call printt1(bla(ispin),icore(i000),nt(1))
c        else
c          call printt1(bla(ispin),icore(i010),nt(2))
c        endif
 100  CONTINUE
C
C STORE THE VALUES OF Z(i;a) AND Z(i^;a^) TO FILE
C
      INCREM=2
      IF(IFLAGS(2).EQ.3.OR.IFLAGS(2).EQ.4) INCREM=0
      CALL PUTLST(ICORE(I000),1,1,1,1+INCREM,LSTZIA)
      IF(IUHF.NE.0) CALL PUTLST(ICORE(I010),1,1,1,2+INCREM,LSTZIA)
      
      RETURN
      END 