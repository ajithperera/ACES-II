      SUBROUTINE DDENOM(ACORE,MAXCOR,IUHF,IRREPX,LISTD0)
C
C Modified by SG (6/95) to use Hartree-Fock orbital energies for
C   Partitioned-EE calculations.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL FIELD,GEOM,ROHF,QRHF,SEMI
      LOGICAL SS, SD, DS, DD
C
      DIMENSION ACORE(MAXCOR/2)
      DIMENSION IOFFOCC(8,2),IOFFVRT(8,2)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/DRVHBAR/SS, SD, DS, DD
C
      DATA ONE/1.D0/
C
C READ IN ORBITAL ENERGIES
C
      IEVALA=1
      IEND=IEVALA+NOCCO(1)+NVRTO(1)
      IF (DD) THEN
        CALL GETREC(20,'JOBARC','FDIAGA  ',IINTFP*(NOCCO(1)+NVRTO(1)),
     &              ACORE(IEVALA)) 
      ELSE
        CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCCO(1)+NVRTO(1)),
     &     ACORE(IEVALA)) 
      ENDIF
      IF(IUHF.EQ.0) THEN
       IEVALB=IEVALA
      ELSE
       IEVALB=IEND
       IEND=IEVALB+NVRTO(2)+NOCCO(2)
       IF (DD) THEN
         CALL GETREC(20,'JOBARC','FDIAGB  ',IINTFP*(NOCCO(2)+NVRTO(2)),
     &               ACORE(IEVALB)) 
       ELSE
         CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCCO(2)+NVRTO(2)),
     &      ACORE(IEVALB)) 
       ENDIF
      ENDIF
      MAXLEN=0
      DO 2 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       MAXLEN=MAX(MAXLEN,IRPDPD(IRREPR,14)*IRPDPD(IRREPL,13))
2     CONTINUE
      IF((MAXLEN+IEND).GE.MAXCOR) CALL INSMEM('DDENOM',MAXLEN+
     &         IEND,MAXCOR)
C
C DETERMINE OFFSETS
C
      IOFFOCC(1,1)=-1
      IOFFOCC(1,2)=-1
      IOFFVRT(1,1)=NOCCO(1)-1
      IOFFVRT(1,2)=NOCCO(2)-1
      DO 1 IRREP=1,NIRREP-1
       IOFFOCC(IRREP+1,1)=IOFFOCC(IRREP,1)+POP(IRREP,1)
       IOFFOCC(IRREP+1,2)=IOFFOCC(IRREP,2)+POP(IRREP,2)
       IOFFVRT(IRREP+1,1)=IOFFVRT(IRREP,1)+VRT(IRREP,1)
       IOFFVRT(IRREP+1,2)=IOFFVRT(IRREP,2)+VRT(IRREP,2)
1     CONTINUE
C
C LOOP OVER IRREPS
C
      DO 1000 IRREPR=1,NIRREP
C
       IRREPL=DIRPRD(IRREPX,IRREPR)
C
C LOOP OVER IRREP RIGHT HAND SIDE
C
       IND=0
C
       DO 900 IRREPRR=1,NIRREP
C
        IRREPRL=DIRPRD(IRREPR,IRREPRR)
        NOCCR=POP(IRREPRR,2)
        NOCCL=POP(IRREPRL,1)
        IOFFJ=IOFFOCC(IRREPRR,2)
        IOFFI=IOFFOCC(IRREPRL,1)
C
        DO 100 IJ=1,NOCCR
        EVALJ=ACORE(IEVALB+IOFFJ+IJ)
        DO 100 II=1,NOCCL
        EVALIJ=EVALJ+ACORE(IEVALA+IOFFI+II)
C
C LOOP OVER IRREP LEFT HAND SIDE
C
        DO 800 IRREPLR=1,NIRREP
C
         IRREPLL=DIRPRD(IRREPLR,IRREPL)
C
         NVRTR=VRT(IRREPLR,2)
         NVRTL=VRT(IRREPLL,1)
         IOFFB=IOFFVRT(IRREPLR,2)
         IOFFA=IOFFVRT(IRREPLL,1)
C
         DO 200 IB=1,NVRTR
         EVALIJB=EVALIJ-ACORE(IEVALB+IOFFB+IB)
C
          DO 200 IA=1,NVRTL
          IND=IND+1
          ACORE(IEND+IND)=(EVALIJB-ACORE(IEVALA+IOFFA+IA))

200      CONTINUE
800     CONTINUE
100     CONTINUE
900    CONTINUE
       CALL PUTLST(ACORE(IEND+1),1,IRPDPD(IRREPR,14),1,IRREPR,LISTD0+3) 
1000  CONTINUE
C
      IF(IUHF.NE.0)THEN
C
C FOR UHF AND ROHF, DO ALSO AAAA AND BBBB CASES
C
       DO 2000 ISPIN=1,2
C
       IF(ISPIN.EQ.1) THEN
        IEVAL=IEVALA
       ELSE IF (ISPIN.EQ.2) THEN
        IEVAL=IEVALB
       ENDIF
C
C LOOP OVER IRREPS
C
        DO 3000 IRREPR=1,NIRREP
C
         IRREPL=DIRPRD(IRREPX,IRREPR)
C
         IND=0
C
C LOOP OVER IRREP RIGHT HAND SIDE
C
         DO 2900 IRREPRR=1,NIRREP
          IRREPRL=DIRPRD(IRREPR,IRREPRR)
          IF(IRREPRR.LT.IRREPRL) GO TO 2900
          NOCCR=POP(IRREPRR,ISPIN)
          NOCCL=POP(IRREPRL,ISPIN)
          IOFFJ=IOFFOCC(IRREPRR,ISPIN)
          IOFFI=IOFFOCC(IRREPRL,ISPIN)
          DO 2100 IJ=1,NOCCR
           IEND1=NOCCL
           IF(IRREPRR.EQ.IRREPRL) IEND1=IJ-1
           EVALJ=ACORE(IEVAL+IOFFJ+IJ)
           DO 2100 II=1,IEND1
           EVALIJ=EVALJ+ACORE(IEVAL+IOFFI+II)
C
C LOOP OVER IRREP LEFT HAND SIDE
C
           DO 2800 IRREPLR=1,NIRREP
            IRREPLL=DIRPRD(IRREPLR,IRREPL)
            IF(IRREPLR.LT.IRREPLL) GO TO 2800
            NVRTL=VRT(IRREPLL,ISPIN)
            NVRTR=VRT(IRREPLR,ISPIN)
            IOFFA=IOFFVRT(IRREPLL,ISPIN)
            IOFFB=IOFFVRT(IRREPLR,ISPIN)
            DO 2200 IB=1,NVRTR
             IEND2=NVRTL
             IF(IRREPLL.EQ.IRREPLR) IEND2=IB-1
             EVALIJB=EVALIJ-ACORE(IEVAL+IOFFB+IB)
             DO 2200 IA=1,IEND2
             IND=IND+1
             ACORE(IEND+IND)=(EVALIJB-ACORE(IEVAL+IOFFA+IA))
2200        CONTINUE
2800       CONTINUE
2100      CONTINUE
2900     CONTINUE
         CALL PUTLST(ACORE(IEND+1),1,IRPDPD(IRREPR,2+ISPIN),1,
     &               IRREPR,LISTD0+ISPIN)
3000    CONTINUE
2000   CONTINUE
      ENDIF
C
      DO 4000 ISPIN=1,1+IUHF
C
       IF(ISPIN.EQ.1) THEN
        IEVAL=IEVALA
       ELSE IF(ISPIN.EQ.2) THEN
        IEVAL=IEVALB
       ENDIF
C
C LOOP OVER IRREPS
C
       IND=0
       DO 5000 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPX,IRREPR)
        NOCCR=POP(IRREPR,ISPIN)
        NVRTL=VRT(IRREPL,ISPIN)
        IOFFI=IOFFOCC(IRREPR,ISPIN)
        IOFFA=IOFFVRT(IRREPL,ISPIN)
        DO 5100 I=1,NOCCR
         EVALI=ACORE(IEVAL+IOFFI+I)
         DO 5200 IA=1,NVRTL
          EVALIA=EVALI-ACORE(IEVAL+IOFFA+IA)
          IND=IND+1
          ACORE(IEND+IND)=EVALIA
5200     CONTINUE
5100    CONTINUE
5000   CONTINUE
C
       CALL PUTLST(ACORE(IEND+1),1,1,1,9,LISTD0+ISPIN)
C
4000  CONTINUE
C
      RETURN
      END
