      SUBROUTINE ROHFU(IRREP,NPERT,NA,NB,UAIA,UAIB,EVALA,EVALB,POPA,
     &                 POPB,VRTA,VRTB,NOCCA,NOCCB)
C
C  THIS ROUTINE FORMS THE ELEMENTS OF THE OCCUPIED-VIRTUAL BLOCK
C  OF THE DENSITY MATRIX
C 
CEND
C
C CODED MARCH/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POPA,POPB,VRTA,VRTB,DIRPRD
      DIMENSION UAIA(NA,NPERT),UAIB(NB,NPERT),EVALA(1),EVALB(1),POPA(8),
     &          POPB(8),VRTA(8),VRTB(8),IV(8,2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      IV(1,1)=NOCCA
      IV(1,2)=NOCCB
      DO 1 IRREPJ=1,NIRREP-1
       IV(IRREPJ+1,1)=IV(IRREPJ,1)+VRTA(IRREPJ)
       IV(IRREPJ+1,2)=IV(IRREPJ,2)+VRTB(IRREPJ)
1     CONTINUE
C
      INDA0=0
      INDB0=0
      IOFFPA=0
      IOFFPB=0
      DO 10 IRREPJ=1,NIRREP
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NOA=POPA(IRREPJ)
       NVA=VRTA(IRREPI)
       NOB=POPB(IRREPJ)
       NVB=VRTB(IRREPI)
       IOFFVA=IV(IRREPI,1)
       IOFFVB=IV(IRREPI,2)
       IF(NOA.EQ.NOB.AND.NVA.EQ.NVB) THEN
C
C NO OPEN-SHELL ORBITALS IN THIS IRREP
C
       DO 15 I=1,NOA
        INDIA=I+IOFFPA
        INDIB=I+IOFFPB
        DO 5 IPERT=1,NPERT
        INDA=INDA0
        INDB=INDB0
        DO 5 IA=1,NVA
        INDAA=IA+IOFFVA
        INDAB=IA+IOFFVB
        INDA=INDA+1
        INDB=INDB+1
        XIA=(UAIA(INDA,IPERT)+UAIB(INDB,IPERT))/
     &                           (EVALA(INDIA)+EVALB(INDIB)
     &                          -EVALA(INDAA)-EVALB(INDAB))
        UAIA(INDA,IPERT)=XIA
        UAIB(INDB,IPERT)=XIA
5      CONTINUE
       INDA0=INDA0+NVA
       INDB0=INDB0+NVB
15     CONTINUE
       ELSE
C
C THERE ARE OPEN-SHELL ORBITALS :
C
C  NOPOO : NUMBER OF OPEN-SHELL ORBITALS WITHIN THE OCCUPIED ORBITALS.
C  NOPVV : NUMBER OF OPEN-SHELL ORBITALS WITHIN THE VIRTUAL ORBITALS.
C
       NOPOO=NOA-NOB
       NOPVV=NVB-NVA
       DO 6 I=1,NOB
        INDIA=I+IOFFPA
        INDIB=I+IOFFPB
        DO 7 IPERT=1,NPERT
        INDB=INDB0
        DO 7 IOPEN=1,NOPVV
         INDAB=IOPEN+IOFFVB
         INDB=INDB+1
         UAIB(INDB,IPERT)=UAIB(INDB,IPERT)/
     &                        (EVALB(INDIB)-EVALB(INDAB))
7       CONTINUE
        INDB0=INDB0+NOPVV
        DO 8 IPERT=1,NPERT
        INDA=INDA0
        INDB=INDB0
        DO 8 IA=1,NVA
         INDAA=IA+IOFFVA
         INDAB=IA+NOPVV+IOFFVB
         INDA=INDA+1
         INDB=INDB+1
         XIA=(UAIA(INDA,IPERT)+UAIB(INDB,IPERT))
     &                         /(EVALA(INDIA)+EVALB(INDIB)
     &                           -EVALA(INDAA)-EVALB(INDAB))
         UAIA(INDA,IPERT)=XIA
         UAIB(INDB,IPERT)=XIA
8       CONTINUE
        INDA0=INDA0+NVA
        INDB0=INDB0+NVA
6      CONTINUE
       DO 9 IOPEN=1,NOPOO
        INDIA=IOPEN+NOB+IOFFPA
        DO 19 IPERT=1,NPERT
        INDA=INDA0
        DO 19 IA=1,NVA
         INDAA=IA+IOFFVA
         INDA=INDA+1
         UAIA(INDA,IPERT)=UAIA(INDA,IPERT)
     &                   /(EVALA(INDIA)-EVALA(INDAA))
19      CONTINUE
        INDA0=INDA0+NVA
9      CONTINUE
       ENDIF
       IOFFPA=IOFFPA+NOA
       IOFFPB=IOFFPB+NOB
10    CONTINUE
      RETURN
      END
