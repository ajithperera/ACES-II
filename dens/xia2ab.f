      SUBROUTINE XIA2AB(W,DVV,XIA,ISPIN,POP1,VRT1,POP2,VRT2,
     &                  DISSYW,NUMSYW,LISTW,IRREP,IUHF,TMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYW,DIRPRD,POP1,VRT1,POP2,VRT2
      DIMENSION W(DISSYW,1),DVV(1),XIA(1),POP1(8),VRT1(8),
     &          POP2(8),VRT2(8),IPD(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C      PICK UP THE INTEGRALS 
C
C
      IF(ISPIN.EQ.2.OR.IUHF.EQ.0) THEN
C
C   OFFSETS FOR DENSITY MATRIX
C
      IPD(1)=0
      DO 1 IRREPJ=1,NIRREP-1
       IPD(IRREPJ+1)=IPD(IRREPJ)+VRT2(IRREPJ)*VRT2(IRREPJ)
1     CONTINUE
C
C  LOOP OVER IRREPS OF VIRTUAL ORBITAL A AND OCCUPIED ORBITAL I
C
C  OFFSETS
C
      IND=0
      IOFFWA=0
      IOFFSET=1
C
      DO 100 IRREP1=1,NIRREP
C
C  POPULATION CORRESPONDING TO A AND I
C
       NVRT1=VRT1(IRREP1)
       NOCC1=POP1(IRREP1)
C
C  DETERMINE IRREP OF M AND N
C
       IRREP2=DIRPRD(IRREP,IRREP1)
C
C  POPULATION WHICH CORRESPONDS TO M AND N
C
       NVRT2=VRT2(IRREP2)
C
C  OFFSET IN D
C
       IOFFD1=IPD(IRREP2)
C
       DO 10 I=1,NOCC1
C
       CALL GETLST(W,IOFFSET,NVRT2,1,IRREP,LISTW)
C
C FOR RHF, SPIN ADAPT THE INTEGRALS
C
       IF(IUHF.EQ.0) THEN
        CALL SPINAD3(IRREP,VRT1,DISSYW,NVRT2,W,TMP,TMP(1+NVRT2))
       ENDIF
C
       IOFFSET=IOFFSET+NVRT2
       IOFFW1=0
       DO 10 IA=1,NVRT1
       IOFFW2=IOFFWA+(IA-1)*NVRT2
       IND=IND+1
       DO 10 N=1,NVRT2
       IOFFW3=IOFFW1+N
       IOFFD2=IOFFD1+(N-1)*NVRT2
       XIA(IND)=XIA(IND)
     &      +SDOT(NVRT2,DVV(IOFFD2+1),1,W(IOFFW2+1,IOFFW3),1)
10     CONTINUE
       IOFFWA=IOFFWA+NVRT2*NVRT1
100    CONTINUE
C
      ELSE IF(ISPIN.EQ.1.AND.IUHF.NE.0) THEN
C
      CALL GETLST(W,1,NUMSYW,1,IRREP,LISTW)
C
C   OFFSETS FOR XIA MATRIX
C
      IPD(1)=0
      DO 2 IRREPJ=1,NIRREP-1
       IPD(IRREPJ+1)=IPD(IRREPJ)+POP1(IRREPJ)*VRT1(IRREPJ)
2     CONTINUE
C
C  LOOP OVER IRREPS OF VIRTUAL ORBITAL A AND OCCUPIED ORBITAL I
C
C  OFFSETS
C
      IOFFD=0
      IOFFWA=0
      IOFFWB=0
C
      DO 200 IRREP2=1,NIRREP
C
C  POPULATION WHICH CORRESPONDS TO M AND N
C
       NVRT2=VRT2(IRREP2)
C
C
C  DETERMINE IRREP OF M AND N
C
       IRREP1=DIRPRD(IRREP,IRREP2)
C
C  POPULATION CORRESPONDING TO A AND I
C
       NVRT1=VRT1(IRREP1)
       NOCC1=POP1(IRREP1)
C  OFFSET IN D
C
       IND=IPD(IRREP1)
C
       DO 20 I=1,NOCC1
       IOFFW2=IOFFWB+I
       DO 20 IA=1,NVRT1
       IOFFW1=IOFFWA+IA
       IND=IND+1
       DO 20 N=1,NVRT2
       IOFFW3=IOFFW2+(N-1)*NOCC1
       IOFFD2=IOFFD+(N-1)*NVRT2
       XIA(IND)=XIA(IND)
     &          +SDOT(NVRT2,DVV(IOFFD2+1),1,W(IOFFW1,IOFFW3),NVRT1)
20     CONTINUE
       IOFFD=IOFFD+NVRT2*NVRT2 
       IOFFWA=IOFFWA+NVRT2*NVRT1
       IOFFWB=IOFFWB+NVRT2*NOCC1
200    CONTINUE
C
      ENDIF
      RETURN
      END
