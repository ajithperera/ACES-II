      SUBROUTINE DIVO21AB(W,DDOO,DIVO,ISPIN,POP1,VRT1,POP2,VRT2,
     &                    DISSYW,NUMSYW,LISTW,IRREP,IRREPX,IUHF,
     &                    TMP) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYW,DIRPRD,POP1,VRT1,POP2,VRT2
      DIMENSION W(DISSYW,1),DDOO(1),DIVO(1),POP1(8),VRT1(8),
     &          POP2(8),VRT2(8),IPD(8),IPW1(8),IPW2(8),TMP(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C PICK UP THE INTEGRALS 
C
      CALL GETLST(W,1,NUMSYW,1,IRREP,LISTW)
C
C FOR RHF SPIN ADAPT THE INTEGRALS
C
      IF(IUHF.EQ.0) THEN
C
       CALL SPINAD3(IRREP,POP1,DISSYW,NUMSYW,W,TMP,TMP(1+NUMSYW))
C
      ENDIF
C
      IF(ISPIN.EQ.1) THEN
C
C OFFSETS FOR DERIVATIVES OF DENSITY MATRIX
C
      IPD(1)=0
      DO 1 IRREPJ=1,NIRREP-1
       IPD(IRREPJ+1)=IPD(IRREPJ)+POP2(IRREPJ)*
     &                           POP2(DIRPRD(IRREPJ,IRREPX))
1     CONTINUE
C
C LOOP OVER IRREPS OF OCCUPIED ORBITAL I
C
C OFFSETS
C
      IND=0
      IOFFWB=0
C
      DO 100 IRREP1L=1,NIRREP
C
C POPULATION WHICH CORRESPONDS TO I
C
       NOCC1=POP1(IRREP1L)
C
C DETERMINE IRREP OF M AND N
C
       IRREP2L=DIRPRD(IRREP,IRREP1L)
       IRREP2R=DIRPRD(IRREP2L,IRREPX)
C
C POPULATION WHICH CORRESPONDS TO M AND N
C
       NOCC2R=POP2(IRREP2R)
       NOCC2L=POP2(IRREP2L)
C
C IRREP OF VIRTUAL ORBITAL A
C
       IRREP1R=DIRPRD(IRREP,IRREP2R)
C
C POPULATION WHICH CORRESPONDS TO A
C
       NVRT1=VRT1(IRREP1R)
C
C DETERMINE OFFSET OF RIGHT HAND SIDE
C
       IOFFWA=0
       DO 11 IRREP1=1,IRREP1R-1
        IOFFWA=IOFFWA+VRT1(IRREP1)*POP2(DIRPRD(IRREP1,IRREP))
11     CONTINUE
C
C OFFSET IN D
C
       IOFFD1=IPD(IRREP2R)
C
       DO 10 I=1,NOCC1
       IOFFW1=IOFFWB+(I-1)*NOCC2L
       DO 10 IA=1,NVRT1
       IOFFW2=IOFFWA+(IA-1)*NOCC2R
       IND=IND+1
       DO 10 N=1,NOCC2R
       IOFFW3=IOFFW2+N
       IOFFD2=IOFFD1+(N-1)*NOCC2L
       DIVO(IND)=DIVO(IND)
     &       -SDOT(NOCC2L,DDOO(IOFFD2+1),1,W(IOFFW1+1,IOFFW3),1)
10     CONTINUE
       IOFFWB=IOFFWB+NOCC2L*NOCC1
100    CONTINUE
C
      ELSE IF(ISPIN.EQ.1.AND.IUHF.NE.0) THEN
C
C  LOOP OVER IRREPS OF VIRTUAL ORBITAL A AND OCCUPIED ORBITAL I
C
C  OFFSETS
C
      IPD(1)=0
      DO 2 IRREPJ=1,NIRREP-1
       IRREPI=DIRPRD(IRREPX,IRREPJ)
       IPD(IRREPJ+1)=IPD(IRREPJ)+POP1(IRREPJ)*VRT1(IRREPI)
2     CONTINUE
C
      IOFFD=0
      IOFFWA=0
C
      DO 200 IRREP2R=1,NIRREP
C
C  POPULATION CORRESPONDING TO N AND M
C
       IRREP2L=DIRPRD(IRREP2R,IRREPX)
       NOCC2L=POP2(IRREP2L)
       NOCC2R=POP2(IRREP2R)
C
C  DETERMINE IRREP OF M AND N
C
       IRREP1R=DIRPRD(IRREP,IRREP2R)
       IRREP1L=DIRPRD(IRREP,IRREP2L) 
C
       IOFFWB=0
       DO 4 IRREP1=1,IRREP2L-1
        IOFFWB=IOFFWB+POP2(IRREP1)*POP1(DIRPRD(IRREP1,IRREP))
4      CONTINUE
C
C  POPULATION WHICH CORRESPONDS TO A AND I
C
       NOCC1=POP1(IRREP1L)
       NVRT1=VRT1(IRREP1R)
C
C  OFFSET IN DIVO
C
       IND=IPD(IRREP1R)
C
C
       DO 20 I=1,NOCC1
       IOFFW1=IOFFWB+I
       DO 20 IA=1,NVRT1
       IOFFW2=IOFFWA+IA
       IND=IND+1
       DO 20 N=1,NOCC2R
       IOFFW3=IOFFW2+(N-1)*NVRT1
       IOFFD2=IOFFD+(N-1)*NOCC2L
       DIVO(IND)=DIVO(IND)
     &      -SDOT(NOCC2L,DDOO(IOFFD2+1),1,W(IOFFW1,IOFFW3),NOCC1)
20     CONTINUE
       IOFFWA=IOFFWA+NVRT1*NOCC2R
       IOFFD=IOFFD+NOCC2L*NOCC2R
200    CONTINUE
C
      ENDIF
      RETURN
      END
