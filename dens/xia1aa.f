         
 

      
      SUBROUTINE XIA1AA(W,DOO,XIA,ISPIN,POP,VRT,
     &                  DISSYW,NUMSYW,NOCCSQ,LISTW,IRREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYW,DIRPRD,POP,VRT
      DIMENSION W(NOCCSQ,1),DOO(1),XIA(1),POP(8),VRT(8),IPD(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C   OFFSETS FOR DENSITY MATRIX
C
      IPD(1)=0
      DO 1 IRREPJ=1,NIRREP-1
       IPD(IRREPJ+1)=IPD(IRREPJ)+POP(IRREPJ)*POP(IRREPJ)
1     CONTINUE
C
C      PICK UP THE INTEGRALS 
C
      CALL GETLST(W,1,NUMSYW,1,IRREP,LISTW)
      CALL SYMEXP2(IRREP,POP,NOCCSQ,DISSYW,NUMSYW,W,W)
C
C  LOOP OVER IRREPS OF VIRTUAL ORBITAL A AND OCCUPIED ORBITAL I
C
C  OFFSETS
C
      IND=0
      IOFFWA=0
      IOFFWB=0
C
      DO 100 IRREP1=1,NIRREP
C
C  POPULATION CORRESPONDING TO A AND I
C
       NVRT1=VRT(IRREP1)
       NOCC1=POP(IRREP1)
C
C  DETERMINE IRREP OF M AND N
C
       IRREP2=DIRPRD(IRREP,IRREP1)
C
C  POPULATION WHICH CORRESPONDS TO M AND N
C
       NOCC2=POP(IRREP2)
C
C  OFFSET IN D
C
       IOFFD1=IPD(IRREP2)
C
       DO 10 I=1,NOCC1
       IOFFW1=IOFFWB+(I-1)*NOCC2
       DO 10 IA=1,NVRT1
       IOFFW2=IOFFWA+(IA-1)*NOCC2
       IND=IND+1
       DO 10 N=1,NOCC2
       IOFFW3=IOFFW2+N
       IOFFD2=IOFFD1+(N-1)*NOCC2
       XIA(IND)=XIA(IND)
     &            +SDOT(NOCC2,DOO(IOFFD2+1),1,W(IOFFW1+1,IOFFW3),1)
10     CONTINUE
       IOFFWA=IOFFWA+NOCC2*NVRT1
       IOFFWB=IOFFWB+NOCC2*NOCC1
100    CONTINUE
C
      RETURN
      END
