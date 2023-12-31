      SUBROUTINE SYPKEV(EVEC,SCR,NBASIS,TOTLEN,ISTART,ISPIN)
C
C THIS ROUTINE ACCEPTS THE FULL SQUARE EIGENVECTOR LIST AND RETURNS
C  THE SAME LIST IN THE FOLLOWING SYMMETRY-ORDERED FORM.
C
C    AO-OCCUPIED MO (ORDERED BY IRREP)
C    AO-VIRTUAL  MO (ORDERED BY IRREP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,MOPOP,AOPOP,TOTLEN
      DIMENSION EVEC(1),SCR(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON /AOOFST/ INDOCC(8,2),INDVRT(8,2)
C
C FIRST FORM THE AO-OCCUPIED BLOCK OF THE EIGENVECTOR MATRIX
C
      IOFF1=1
      IOFF2=1
      DO 10 IRREP=1,NIRREP
       INDOCC(IRREP,ISPIN)=IOFF2+ISTART
       DO 20 IOCC=1,POP(IRREP,ISPIN)
        LENGTH=AOPOP(IRREP)
        CALL SCOPY(LENGTH,EVEC(IOFF1),1,SCR(IOFF2),1)
        IOFF2=IOFF2+LENGTH
        IOFF1=IOFF1+NBASIS
20     CONTINUE
       IOFF1=IOFF1+AOPOP(IRREP)
10    CONTINUE
C
C NOW FORM THE AO-VIRTUAL BLOCK OF THE EIGENVECTOR MATRIX
C
      IOFF1=NBASIS*NOCCO(ISPIN)+1
      DO 110 IRREP=1,NIRREP
       INDVRT(IRREP,ISPIN)=IOFF2+ISTART
       DO 120 IVRT=1,VRT(IRREP,ISPIN)
        LENGTH=AOPOP(IRREP)
        CALL SCOPY(LENGTH,EVEC(IOFF1),1,SCR(IOFF2),1)
        IOFF2=IOFF2+LENGTH
        IOFF1=IOFF1+NBASIS
120    CONTINUE
       IOFF1=IOFF1+AOPOP(IRREP)
110   CONTINUE
C
C NOW OVERWRITE INPUT ARRAY WITH EIGENVECTOR ARRAY AND RETURN
C
      TOTLEN=IOFF2-1
      CALL SCOPY(TOTLEN,SCR,1,EVEC,1)
      RETURN
      END
