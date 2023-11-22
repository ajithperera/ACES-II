
      SUBROUTINE SPINMOD(IRREP,NUM1,DIS,NUM,A,SCR,ISCR, FD, FX)
C
C THIS ROUTINE SPIN ADAPTS A MATRIX A BY
C
C  A(J,I,K,L) ----> FD*A(J,I,K,L) + FX A(I,J,K,L) 
C
C NOTE THAT THE SPIN ADAPTION IS CARRIED OUT HERE IN PLACE AND
C ONLY TWO ADDITIONAL ARRAYS OF SIZE NUM ARE REQUIRED
C
C INPUT : IRREP .......... IRREP OF THE GIVEN PAR OF A
C         NUM1 ........... POPULATION VECTOR OF I AND J
C         DIS ............ DISTRIBUTION SIZE OF A
C         NUM ............ NUMBER OF DISTRIBUTIONS IN A
C         A .............. INPUT MATRIX A
C         SCR,ISCR ....... TWO SCRATCH ARRAYS OF DIMENSION NUM
C
C  OUTPUT : A ............ SPIN ADAPTED MATRIX A
C 
CEND
C
C CODED SEPTEMBER/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIS,DIRPRD
      DIMENSION A(DIS,NUM),SCR(1),ISCR(2),NUM1(8),IP(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
C
      DATA TWO /2.0D+00/
C
      IP(1)=0
      DO 10 IRREPJ=1,NIRREP-1
       IRREPI=DIRPRD(IRREP,IRREPJ)
       IP(IRREPJ+1)=IP(IRREPJ)+NUM1(IRREPJ)*NUM1(IRREPI)
10    CONTINUE
       
C
C  GET FIRST THE NEW ADDRESSES AND STORE THEM IN ISCR
C
      NTOTAL=0
      DO 20  IRREPJ=1,NIRREP
       NUMJ=NUM1(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NUMI=NUM1(IRREPI)
       NTOTAL=NTOTAL+NUMI*NUMJ
       DO 15 J=1,NUMJ
       DO 15 I=1,NUMI
        IND1=IP(IRREPJ)+(J-1)*NUMI+I
        IND2=IP(IRREPI)+(I-1)*NUMJ+J
        ISCR(IND1)=IND2
15     CONTINUE
20    CONTINUE
C
C  NOW SPIN ADAPT A
C
      DO 100 IJ=1,NTOTAL
C
C  IF THIS ELEMENT HAS BEEN ALREADY TREATED SKIP
C
       IF(ISCR(IJ).NE.0) THEN
C
C  GET ADDRESS OF TRANPOSED ELEMENT AND SET IT TO ZERO
C
        IJTR=ISCR(IJ)  
        ISCR(IJTR)=0
        IF(IJ.NE.IJTR) THEN
C
C  COPY ELEMENTS OF A(IJ,..) TO SCR
C
         DO 4 L=1,NUM
          SCR(L)=A(IJ,L)
4        CONTINUE
C
C  SPIN ADAPT A(IJ,L)
C
CDIR$ IVDEP
*VOCL LOOP,NOVREC
         DO 5 L=1,NUM
          A(IJ,L)=FD*A(IJ,L)+FX*A(IJTR,L)
5        CONTINUE
C
C  SPIN ADAPT A(IJTR,L)
C
         DO 6 L=1,NUM
          A(IJTR,L)=FD*A(IJTR,L)+FD*SCR(L)
6        CONTINUE
C
        ENDIF
       ENDIF
C ALL DONE FOR IJ AND IJTR
C
100   CONTINUE 
C
      RETURN
      END