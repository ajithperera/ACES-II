
      
      SUBROUTINE SYMRHF(IRREP,NUML,NUMR,DISSIZE,A,SCR,ISCRL,ISCRR)
C
C  THIS ROUTINE FORMS THE TERM 
C
C   A(I,J,K,L) --> A(I,J,K,L) + A(J,I,L,K)  
C
C  WHICH IS REQUIRED IN RHF CCSD CALCULATIONS
C
C  INPUT :  IRREP  ...   IRREP OF   G(K)*G(L) = G(I)*G(J)
C           NUMR   ....  POPULATION VECTOR FOR K,L
C           NUML   ....  POPULATION VECTOR FOR I,J
C           DISSIZE ...  DISTRIBUTION SIZE OF A
C           A      ....  HOLDS THE MATRIX A         
C           SCR,ISCRL,ISCRR ... THREE SCRATCH ARRAYS OF SIZE DISSIZE
C
C  OUTPUT :  A     ....  THE RHF SYMMETRIZED MATRIX
C
CEND
C
C  WRITTEN IN SEPTEMBER/90   JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZE,DIRPRD
      DIMENSION A(DISSIZE,1),SCR(1),ISCRL(1),ISCRR(2),NUMR(8),NUML(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      DIMENSION IPL(8),IPR(8)
C
      IPR(1)=0
      IPL(1)=0          
      DO 10 IRREPJ=1,NIRREP-1
       IRREPI=DIRPRD(IRREP,IRREPJ)
       IPR(IRREPJ+1)=IPR(IRREPJ)+NUMR(IRREPJ)*NUMR(IRREPI)
       IPL(IRREPJ+1)=IPL(IRREPJ)+NUML(IRREPJ)*NUML(IRREPI)
10    CONTINUE
       
C
C  GET FIRST THE NEW ADDRESSES AND STORE THEM IN ISCR
C
C ADDRESSES FOR THE RIGHT SIDE
C
      NTOTAL=0
      DO 20  IRREPJ=1,NIRREP
       NUMJ=NUMR(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NUMI=NUMR(IRREPI)
       NTOTAL=NTOTAL+NUMI*NUMJ
       DO 15 J=1,NUMJ
       DO 15 I=1,NUMI
        IND1=IPR(IRREPJ)+(J-1)*NUMI+I
        IND2=IPR(IRREPI)+(I-1)*NUMJ+J
        ISCRR(IND1)=IND2
15     CONTINUE
20    CONTINUE
C
C  ADDRESSES FOR THE LEFT SIDE
C
      DO 120  IRREPJ=1,NIRREP
       NUMJ=NUML(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NUMI=NUML(IRREPI)
        DO 115 J=1,NUMJ
        DO 115 I=1,NUMI
         IND1=IPL(IRREPJ)+(J-1)*NUMI+I
         IND2=IPL(IRREPI)+(I-1)*NUMJ+J
         ISCRL(IND1)=IND2
115     CONTINUE
120    CONTINUE
C
C  NOW FORM A(IJ,KL) + A(JI,LK)
C
      DO 100 IJ=1,NTOTAL
       IF(ISCRR(IJ).NE.0) THEN
        IJTR=ISCRR(IJ) 
        ISCRR(IJTR)=0
        DO 3 L=1,DISSIZE
         SCR(L)=A(L,IJ)+A(ISCRL(L),IJTR)
3       CONTINUE
        DO 4 L=1,DISSIZE
         A(L,IJ)=SCR(L)
4       CONTINUE
CDIR$ IVDEP
        DO 5 L=1,DISSIZE
         A(ISCRL(L),IJTR)=SCR(L)
5       CONTINUE
       ENDIF
100   CONTINUE
      RETURN
      END
