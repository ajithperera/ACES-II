C
      SUBROUTINE QSYMRHF(IRREPL, IRREPR, NUML, NUMR, DISSIZE, A,
     &                   SCR, ISCRL, ISCRR)
C
C This routine forms the term 
C    A(IJ,KL) --> A(IJ,KL) + A(JI,LK)  
C which is required in RHF CCSD calculations
C
C Input arguments
C
C  IRREPL           : Irrep of the left hand side
C  IRREPR           : Irrep of the right hand side
C  NUMR             : Population vector for K, L
C  NUML             : Population vector for I,J
C  DISSIZE          : Distribution size of A
C  A                : Holds the matrix A         
C  SCR,ISCRL,ISCRR  : These scratch arrays of size DISSIZE
C
C Output argument
C
C  A               : The RHF symmetrized matrix
C
C  Written in SEPTEMBER/90 by JG Modified JULY/1994 Ajith
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZE, DIRPRD
      DIMENSION A(DISSIZE,1),SCR(1),ISCRL(1),ISCRR(2),NUMR(8),NUML(8)
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &               DIRPRD(8,8)
      DIMENSION IPL(8),IPR(8)
C     
      IPR(1) = 0
      IPL(1) = 0          
C
      DO 1000 IRREPJR = 1, NIRREP - 1
         IRREPIR = DIRPRD(IRREPR, IRREPJR)
         IPR(IRREPJR + 1) = IPR(IRREPJR) + NUMR(IRREPJR)*NUMR(IRREPIR)
 1000 CONTINUE
C
       DO 2000 IRREPJL = 1, NIRREP - 1
         IRREPIL = DIRPRD(IRREPL, IRREPJL)
         IPL(IRREPJL + 1) = IPL(IRREPJL) + NUML(IRREPJL)*NUML(IRREPIL)
 2000 CONTINUE
C     
C First get the new addresses and store in ISCR, Address for the
C right hand side.
C     
      NTOTAL = 0
      DO 20  IRREPJR = 1, NIRREP
         IRREPIR = DIRPRD(IRREPR, IRREPJR)
C
         NUMJR  = NUMR(IRREPJR)
         NUMIR  = NUMR(IRREPIR)
         NTOTAL = NTOTAL + NUMIR*NUMJR
C
         DO 15 J = 1, NUMJR
            DO 16 I = 1, NUMIR
               IND1 = IPR(IRREPJR) + (J - 1)*NUMIR + I
               IND2 = IPR(IRREPIR) + (I - 1)*NUMJR + J
               ISCRR(IND1) = IND2
 16         CONTINUE
 15      CONTINUE
C
 20   CONTINUE   
C     
C Addresses for the left side
C     
      DO 120  IRREPJL = 1, NIRREP
         IRREPIL = DIRPRD(IRREPL, IRREPJL)
C     
         NUMJL = NUML(IRREPJL)
         NUMIL = NUML(IRREPIL)
C
         DO 115 J = 1, NUMJL
            DO 116 I = 1, NUMIL
               IND1 = IPL(IRREPJL) + (J - 1)*NUMIL + I
               IND2 = IPL(IRREPIL) + (I - 1)*NUMJL + J
               ISCRL(IND1) = IND2
 116        CONTINUE
 115     CONTINUE
C
 120  CONTINUE        
C     
C Now form A(IJ, KL) + A(JI, LK)
C     
      DO 100 IJ = 1, NTOTAL
C
         IF(ISCRR(IJ) .NE. 0) THEN
            IJTR = ISCRR(IJ) 
            ISCRR(IJTR) = 0
C
            DO 3 L = 1, DISSIZE
               SCR(L) = A(L, IJ) + A(ISCRL(L), IJTR)
 3          CONTINUE
C    
            DO 4 L = 1, DISSIZE
               A(L, IJ) = SCR(L)
 4          CONTINUE
C
            DO 5 L = 1, DISSIZE
               A(ISCRL(L), IJTR) = SCR(L)
 5          CONTINUE
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END
