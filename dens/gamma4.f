      SUBROUTINE GAMMA4(ICORE,MAXCOR,IUHF)
C
C   THIS ROUTINE CALCULATES THE FOURTH GAMMA INTERMEDIATE
C
C  MBPT(3)
C
C   G(IA,JB) = - SUM M,E  T[1](IM,BE) T[1](MJ,EA)
C
C  ROHF-MBPT(3)
C
C   G(IA,JB) = - SUM M,E  T[1](IM,BE) T[1](MJ,EA)
C
C              + T[1](I,B) T[1](J,A)
C
C  MBPT(4) 
C
C   G(IA,JB) =  - 1/2 P((IA),(JB)) SUM M,E [2 T(IM,BE) - T1(IM,BE)] T1(MJ,EA)
C
C  CCD
C   
C   G(IA,JB) = - 1/2 P((IA),(JB)) SUM M,E T(IM,BE) L(MJ,EA)
C
C  QCISD
C
C   G(IA,JB) = - 1/2 P((IA),(JB)) SUM M,E T(IM,BE) L(MJ,EA)
C
C              + 1/2 P((IA),(JB)) T(I,B) L(J,A)
C
C  CCSD
C
C   G(IA,JB) = - 1/2 P((IA),(JB)) SUM M,E [T(IM,BE)- T(I,E) T(M,B)] L(MJ,EA) 
C
C              + 1/2 P((IA),(JB)) T(I,B) L(J,A)
C
C  THIS TERM IS VERY SIMILAR TO THE T1(IJ,AB) CONTRIBUTION TO
C  THE W-RING INTERMEDIATE
C
C  THE SPIN CASES ARE
C
C    AAAA : =  - SUM M,E T(IM,BE) T(MJ,EA) - SUM m,e T(Im,Be) T(Jm,Ae)
C
C    ABAB : =  - SUM m,E T(Im,Eb) T(Jm,Ea)
C
C    ABBA : =  - SUM M,E T(IM,BE) T(Mj,Ea) - SUM m,e T(Im,Be) T(mj,ea)
C
C FOR UHF IN ADDITION THE BBBB BABA AND BABA SPIN CASES HAS TO CALCULATES
C
C  THE GAMMA INTERMEDIATES ARE STORED ON THE 
C
C  AAAA : LIST 123  (CALCULATED FROM ABAB AND ABBA IN QUIKAA)
C  BBBB : LIST 124  (UHF ONLY)
C  ABBA : LIST 125
C  BAAB : LIST 126  (UHF ONLY)
C  ABAB : LIST 118
C  BABA : LIST 117  (UHF ONLY)
C
CEND
C 
C  CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL CIS,EOM,bRedundant
      DIMENSION ICORE(MAXCOR)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/EXCITE/CIS,EOM
      COMMON/FLAGS2/IFLAGS2(500)
C
      DATA ONE/1.D0/
      bRedundant = (iflags2(155).eq.0)
C
C   CASES REQUIRED FOR RHF AND UHF
C
      CALL G4ALL(ICORE,MAXCOR,'ABBA',IUHF,bRedundant)
      CALL G4ALL(ICORE,MAXCOR,'ABAB',IUHF,bRedundant)
      IF(IUHF.EQ.1) THEN
       CALL G4ALL(ICORE,MAXCOR,'AAAA',IUHF,bRedundant)
C
C   BAAB HAS BEEN CALCULATED ALREADY IN THE CALL WITH ABBA
C
       CALL G4ALL(ICORE,MAXCOR,'BBBB',IUHF,bRedundant)
       CALL G4ALL(ICORE,MAXCOR,'BABA',IUHF,bRedundant)
      ENDIF
C
C  FOR CCSD CALCULATE IN ADDITION A T1**2 L2 CONTRIBUTION
C
      IF(CCSD) CALL T12ING4(ICORE,MAXCOR,IUHF,bRedundant)
C
C  FORM AAAA PART FOR RHF REFERENCE FUNCTIONS
C
      IF(IUHF.EQ.0) THEN
       CALL QUIKAA(ICORE,MAXCOR)
      ENDIF
C
      TWO=2.0D0/DFLOAT(1+IUHF)
      CALL CHECKGAM(ICORE,23,123,TWO)
      CALL CHECKGAM(ICORE,25,125,TWO)
      CALL CHECKGAM(ICORE,18,118,TWO)
      IF(IUHF.NE.0) THEN
       CALL CHECKGAM(ICORE,24,124,TWO)
       CALL CHECKGAM(ICORE,26,126,TWO)
       CALL CHECKGAM(ICORE,17,117,TWO)
      ENDIF
C
      RETURN
      END