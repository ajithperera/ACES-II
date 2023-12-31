      SUBROUTINE MOOVLP(EVECA,EVECB,S,SCR,SMOAB,NBAS,NCOMP)
C
C THIS ROUTINE FORMS THE ALPHA-BETA MO OVERLAP MATRIX:
C
C      DELTA(Ij), DELTA(Aj), DELTA(Ib) AND DELTA(Ab)
C
C WHICH ARE NEEDED TO COMPUTE THE CORRELATION CORRECTION TO
C  THE SPIN MULTIPLICITY.  ONLY THE SYMMETRY ALLOWED ELEMENTS
C  ARE COMPUTED.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION EVECA(NBAS,NCOMP),EVECB(NBAS,NCOMP),S(NBAS,NBAS)
      DIMENSION SCR(NBAS*NBAS),SMOAB(NCOMP,NCOMP),NBIRR(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      DATA ONE /1.0/
      DATA ONEM/-1.0/
      DATA ZILCH /0.0/
C
C READ IN AO OVERLAP MATRIX AND TRANSFORM TO A-B MO BASIS
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,S)
      CALL GETREC(20,'JOBARC','SCFEVECA',NBAS*NCOMP*IINTFP,EVECA)
      CALL GETREC(20,'JOBARC','SCFEVECB',NBAS*NCOMP*IINTFP,EVECB)
      CALL XGEMM('N','N',NBAS,NCOMP,NBAS,ONE,S,NBAS,
     &           EVECB,NBAS,ZILCH,SCR,NBAS)
      CALL XGEMM('T','N',NCOMP,NCOMP,NBAS,ONE,EVECA,NBAS,
     &           SCR,NBAS,ZILCH,SMOAB,NCOMP)
C#ifdef _DEBUG_LVLM
C      write(6,"(a)") "The MO(A,b) overlap"
C      call output(SMOAB, 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
C      call checksum("MOOVERLAP", SMOAB, Nbas*Nbas)
C#endif
      
C
C NOW WRITE OUT SYMMETRY BLOCKED MO OVERLAPS
C
C FIRST DO OO BLOCK
C
      ITHRU=0    
      IOFFI=0
      IOFFJ=0
      DO 10 IRREP=1,NIRREP
       DO 100 J=1,POP(IRREP,2)
        DO 200 I=1,POP(IRREP,1)
         ITHRU=ITHRU+1
         SCR(ITHRU)=SMOAB(IOFFI+I,IOFFJ+J)
200     CONTINUE
100    CONTINUE
       IOFFI=IOFFI+POP(IRREP,1)
       IOFFJ=IOFFJ+POP(IRREP,2)
10    CONTINUE
      CALL UPDMOI(1,ITHRU,5,90,0,0)
      CALL PUTLST(SCR,1,1,1,5,90)
C
C NOW DO VO BLOCK
C
      ITHRU=0
      IOFFI=NOCCO(1)
      IOFFJ=0
      DO 11 IRREP=1,NIRREP
       DO 101 J=1,POP(IRREP,2)
        DO 201 I=1,VRT(IRREP,1)
         ITHRU=ITHRU+1
         SCR(ITHRU)=SMOAB(IOFFI+I,IOFFJ+J)
201     CONTINUE
101    CONTINUE
       IOFFI=IOFFI+VRT(IRREP,1)
       IOFFJ=IOFFJ+POP(IRREP,2)
11    CONTINUE
      CALL UPDMOI(1,ITHRU,6,90,0,0)
      CALL PUTLST(SCR,1,1,1,6,90)
C
C NOW DO OV BLOCK
C
      ITHRU=0
      IOFFI=0
      IOFFJ=NOCCO(2)
      DO 12 IRREP=1,NIRREP
       DO 102 J=1,VRT(IRREP,2)
        DO 202 I=1,POP(IRREP,1)
         ITHRU=ITHRU+1
         SCR(ITHRU)=SMOAB(IOFFI+I,IOFFJ+J)
202     CONTINUE
102    CONTINUE
       IOFFI=IOFFI+POP(IRREP,1)
       IOFFJ=IOFFJ+VRT(IRREP,2)
12    CONTINUE
      CALL UPDMOI(1,ITHRU,7,90,0,0)
      CALL PUTLST(SCR,1,1,1,7,90)
C
C NOW DO VV BLOCK
C
      ITHRU=0
      IOFFI=NOCCO(1)
      IOFFJ=NOCCO(2)
      DO 13 IRREP=1,NIRREP
       DO 103 J=1,VRT(IRREP,2)
        DO 203 I=1,VRT(IRREP,1)
         ITHRU=ITHRU+1
         SCR(ITHRU)=SMOAB(IOFFI+I,IOFFJ+J)
203     CONTINUE
103    CONTINUE
       IOFFI=IOFFI+VRT(IRREP,1)
       IOFFJ=IOFFJ+VRT(IRREP,2)
13    CONTINUE
      CALL UPDMOI(1,ITHRU,8,90,0,0)
      CALL PUTLST(SCR,1,1,1,8,90)
      RETURN
      END
