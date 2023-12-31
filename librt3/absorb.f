C
C     trplib.f is a collection of well over 100 subroutines which can be
C     used with vcc.f, appropriate ACES II (CRAPS) libraries, and other
C     modules for CCSDT energy calculations for arbitrary orthonormal
C     reference determinants and CCSD(T) and QCISD(T) energy derivative
C     calculations.
C
      SUBROUTINE ABSORB(ABSVRT,ABSOCC,MAXNBF)
      IMPLICIT INTEGER (A-Z)
C
C     NOTE THAT USE OF MAXNBF IS A BIT WASTEFUL
C
      DIMENSION ABSVRT(MAXNBF,8,2),ABSOCC(MAXNBF,8,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
C
C     SUBROUTINE TO COMPUTE TABLES SUCH THAT GIVEN AN OCCUPIED OR
C     VIRTUAL ORBITAL IN A GIVEN SYMMETRY BLOCK, WE CAN GET IT'S
C     "ABSOLUTE" POSITION, IE THE POSITION IN THE EVAL ARRAY.
C
      CALL IZERO(ABSVRT,MAXNBF*16)
      CALL IZERO(ABSOCC,MAXNBF*16)
C
      DO   30 ISPIN=1,2
      ABSVAL = 0
      DO   20 IRPA=1,NIRREP
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO 20
      DO   10    A=1,VRT(IRPA,ISPIN)
      ABSVAL = ABSVAL + 1
      ABSVRT(A,IRPA,ISPIN) = ABSVAL
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
C
      DO  130 ISPIN=1,2
      ABSVAL = 0
      DO  120  IRPI=1,NIRREP
      IF(POP(IRPI,ISPIN).EQ.0) GOTO 120
      DO  110     I=1,POP(IRPI,ISPIN)
      ABSVAL = ABSVAL + 1
      ABSOCC(I,IRPI,ISPIN) = ABSVAL
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
      RETURN
      END
