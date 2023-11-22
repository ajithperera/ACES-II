      SUBROUTINE DIMT3(DISSIZ,NDIS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZ,A,B,C,DIRPRD,POP,VRT
C
C     This subroutine calculates the distribution sizes and the numbers
C     of distributions of the different symmetries and spin cases of
C     T3 amplitudes.
C
      DIMENSION DISSIZ(8,4),NDIS(8,4)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      EQUIVALENCE(IPRINT,IFLAGS(1))
C
      CALL IZERO(DISSIZ,8*4)
      CALL IZERO(NDIS  ,8*4)
C
C     AAA/BBB
C
      DO  30 ISPIN=1,2
C
      DO 233 IRPC=1,NIRREP
      DO 232 IRPB=1,IRPC
      DO 231 IRPA=1,IRPB
      IRPBC  = DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPA.EQ.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) * (VRT(IRPC,ISPIN) - 1) *
     1                            (VRT(IRPC,ISPIN) - 2)) / 6
      ENDIF
C
      IF(IRPA.EQ.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) *  VRT(IRPB,ISPIN)      *
     1                            (VRT(IRPB,ISPIN) - 1)) / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) * (VRT(IRPC,ISPIN) - 1) *
     1                             VRT(IRPA,ISPIN)    )  / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH =  VRT(IRPC,ISPIN) *  VRT(IRPB,ISPIN)      *
     1                             VRT(IRPA,ISPIN)
      ENDIF
C
      DISSIZ(IRPABC,1+3*(ISPIN-1)) = DISSIZ(IRPABC,1+3*(ISPIN-1)) +
     1                               LENGTH
  231 CONTINUE
  232 CONTINUE
  233 CONTINUE
C
c      DO  10 C=3,NVRTO(ISPIN)
c      DO  10 B=2,C-1
c      DO  10 A=1,B-1
c      IRPA = IRREPS(NOCCO(ISPIN) + A,ISPIN)
c      IRPB = IRREPS(NOCCO(ISPIN) + B,ISPIN)
c      IRPC = IRREPS(NOCCO(ISPIN) + C,ISPIN)
c      IRPBC = DIRPRD(IRPB,IRPC)
c      IRPABC = DIRPRD(IRPBC,IRPA)
c      DISSIZ(IRPABC,1+3*(ISPIN-1)) = DISSIZ(IRPABC,1+3*(ISPIN-1)) + 1
c   10 CONTINUE
C
c      DO  20 K=3,NOCCO(ISPIN)
c      DO  20 J=2,K-1
c      DO  20 I=1,J-1
c      IRPI = IRREPS(               I,ISPIN)
c      IRPJ = IRREPS(               J,ISPIN)
c      IRPK = IRREPS(               K,ISPIN)
c      IRPJK = DIRPRD(IRPJ,IRPK)
c      IRPIJK = DIRPRD(IRPJK,IRPI)
c      NDIS(IRPIJK,1+3*(ISPIN-1)) = NDIS(IRPIJK,1+3*(ISPIN-1)) + 1
c   20 CONTINUE
C
      DO 333 IRPC=1,NIRREP
      DO 332 IRPB=1,IRPC
      DO 331 IRPA=1,IRPB
      IRPBC  = DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
      IF(IRPA.EQ.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (POP(IRPC,ISPIN) * (POP(IRPC,ISPIN) - 1) *
     1                            (POP(IRPC,ISPIN) - 2)) / 6
      ENDIF
C
      IF(IRPA.EQ.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH = (POP(IRPC,ISPIN) *  POP(IRPB,ISPIN)      *
     1                            (POP(IRPB,ISPIN) - 1)) / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (POP(IRPC,ISPIN) * (POP(IRPC,ISPIN) - 1) *
     1                             POP(IRPA,ISPIN)    )  / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH =  POP(IRPC,ISPIN) *  POP(IRPB,ISPIN)      *
     1                             POP(IRPA,ISPIN)
      ENDIF
C
      NDIS(IRPABC,1+3*(ISPIN-1)) = NDIS(IRPABC,1+3*(ISPIN-1)) +
     1                             LENGTH
  331 CONTINUE
  332 CONTINUE
  333 CONTINUE
   30 CONTINUE
C
C     AAB/BBA
C
      DO  60 IPASS=1,2
C
      IF(IPASS.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      DO  50 IRPC =1,NIRREP
      DO  50 IRPAB=1,NIRREP
      IRPABC=DIRPRD(IRPAB,IRPC)
      LEN   =VRT(IRPC,ISPIN2) * IRPDPD(IRPAB,ISPIN1)
      DISSIZ(IRPABC,ISPIN1 + 1) = DISSIZ(IRPABC,ISPIN1 + 1) + LEN
C
      LEN   =POP(IRPC,ISPIN2) * IRPDPD(IRPAB,ISPIN1 + 2)
      NDIS(IRPABC,ISPIN1 + 1) = NDIS(IRPABC,ISPIN1 + 1) + LEN
   50 CONTINUE
   60 CONTINUE

      IF(IPRINT.GT.10)THEN
      WRITE(LUOUT,1010)
      DO  70 ISPIN=1,4
      DO  70 IRREP=1,NIRREP
      WRITE(LUOUT,1020) DISSIZ(IRREP,ISPIN),NDIS(IRREP,ISPIN),
     1                  IRREP,ISPIN
   70 CONTINUE
      ENDIF
C
      RETURN
 1010 FORMAT(' @DIMT3-I, Dimensions of T3 amplitudes : ',/,
     1 '               # abc          # ijk       Symmetry    Spin ')
 1020 FORMAT(10X,2I15,5X,I1,10X,I1)
      END
