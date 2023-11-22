      SUBROUTINE OFFT3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      COMMON /FLAGS/ IFLAGS(100)
      EQUIVALENCE(IPRINT,IFLAGS(1))
C
C     This subroutine calculates some offsets for T3 i/o.
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT,MOINTS
C
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
C     IJKPOS - For a given ordering of IRPI,IRPJ,IRPK (i.e. either
C              IRPI LE IRPJ ; IRPK or IRPI LE IRPJ LE IRPK), this specifies
C              the position of this triplet within the group of IRPIJK
C              triplets. For the first ordering, there appear to be
C              h * (h+1)/2 triplets for each IRPIJK, while for the second
C              ordering, there are apparently (h+1)*(h+2)/6 triplets.
C              The first is the origin of the 36 in the dimensions of
C              IJKLEN and IJKOFF. Note that h must be even for these
C              dimensions to work (in fact it must 2**n).
C
C     IJKLEN - The length of the IJKth triplet of the IRPIJKth irrep for
C              spin case ISPIN.
C
C     IJKOFF - The offset of the IJKth triplet of the IRPIJKth irrep for
C              spin case ISPIN.
C
      CALL IZERO(IJKPOS, 8*8*8*2)
      CALL IZERO(IJKLEN,36*8*4)
      CALL IZERO(IJKOFF,36*8*4)
C
C     IRPI LE IRPJ LE IRPK
C
      IF(IPRINT.GT.10) WRITE(LUOUT,1010)
C
      IORDER = 1
      DO  50 IRPIJK=1,NIRREP
C
      IJK = 0
      DO  40   IRPK=1,NIRREP
      IRPIJ = DIRPRD(IRPK,IRPIJK)
C
      DO  30   IRPJ=1,IRPK
      IRPI  = DIRPRD(IRPJ,IRPIJ)
C
      IF(IRPI.GT.IRPJ) GOTO 30
C
      IJK = IJK + 1
      IJKPOS(IRPI,IRPJ,IRPK,IORDER) = IJK
      IF(IPRINT.GT.10)THEN
      WRITE(LUOUT,1020) IORDER,IRPIJK,IRPI,IRPJ,IRPK,
     1                         IJKPOS(IRPI,IRPJ,IRPK,IORDER)
      ENDIF
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
C
C     IRPI LE IRPJ; IRPK
C
      IF(IPRINT.GT.10) WRITE(LUOUT,1010)
C
      IORDER = 2
      DO  80 IRPIJK=1,NIRREP
C
      IJK = 0
      DO  70   IRPK=1,NIRREP
      IRPIJ = DIRPRD(IRPK,IRPIJK)
C
      DO  60   IRPJ=1,NIRREP
      IRPI  = DIRPRD(IRPJ,IRPIJ)
C
      IF(IRPI.GT.IRPJ) GOTO 60
C
      IJK = IJK + 1
      IJKPOS(IRPI,IRPJ,IRPK,IORDER) = IJK
      IF(IPRINT.GT.10)THEN
      WRITE(LUOUT,1020) IORDER,IRPIJK,IRPI,IRPJ,IRPK,
     1                         IJKPOS(IRPI,IRPJ,IRPK,IORDER)
      ENDIF
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
C     Compute the lengths of the various blocks.
C
C     IRPI LE IRPJ LE IRPK
C
      DO 140  ISPIN=1,2
      DO 130 IRPIJK=1,NIRREP
C
      IJK = 0
      DO 120   IRPK=1,NIRREP
      IRPIJ = DIRPRD(IRPK,IRPIJK)
C
      DO 110   IRPJ=1,IRPK
      IRPI  = DIRPRD(IRPJ,IRPIJ)
C
      IF(IRPI.GT.IRPJ) GOTO 110
C
      IJK = IJK + 1
C
      IF(IRPI.EQ.IRPJ.AND.IRPI.EQ.IRPK)THEN
      IJKLEN(IJK,IRPIJK,1 + 3*(ISPIN-1)) = (POP(IRPK,ISPIN) *
     1 (POP(IRPK,ISPIN)-1) * (POP(IRPK,ISPIN)-2)) / 6
      ENDIF
C
      IF(IRPI.EQ.IRPJ.AND.IRPI.LT.IRPK)THEN
      IJKLEN(IJK,IRPIJK,1 + 3*(ISPIN-1)) = (POP(IRPK,ISPIN) *
     1  POP(IRPJ,ISPIN)    * (POP(IRPJ,ISPIN)-1)) / 2
      ENDIF
C
      IF(IRPI.LT.IRPJ.AND.IRPJ.EQ.IRPK)THEN
      IJKLEN(IJK,IRPIJK,1 + 3*(ISPIN-1)) = (POP(IRPK,ISPIN) *
     1 (POP(IRPK,ISPIN)-1) *  POP(IRPI,ISPIN))    / 2
      ENDIF
C
      IF(IRPI.LT.IRPJ.AND.IRPJ.LT.IRPK)THEN
      IJKLEN(IJK,IRPIJK,1 + 3*(ISPIN-1)) =  POP(IRPK,ISPIN) *
     1  POP(IRPJ,ISPIN)    *  POP(IRPI,ISPIN)
      ENDIF
C
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
C     IRPI LE IRPJ; IRPK
C
      DO 190 ISPIN=2,3
C
      IF(ISPIN.EQ.2)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      DO 180 IRPIJK=1,NIRREP
C
      IJK = 0
      DO 170   IRPK=1,NIRREP
      IRPIJ = DIRPRD(IRPK,IRPIJK)
C
      DO 160   IRPJ=1,NIRREP
      IRPI  = DIRPRD(IRPJ,IRPIJ)
C
      IF(IRPI.GT.IRPJ) GOTO 160
C
      IJK = IJK + 1
C
      IF(IRPI.EQ.IRPJ)THEN
      IJKLEN(IJK,IRPIJK,ISPIN) = (POP(IRPK,ISPIN2) *
     1  POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1)) /2
      ENDIF
C
      IF(IRPI.LT.IRPJ)THEN
      IJKLEN(IJK,IRPIJK,ISPIN) =  POP(IRPK,ISPIN2) *
     1  POP(IRPJ,ISPIN1) *  POP(IRPI,ISPIN1)
      ENDIF
C
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
C
      IF(IPRINT.GT.10) WRITE(LUOUT,1030)
C
C     Compute the offsets of the blocks.
C
      NCOMB(1) = ((NIRREP + 1) * (NIRREP + 2))/6
      NCOMB(4) = ((NIRREP + 1) * (NIRREP + 2))/6
      NCOMB(2) = ((NIRREP + 1) *  NIRREP     )/2
      NCOMB(3) = ((NIRREP + 1) *  NIRREP     )/2
C
      DO 230  ISPIN=1,4
      DO 220 IRPIJK=1,NIRREP
      DO 210    IJK=1,NCOMB(ISPIN)
C
      IF(IJK.EQ.1)THEN
      IJKOFF(IJK,IRPIJK,ISPIN) = 0
      ELSE
      IJKOFF(IJK,IRPIJK,ISPIN) = IJKOFF(IJK-1,IRPIJK,ISPIN) +
     1                           IJKLEN(IJK-1,IRPIJK,ISPIN)
      ENDIF
      IF(IPRINT.GT.10)THEN
      WRITE(LUOUT,1040) IRPIJK,ISPIN,IJK,IJKLEN(IJK,IRPIJK,ISPIN),
     1                                   IJKOFF(IJK,IRPIJK,ISPIN)
      ENDIF
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
C
      RETURN
 1010 FORMAT(' @OFFT3-I, Order Symmetry IRPI IRPJ IRPK Position ',/,
     1       '           -------------------------------------- ')
 1020 FORMAT(13X,I1,7X,I1,6X,I1,4X,I1,4X,I1,5X,I2)
 1030 FORMAT(' @OFFT3-I, Symmetry Spin IJK   Length   Offset ',/,
     1       '           ----------------------------------- ')
 1040 FORMAT(15X,I1,6X,I1,3X,I2,3X,I6,3X,I6)
      END
