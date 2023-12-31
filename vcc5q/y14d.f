      SUBROUTINE Y14D(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZY,B,E,BE,BELM
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      WRITE(6,1000)
 1000 FORMAT(' @Y14D-I, Y1 * Y2 contributions. ')
C
C     Y(B,E,L,M) = - D(B,E) * D(L,M)
C
C     BELM = BE * LM
C     belm = be * lm
C
      DO  100 ISPIN=1,IUHF+1
C
      I000 = 1
      I010 = I000 + NFEA(ISPIN)
      I020 = I010 + NFMI(ISPIN)
      I030 = I020 + IRPDPD(1,18+ISPIN) * IRPDPD(1,20+ISPIN)
      I040 = I030 + IRPDPD(1,18+ISPIN) * IRPDPD(1,20+ISPIN)
C
      NEED = I040 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
 1020 FORMAT(' @Y13D-I, Insufficient memory. Need ',I15,' . Got ',I15)
      STOP 'Y14D'
      ENDIF
C
      CALL ZERO(CORE(I020),IRPDPD(1,18+ISPIN) * IRPDPD(1,20+ISPIN))
C
      CALL GETLST(CORE(I000),1,1,2,ISPIN,92)
      CALL GETLST(CORE(I010),1,1,2,ISPIN,91)
C
      IOFFY  = I020 - 1
      DO   90 IRPM=1,NIRREP
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 90
      IRPL = IRPM
      DO   80 IRPE=1,NIRREP
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 80
      IRPB = IRPE
C
      IOFFBE = IOFFVV(IRPE,1,2+ISPIN)
      IOFFLM = IOFFOO(IRPM,1,2+ISPIN)

      DO   70 M=1,POP(IRPM,ISPIN)
      DO   60 L=1,POP(IRPL,ISPIN)
      DO   50 E=1,VRT(IRPE,ISPIN)
      DO   40 B=1,VRT(IRPB,ISPIN)
C
      BE = IOFFBE + (E-1)*VRT(IRPB,ISPIN) + B
C
      LM = IOFFLM + (M-1)*POP(IRPL,ISPIN) + L

C
      BELM = (LM-1) * IRPDPD(1,18+ISPIN) + BE
C
      CORE(IOFFY + BELM) = CORE(I000 - 1 + BE) *
     1                     CORE(I010 - 1 + LM)
   40 CONTINUE
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
C
   80 CONTINUE
   90 CONTINUE
C
      DISSIZY = IRPDPD(1,18+ISPIN)
Csurely next line is right, not :      NDISY   = IRPDPD(1,23-ISPIN)
      NDISY   = IRPDPD(1,20+ISPIN)
      CALL GETLIST(CORE(I030),1,NDISY,2,1,24+ISPIN)
      NSIZ = DISSIZY * NDISY
      CALL    VADD(CORE(I030),CORE(I030),CORE(I020),NSIZ, 1.0D+00)
c     CALL    VADD(CORE(I030),CORE(I030),CORE(I020),NSIZ,-1.0D+00)
      CALL PUTLIST(CORE(I030),1,NDISY,2,1,24+ISPIN)
  100 CONTINUE
C
C     BElm = BE * lm
C     beLM = be * LM
C
      DO  200 ISPIN=1,IUHF+1
C
      I000 = 1
      I010 = I000 + NFEA(  ISPIN)
      I020 = I010 + NFMI(3-ISPIN)
      I030 = I020 + IRPDPD(1,18+ISPIN) * IRPDPD(1,23-ISPIN)
      I040 = I030 + IRPDPD(1,18+ISPIN) * IRPDPD(1,23-ISPIN)
C
      NEED = I040 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y14D'
      ENDIF
C
      CALL ZERO(CORE(I020),IRPDPD(1,18+ISPIN) * IRPDPD(1,23-ISPIN))
C
      CALL GETLST(CORE(I000),1,1,2,  ISPIN       ,92)
      CALL GETLST(CORE(I010),1,1,2,3-ISPIN-1+IUHF,91)
C
      IOFFY  = I020 - 1
      DO  190 IRPM=1,NIRREP
      IF(POP(IRPM,3-ISPIN).EQ.0) GOTO 190
      IRPL = IRPM
      DO  180 IRPE=1,NIRREP
      IF(VRT(IRPE,  ISPIN).EQ.0) GOTO 180
      IRPB = IRPE
C
      IOFFBE = IOFFVV(IRPE,1,2+ISPIN)
      IOFFLM = IOFFOO(IRPM,1,5-ISPIN)

      DO  170 M=1,POP(IRPM,3-ISPIN)
      DO  160 L=1,POP(IRPL,3-ISPIN)
      DO  150 E=1,VRT(IRPE,ISPIN)
      DO  140 B=1,VRT(IRPB,ISPIN)
C
      BE = IOFFBE + (E-1)*VRT(IRPB,  ISPIN) + B
C
      LM = IOFFLM + (M-1)*POP(IRPL,3-ISPIN) + L
C
      BELM = (LM-1) * IRPDPD(1,18+ISPIN) + BE
C
      CORE(IOFFY + BELM) = CORE(I000 - 1 + BE) *
     1                     CORE(I010 - 1 + LM)
  140 CONTINUE
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
C
  180 CONTINUE
  190 CONTINUE
C
      DISSIZY = IRPDPD(1,18+ISPIN)
      NDISY   = IRPDPD(1,23-ISPIN)
      CALL GETLIST(CORE(I030),1,NDISY,2,1,28+ISPIN)
      NSIZ = DISSIZY * NDISY
      CALL    VADD(CORE(I030),CORE(I030),CORE(I020),NSIZ, 1.0D+00)
C     CALL    VADD(CORE(I030),CORE(I030),CORE(I020),NSIZ,-1.0D+00)
      CALL PUTLIST(CORE(I030),1,NDISY,2,1,28+ISPIN)
  200 CONTINUE
      RETURN
      END
