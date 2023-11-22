      SUBROUTINE T3WT314P(D3T3EXP,CORE,MAXCOR,ISPIN,LEN,LENEXP,
     1                    IRPIJK,IJKVAL,IADBLK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD,DISSIZ
      DIMENSION D3T3EXP(LENEXP),CORE(1)
      DIMENSION IADBLK(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT1(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      COMMON /FILES/  LUOUT,MOINTS
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      INDEX(I) = I*(I-1)/2
C
      I000  = 1
      I010  = I000 + LEN
      I020  = I010 + LENEXP
      NEED  = I020 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3WT314P-F, Insufficient memory. Need ',I10,' Got ',I10)
      CALL INSMEM('T3WT314P',NEED,MAXCOR)
      ENDIF
C
      CALL GETLIST(CORE(I000),IJKVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)
C
      CALL ZERO(CORE(I010),LENEXP)
      CALL SYMEXPT3(CORE(I000),CORE(I010),IADBLK,ISPIN,IRPIJK)
C
      IOFF1 = I010
      IOFF2 = 1
      DO   10 IRPC=1,NIRREP
      IRPEF = DIRPRD(IRPIJK,IRPC)
      IRPAB = IRPEF
C
      I030 = I020 + IRPDPD(IRPAB,ISPIN) * IRPDPD(IRPEF,ISPIN)
      NEED  = I030 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
      CALL INSMEM('T3WT314P',NEED,MAXCOR)
      ENDIF
C
      IF(VRT(IRPC,ISPIN).EQ.0.OR.IRPDPD(IRPAB,ISPIN).EQ.0) GOTO 10
C
      CALL GETLST(CORE(I020),1,IRPDPD(IRPEF,ISPIN),2,IRPEF,230+ISPIN)
C
C     D3T3EXP(A<B,C) = SUM_{EF} W(A<B,E<F) * T3EXPOLD(E<F,C)
C
      DISSIZ = IRPDPD(IRPAB,ISPIN)
      NDIS   = VRT(IRPC,ISPIN)
      NSUM   = DISSIZ
      CALL XGEMM('N','N',DISSIZ,NDIS,NSUM,1.0D+00,
     1           CORE(I020),DISSIZ,CORE(IOFF1),NSUM,1.0D+00,
     1           D3T3EXP(IOFF2),DISSIZ)
C
      IOFF1 = IOFF1 + NSUM   * NDIS
      IOFF2 = IOFF2 + DISSIZ * NDIS
   10 CONTINUE
C
      RETURN
      END
