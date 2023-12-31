      SUBROUTINE T3FT323O(T3,CORE,MAXCOR,IUHF,ISPIN1,ISPIN2,LEN,
     1                    IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD
      LOGICAL NONEQL,IJEQL,MJEQL,IMEQL
      DIMENSION T3(LEN),CORE(1)
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
      I010  = I000 + MAX(NFMI(1),NFMI(2))
      I020  = I010 + LEN
      NEED  = I020 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3FT323O-F, Insufficient memory. Need ',I10,' Got ',I10)
      CALL INSMEM('T3FT323O',NEED,MAXCOR)
      ENDIF
C
C     T3(ABc,IJm) * F(m,k)
C     T3(abC,ijM) * F(M,K)
C
C     If this is a CCSDT calculation, we read the FMI intermediate.
C     Otherwise, we read the OCC-OCC block of the Fock matrix.
C
      IF(ICLLVL.EQ.18)THEN
      LISTF = 91
      LPART = ISPIN2
      IF(IUHF.EQ.0) LPART = 1
      ELSE
      LISTF = 91
      LPART = ISPIN2 + 2
      IF(IUHF.EQ.0) LPART = 3
      ENDIF
C
      CALL GETLST(CORE(I000),1,1,2,LPART,LISTF)
C     CALL RMDIAG(CORE(I000),POP(1,ISPIN2),NIRREP)
C
      IRPM = IRPK
C
      IJM  = IJKPOS(IRPI,IRPJ,IRPM,2)
      IOFF = IJKOFF(IJM,IRPIJK,ISPIN1 + 1)
C
      IJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ) IJEQL   = .TRUE.
      IF(IRPI.NE.IRPJ) NONEQL  = .TRUE.
C
      IF(IJEQL)THEN
      NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1) - 1))/2
      IJ  = INDEX(J-1) + I
      ELSE
      NIJ = POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
      IJ  = (J-1) * POP(IRPI,ISPIN1) + I
      ENDIF
C
      IF(POP(IRPM,ISPIN2).GT.0)THEN
C
      DO   10 M=1,POP(IRPM,ISPIN2)
C
C no, silly      IF(M.EQ.K) GOTO 10
C
      IJMVAL = IOFF + (M-1)*NIJ + IJ
C
      FMK = CORE(I000 - 1 + IOFFOO(IRPM,1,2+ISPIN2) + 
     1                      (K-1)*POP(IRPM,ISPIN2) + M)
C
      CALL GETLIST(CORE(I010),IJMVAL,1,1,IRPIJK,ISPIN1 + 1 + 4)
      CALL VADD(T3,T3,CORE(I010),LEN,-FMK)
   10 CONTINUE
C
      ENDIF
C
C
C     T3(ABc,MJk) * F(M,I)    -    T3(ABc,MIk) * F(M,J)
C     T3(abC,mjK) * F(m,i)    -    T3(abC,miK) * F(m,j)
C
C     If this is a CCSDT calculation, we read the FMI intermediate.
C     Otherwise, we read the OCC-OCC block of the Fock matrix.
C
      IF(ICLLVL.EQ.18)THEN
      LISTF = 91
      LPART = ISPIN1
      ELSE
      LISTF = 91
      LPART = ISPIN1 + 2
      ENDIF
C
      CALL GETLST(CORE(I000),1,1,2,LPART,LISTF)
C     CALL RMDIAG(CORE(I000),POP(1,ISPIN1),NIRREP)
      DO   50 IRPM=1,NIRREP
C
      IF(IRPM.NE.IRPI.AND.IRPM.NE.IRPJ) GOTO 50
      IF(POP(IRPM,ISPIN1).LE.0)         GOTO 50
C
      IF(IRPM.EQ.IRPI)THEN
C
C      T(ABC,MJK) * F(M,I)
C
      MJK  = IJKPOS(IRPM,IRPJ,IRPK,2)
      IOFF = IJKOFF(MJK,IRPIJK,ISPIN1 + 1)
C
      MJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPJ) MJEQL   = .TRUE.
      IF(IRPM.NE.IRPJ) NONEQL  = .TRUE.
C
      IF(MJEQL)THEN
      NMJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1) - 1))/2
      ELSE
      NMJ = POP(IRPM,ISPIN1) * POP(IRPJ,ISPIN1)
      ENDIF
C
      DO   20 M=1,POP(IRPM,ISPIN1)
      SIGN = 1.0D+00
C
      IF(IRPJ.EQ.IRPM.AND.J.EQ.M) GOTO 20
C
      IF(NONEQL) MJKVAL = IOFF + (K-1)*NMJ + (J-1)*POP(IRPM,ISPIN1) + M
C
      IF(MJEQL)THEN
      MJ = INDEX(MAX(J,M)-1) + MIN(J,M)
      MJKVAL = IOFF + (K-1) * NMJ + MJ
      IF(J.LT.M) SIGN = -1.0D+00
      ENDIF
C
      FMI = CORE(I000 - 1 + IOFFOO(IRPM,1,2+ISPIN1) + 
     1                      (I-1)*POP(IRPM,ISPIN1) + M) * SIGN
C
      CALL GETLIST(CORE(I010),MJKVAL,1,1,IRPIJK,ISPIN1 + 1 + 4)
C
      CALL VADD(T3,T3,CORE(I010),LEN,-FMI)
   20 CONTINUE
C
      ENDIF
C
      IF(IRPM.EQ.IRPJ)THEN
C
C     T(ABC,IMK) * F(M,J)
C
      IMK  = IJKPOS(IRPI,IRPM,IRPK,2)
      IOFF = IJKOFF(IMK,IRPIJK,ISPIN1 + 1)
C
      IMEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPI) IMEQL   = .TRUE.
      IF(IRPM.NE.IRPI) NONEQL  = .TRUE.
C
      IF(IMEQL)THEN
      NIM = (POP(IRPI,ISPIN1) * (POP(IRPI,ISPIN1) - 1))/2
      ELSE
      NIM = POP(IRPM,ISPIN1) * POP(IRPI,ISPIN1)
      ENDIF
C
      DO   30 M=1,POP(IRPM,ISPIN1)
      SIGN = 1.0D+00
C
      IF(IRPI.EQ.IRPM.AND.I.EQ.M) GOTO 30
C
      IF(NONEQL) IMKVAL = IOFF + (K-1)*NIM + (M-1)*POP(IRPI,ISPIN1) + I
C
      IF(IMEQL)THEN
      IM = INDEX(MAX(I,M)-1) + MIN(I,M)
      IMKVAL = IOFF + (K-1) * NIM + IM
      IF(I.GT.M) SIGN = -1.0D+00
      ENDIF
C
      FMJ = CORE(I000 - 1 + IOFFOO(IRPM,1,2+ISPIN1) + 
     1                      (J-1)*POP(IRPM,ISPIN1) + M) * SIGN
C
      CALL GETLIST(CORE(I010),IMKVAL,1,1,IRPIJK,ISPIN1 + 1 + 4)
      CALL VADD(T3,T3,CORE(I010),LEN,-FMJ)
   30 CONTINUE
C
      ENDIF
C
   50 CONTINUE
      RETURN
      END
