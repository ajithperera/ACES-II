      SUBROUTINE T3WT323HH(T3,CORE,MAXCOR,IUHF,ISPIN1,ISPIN2,LEN,
     1                     IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD
      LOGICAL MNEQL,JMEQL,IMEQL,NONEQL
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
C     T(ABc,MNk) * W(MN,IJ) - T(ABc,JMn) * W(Mn,Ik) + T(ABc,IMn) * W(Mn,Jk)
C     T(abC,mnK) * W(mn,ij) - T(abC,jmN) * W(Nm,Ki) + T(abC,imN) * W(Nm,Kj)
C
      INDEX(I) = I*(I-1)/2
C
      IRPIJ = DIRPRD(IRPI,IRPJ)
      IRPIK = DIRPRD(IRPI,IRPK)
      IRPJK = DIRPRD(IRPJ,IRPK)
C
      I000  = 1
      I001  = I000 + MAX(IRPDPD(IRPIJ,ISPIN1+2),
     1                   IRPDPD(IRPIJ,14),
     1                   IRPDPD(IRPIK,14),IRPDPD(IRPJK,14))
      I002  = I001 + MAX(IRPDPD(IRPIJ,ISPIN1+2),
     1                   IRPDPD(IRPIJ,14),
     1                   IRPDPD(IRPIK,14),IRPDPD(IRPJK,14))
      I010  = I002 + MAX(IRPDPD(IRPIJ,ISPIN1+2),
     1                   IRPDPD(IRPIJ,14),
     1                   IRPDPD(IRPIK,14),IRPDPD(IRPJK,14))
      I020  = I010 + LEN
      NEED  = I020 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3WT323HH-F, Insufficient memory. Need ',I10,' Got ',
     1                                                   I10)
      CALL INSMEM('T3WT323HH',NEED,MAXCOR)
      ENDIF
C
      DO  100 IRPMN=1,NIRREP
C
      IF(IRPMN.NE.IRPIJ)              GOTO 100
      IF(IRPDPD(IRPMN,ISPIN1+2).EQ.0) GOTO 100
C
C     T(ABc,MNk) * W(MN,IJ)
C     T(abC,mnK) * W(mn,ij)
C
      IF(IUHF.EQ.0)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,5)      + (J-1)*POP(IRPI,1) + I
      JI = IOFFOO(IRPI,IRPIJ,5)      + (I-1)*POP(IRPJ,1) + J
      IF(ICLLVL.EQ.18)THEN
      LIST = 53
      ELSE
      LIST = 13
      ENDIF
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LIST)
      CALL GETLST(CORE(I001),JI,1,2,IRPIJ,LIST)
      CALL   VADD(CORE(I000),CORE(I000),CORE(I001),IRPDPD(IRPIJ,14),
     1            -1.0D+00)
      ELSE
C
      IF(IRPIJ.EQ.1)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 50 + ISPIN1
      ELSE
      LIST = 10 + ISPIN1
      ENDIF
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LIST)
      ENDIF
C
      DO   40 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(IRPM.GT.IRPN)                                   GOTO 40
      IF(POP(IRPM,ISPIN1).EQ.0.OR.POP(IRPN,ISPIN1).EQ.0) GOTO 40
      IF(IRPM.EQ.IRPN.AND.POP(IRPM,ISPIN1).LT.2)         GOTO 40
C
      MNK  = IJKPOS(IRPM,IRPN,IRPK,2)
      IOFF = IJKOFF(MNK,IRPIJK,ISPIN1 + 1)
C
      MNEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPN) MNEQL  = .TRUE.
      IF(IRPM.NE.IRPN) NONEQL = .TRUE.
C
      IF(MNEQL)THEN
      NLOW  = 2
      NHIGH = POP(IRPN,ISPIN1)
      MLOW  = 1
      ELSE
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN1)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN1)
      ENDIF
C
      DO   30 N=NLOW,NHIGH
      IF(IRPM.EQ.IRPN) MHIGH = N-1
      DO   20 M=MLOW,MHIGH
C
      IF(IUHF.EQ.0)THEN
      MN = IOFFOO(IRPN,IRPMN,5)      + (N-1)*POP(IRPM,1) + M
      ELSE
      IF(IRPMN.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,ISPIN1) + INDEX(N-1) + M
      ELSE
      MN = IOFFOO(IRPN,IRPMN,ISPIN1) + (N-1)*POP(IRPM,ISPIN1) + M
      ENDIF
      ENDIF
C
      IF(NONEQL)THEN
      MNKVAL = IOFF + (K-1) * POP(IRPM,ISPIN1) * POP(IRPN,ISPIN1) +
     1                (N-1) * POP(IRPM,ISPIN1) + M
      ENDIF
C
      IF(MNEQL)THEN
      MNKVAL = IOFF + (K-1)*(POP(IRPN,ISPIN1)*(POP(IRPN,ISPIN1)-1)/2) +
     1                INDEX(N-1) + M
      ENDIF
C
      CALL GETLIST(CORE(I010),MNKVAL,1,1,IRPIJK,ISPIN1 + 1 + 4)
C
      CALL VADD(T3,T3,CORE(I010),LEN,CORE(I000 - 1 + MN))
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
C
  100 CONTINUE
C
      DO  300 IRPMN=1,NIRREP
C
      IF(IRPMN.NE.IRPIK.AND.IRPMN.NE.IRPJK) GOTO 300
      IF(IRPDPD(IRPMN,14).EQ.0)             GOTO 300
C
      IF(IRPMN.EQ.IRPIK)THEN
C
C     -T(ABc,JMn) * W(Mn,Ik)
C     -T(abC,jmN) * W(Nm,Ki)
C
      IF(ISPIN1.EQ.1)THEN
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,ISPIN1) + I
      ELSE
      IK = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,ISPIN2) + K
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 53
      ELSE
      LIST = 13
      ENDIF
      CALL GETLST(CORE(I000),IK,1,2,IRPIK,LIST)
C
      DO  140 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.POP(IRPN,ISPIN2).EQ.0) GOTO 140
      IF(IRPM.EQ.IRPJ.AND.POP(IRPM,ISPIN1).LT.2)         GOTO 140
C
      IF(IRPJ.LE.IRPM) JMN = IJKPOS(IRPJ,IRPM,IRPN,2)
      IF(IRPJ.GT.IRPM) JMN = IJKPOS(IRPM,IRPJ,IRPN,2)
      IOFF = IJKOFF(JMN,IRPIJK,1 + ISPIN1)
C
      JMEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPJ) JMEQL  = .TRUE.
      IF(IRPM.NE.IRPJ) NONEQL = .TRUE.
C
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN2)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN1)
C
      DO  130 N=NLOW,NHIGH
      DO  120 M=MLOW,MHIGH
      SIGN = -1.0D+00
C
      IF(ISPIN1.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,5) + (N-1)*POP(IRPM,ISPIN1) + M
      ELSE
      MN = IOFFOO(IRPM,IRPMN,5) + (M-1)*POP(IRPN,ISPIN2) + N
      ENDIF
C
      IF(IRPM.EQ.IRPJ.AND.M.EQ.J) GOTO 120
C
      IF(NONEQL)THEN
      IF(IRPM.GT.IRPJ) JMNVAL = IOFF +
     1                 (N-1) * POP(IRPM,ISPIN1) * POP(IRPJ,ISPIN1) +
     1                 (M-1) * POP(IRPJ,ISPIN1) + J
      IF(IRPM.LT.IRPJ)THEN
                       JMNVAL = IOFF +
     1                 (N-1) * POP(IRPM,ISPIN1) * POP(IRPJ,ISPIN1) +
     1                 (J-1) * POP(IRPM,ISPIN1) + M
      SIGN =  1.0D+00
      ENDIF
      ENDIF
C
      IF(JMEQL)THEN
      NJM = (POP(IRPM,ISPIN1) * (POP(IRPM,ISPIN1) - 1))/2
      JM  = INDEX(MAX(M,J) - 1) + MIN(M,J)
      JMNVAL = IOFF + (N-1) * NJM + JM
      IF(M.LT.J) SIGN =  1.0D+00
      ENDIF
C
      CALL GETLIST(CORE(I010),JMNVAL,1,1,IRPIJK,1 + ISPIN1 + 4)
C
      CALL VADD(T3,T3,CORE(I010),LEN,CORE(I000 - 1 + MN) * SIGN)
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C
      ENDIF
C
      IF(IRPMN.EQ.IRPJK)THEN
C
C      T(ABc,IMn) * W(Mn,Jk)
C      T(abC,imN) * W(Nm,Kj)
C
      IF(ISPIN1.EQ.1)THEN
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,ISPIN1) + J
      ELSE
      JK = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,ISPIN2) + K
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 53
      ELSE
      LIST = 13
      ENDIF
      CALL GETLST(CORE(I000),JK,1,2,IRPJK,LIST)
C
      DO  240 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.POP(IRPN,ISPIN2).EQ.0) GOTO 240
      IF(IRPM.EQ.IRPI.AND.POP(IRPM,ISPIN1).LT.2)         GOTO 240
C
      IF(IRPI.LE.IRPM) IMN = IJKPOS(IRPI,IRPM,IRPN,2)
      IF(IRPI.GT.IRPM) IMN = IJKPOS(IRPM,IRPI,IRPN,2)
      IOFF = IJKOFF(IMN,IRPIJK,1 + ISPIN1)
C
      IMEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPI) IMEQL  = .TRUE.
      IF(IRPM.NE.IRPI) NONEQL = .TRUE.
C
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN2)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN1)
C
      DO  230 N=NLOW,NHIGH
      DO  220 M=MLOW,MHIGH
      SIGN =  1.0D+00
C
      IF(ISPIN1.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,5) + (N-1)*POP(IRPM,ISPIN1) + M
      ELSE
      MN = IOFFOO(IRPM,IRPMN,5) + (M-1)*POP(IRPN,ISPIN2) + N
      ENDIF
C
      IF(IRPM.EQ.IRPI.AND.M.EQ.I) GOTO 220
C
      IF(NONEQL)THEN
      IF(IRPM.GT.IRPI) IMNVAL = IOFF +
     1                 (N-1) * POP(IRPM,ISPIN1) * POP(IRPI,ISPIN1) +
     1                 (M-1) * POP(IRPI,ISPIN1) + I
      IF(IRPM.LT.IRPI)THEN
                       IMNVAL = IOFF +
     1                 (N-1) * POP(IRPM,ISPIN1) * POP(IRPI,ISPIN1) +
     1                 (I-1) * POP(IRPM,ISPIN1) + M
      SIGN = -1.0D+00
      ENDIF
      ENDIF
C
      IF(IMEQL)THEN
      NIM = (POP(IRPM,ISPIN1) * (POP(IRPM,ISPIN1) - 1))/2
      IM  = INDEX(MAX(M,I) - 1) + MIN(M,I)
      IMNVAL = IOFF + (N-1) * NIM + IM
      IF(M.LT.I) SIGN = -1.0D+00
      ENDIF
C
      CALL GETLIST(CORE(I010),IMNVAL,1,1,IRPIJK,1 + ISPIN1 + 4)
C
      CALL VADD(T3,T3,CORE(I010),LEN,CORE(I000 - 1 + MN) * SIGN)
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
C
      ENDIF
C
  300 CONTINUE
      RETURN
      END
