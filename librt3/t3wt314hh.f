      SUBROUTINE T3WT314HH(T3,CORE,MAXCOR,ISPIN,LEN,IRPI,IRPJ,IRPK,
     1                    IRPIJK,I,J,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER VRT,POP,DIRPRD
      LOGICAL MNKEQL,MKEQL,NKEQL,MNEQL,NONEQL,
     1        MNJEQL,MJEQL,NJEQL,
     1        MNIEQL,MIEQL,NIEQL
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
C     T(abc,mnk) * W(mn,ij) - T(abc,mnj) * W(mn,ik) + T(abc,imn) * W(mn,jk)
C
      INDEX(I) = I*(I-1)/2
C
      IRPIJ = DIRPRD(IRPI,IRPJ)
      IRPIK = DIRPRD(IRPI,IRPK)
      IRPJK = DIRPRD(IRPJ,IRPK)
C
      I000  = 1
      I010  = I000 + MAX(IRPDPD(IRPIJ,ISPIN+2),
     1                   IRPDPD(IRPIK,ISPIN+2),
     1                   IRPDPD(IRPJK,ISPIN+2))
      I020  = I010 + LEN
      NEED  = I020 * IINTFP
      IF(NEED.GT.MAXCOR)THEN
      WRITE(LUOUT,1010) NEED,MAXCOR
 1010 FORMAT(' @T3WT314HH-F, Insufficient memory. Need ',I10,' Got ',
     1                                                   I10)
      CALL INSMEM('T3WT314HH',NEED,MAXCOR)
      ENDIF
C
      DO  300 IRPMN=1,NIRREP
C
      IF(IRPMN.NE.IRPIJ.AND.IRPMN.NE.IRPIK.AND.IRPMN.NE.IRPJK) GOTO 300
      IF(IRPDPD(IRPMN,ISPIN+2).EQ.0)                           GOTO 300
C
      IF(IRPMN.EQ.IRPIJ)THEN
C
C     T(ABC,MNK) * W(MN,IJ)
C
      IF(IRPIJ.EQ.1)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 50 + ISPIN
      ELSE
      LIST = 10 + ISPIN
      ENDIF
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LIST)
C
      DO   40 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(IRPM.GT.IRPN) GO TO 40
      IF(POP(IRPM,ISPIN).EQ.0.OR.POP(IRPN,ISPIN).EQ.0) GOTO 40
      IF(IRPM.EQ.IRPN.AND.POP(IRPM,ISPIN).LT.2)        GOTO 40
      IF(IRPM.EQ.IRPN.AND.IRPM.EQ.IRPK.AND.POP(IRPM,ISPIN).LT.3) GOTO 40
      IF(IRPM.EQ.IRPK.AND.POP(IRPM,ISPIN).LT.2)        GOTO 40
      IF(IRPN.EQ.IRPK.AND.POP(IRPN,ISPIN).LT.2)        GOTO 40
C
      IF(IRPN.LE.IRPK)                  MNK = IJKPOS(IRPM,IRPN,IRPK,1)
      IF(IRPM.GE.IRPK)                  MNK = IJKPOS(IRPK,IRPM,IRPN,1)
      IF(IRPM.LT.IRPK.AND.IRPK.LT.IRPN) MNK = IJKPOS(IRPM,IRPK,IRPN,1)
      IOFF = IJKOFF(MNK,IRPIJK,1 + 3*(ISPIN-1))
C
      MNKEQL = .FALSE.
      MNEQL  = .FALSE.
      MKEQL  = .FALSE.
      NKEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPK.AND.IRPN.EQ.IRPK) MNKEQL = .TRUE.
      IF(IRPM.EQ.IRPN.AND.IRPM.NE.IRPK) MNEQL  = .TRUE.
      IF(IRPM.EQ.IRPK.AND.IRPN.NE.IRPK) MKEQL  = .TRUE.
      IF(IRPN.EQ.IRPK.AND.IRPM.NE.IRPK) NKEQL  = .TRUE.
      IF(IRPM.NE.IRPK.AND.IRPN.NE.IRPK.AND.IRPM.NE.IRPN) NONEQL = .TRUE.
C
      IF(IRPM.EQ.IRPN)THEN
      NLOW  = 2
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      ELSE
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN)
      ENDIF
C
      DO   30 N=NLOW,NHIGH
      IF(IRPM.EQ.IRPN) MHIGH = N-1
      DO   20 M=MLOW,MHIGH
      SIGN = 1.0D+00
C
      IF(IRPMN.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + INDEX(N-1) + M
      ELSE
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + (N-1)*POP(IRPM,ISPIN) + M
      ENDIF
C
      IF((IRPM.EQ.IRPK.AND.M.EQ.K).OR.(IRPN.EQ.IRPK.AND.N.EQ.K)) GOTO 20
C
C     write(6,*) ' n,m,k ',n,m,k
C
      IF(NONEQL)THEN
      IF(IRPK.GT.IRPN) MNKVAL = IOFF + 
     1                 (K-1) * POP(IRPM,ISPIN) * POP(IRPN,ISPIN) +
     1                 (N-1) * POP(IRPM,ISPIN) + M
      IF(IRPK.LT.IRPM) MNKVAL = IOFF + 
     1                 (N-1) * POP(IRPK,ISPIN) * POP(IRPM,ISPIN) +
     1                 (M-1) * POP(IRPK,ISPIN) + K
      IF(IRPK.LT.IRPN.AND.IRPK.GT.IRPM)THEN
      SIGN = -1.0D+00
      MNKVAL = IOFF +  (N-1) * POP(IRPM,ISPIN) * POP(IRPK,ISPIN) +
     1                 (K-1) * POP(IRPM,ISPIN) + M
      ENDIF
      ENDIF
C
      IF(MNEQL)THEN
      IF(IRPK.GT.IRPN) MNKVAL = IOFF +
     1                 (K-1)*(POP(IRPN,ISPIN)*(POP(IRPN,ISPIN)-1)/2) +
     1                 INDEX(N-1) + M
      IF(IRPK.LT.IRPN) MNKVAL = IOFF +
     1                 (INDEX(N-1) + M - 1) * POP(IRPK,ISPIN) + K
      ENDIF
C
      IF(MKEQL)THEN
      NMK = (POP(IRPM,ISPIN) * (POP(IRPM,ISPIN) - 1))/2
      MK  = INDEX(MAX(M,K) - 1) + MIN(M,K)
      MNKVAL = IOFF + (N-1) * NMK + MK
      IF(M.LT.K) SIGN = -1.0D+00
      ENDIF
C
      IF(NKEQL)THEN
      NK  = INDEX(MAX(N,K) - 1) + MIN(N,K)
      MNKVAL = IOFF + (NK-1) * POP(IRPM,ISPIN) + M
      IF(N.GT.K) SIGN = -1.0D+00
      ENDIF
C
      IF(MNKEQL)THEN
      MAXVAL = MAX(M,N,K)
      MINVAL = MIN(M,N,K)
      IF(K.GT.N) MIDVAL = N
      IF(K.LT.M) MIDVAL = M
      IF(M.LT.K.AND.N.GT.K)THEN
      MIDVAL = K
      SIGN = -1.0D+00
      ENDIF
      MNKVAL = IOFF + ((MAXVAL-1)*(MAXVAL-2)*(MAXVAL-3))/6 + 
     1         INDEX(MIDVAL-1) + MINVAL
C     write(6,*) ' mnkval ',mnkval,ioff,ispin,maxval,minval,midval
      ENDIF
C
      CALL GETLIST(CORE(I010),MNKVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)
C
      CALL VADD(T3,T3,CORE(I010),LEN,CORE(I000 - 1 + MN) * SIGN)
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
C
      ENDIF
C
      IF(IRPMN.EQ.IRPIK)THEN
C
C     -T(ABC,MNJ) * W(MN,IK)
C
      IF(IRPIK.EQ.1)THEN
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 50 + ISPIN
      ELSE
      LIST = 10 + ISPIN
      ENDIF
      CALL GETLST(CORE(I000),IK,1,2,IRPIK,LIST)
C
      DO  140 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(IRPM.GT.IRPN) GO TO 140
      IF(POP(IRPM,ISPIN).EQ.0.OR.POP(IRPN,ISPIN).EQ.0) GOTO 140
      IF(IRPM.EQ.IRPN.AND.POP(IRPM,ISPIN).LT.2)        GOTO 140
      IF(IRPM.EQ.IRPN.AND.IRPM.EQ.IRPJ.AND.POP(IRPM,ISPIN).LT.3)
     1                                                 GOTO 140
      IF(IRPM.EQ.IRPJ.AND.POP(IRPM,ISPIN).LT.2)        GOTO 140
      IF(IRPN.EQ.IRPJ.AND.POP(IRPN,ISPIN).LT.2)        GOTO 140
C
      IF(IRPN.LE.IRPJ)                  MNJ = IJKPOS(IRPM,IRPN,IRPJ,1)
      IF(IRPM.GE.IRPJ)                  MNJ = IJKPOS(IRPJ,IRPM,IRPN,1)
      IF(IRPM.LT.IRPJ.AND.IRPJ.LT.IRPN) MNJ = IJKPOS(IRPM,IRPJ,IRPN,1)
      IOFF = IJKOFF(MNJ,IRPIJK,1 + 3*(ISPIN-1))
C
      MNJEQL = .FALSE.
      MNEQL  = .FALSE.
      MJEQL  = .FALSE.
      NJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPJ.AND.IRPN.EQ.IRPJ) MNJEQL = .TRUE.
      IF(IRPM.EQ.IRPN.AND.IRPM.NE.IRPJ) MNEQL  = .TRUE.
      IF(IRPM.EQ.IRPJ.AND.IRPN.NE.IRPJ) MJEQL  = .TRUE.
      IF(IRPN.EQ.IRPJ.AND.IRPM.NE.IRPJ) NJEQL  = .TRUE.
      IF(IRPM.NE.IRPJ.AND.IRPN.NE.IRPJ.AND.IRPM.NE.IRPN) NONEQL = .TRUE.
C
      IF(IRPM.EQ.IRPN)THEN
      NLOW  = 2
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      ELSE
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN)
      ENDIF
C
      DO  130 N=NLOW,NHIGH
      IF(IRPM.EQ.IRPN) MHIGH = N-1
      DO  120 M=MLOW,MHIGH
      SIGN = -1.0D+00
C
      IF(IRPMN.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + INDEX(N-1) + M
      ELSE
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + (N-1)*POP(IRPM,ISPIN) + M
      ENDIF
C
      IF((IRPM.EQ.IRPJ.AND.M.EQ.J).OR.(IRPN.EQ.IRPJ.AND.N.EQ.J))
     1                                                  GOTO 120
C
      IF(NONEQL)THEN
      IF(IRPJ.GT.IRPN) MNJVAL = IOFF + 
     1                 (J-1) * POP(IRPM,ISPIN) * POP(IRPN,ISPIN) +
     1                 (N-1) * POP(IRPM,ISPIN) + M
      IF(IRPJ.LT.IRPM) MNJVAL = IOFF + 
     1                 (N-1) * POP(IRPJ,ISPIN) * POP(IRPM,ISPIN) +
     1                 (M-1) * POP(IRPJ,ISPIN) + J
      IF(IRPJ.LT.IRPN.AND.IRPJ.GT.IRPM)THEN
      SIGN =  1.0D+00
      MNJVAL = IOFF +  (N-1) * POP(IRPM,ISPIN) * POP(IRPJ,ISPIN) +
     1                 (J-1) * POP(IRPM,ISPIN) + M
      ENDIF
      ENDIF
C
      IF(MNEQL)THEN
      IF(IRPJ.GT.IRPN) MNJVAL = IOFF +
     1                 (J-1)*(POP(IRPN,ISPIN)*(POP(IRPN,ISPIN)-1)/2) +
     1                 INDEX(N-1) + M
      IF(IRPJ.LT.IRPN) MNJVAL = IOFF +
     1                 (INDEX(N-1) + M - 1) * POP(IRPJ,ISPIN) + J
      ENDIF
C
      IF(MJEQL)THEN
      NMJ = (POP(IRPM,ISPIN) * (POP(IRPM,ISPIN) - 1))/2
      MJ  = INDEX(MAX(M,J) - 1) + MIN(M,J)
      MNJVAL = IOFF + (N-1) * NMJ + MJ
      IF(M.LT.J) SIGN =  1.0D+00
      ENDIF
C
      IF(NJEQL)THEN
      NJ  = INDEX(MAX(N,J) - 1) + MIN(N,J)
      MNJVAL = IOFF + (NJ-1) * POP(IRPM,ISPIN) + M
      IF(N.GT.J) SIGN =  1.0D+00
      ENDIF
C
      IF(MNJEQL)THEN
      MAXVAL = MAX(M,N,J)
      MINVAL = MIN(M,N,J)
      IF(J.GT.N) MIDVAL = N
      IF(J.LT.M) MIDVAL = M
      IF(M.LT.J.AND.N.GT.J)THEN
      MIDVAL = J
      SIGN =  1.0D+00
      ENDIF
      MNJVAL = IOFF + ((MAXVAL-1)*(MAXVAL-2)*(MAXVAL-3))/6 + 
     1         INDEX(MIDVAL-1) + MINVAL
      ENDIF
C
      CALL GETLIST(CORE(I010),MNJVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)
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
C     T(ABC,MNI) * W(MN,JK)
C
      IF(IRPJK.EQ.1)THEN
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
C     Have to handle here the question of CCSDT or CCSDT-4 and the
C     possible need for transposition. For now just read an integral.
      IF(ICLLVL.EQ.18)THEN
      LIST = 50 + ISPIN
      ELSE
      LIST = 10 + ISPIN
      ENDIF
      CALL GETLST(CORE(I000),JK,1,2,IRPJK,LIST)
C
      DO  240 IRPN=1,NIRREP
      IRPM = DIRPRD(IRPN,IRPMN)
C
      IF(IRPM.GT.IRPN) GO TO 240
      IF(POP(IRPM,ISPIN).EQ.0.OR.POP(IRPN,ISPIN).EQ.0) GOTO 240
      IF(IRPM.EQ.IRPN.AND.POP(IRPM,ISPIN).LT.2)        GOTO 240
C
      IF(IRPM.EQ.IRPN.AND.IRPM.EQ.IRPI.AND.POP(IRPM,ISPIN).LT.3)
     1                                                 GOTO 240
      IF(IRPM.EQ.IRPI.AND.POP(IRPM,ISPIN).LT.2)        GOTO 240
      IF(IRPN.EQ.IRPI.AND.POP(IRPN,ISPIN).LT.2)        GOTO 240
C
      IF(IRPN.LE.IRPI)                  MNI = IJKPOS(IRPM,IRPN,IRPI,1)
      IF(IRPM.GE.IRPI)                  MNI = IJKPOS(IRPI,IRPM,IRPN,1)
      IF(IRPM.LT.IRPI.AND.IRPI.LT.IRPN) MNI = IJKPOS(IRPM,IRPI,IRPN,1)
      IOFF = IJKOFF(MNI,IRPIJK,1 + 3*(ISPIN-1))
C
      MNIEQL = .FALSE.
      MNEQL  = .FALSE.
      MIEQL  = .FALSE.
      NIEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPM.EQ.IRPI.AND.IRPN.EQ.IRPI) MNIEQL = .TRUE.
      IF(IRPM.EQ.IRPN.AND.IRPM.NE.IRPI) MNEQL  = .TRUE.
      IF(IRPM.EQ.IRPI.AND.IRPN.NE.IRPI) MIEQL  = .TRUE.
      IF(IRPN.EQ.IRPI.AND.IRPM.NE.IRPI) NIEQL  = .TRUE.
      IF(IRPM.NE.IRPI.AND.IRPN.NE.IRPI.AND.IRPM.NE.IRPN) NONEQL = .TRUE.
C
      IF(IRPM.EQ.IRPN)THEN
      NLOW  = 2
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      ELSE
      NLOW  = 1
      NHIGH = POP(IRPN,ISPIN)
      MLOW  = 1
      MHIGH = POP(IRPM,ISPIN)
      ENDIF
C
      DO  230 N=NLOW,NHIGH
      IF(IRPM.EQ.IRPN) MHIGH = N-1
      DO  220 M=MLOW,MHIGH
      SIGN = 1.0D+00
C
      IF(IRPMN.EQ.1)THEN
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + INDEX(N-1) + M
      ELSE
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + (N-1)*POP(IRPM,ISPIN) + M
      ENDIF
C
      IF((IRPM.EQ.IRPI.AND.M.EQ.I).OR.(IRPN.EQ.IRPI.AND.N.EQ.I))
     1                                                  GOTO 220
C
      IF(NONEQL)THEN
      IF(IRPI.GT.IRPN) MNIVAL = IOFF + 
     1                 (I-1) * POP(IRPM,ISPIN) * POP(IRPN,ISPIN) +
     1                 (N-1) * POP(IRPM,ISPIN) + M
      IF(IRPI.LT.IRPM) MNIVAL = IOFF + 
     1                 (N-1) * POP(IRPI,ISPIN) * POP(IRPM,ISPIN) +
     1                 (M-1) * POP(IRPI,ISPIN) + I
      IF(IRPI.LT.IRPN.AND.IRPI.GT.IRPM)THEN
      SIGN = -1.0D+00
      MNIVAL = IOFF +  (N-1) * POP(IRPM,ISPIN) * POP(IRPI,ISPIN) +
     1                 (I-1) * POP(IRPM,ISPIN) + M
      ENDIF
      ENDIF
C
      IF(MNEQL)THEN
      IF(IRPI.GT.IRPN) MNIVAL = IOFF +
     1                 (I-1)*(POP(IRPN,ISPIN)*(POP(IRPN,ISPIN)-1)/2) +
     1                 INDEX(N-1) + M
      IF(IRPI.LT.IRPN) MNIVAL = IOFF +
     1                 (INDEX(N-1) + M - 1) * POP(IRPI,ISPIN) + I
      ENDIF
C
      IF(MIEQL)THEN
      NMI = (POP(IRPM,ISPIN) * (POP(IRPM,ISPIN) - 1))/2
      MI  = INDEX(MAX(M,I) - 1) + MIN(M,I)
      MNIVAL = IOFF + (N-1) * NMI + MI
      IF(M.LT.I) SIGN = -1.0D+00
      ENDIF
C
      IF(NIEQL)THEN
      NI  = INDEX(MAX(N,I) - 1) + MIN(N,I)
      MNIVAL = IOFF + (NI-1) * POP(IRPM,ISPIN) + M
      IF(N.GT.I) SIGN = -1.0D+00
      ENDIF
C
      IF(MNIEQL)THEN
      MAXVAL = MAX(M,N,I)
      MINVAL = MIN(M,N,I)
      IF(I.GT.N) MIDVAL = N
      IF(I.LT.M) MIDVAL = M
      IF(M.LT.I.AND.N.GT.I)THEN
      MIDVAL = I
      SIGN = -1.0D+00
      ENDIF
      MNIVAL = IOFF + ((MAXVAL-1)*(MAXVAL-2)*(MAXVAL-3))/6 + 
     1         INDEX(MIDVAL-1) + MINVAL
      ENDIF
C
      CALL GETLIST(CORE(I010),MNIVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)
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
