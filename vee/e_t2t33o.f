      SUBROUTINE E_T2T33O(T3,T3EXP,CORE,IADT3,IADT3EXP,
     1                    IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &                    I,J,K,IRREPX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION T3(1),T3EXP(1),CORE(1)
      INTEGER DIRPRD,POP,VRT
      DIMENSION IADT3EXP(8),IADT3(8)
      DIMENSION IADW(8),LENW(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
C
      INDEX(I) = I*(I-1)/2
C
C     Z2 lists : LIST2I1, LIST2I2, LIST2I3
C     W  lists : LWIC21, LWIC22, LWIC23, LWIC24
C     T3       : (a<b;C)
C     T3EXP    : (bC ,a)
C
C     Expected list ordering :
C
C              a<b ; i<j (LIST2I1, LIST2I2)
C              b,A ; I,j (LIST2I3)
C              k,a ; i<j (LWIC21, LWIC22)
C              A,k ; I,j (LWIC23)
C              K,a ; I,j (LWIC24)
C
      I000 = 1
C
C-----------------------------------------------------------------------
C     Z(ab,im) = - W(Kj,Cm) * T(abC,ijK)   [ T(a<b,C) * W(C) ]
C-----------------------------------------------------------------------
C
      DO   10 IRPM=1,NIRREP
      IRPC = DIRPRD(IRPM,IRPJK)
      LENW(IRPM) = VRT(IRPC,1)*POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADW(IRPM) = I000
      ELSE
      IADW(IRPM) = IADW(IRPM-1) + LENW(IRPM-1)
      ENDIF
   10 CONTINUE
C
      KJ = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,1) + K
      CALL GETLST(CORE(I000),KJ,1,2,IRPJK,LWIC23)
C
C-----------------------------------------------------------------------
C     Address for Z2.
C-----------------------------------------------------------------------
C
      I010 = IADW(NIRREP) + LENW(NIRREP)
C
C
      DO   50 IRPM=1,NIRREP
      IF(POP(IRPM,2).EQ.0) GOTO 50
C
      IRPIM = DIRPRD(IRPI,IRPM)
      IRPC  = DIRPRD(IRPM,IRPJK)
      IRPAB = DIRPRD(IRREPX,IRPIM)
C
      IF(IRPM.GT.IRPI)THEN
      DO   20 M=1,POP(IRPM,2)
      IM = IOFFOO(IRPM,IRPIM,2) + (M-1)*POP(IRPI,2) + I
      CALL GETLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,1)
      CALL PUTLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
   20 CONTINUE
      ENDIF
C
      IF(IRPM.LT.IRPI)THEN
      DO   30 M=1,POP(IRPM,2)
      IM = IOFFOO(IRPI,IRPIM,2) + (I-1)*POP(IRPM,2) + M
      CALL GETLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,0)
      CALL PUTLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
   30 CONTINUE
      ENDIF
C
      IF(IRPM.EQ.IRPI)THEN
C
      DO   40 M=1,POP(IRPM,2)
      IF(M.GT.I)THEN
      IM = IOFFOO(IRPM,IRPIM,2) + INDEX(M-1) + I
      CALL GETLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,1)
      CALL PUTLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      ENDIF
C
      IF(M.LT.I)THEN
      IM = IOFFOO(IRPI,IRPIM,2) + INDEX(I-1) + M
      CALL GETLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,0)
      CALL PUTLST(CORE(I010),IM,1,1,IRPIM,LIST2I2)
      ENDIF
   40 CONTINUE
C
      ENDIF
C
   50 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(ab,jm) =   W(Ki,Cm) * T(abC,ijK)      [ T(a<b,C) * W(C) ]
C-----------------------------------------------------------------------
C
      DO  110 IRPM=1,NIRREP
      IRPC = DIRPRD(IRPM,IRPIK)
      LENW(IRPM) = VRT(IRPC,1)*POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADW(IRPM) = I000
      ELSE
      IADW(IRPM) = IADW(IRPM-1) + LENW(IRPM-1)
      ENDIF
  110 CONTINUE
C
      KI = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,1) + K
      CALL GETLST(CORE(I000),KI,1,2,IRPIK,LWIC23)
C
C-----------------------------------------------------------------------
C     Address for Z2.
C-----------------------------------------------------------------------
C
      I010 = IADW(NIRREP) + LENW(NIRREP)
C
C
      DO  150 IRPM=1,NIRREP
      IF(POP(IRPM,2).EQ.0) GOTO 150
C
      IRPJM = DIRPRD(IRPJ,IRPM)
      IRPC  = DIRPRD(IRPM,IRPIK)
      IRPAB = DIRPRD(IRREPX,IRPJM)
C
      IF(IRPM.GT.IRPJ)THEN
      DO  120 M=1,POP(IRPM,2)
      JM = IOFFOO(IRPM,IRPJM,2) + (M-1)*POP(IRPJ,2) + J
      CALL GETLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,0)
      CALL PUTLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
  120 CONTINUE
      ENDIF
C
      IF(IRPM.LT.IRPJ)THEN
      DO  130 M=1,POP(IRPM,2)
      JM = IOFFOO(IRPJ,IRPJM,2) + (J-1)*POP(IRPM,2) + M
      CALL GETLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,1)
      CALL PUTLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
  130 CONTINUE
      ENDIF
C
      IF(IRPM.EQ.IRPJ)THEN
C
      DO  140 M=1,POP(IRPM,2)
      IF(M.GT.J)THEN
      JM = IOFFOO(IRPM,IRPJM,2) + INDEX(M-1) + J
      CALL GETLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,0)
      CALL PUTLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      ENDIF
C
      IF(M.LT.J)THEN
      JM = IOFFOO(IRPJ,IRPJM,2) + INDEX(J-1) + M
      CALL GETLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPC,1)),
     1            CORE(I010),IRPDPD(IRPAB,2),VRT(IRPC,1),0,1)
      CALL PUTLST(CORE(I010),JM,1,1,IRPJM,LIST2I2)
      ENDIF
  140 CONTINUE
C
      ENDIF
C
  150 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(Cb,Km) =   W(ij,ma) T(abC,ijK)    [ T3EXP(bC,a) * W(a) ]
C-----------------------------------------------------------------------
C
      DO  210 IRPM=1,NIRREP
      IRPA = DIRPRD(IRPM,IRPIJ)
      LENW(IRPM) = VRT(IRPA,2)*POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADW(IRPM) = I000
      ELSE
      IADW(IRPM) = IADW(IRPM-1) + LENW(IRPM-1)
      ENDIF
  210 CONTINUE
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,2) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,2) + (J-1)*POP(IRPI,2) + I
      ENDIF
C
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LWIC22)
C
C     Convert W(m,a) to W(a,m).
C
      I010 = I000 + IRPDPD(IRPIJ,10)
      I020 = I010 + IRPDPD(IRPIJ,10)
      I030 = I020 + IRPDPD(IRPIJ,10)
      CALL SYMTR3(IRPIJ,POP(1,2),VRT(1,2),IRPDPD(IRPIJ,10),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
C     Set address for Z2.
C
      I010 = IADW(NIRREP) + LENW(NIRREP)
C
      DO  250 IRPM=1,NIRREP
      IF(POP(IRPM,2).EQ.0) GOTO 250
C
      IRPMK = DIRPRD(IRPM,IRPK)
      IRPA =  DIRPRD(IRPM,IRPIJ)
      IRPBC = DIRPRD(IRREPX,IRPMK)
      DO  220 M=1,POP(IRPM,2)
      MK = IOFFOO(IRPM,IRPMK,5) + (M-1)*POP(IRPK,1) + K
      CALL GETLST(CORE(I010),MK,1,1,IRPMK,LIST2I3)
      CALL MATVEC(T3EXP(IADT3EXP(IRPA)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPA,2)),
     1            CORE(I010),
     1            IRPDPD(IRPBC,13),VRT(IRPA,2),0,0)
      CALL PUTLST(CORE(I010),MK,1,1,IRPMK,LIST2I3)
  220 CONTINUE
  250 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(Cb,Mi) =    W(Kj,Ma) T(abC,ijK)    [ T3EXP(bC,a) W(a) ]
C-----------------------------------------------------------------------
C
      DO  310 IRPM=1,NIRREP
      IRPA = DIRPRD(IRPM,IRPJK)
      LENW(IRPM) = VRT(IRPA,2)*POP(IRPM,1)
      IF(IRPM.EQ.1)THEN
      IADW(IRPM) = I000
      ELSE
      IADW(IRPM) = IADW(IRPM-1) + LENW(IRPM-1)
      ENDIF
  310 CONTINUE
C
      JK = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,1) + K
      CALL GETLST(CORE(I000),JK,1,1,IRPJK,LWIC24)
C
C     Convert W(M,a) to W(a,M).
C
      I010 = I000 + IRPDPD(IRPJK,12)
      I020 = I010 + IRPDPD(IRPJK,12)
      I030 = I020 + IRPDPD(IRPJK,12)
      CALL SYMTR3(IRPJK,POP(1,1),VRT(1,2),IRPDPD(IRPJK,12),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
      I010 = IADW(NIRREP) + LENW(NIRREP)
C
      DO  350 IRPM=1,NIRREP
      IF(POP(IRPM,1).EQ.0) GOTO 350
C
      IRPIM = DIRPRD(IRPI,IRPM)
      IRPA  = DIRPRD(IRPM,IRPJK)
      IRPBC = DIRPRD(IRREPX,IRPIM)
      DO  320 M=1,POP(IRPM,1)
      IM = IOFFOO(IRPI,IRPIM,5) + (I-1)*POP(IRPM,1) + M
      CALL GETLST(CORE(I010),IM,1,1,IRPIM,LIST2I3)
C
      CALL MATVEC(T3EXP(IADT3EXP(IRPA)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPA,2)),
     1            CORE(I010),
     1            IRPDPD(IRPBC,13),VRT(IRPA,2),0,0)
      CALL PUTLST(CORE(I010),IM,1,1,IRPIM,LIST2I3)
  320 CONTINUE
  350 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(Cb,Mj) =  - W(Ki,Ma) T(abC,ijK)    [ T3EXP(bC,a) W(a) ]
C-----------------------------------------------------------------------
C
      DO  410 IRPM=1,NIRREP
      IRPA = DIRPRD(IRPM,IRPIK)
      LENW(IRPM) = VRT(IRPA,2)*POP(IRPM,1)
      IF(IRPM.EQ.1)THEN
      IADW(IRPM) = I000
      ELSE
      IADW(IRPM) = IADW(IRPM-1) + LENW(IRPM-1)
      ENDIF
  410 CONTINUE
C
      IK = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,1) + K
      CALL GETLST(CORE(I000),IK,1,1,IRPIK,LWIC24)
C
C     Convert W(M,a) to W(a,M).
C
      I010 = I000 + IRPDPD(IRPIK,12)
      I020 = I010 + IRPDPD(IRPIK,12)
      I030 = I020 + IRPDPD(IRPIK,12)
      CALL SYMTR3(IRPIK,POP(1,1),VRT(1,2),IRPDPD(IRPIK,12),1,CORE(I000),
     &            CORE(I010),CORE(I020),CORE(I030))
C
      I010 = IADW(NIRREP) + LENW(NIRREP)
C
      DO  450 IRPM=1,NIRREP
      IF(POP(IRPM,1).EQ.0) GOTO 450
C
      IRPJM = DIRPRD(IRPJ,IRPM)
      IRPA  = DIRPRD(IRPM,IRPIK)
      IRPBC = DIRPRD(IRREPX,IRPJM)
      DO  420 M=1,POP(IRPM,1)
      JM = IOFFOO(IRPJ,IRPJM,5) + (J-1)*POP(IRPM,1) + M
      CALL GETLST(CORE(I010),JM,1,1,IRPJM,LIST2I3)
C
      CALL MATVEC(T3EXP(IADT3EXP(IRPA)),
     1            CORE(IADW(IRPM) + (M-1)*VRT(IRPA,2)),
     1            CORE(I010),
     1            IRPDPD(IRPBC,13),VRT(IRPA,2),0,1)
      CALL PUTLST(CORE(I010),JM,1,1,IRPJM,LIST2I3)
  420 CONTINUE
  450 CONTINUE
      RETURN
      END
