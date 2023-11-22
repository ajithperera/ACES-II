      SUBROUTINE Y12B(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZ,DISSIZA,DISSIZB
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1),DISSIZ(8),DISSIZA(8),DISSIZB(8),
     1          IOFFY4(8),IOFFY4R(8),IOFFW(8),IOFFWR(8)
      DIMENSION IOFFY4A(8),IOFFY4AR(8),IOFFY4B(8),IOFFY4BR(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I * (I-1)/2
C
      WRITE(6,1000)
 1000 FORMAT(' @Y12B-I, Y4 * Y4 contributions. ')
C
C     AAAA/BBBB
C
      IF(IUHF.GT.0)THEN
C
      DO  400 ISPIN=1,2
C
      LEN = 0
      DO   20 IRREP=1,NIRREP
      DISSIZ(IRREP) = 0
      DO   10 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZ(IRREP)  =  DISSIZ(IRREP) +
     1                  POP(JRREP,ISPIN) * POP(KRREP,ISPIN)
   10 CONTINUE
      LEN = LEN + DISSIZ(IRREP) * DISSIZ(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY4(IRREP) = 0
      ELSE
      IOFFY4(IRREP) = IOFFY4(IRREP-1) + DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      ENDIF
   20 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NOCA * NOCB + NOCA * NOCB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
 1020 FORMAT(' @Y12B-I, Insufficient memory. Need ',I15,' . Got ',I15)
      STOP 'Y12B'
      ENDIF
C
C     Read ijkn intermediates and expand each symmetry block in turn to
C     (i,j,k,n).
C
      DO   30 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY4(IRREP)),
     1             1,IRPDPD(IRREP,ISPIN+2),2,IRREP,6+ISPIN)
C
C     (i<j,k<n) ---> (i,j,k<n)
      CALL SYMEXP2(IRREP,POP(1,ISPIN),DISSIZ(IRREP),
     1             IRPDPD(IRREP,2+ISPIN),IRPDPD(IRREP,2+ISPIN),
     1             CORE(I010+IOFFY4(IRREP)),CORE(I010+IOFFY4(IRREP)))
C     (i,j,k<n) ---> (i,j,k,n)
      CALL  SYMEXP(IRREP,POP(1,ISPIN),DISSIZ(IRREP),
     1             CORE(I010+IOFFY4(IRREP)))
C
   30 CONTINUE
C
C     Now reorder (i,j,k,n) to (i,k,j,n). This is also effectively
C     reordering (i,m,k,l) to (i,k,m,l).
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),
     1            CORE(I020),1,'1324')
C
C     Now perform massive matrix multiply (ik,ml) * (ik,jn) = (ml,jn)
C
C     Y4 intermediates are now at I000. Target can be at I010. Allocation
C     is same.
C
      DO   40 IRPIK=1,NIRREP
C
      CALL XGEMM('T','N',DISSIZ(IRPIK),DISSIZ(IRPIK),DISSIZ(IRPIK),
     1            1.0D+00,
     1           CORE(I000+IOFFY4(IRPIK)),DISSIZ(IRPIK),
     1           CORE(I000+IOFFY4(IRPIK)),DISSIZ(IRPIK),
     1            0.0D+00,
     1           CORE(I010+IOFFY4(IRPIK)),DISSIZ(IRPIK))
C
   40 CONTINUE
C
C     Now loop over symmetry blocks of Y12 and insert appropriate pieces.
C
      DO  190 IRPJL=1,NIRREP
C
      I030 = I020 + IRPDPD(IRPJL,2+ISPIN) * IRPDPD(IRPJL,2+ISPIN)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPJL,2+ISPIN),2,IRPJL,9+ISPIN)
C
      IRPMN = IRPJL
C
      DO  180 IRPL =1,NIRREP
C
      IRPJ = DIRPRD(IRPL,IRPJL)
      IF(IRPJ.GT.IRPL) GOTO 180
C
      DO  170 IRPN =1,NIRREP
C
      IRPM = DIRPRD(IRPN,IRPMN)
      IF(IRPM.GT.IRPN) GOTO 170
C
      IRPML = DIRPRD(IRPM,IRPL)
      IRPJN = DIRPRD(IRPJ,IRPN)
      IRPNL = DIRPRD(IRPN,IRPL)
      IRPJM = DIRPRD(IRPJ,IRPM)
C
      IF(IRPJL.EQ.1)THEN
C
      IF(POP(IRPL,ISPIN).GE.2.AND.POP(IRPN,ISPIN).GE.2)THEN
C
      DO   80 L=2,POP(IRPL,ISPIN)
      DO   70 J=1,L-1
      DO   60 N=2,POP(IRPN,ISPIN)
      DO   50 M=1,N-1
C
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + INDEX(L-1) + J
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + INDEX(N-1) + M
C
      MNJL = (JL-1) * IRPDPD(IRPMN,ISPIN+2) + MN
C
      ML = IOFFOO(IRPL,IRPML,2+ISPIN) + (L-1)*POP(IRPM,ISPIN) + M
      JN = IOFFOO(IRPN,IRPJN,2+ISPIN) + (N-1)*POP(IRPJ,ISPIN) + J
      NL = IOFFOO(IRPL,IRPNL,2+ISPIN) + (L-1)*POP(IRPN,ISPIN) + N
      JM = IOFFOO(IRPM,IRPJM,2+ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
C
      MLJN = I010 + IOFFY4(IRPJN) + (JN-1)*DISSIZ(IRPML) + ML - 1
      NLJM = I010 + IOFFY4(IRPJM) + (JM-1)*DISSIZ(IRPNL) + NL - 1
C
      CORE(I020 + MNJL - 1) = CORE(I020 + MNJL - 1)
     1                      + CORE(MLJN) - CORE(NLJM)
C
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
      ENDIF
C
      ELSE
C
      IF(POP(IRPL,ISPIN).GE.1.AND.POP(IRPN,ISPIN).GE.1.AND.
     1   POP(IRPJ,ISPIN).GE.1.AND.POP(IRPM,ISPIN).GE.1)THEN
C
      DO  120 L=1,POP(IRPL,ISPIN)
      DO  110 J=1,POP(IRPJ,ISPIN)
      DO  100 N=1,POP(IRPN,ISPIN)
      DO   90 M=1,POP(IRPM,ISPIN)
C
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + (L-1)*VRT(IRPJ,ISPIN) + J
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + (N-1)*VRT(IRPM,ISPIN) + M
C
      MNJL = (JL-1) * IRPDPD(IRPMN,ISPIN+2) + MN
C
      ML = IOFFOO(IRPL,IRPML,2+ISPIN) + (L-1)*POP(IRPM,ISPIN) + M
      JN = IOFFOO(IRPN,IRPJN,2+ISPIN) + (N-1)*POP(IRPJ,ISPIN) + J
      NL = IOFFOO(IRPL,IRPNL,2+ISPIN) + (L-1)*POP(IRPN,ISPIN) + N
      JM = IOFFOO(IRPM,IRPJM,2+ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
C
      MLJN = I010 + IOFFY4(IRPJN) + (JN-1)*DISSIZ(IRPML) + ML - 1
      NLJM = I010 + IOFFY4(IRPJM) + (JM-1)*DISSIZ(IRPNL) + NL - 1
C
      CORE(I020 + MNJL - 1) = CORE(I020 + MNJL - 1)
     1                      + CORE(MLJN) - CORE(NLJM)
C
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
C
      ENDIF
C
      ENDIF
C
  170 CONTINUE
  180 CONTINUE
C
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPJL,2+ISPIN),2,IRPJL,9+ISPIN)
C
  190 CONTINUE
C
C     alpha-beta contribution.
C
      IF(ISPIN.EQ.1)THEN
      ISPIN1 = 1
      ISPIN2 = 2
      ELSE
      ISPIN1 = 2
      ISPIN2 = 1
      ENDIF
C
      LEN  = 0
      LEN2 = 0
      LENW = 0
      DO  220 IRREP=1,NIRREP
      DISSIZ(IRREP) = 0
      DISSIZA(IRREP) = 0
      DISSIZB(IRREP) = 0
      DO  210 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZ(IRREP)   =  DISSIZ(IRREP) +
     1                   POP(JRREP,ISPIN1) * POP(KRREP,ISPIN2)
      DISSIZA(IRREP)  =  DISSIZA(IRREP) +
     1                   POP(JRREP,ISPIN1) * POP(KRREP,ISPIN1)
      DISSIZB(IRREP)  =  DISSIZB(IRREP) +
     1                   POP(JRREP,ISPIN2) * POP(KRREP,ISPIN2)
  210 CONTINUE
      LEN  = LEN  + DISSIZA(IRREP) * DISSIZB(IRREP)
      LEN2 = LEN2 + DISSIZ(IRREP)  * DISSIZ(IRREP)
      LENW = LENW + DISSIZA(IRREP) * DISSIZA(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY4(IRREP)  = 0
      IOFFY4R(IRREP) = 0
      IOFFW(IRREP)   = 0
      ELSE
      IOFFY4(IRREP)  = IOFFY4(IRREP-1) + 
     1                 DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      IOFFY4R(IRREP) = IOFFY4R(IRREP-1) + 
     1                 DISSIZA(IRREP-1)*DISSIZB(IRREP-1)
      IOFFW(IRREP)   = IOFFW(IRREP-1) +
     1                 DISSIZA(IRREP-1)*DISSIZA(IRREP-1)
      ENDIF
  220 CONTINUE
C
      IF(LEN.NE.LEN2)THEN
      WRITE(6,1030)
 1030 FORMAT(' @Y12B-I, A screwup... No way anything can be right. ')
      ENDIF
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NOCA * NOCB + NOCA * NOCB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
C     Read ijkn intermediates.
C
      DO  230 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY4(IRREP)),
     1             1,IRPDPD(IRREP,14),2,IRREP,9)
  230 CONTINUE
C
C     ISPIN = 1 : (J,i,N,k), (M,i,L,k) ---> (J,N,i,k), (M,L,i,k)
C     ISPIN = 2 : (I,j,K,n), (I,m,K,l) ---> (I,K,j,n), (I,K,m,l)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I020),1,'1324')
C
C     Now perform massive matrix multiply :
C
C     ISPIN = 1 : (ML,ik) * (JN,ik) = (ML,JN)
C     ISPIN = 2 : (IK,ml) * (IK,jn) = (ml,jn)
C
C     Y4 intermediates are now at I000. Target can be at I010.
C
      I020 = I010 + LENW
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  240 IRPIK=1,NIRREP
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',DISSIZA(IRPIK),DISSIZA(IRPIK),DISSIZB(IRPIK),
     1            1.0D+00,
     1           CORE(I000+IOFFY4R(IRPIK)),DISSIZA(IRPIK),
     1           CORE(I000+IOFFY4R(IRPIK)),DISSIZA(IRPIK),
     1            0.0D+00,
     1           CORE(I010+IOFFW(IRPIK)),DISSIZA(IRPIK))
      ELSE
      CALL XGEMM('T','N',DISSIZA(IRPIK),DISSIZA(IRPIK),DISSIZB(IRPIK),
     1            1.0D+00,
     1           CORE(I000+IOFFY4R(IRPIK)),DISSIZB(IRPIK),
     1           CORE(I000+IOFFY4R(IRPIK)),DISSIZB(IRPIK),
     1            0.0D+00,
     1           CORE(I010+IOFFW(IRPIK)),DISSIZA(IRPIK))
      ENDIF
C
  240 CONTINUE
C
C     Now loop over symmetry blocks of Y12 and insert appropriate pieces.
C
      DO  390 IRPJL=1,NIRREP
C
      I030 = I020 + IRPDPD(IRPJL,2+ISPIN) * IRPDPD(IRPJL,2+ISPIN)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPJL,2+ISPIN),2,IRPJL,9+ISPIN)
C
      IRPMN = IRPJL
C
      DO  380 IRPL =1,NIRREP
C
      IRPJ = DIRPRD(IRPL,IRPJL)
      IF(IRPJ.GT.IRPL) GOTO 380
C
      DO  370 IRPN =1,NIRREP
C
      IRPM = DIRPRD(IRPN,IRPMN)
      IF(IRPM.GT.IRPN) GOTO 370
C
      IRPML = DIRPRD(IRPM,IRPL)
      IRPJN = DIRPRD(IRPJ,IRPN)
      IRPNL = DIRPRD(IRPN,IRPL)
      IRPJM = DIRPRD(IRPJ,IRPM)
C
      IF(IRPJL.EQ.1)THEN
C
      IF(POP(IRPL,ISPIN).GE.2.AND.POP(IRPN,ISPIN).GE.2)THEN
C
      DO  280 L=2,POP(IRPL,ISPIN)
      DO  270 J=1,L-1
      DO  260 N=2,POP(IRPN,ISPIN)
      DO  250 M=1,N-1
C
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + INDEX(L-1) + J
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + INDEX(N-1) + M
C
      MNJL = (JL-1) * IRPDPD(IRPMN,ISPIN+2) + MN
C
      ML = IOFFOO(IRPL,IRPML,2+ISPIN) + (L-1)*POP(IRPM,ISPIN) + M
      JN = IOFFOO(IRPN,IRPJN,2+ISPIN) + (N-1)*POP(IRPJ,ISPIN) + J
      NL = IOFFOO(IRPL,IRPNL,2+ISPIN) + (L-1)*POP(IRPN,ISPIN) + N
      JM = IOFFOO(IRPM,IRPJM,2+ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
C
      MLJN = I010 + IOFFW(IRPJN) + (JN-1)*DISSIZA(IRPML) + ML - 1
      NLJM = I010 + IOFFW(IRPJM) + (JM-1)*DISSIZA(IRPNL) + NL - 1
C
      CORE(I020 + MNJL - 1) = CORE(I020 + MNJL - 1)
     1                      + CORE(MLJN) - CORE(NLJM)
C
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
C
      ENDIF
C
      ELSE
C
      IF(POP(IRPL,ISPIN).GE.1.AND.POP(IRPN,ISPIN).GE.1.AND.
     1   POP(IRPJ,ISPIN).GE.1.AND.POP(IRPM,ISPIN).GE.1)THEN
C
      DO  320 L=1,POP(IRPL,ISPIN)
      DO  310 J=1,POP(IRPJ,ISPIN)
      DO  300 N=1,POP(IRPN,ISPIN)
      DO  290 M=1,POP(IRPM,ISPIN)
C
      JL = IOFFOO(IRPL,IRPJL,ISPIN) + (L-1)*VRT(IRPJ,ISPIN) + J
      MN = IOFFOO(IRPN,IRPMN,ISPIN) + (N-1)*VRT(IRPM,ISPIN) + M
C
      MNJL = (JL-1) * IRPDPD(IRPMN,ISPIN+2) + MN
C
      ML = IOFFOO(IRPL,IRPML,2+ISPIN) + (L-1)*POP(IRPM,ISPIN) + M
      JN = IOFFOO(IRPN,IRPJN,2+ISPIN) + (N-1)*POP(IRPJ,ISPIN) + J
      NL = IOFFOO(IRPL,IRPNL,2+ISPIN) + (L-1)*POP(IRPN,ISPIN) + N
      JM = IOFFOO(IRPM,IRPJM,2+ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
C
      MLJN = I010 + IOFFW(IRPJN) + (JN-1)*DISSIZA(IRPML) + ML - 1
      NLJM = I010 + IOFFW(IRPJM) + (JM-1)*DISSIZA(IRPNL) + NL - 1
C
      CORE(I020 + MNJL - 1) = CORE(I020 + MNJL - 1)
     1                      + CORE(MLJN) - CORE(NLJM)
C
  290 CONTINUE
  300 CONTINUE
  310 CONTINUE
  320 CONTINUE
C
      ENDIF
C
      ENDIF
C
  370 CONTINUE
  380 CONTINUE
C
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPJL,2+ISPIN),2,IRPJL,9+ISPIN)
C
  390 CONTINUE
C
  400 CONTINUE
C
      ENDIF
C
C     ABAB
C
C                                    SSTGEN
C     MnJl = iJKn * iMKl (JiKn * MiKl ---> JnKi * MlKi; 
C     target : MlJn or JnMl)
C
      LEN  = 0
      LEN2 = 0
      LENW = 0
      DO  420 IRREP=1,NIRREP
      DISSIZ(IRREP) = 0
      DO  410 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZ(IRREP)   =  DISSIZ(IRREP) +
     1                   POP(JRREP,1) * POP(KRREP,2)
  410 CONTINUE
      LEN  = LEN  + DISSIZ(IRREP) * DISSIZ(IRREP)
      LENW = LENW + DISSIZ(IRREP) * DISSIZ(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY4(IRREP)  = 0
      IOFFW(IRREP)   = 0
      ELSE
      IOFFY4(IRREP)  = IOFFY4(IRREP-1) + 
     1                 DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      IOFFW(IRREP)   = IOFFW(IRREP-1) +
     1                 DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      ENDIF
  420 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + 2 * (NOCA * NOCB + NOCA * NOCB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
C     Read ijkn intermediates.
C
      DO  430 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY4(IRREP)),
     1             1,IRPDPD(IRREP,14),2,IRREP,9)
  430 CONTINUE
C
C     (J,i,K,n), (M,i,K,l) ---> (J,n,i,K), (M,l,i,K)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I020),1,'1432')
C
C     Now perform massive matrix multiply :
C
C     ISPIN = 1 : (Ml,iK) * (Jn,iK) = (Ml,Jn)
C
C     Y4 intermediates are now at I000. Target can be at I010.
C
      I020 = I010 + LENW
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  440 IRPIK=1,NIRREP
C
      CALL XGEMM('N','T',DISSIZ(IRPIK),DISSIZ(IRPIK),DISSIZ(IRPIK),
     1            1.0D+00,
     1           CORE(I000+IOFFY4(IRPIK)),DISSIZ(IRPIK),
     1           CORE(I000+IOFFY4(IRPIK)),DISSIZ(IRPIK),
     1            0.0D+00,
     1           CORE(I010+IOFFW(IRPIK)),DISSIZ(IRPIK))
C
  440 CONTINUE
C
C     We have matrix (Ml,Jn). We want (Mn,Jl). Apply SSTGEN and sum
C     directly to symmetry blocks of Y12.
C
C     (M,l,J,n) ---> (M,n,J,l)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I020),1,'1432')
C
C     Now loop over symmetry blocks of Y12 and insert appropriate pieces.
C
      DO  590 IRPJL=1,NIRREP
C
      I020 = I010 + IRPDPD(IRPJL,14) * IRPDPD(IRPJL,14)
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      CALL GETLIST(CORE(I010),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(I000 + IOFFW(IRPJL)),
     1          IRPDPD(IRPJL,14)*IRPDPD(IRPJL,14), 1.0D+00)
C
      CALL PUTLIST(CORE(I010),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
  590 CONTINUE
C
C                                        SSTGEN
C     MnJl = - iJkM * inkl   (JiMk * inkl ---> JMik * iknl; 
C     target : JMnl or nlJM; JMnl ---> JnMl ---> JlMn)
C
      LEN  = 0
      LEN2 = 0
      LENA = 0
      LENB = 0
      LENW = 0
      DO  620 IRREP=1,NIRREP
      DISSIZ(IRREP) = 0
      DISSIZA(IRREP) = 0
      DISSIZB(IRREP) = 0
      DO  610 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZ(IRREP)   =  DISSIZ(IRREP) +
     1                   POP(JRREP,1) * POP(KRREP,2)
      DISSIZA(IRREP)   =  DISSIZA(IRREP) +
     1                   POP(JRREP,1) * POP(KRREP,1)
      DISSIZB(IRREP)   = DISSIZB(IRREP) +
     1                   POP(JRREP,2) * POP(KRREP,2)
  610 CONTINUE
      LENA = LENA + DISSIZ(IRREP)  * DISSIZ(IRREP)
      LENB = LENB + DISSIZB(IRREP) * DISSIZB(IRREP)
      LENW = LENW + DISSIZA(IRREP) * DISSIZB(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY4A(IRREP)  = 0
      IOFFY4AR(IRREP)  = 0
      IOFFY4B(IRREP)  = 0
      IOFFW(IRREP)   = 0
      IOFFWR(IRREP) = 0
      ELSE
      IOFFY4A(IRREP)  = IOFFY4A(IRREP-1) + 
     1                 DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      IOFFY4AR(IRREP)  = IOFFY4AR(IRREP-1) + 
     1                 DISSIZA(IRREP-1)*DISSIZB(IRREP-1)
      IOFFY4B(IRREP)  = IOFFY4B(IRREP-1) + 
     1                 DISSIZB(IRREP-1)*DISSIZB(IRREP-1)
      IOFFW(IRREP)   = IOFFW(IRREP-1) +
     1                 DISSIZA(IRREP-1)*DISSIZB(IRREP-1)
      IOFFWR(IRREP) = IOFFWR(IRREP-1) +
     1                DISSIZ(IRREP-1) * DISSIZ(IRREP-1)
      ENDIF
  620 CONTINUE
C
      I000 = 1
      I010 = I000 + LENA
      I020 = I010 + LENA
      I030 = I020 + 2 * (NOCA * NOCB + NOCA * NOCB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
C     Read JiMk intermediates.
C
      DO  630 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY4A(IRREP)),
     1             1,IRPDPD(IRREP,14),2,IRREP,9)
  630 CONTINUE
C
C     (J,i,M,k) ---> (J,M,i,K)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LENA,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I020),1,'1324')
C
C     (JM,ik) at I000. Read inkl at I020 and form iknl at I010.
C
      I020 = I010 + LENB
      I030 = I020 + LENB
      I040 = I030 + 2 * (NOCA*NOCB + NOCA*NOCB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  640 IRREP=1,NIRREP
      CALL GETLIST(CORE(I020+IOFFY4B(IRREP)),
     1             1,IRPDPD(IRREP,4),2,IRREP,8+IUHF-1)
C
C     (i<n,k<l) ---> (i,n,k<l)
      CALL SYMEXP2(IRREP,POP(1,2),DISSIZB(IRREP),
     1             IRPDPD(IRREP,4),IRPDPD(IRREP,4),
     1             CORE(I020+IOFFY4B(IRREP)),CORE(I020+IOFFY4B(IRREP)))
C     (i,n,k<l) ---> (i,n,k,l)
      CALL  SYMEXP(IRREP,POP(1,2),DISSIZB(IRREP),
     1             CORE(I020+IOFFY4B(IRREP)))
C
  640 CONTINUE
C
C     (i,n,k,l) ---> (i,k,n,l)
C
      CALL SSTGEN(CORE(I020),CORE(I010),LENB,
     1            POP(1,2),POP(1,2),POP(1,2),POP(1,2),
     1            CORE(I030),1,'1324')
C
C     Now perform massive matrix multiply :
C
C     (JM,iK) * (ik,nl) = (JM,nl)
C
C     Y4 intermediates are now at I000 and I010. Target can be at I020.

      I030 = I020 + LENW
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  650 IRPIK=1,NIRREP
C
      CALL XGEMM('N','N',DISSIZA(IRPIK),DISSIZB(IRPIK),DISSIZB(IRPIK),
     1           -1.0D+00,
     1           CORE(I000+IOFFY4AR(IRPIK)),DISSIZA(IRPIK),
     1           CORE(I010+IOFFY4B(IRPIK)),DISSIZB(IRPIK),
     1            0.0D+00,
     1           CORE(I020+IOFFW(IRPIK)),DISSIZA(IRPIK))
C
  650 CONTINUE
C
C     We have matrix (JM,nl). Form (Jn,Ml) and then (Jl,Mn).
C
C     (JM,nl) ---> (Jn,Ml)
C
      CALL SSTGEN(CORE(I020),CORE(I000),LENW,
     1            POP(1,1),POP(1,1),POP(1,2),POP(1,2),
     1            CORE(I030),1,'1324')
C
C     (Jn,Ml) ---> (Jl,Mn)
C
      CALL SSTGEN(CORE(I000),CORE(I020),LENW,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I030),1,'1432')
C
C     Now loop over symmetry blocks of Y12 and insert appropriate pieces.
C
      DO  790 IRPJL=1,NIRREP
C
      I040 = I030 + IRPDPD(IRPJL,14) * IRPDPD(IRPJL,14)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      CALL GETLIST(CORE(I030),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
      CALL VADD(CORE(I030),CORE(I030),CORE(I020 + IOFFWR(IRPJL)),
     1          IRPDPD(IRPJL,14)*IRPDPD(IRPJL,14), 1.0D+00)
C
      CALL PUTLIST(CORE(I030),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
  790 CONTINUE
C
C      RETURN
C                                        SSTGEN
C     MnJl = - IJKM * InKl   (IJKM * InKl ---> IKJM * IKnl; 
C     target : JMnl; JMnl ---> JnMl ---> JlMn)
C
      LEN  = 0
      LEN2 = 0
      LENA = 0
      LENB = 0
      LENW = 0
      DO  820 IRREP=1,NIRREP
      DISSIZ(IRREP) = 0
      DISSIZA(IRREP) = 0
      DISSIZB(IRREP) = 0
      DO  810 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZ(IRREP)   =  DISSIZ(IRREP) +
     1                   POP(JRREP,1) * POP(KRREP,2)
      DISSIZA(IRREP)   =  DISSIZA(IRREP) +
     1                   POP(JRREP,1) * POP(KRREP,1)
      DISSIZB(IRREP)   = DISSIZB(IRREP) +
     1                   POP(JRREP,2) * POP(KRREP,2)
  810 CONTINUE
      LENA = LENA + DISSIZA(IRREP) * DISSIZA(IRREP)
      LENB = LENB + DISSIZA(IRREP) * DISSIZB(IRREP)
      LENW = LENW + DISSIZA(IRREP) * DISSIZB(IRREP)
      IF(IRREP.EQ.1)THEN
      IOFFY4A(IRREP)  = 0
      IOFFY4BR(IRREP)  = 0
      IOFFY4B(IRREP)  = 0
      IOFFW(IRREP)   = 0
      IOFFWR(IRREP) = 0
      ELSE
      IOFFY4A(IRREP)  = IOFFY4A(IRREP-1) + 
     1                 DISSIZA(IRREP-1)*DISSIZA(IRREP-1)
      IOFFY4BR(IRREP)  = IOFFY4BR(IRREP-1) + 
     1                 DISSIZA(IRREP-1)*DISSIZB(IRREP-1)
      IOFFY4B(IRREP)  = IOFFY4B(IRREP-1) + 
     1                 DISSIZ(IRREP-1)*DISSIZ(IRREP-1)
      IOFFW(IRREP)   = IOFFW(IRREP-1) +
     1                 DISSIZA(IRREP-1)*DISSIZB(IRREP-1)
      IOFFWR(IRREP) = IOFFWR(IRREP-1) +
     1                DISSIZ(IRREP-1) * DISSIZ(IRREP-1)
      ENDIF
  820 CONTINUE
C
      I000 = 1
      I010 = I000 + LENA
      I020 = I010 + LENA
      I030 = I020 + 2 * (NOCA * NOCB + NOCA * NOCB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
C     Read IJKM intermediates.
C
      DO  830 IRREP=1,NIRREP
      CALL GETLIST(CORE(I010+IOFFY4A(IRREP)),
     1             1,IRPDPD(IRREP,3),2,IRREP,7)
C
C     (I<J,K<M) ---> (I,J,K<M)
      CALL SYMEXP2(IRREP,POP(1,1),DISSIZA(IRREP),
     1             IRPDPD(IRREP,3),IRPDPD(IRREP,3),
     1             CORE(I010+IOFFY4A(IRREP)),CORE(I010+IOFFY4A(IRREP)))
C     (I,J,K<M) ---> (I,J,K,M)
      CALL  SYMEXP(IRREP,POP(1,1),DISSIZA(IRREP),
     1             CORE(I010+IOFFY4A(IRREP)))
  830 CONTINUE
C
C     (I,J,K,M) ---> (I,K,J,M)
C
      CALL SSTGEN(CORE(I010),CORE(I000),LENA,
     1            POP(1,1),POP(1,1),POP(1,1),POP(1,1),
     1            CORE(I020),1,'1324')
C
C     (IK,JM) at I000. Read InKl at I020 and form IKnl at I010.
C
      I020 = I010 + LENB
      I030 = I020 + LENB
      I040 = I030 + 2 * (NOCA*NOCB + NOCA*NOCB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  840 IRREP=1,NIRREP
      CALL GETLIST(CORE(I020+IOFFY4B(IRREP)),
     1             1,IRPDPD(IRREP,14),2,IRREP,9)
C
  840 CONTINUE
C
C     (InKl) ---> (I,K,n,l)
C
      CALL SSTGEN(CORE(I020),CORE(I010),LENB,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I030),1,'1324')
C
C     Now perform massive matrix multiply :
C
C     (IK,JM) * (IK,nl) = (JM,nl)
C
C     Y4 intermediates are now at I000 and I010. Target can be at I020.
C
      I030 = I020 + LENW
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      DO  850 IRPIK=1,NIRREP
C
      CALL XGEMM('T','N',DISSIZA(IRPIK),DISSIZB(IRPIK),DISSIZA(IRPIK),
     1           -1.0D+00,
     1           CORE(I000+IOFFY4A(IRPIK)),DISSIZA(IRPIK),
     1           CORE(I010+IOFFY4BR(IRPIK)),DISSIZA(IRPIK),
     1            0.0D+00,
     1           CORE(I020+IOFFW(IRPIK)),DISSIZA(IRPIK))
C
  850 CONTINUE
C
C     We have matrix (JM,nl). Form (Jn,Ml) and then (Jl,Mn).
C
C     (JM,nl) ---> (Jn,Ml)
C
      CALL SSTGEN(CORE(I020),CORE(I000),LENW,
     1            POP(1,1),POP(1,1),POP(1,2),POP(1,2),
     1            CORE(I030),1,'1324')
C
C     (Jn,Ml) ---> (Jl,Mn)
C
      CALL SSTGEN(CORE(I000),CORE(I020),LENW,
     1            POP(1,1),POP(1,2),POP(1,1),POP(1,2),
     1            CORE(I030),1,'1432')
C
C     Now loop over symmetry blocks of Y12 and insert appropriate pieces.
C
      DO  990 IRPJL=1,NIRREP
C
      I040 = I030 + IRPDPD(IRPJL,14) * IRPDPD(IRPJL,14)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y12B'
      ENDIF
C
      CALL GETLIST(CORE(I030),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
      CALL VADD(CORE(I030),CORE(I030),CORE(I020 + IOFFWR(IRPJL)),
     1          IRPDPD(IRPJL,14)*IRPDPD(IRPJL,14), 1.0D+00)
C
      CALL PUTLIST(CORE(I030),1,IRPDPD(IRPJL,14),2,IRPJL,12)
C
  990 CONTINUE
      RETURN
      END
