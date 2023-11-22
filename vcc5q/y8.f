      SUBROUTINE Y8(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZV,DISSIZT,DISSIZW,SIZAAO,SIZAAV,SIZBBO,SIZBBV
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1)
      DIMENSION SIZAAO(8),SIZAAV(8),SIZBBO(8),SIZBBV(8)
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
 1000 FORMAT(' @Y8-I,   <md||le> * T2 contribution. ')
C
      DO   20 IRREP=1,NIRREP
      SIZAAO(IRREP) = 0
      SIZAAV(IRREP) = 0
      SIZBBO(IRREP) = 0
      SIZBBV(IRREP) = 0
C
      DO   10 JRREP=1,NIRREP
      KRREP=DIRPRD(IRREP,JRREP)
C
      SIZAAO(IRREP) = SIZAAO(IRREP) + POP(JRREP,1)*POP(KRREP,1)
      SIZAAV(IRREP) = SIZAAV(IRREP) + VRT(JRREP,1)*VRT(KRREP,1)
      SIZBBO(IRREP) = SIZBBO(IRREP) + POP(JRREP,2)*POP(KRREP,2)
      SIZBBV(IRREP) = SIZBBV(IRREP) + VRT(JRREP,2)*VRT(KRREP,2)
C
   10 CONTINUE
   20 CONTINUE
C
C
C     CEKM = MDLE * CDKL
C     cekm = mdle * cdkl
C
      DO  100 ISPIN=1,IUHF+1
C
C     Read <MD//LE) or <md//le>. These are stored as (E,M,D,L).
C     Read t(CD,KL) or t(cd,kl). These are stored as (C,L,D,K) and
C     will be resorted to (C,K,D,L).
C
      LEN = 0
      DO   30 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,8+ISPIN) * IRPDPD(IRREP,8+ISPIN)
   30 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
 1020 FORMAT(' @Y8-I, Insufficient memory. Need ',I15,' . Got ',I15)
      STOP 'Y8'
      ENDIF
C
      CALL GETALL(CORE(I000),LEN,1,22+ISPIN)
      CALL GETALL(CORE(I020),LEN,1,33+ISPIN)
C
C     Sort t(c,l,d,k) to t(c,k,d,l).
C
      CALL SSTGEN(CORE(I020),CORE(I010),LEN,
     1            VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     1            CORE(I030),1,'1432')
C
      IOFFV = I000
      IOFFT = I010
      IOFFW = I020
      DO   40 IRPDL=1,NIRREP
C
      DISSIZV = IRPDPD(IRPDL,8+ISPIN)
      NDISV   = IRPDPD(IRPDL,8+ISPIN)
      DISSIZT = IRPDPD(IRPDL,8+ISPIN)
      NDIST   = IRPDPD(IRPDL,8+ISPIN)
C
      CALL XGEMM('N','T',DISSIZT,DISSIZV,NDISV,
     1            1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
C
      IOFFT = IOFFT + DISSIZT *   NDIST
      IOFFV = IOFFV + DISSIZV *   NDISV
      IOFFW = IOFFW + DISSIZT * DISSIZV
C
   40 CONTINUE
C
C     We have W(C,K,E,M) at I020. Reorder to (C,E,K,M) at I000.
C
      CALL SSTGEN(CORE(I020),CORE(I000),LEN,
     1            VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     1            CORE(I030),1,'1324')
C
      IOFFW = I000
      DO 50 IRPKM=1,NIRREP
C
      IF(ISPIN.EQ.1)THEN
      DISSIZW = SIZAAV(IRPKM)
      NDISW   = SIZAAO(IRPKM)
      ELSE
      DISSIZW = SIZBBV(IRPKM)
      NDISW   = SIZBBO(IRPKM)
      ENDIF
      CALL GETLIST(CORE(I010),1,NDISW,2,IRPKM,12+ISPIN)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(IOFFW),DISSIZW*NDISW,
     1           1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISW,2,IRPKM,12+ISPIN)
C
C
      IOFFW = IOFFW + DISSIZW * NDISW
C
   50 CONTINUE
  100 CONTINUE
C
C     CEKM = MdlE * CdKl
C     cekm = mDLe * cDkL
C
      DO  200 ISPIN=1,IUHF+1
C
C     Read <MD//LE) or <md//le>. These are stored as (E,M,D,L).
C     Read t(CD,KL) or t(cd,kl). These are stored as (C,L,D,K) and
C     will be resorted to (C,K,D,L).
C
      LENV = 0
      LENT = 0
      LENW = 0
      DO  130 IRREP=1,NIRREP
      LENV = LENV + IRPDPD(IRREP,11)      * IRPDPD(IRREP,12)
      LENT = LENT + IRPDPD(IRREP,13)      * IRPDPD(IRREP,14)
      LENW = LENW + IRPDPD(IRREP,8+ISPIN) * IRPDPD(IRREP,8+ISPIN)
  130 CONTINUE
      LEN = MAX(LENV,LENT,LENW)
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y8'
      ENDIF
C
C     ISPIN = 1 : -V(E,M,d,l) at I000
C     ISPIN = 2 : -V(D,L,e,m) at I000
C
      CALL GETALL(CORE(I000),LEN,1,18)
C
C     ISPIN = 1 :  T(C,K,d,l) at I010
C     ISPIN = 2 :  T(D,L,c,k) at I010
      CALL GETALL(CORE(I010),LEN,1,37)
C
      IOFFV = I000
      IOFFT = I010
      IOFFW = I020
      DO  140 IRPDL=1,NIRREP
C
      DISSIZV = IRPDPD(IRPDL, 9)
      NDISV   = IRPDPD(IRPDL,10)
      DISSIZT = IRPDPD(IRPDL, 9)
      NDIST   = IRPDPD(IRPDL,10)
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',DISSIZT,DISSIZV,NDIST,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
      ELSE
      CALL XGEMM('T','N',NDIST,NDISV,DISSIZT,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),NDIST)
      ENDIF
C
      IOFFT = IOFFT + DISSIZT *   NDIST
      IOFFV = IOFFV + DISSIZV *   NDISV
      IF(ISPIN.EQ.1)THEN
      IOFFW = IOFFW + DISSIZT * DISSIZV
      ELSE
      IOFFW = IOFFW + NDIST   * NDISV
      ENDIF
C
  140 CONTINUE
C
C     We have W(C,K,E,M) at I020. Reorder to (C,E,K,M) at I000.
C
      CALL SSTGEN(CORE(I020),CORE(I000),LEN,
     1            VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     1            CORE(I030),1,'1324')
C
      IOFFW = I000
      DO 150 IRPKM=1,NIRREP
C
      IF(ISPIN.EQ.1)THEN
      DISSIZW = SIZAAV(IRPKM)
      NDISW   = SIZAAO(IRPKM)
      ELSE
      DISSIZW = SIZBBV(IRPKM)
      NDISW   = SIZBBO(IRPKM)
      ENDIF
      CALL GETLIST(CORE(I010),1,NDISW,2,IRPKM,12+ISPIN)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(IOFFW),DISSIZW*NDISW,
     1           1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISW,2,IRPKM,12+ISPIN)
C
C
      IOFFW = IOFFW + DISSIZW * NDISW
C
  150 CONTINUE
  200 CONTINUE
C
C     CeKm = mDLe * CDKL
C     cEkM = MdlE * cdkl
C
      DO  300 ISPIN=1,IUHF+1
C
C     Read <mD//Le) or <Md//lE>. These are stored as (D,L,e,m) or
C     (E,M,d,l). Read t(CD,KL) or t(cd,kl). These are stored as (C,L,D,K) 
C     and will be resorted to (C,K,D,L).
C
      LENV = 0
      LENT = 0
      LENW = 0
      DO  230 IRREP=1,NIRREP
      LENV = LENV + IRPDPD(IRREP,11)      * IRPDPD(IRREP,12)
      LENT = LENT + IRPDPD(IRREP,8+ISPIN) * IRPDPD(IRREP,8+ISPIN)
      LENW = LENW + IRPDPD(IRREP,13)      * IRPDPD(IRREP,14)
  230 CONTINUE
      LEN = MAX(LENV,LENT,LENW)
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y8'
      ENDIF
C
C     ISPIN = 1 : -V(D,L,e,m) at I000
C     ISPIN = 2 : -V(E,M,d,l) at I000
C
      CALL GETALL(CORE(I000),LEN,1,18)
C
C     ISPIN = 1 :  T(C,L,D,K) at I020
C     ISPIN = 2 :  T(c,l,d,k) at I020
      CALL GETALL(CORE(I020),LEN,1,33+ISPIN)
C
C     Sort t(c,l,d,k) to t(c,k,d,l).
C
      CALL SSTGEN(CORE(I020),CORE(I010),LEN,
     1            VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     1            CORE(I030),1,'1432')
C
      IOFFV = I000
      IOFFT = I010
      IOFFW = I020
      DO  240 IRPDL=1,NIRREP
C
      DISSIZV = IRPDPD(IRPDL, 9)
      NDISV   = IRPDPD(IRPDL,10)
      DISSIZT = IRPDPD(IRPDL, 8+ISPIN)
      NDIST   = IRPDPD(IRPDL, 8+ISPIN)
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','N',DISSIZT,NDISV,NDIST,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
      ELSE
      CALL XGEMM('N','T',DISSIZT,DISSIZV,NDISV,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
      ENDIF
C
      IOFFT = IOFFT + DISSIZT *   NDIST
      IOFFV = IOFFV + DISSIZV *   NDISV
      IF(ISPIN.EQ.1)THEN
      IOFFW = IOFFW + DISSIZT * NDISV
      ELSE
      IOFFW = IOFFW + DISSIZT * DISSIZV
      ENDIF
C
  240 CONTINUE
C
C     We have W(C,K,e,m) or W(c,k,E,M) at I020. Reorder to (C,e,K,m) or
C     (c,E,k,M) at I000.
C
      CALL SSTGEN(CORE(I020),CORE(I000),LEN,
     1        VRT(1,ISPIN),POP(1,ISPIN),VRT(1,3-ISPIN),POP(1,3-ISPIN),
     1            CORE(I030),1,'1324')
C
      IOFFW = I000
      DO 250 IRPKM=1,NIRREP
C
      DISSIZW = IRPDPD(IRPKM,13)
      NDISW   = IRPDPD(IRPKM,14)
      CALL GETLIST(CORE(I010),1,NDISW,2,IRPKM,14+ISPIN)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(IOFFW),DISSIZW*NDISW,
     1           1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISW,2,IRPKM,14+ISPIN)
C
      IOFFW = IOFFW + DISSIZW * NDISW
  250 CONTINUE
  300 CONTINUE
C
C     CeKm = mdle * CdKl
C     cEkM = MDLE * cDkL
C
      DO  400 ISPIN=1,IUHF+1
C
C     Read <md//le) or <MD//LE>. These are stored as (e,m,d,l) or
C     (E,M,D,L). Read t(Cd,Kl) or t(cD,kL). These are stored as (C,K,d,l) 
C     or (D,L,c,k).
C
      LENV = 0
      LENT = 0
      LENW = 0
      DO  330 IRREP=1,NIRREP
      LENV = LENV + IRPDPD(IRREP,11-ISPIN) * IRPDPD(IRREP,11-ISPIN)
      LENT = LENT + IRPDPD(IRREP,13)       * IRPDPD(IRREP,14)
      LENW = LENW + IRPDPD(IRREP,13)       * IRPDPD(IRREP,14)
  330 CONTINUE
      LEN = MAX(LENV,LENT,LENW)
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y8'
      ENDIF
C
C     ISPIN = 1 :  V(e,m,d,l) at I000
C     ISPIN = 2 :  V(E,M,D,L) at I000
C
      IF(IUHF.GT.0)THEN
      CALL GETALL(CORE(I000),LEN,1,25-ISPIN)
      ELSE
      CALL GETALL(CORE(I000),LEN,1,23)
      ENDIF
C
C     ISPIN = 1 :  T(C,K,d,l) at I010
C     ISPIN = 2 :  T(D,L,c,k) at I010
      CALL GETALL(CORE(I010),LEN,1,37)
C
      IOFFV = I000
      IOFFT = I010
      IOFFW = I020
      DO  340 IRPDL=1,NIRREP
C
      DISSIZV = IRPDPD(IRPDL,11-ISPIN)
      NDISV   = IRPDPD(IRPDL,11-ISPIN)
      DISSIZT = IRPDPD(IRPDL, 9)
      NDIST   = IRPDPD(IRPDL,10)
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',DISSIZT,DISSIZV,NDIST,
     1             1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
      ELSE
      CALL XGEMM('T','T',NDIST,DISSIZV,NDISV,
     1             1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),NDIST)
      ENDIF
C
      IOFFT = IOFFT + DISSIZT *   NDIST
      IOFFV = IOFFV + DISSIZV *   NDISV
      IF(ISPIN.EQ.1)THEN
      IOFFW = IOFFW + DISSIZT * DISSIZV
      ELSE
      IOFFW = IOFFW + NDIST   * DISSIZV
      ENDIF
C
  340 CONTINUE
C
C     We have W(C,K,e,m) or W(c,k,E,M) at I020. Reorder to (C,e,K,m) or
C     (c,E,k,M) at I000.
C
      CALL SSTGEN(CORE(I020),CORE(I000),LEN,
     1        VRT(1,ISPIN),POP(1,ISPIN),VRT(1,3-ISPIN),POP(1,3-ISPIN),
     1            CORE(I030),1,'1324')
C
      IOFFW = I000
      DO 350 IRPKM=1,NIRREP
C
      DISSIZW = IRPDPD(IRPKM,13)
      NDISW   = IRPDPD(IRPKM,14)
      CALL GETLIST(CORE(I010),1,NDISW,2,IRPKM,14+ISPIN)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(IOFFW),DISSIZW*NDISW,
     1           1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISW,2,IRPKM,14+ISPIN)
C
      IOFFW = IOFFW + DISSIZW * NDISW
  350 CONTINUE
  400 CONTINUE
C
C     CekM = MdLe * CdkL
C     cEKm = mDlE * cDKl
C
      DO  500 ISPIN=1,IUHF+1
C
C     Read <Md//Le) or <mD//lE>. These are stored as (e,M,d,L) or
C     (E,m,D,l). Read t(Cd,kL) or t(cD,Kl). These are stored as (C,k,d,L) 
C     or (D,l,c,K).
C
      LENV = 0
      LENT = 0
      LENW = 0
      DO  430 IRREP=1,NIRREP
      LENV = LENV + IRPDPD(IRREP,13-ISPIN) * IRPDPD(IRREP,13-ISPIN)
      LENT = LENT + IRPDPD(IRREP,13)       * IRPDPD(IRREP,14)
      LENW = LENW + IRPDPD(IRREP,13)       * IRPDPD(IRREP,14)
  430 CONTINUE
      LEN = MAX(LENV,LENT,LENW)
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y8'
      ENDIF
C
C     ISPIN = 1 :  V(e,M,d,L) at I000
C     ISPIN = 2 :  V(E,m,D,l) at I000
C
      IF(IUHF.GT.0)THEN
      CALL GETALL(CORE(I000),LEN,1,27-ISPIN)
      ELSE
      CALL GETALL(CORE(I000),LEN,1,25)
      ENDIF
C
C     ISPIN = 1 :  T(C,k,d,L) at I010
C     ISPIN = 2 :  T(D,l,c,K) at I010
      CALL GETALL(CORE(I010),LEN,1,39)
C
      IOFFV = I000
      IOFFT = I010
      IOFFW = I020
      DO  440 IRPDL=1,NIRREP
C
      DISSIZV = IRPDPD(IRPDL,13-ISPIN)
      NDISV   = IRPDPD(IRPDL,13-ISPIN)
      DISSIZT = IRPDPD(IRPDL,11)
      NDIST   = IRPDPD(IRPDL,12)
C
      IF(ISPIN.EQ.1)THEN
      CALL XGEMM('N','T',DISSIZT,DISSIZV,NDIST,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),DISSIZT)
      ELSE
      CALL XGEMM('T','T',NDIST,DISSIZV,NDISV,
     1            -1.0D+00,
     1            CORE(IOFFT),DISSIZT,
     1            CORE(IOFFV),DISSIZV,0.0D+00,
     1            CORE(IOFFW),NDIST)
      ENDIF
C
      IOFFT = IOFFT + DISSIZT *   NDIST
      IOFFV = IOFFV + DISSIZV *   NDISV
      IF(ISPIN.EQ.1)THEN
      IOFFW = IOFFW + DISSIZT * DISSIZV
      ELSE
      IOFFW = IOFFW + NDIST   * DISSIZV
      ENDIF
C
  440 CONTINUE
C
C     We have W(C,k,e,M) or W(c,K,E,m) at I020. Reorder to (C,e,k,M) or
C     (c,E,K,m) at I000.
C
      CALL SSTGEN(CORE(I020),CORE(I000),LEN,
     1        VRT(1,ISPIN),POP(1,3-ISPIN),VRT(1,3-ISPIN),POP(1,ISPIN),
     1            CORE(I030),1,'1324')
C
      IOFFW = I000
      DO 450 IRPKM=1,NIRREP
C
      DISSIZW = IRPDPD(IRPKM,13)
      NDISW   = IRPDPD(IRPKM,14)
      CALL GETLIST(CORE(I010),1,NDISW,2,IRPKM,16+ISPIN)
C
      CALL VADD(CORE(I010),CORE(I010),CORE(IOFFW),DISSIZW*NDISW,
     1           1.0D+00)
      CALL PUTLIST(CORE(I010),1,NDISW,2,IRPKM,16+ISPIN)
C
      IOFFW = IOFFW + DISSIZW * NDISW
  450 CONTINUE
  500 CONTINUE
      RETURN
      END
