      SUBROUTINE NT3AAB(D1T1A,D1T1B,S1A,S1B,SCR1,SCR2,SCR3,
     &                  EVAL,CORE,MAXCOR,
     &                  INT1,INT2,NONHF,IUHF,IRREPX,ISIDE,ROOT,
     &                  NOCA,NOCB,NVRTA,NVRTB,EL3R3,OVRLAP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ABSOCC,ABSVRT,DISTSZ,DIRPRD,POP,VRT
      LOGICAL INT1,INT2,NONHF
      LOGICAL IJEQL,NONEQL
      LOGICAL IJKEQL2,JKEQL2,IJEQL2,NONEQL2
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      logical threebod
      DOUBLE PRECISION ROOT
      DOUBLE PRECISION EL3R3,ONETHD,OVRLAP
      DOUBLE PRECISION D1T1A(1),D1T1B(1),S1A(1),S1B(1),
     1                 SCR1(1),SCR2(1),SCR3(1)
      DOUBLE PRECISION CORE(1)
      DOUBLE PRECISION EVAL(NOCA + NVRTA,2)
      DOUBLE PRECISION AABL3R3
      DOUBLE PRECISION DIJK
      PARAMETER(MAXNBF=1000)
C
      DIMENSION LEN(8,8),IADZ3(8),LENZ3(9),IADZ3EXP(8),LENZ3EXP(9)
      DIMENSION IADT2(8),LENT2(8),IADW2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(INONHF,IFLAGS(38))
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXNBF,8,2),ABSOCC(MAXNBF,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      INDEX(I) = I*(I-1)/2
c
cdebug aid
      threebod = .true.
C
      ISPIN1 = 1
      ISPIN2 = 2
C
      ONETHD = 1.0D+00/3.0D+00
      AABL3R3 = 0.0D+00
c      OVRLAP = 0.0D+00
      total = 0.0D0
C
      WRITE(6,1015)
 1015 FORMAT(' @NT3AAB-I, Spin case AAB ')
C
C-----------------------------------------------------------------------
C     Initialize D1T1 and "S1" increments.
C-----------------------------------------------------------------------
C
      CALL ZERO(D1T1A,IRPDPD(IRREPX, 9))
      CALL ZERO(D1T1B,IRPDPD(IRREPX,10))
      CALL ZERO(S1A  ,IRPDPD(IRREPX, 9))
      CALL ZERO(S1B  ,IRPDPD(IRREPX,10))
C
C-----------------------------------------------------------------------
C     Get diagonal Fock matrix elements.
C-----------------------------------------------------------------------
C
      CALL GETREC(20,'JOBARC','SCFEVALA',(NOCA+NVRTA)*IINTFP,EVAL(1,1))
      CALL GETREC(20,'JOBARC','SCFEVALB',(NOCB+NVRTB)*IINTFP,EVAL(1,2))
C      
C-----------------------------------------------------------------------
C     Determine the lengths and dimensions of the NIRREP**2 T3 AAB arrays.
C-----------------------------------------------------------------------
C
C      WRITE(6,1070)
C 1070 FORMAT(' @TRPS2-I, Triples symmetry block information ')
      DO   20 IRPIJK=1,NIRREP
      IRPABC = DIRPRD(IRREPX,IRPIJK)
      DO   10 IRPC  =1,NIRREP
      IRPAB = DIRPRD(IRPABC,IRPC)
      LEN(IRPC,IRPIJK) = IRPDPD(IRPAB,1) * VRT(IRPC,2)
C      WRITE(6,1080) IRPIJK,IRPC,IRPAB,LEN(IRPC,IRPIJK)
 1080 FORMAT('   FOR IJK BLOCK ',I4,/,
     1       '       C   BLOCK ',I4,' AB BLOCK ',I4,' LENGTH IS ',I6)
   10 CONTINUE
   20 CONTINUE
C
      DO 1000 IRPIJK=1,NIRREP
C
      IRPABC = DIRPRD(IRREPX,IRPIJK)
C
      if(distsz(irpabc,2).eq.0.or.ndists(irpijk,2).eq.0) goto 1000
C
      DO  990 IRPK  =1,NIRREP
      IF(POP(IRPK,ISPIN2).EQ.0) GOTO 990
C
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN2)
C
      DO  980 IRPJ=1,NIRREP
      IF(POP(IRPJ,ISPIN1).EQ.0) GOTO 980
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IF(IRPI.GT.IRPJ) GOTO 980
C
      IF(POP(IRPI,ISPIN1).EQ.0) GOTO 980
C
      IJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ)THEN
      IJEQL  = .TRUE.
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1))/2
      ELSE
      NONEQL = .TRUE.
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN1)
      NIJ =  POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
      ENDIF
C
      NIK = POP(IRPI,ISPIN1) * POP(IRPK,ISPIN2)
      NJK = POP(IRPJ,ISPIN1) * POP(IRPK,ISPIN2)
C
      IF(IJEQL.AND.POP(IRPJ,ISPIN1).LT.2) GOTO 980
C
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
      DO  30 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADZ3(IRREP) = 1
      ELSE
      IADZ3(IRREP) = IADZ3(IRREP-1) + LEN(IRREP-1,IRPIJK)
      ENDIF
      LENZ3(IRREP) = LEN(IRREP,IRPIJK)
   30 CONTINUE
C
      LENZ3(9) = 0
      DO  40 IRREP=1,NIRREP
      LENZ3(9) = LENZ3(9) + LENZ3(IRREP)
   40 CONTINUE
C
      DO  50 IRREP=1,NIRREP
      IRPA = DIRPRD(IRPABC,IRREP)
      LENZ3EXP(IRREP) = VRT(IRPA,1) * IRPDPD(IRREP,13)
   50 CONTINUE
C
      DO  60 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADZ3EXP(IRREP) = 1
      ELSE
      IADZ3EXP(IRREP) = IADZ3EXP(IRREP-1) + LENZ3EXP(IRREP-1)
      ENDIF
   60 CONTINUE
C
      LENZ3EXP(9) = 0
      DO  70 IRREP=1,NIRREP
      LENZ3EXP(9) = LENZ3EXP(9) + LENZ3EXP(IRREP)
   70 CONTINUE
C
      I000 = 1
      I010 = I000 + LENZ3(9)
      I020 = I010 + LENZ3EXP(9)
      I030 = I020 + LENZ3(9)
      I040 = I030 + LENZ3EXP(9)
      I050 = I040 + IRPDPD(1,13)
      I060 = I050 + IRPDPD(1,13)
      I070 = I060 + LENZ3(9)
      I090 = I070
      I100 = I090 + LENZ3(9)
      ISTART = I100
C
      NEED  = IINTFP * ISTART
      MLEFT = MAXCOR - NEED
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,9010)
      CALL INSMEM('NT3AAB',NEED,MAXCOR)
      ENDIF
C
C     COMPUTE DENOMINATOR ARRAY
C
      CALL MKD32(CORE(I060),EVAL,IADZ3,IRPABC,NOCA,NOCB,NVRTA,NVRTB)

C
      DO  900    K=KLOW,KHIGH
C
      DO  890    J=JLOW,JHIGH
C
      IF(IJEQL) IHIGH = J-1
C
      DO  880    I=ILOW,IHIGH
C
      CALL ZERO(CORE(I000),LENZ3(9))
      CALL ZERO(CORE(I010),LENZ3EXP(9))
      CALL ZERO(CORE(I020),LENZ3(9))
      CALL ZERO(CORE(I030),LENZ3EXP(9))
C
      IF(ISIDE.EQ.1)THEN
        LIST2OFF =  400
      ELSE
        LIST2OFF =  403
      ENDIF
C
C-----------------------------------------------------------------------
C     On first pass (ISIDE=1) compute D*R3 at I000 and I010.
C     On second pass (ISIDE=2) compute D*L3 at I000 and I010.
C-----------------------------------------------------------------------
C
      IF(ISIDE.EQ.1)THEN
        LIST2OFF =  400
        LWOOFF = LWIC11 - 1
        LWVOFF = LWIC15 - 1
        IRREPT = IRREPX
        IRREPW = 1
        CALL GT3WT22(CORE(I000),CORE(I010),CORE(ISTART),
     &               IADZ3,IADZ3EXP,IADT2,LENT2,
     &               I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &               IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MLEFT,
     &               total,irrepx)
C      call checksum("T3-IN-NT3AAB",  CORE(I000), LENZ3(9),sum)
C      tot = tot + sum
C      Write(6,*)
C      write(6,*) tot

C
       if(threebod)then
        IF(IFLAGS(2).EQ.22.AND.IFLAGS2(124).GE.5)THEN
          IRREPT = 1
          IRREPW = IRREPX
          LIST2OFF = 43
          LWOOFF = 376
          LWVOFF = 380
          CALL GT3WT22(CORE(I000),CORE(I010),CORE(ISTART),
     &                 IADZ3,IADZ3EXP,IADT2,LENT2,
     &                 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &                 IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MLEFT)
        ENDIF
       endif
      ENDIF

      IF(ISIDE.EQ.2)THEN
C
        LIST2OFF =  403
        LWOOFF = LWIC21 - 1
        LWVOFF = LWIC25 - 1
        IRREPT = IRREPX
        IRREPW = 1
        CALL GT3WT22(CORE(I000),CORE(I010),CORE(ISTART),
     &               IADZ3,IADZ3EXP,IADT2,LENT2,
     &               I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &               IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MLEFT)
c      CALL S1S223(CORE(ISTART),CORE(I000),CORE(I010),
c     &            IADZ3,IADZ3EXP,
c     &            ISPIN1,ISPIN2,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
c     &            IRPABC,IRREPX,1,IUHF,4)
C <ij||ab> L1
        CALL E_S1S223(CORE(ISTART),CORE(I000),CORE(I010),
     &                IADZ3,IADZ3EXP,ISPIN1,ISPIN2,I,J,K,IRPI,IRPJ,
     &                IRPK,IRPIJ,IRPIK,IRPJK,
     &                IRPABC,IRREPX,1,IUHF,410,2,13,.FALSE.)
C
        IF(IFLAGS2(124).GE.3)THEN
C F L2
          CALL E_S1S223(CORE(ISTART),CORE(I000),CORE(I010),
     &                  IADZ3,IADZ3EXP,
     &                  ISPIN1,ISPIN2,I,J,K,IRPI,IRPJ,IRPK,
     &                  IRPIJ,IRPIK,IRPJK,
     &                  IRPABC,1,IRREPX,IUHF, 93,0,403,.FALSE.)
        ENDIF
      ENDIF
C
      CALL EXPSC2(CORE(I000),CORE(I010),IADZ3,IADZ3EXP,IRPABC)
C
C     Remove denominators (D3T3 ---> T3)
C
      DIJK = EVAL(ABSOCC(I,IRPI,1),1) + EVAL(ABSOCC(J,IRPJ,1),1)
     1     + EVAL(ABSOCC(K,IRPK,2),2) + ROOT
      CALL RMD314(CORE(I000),CORE(I060),LENZ3(9),DIJK)
C
C     Generate R3 vector also in second pass.
C
      IF(ISIDE.EQ.2)THEN
        LIST2OFF = 400
        LWOOFF = LWIC11 - 1
        LWVOFF = LWIC15 - 1
        IRREPT = IRREPX
        IRREPW = 1
        CALL GT3WT22(CORE(I020),CORE(I030),CORE(ISTART),
     &               IADZ3,IADZ3EXP,IADT2,LENT2,
     &               I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &               IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MLEFT)
c        CALL T3WT22( CORE(I020),CORE(I030),CORE(ISTART),
c     1               IADZ3,IADZ3EXP,IADT2,LENT2,
c     1               I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
c     1               IRREPX,LIST2OFF,IUHF,MLEFT)
C
       if(threebod)then
        IF(IFLAGS(2).EQ.22.AND.IFLAGS2(124).GE.5)THEN
          IRREPT = 1
          IRREPW = IRREPX
          LIST2OFF = 43
          LWOOFF = 376
          LWVOFF = 380
          CALL GT3WT22(CORE(I020),CORE(I030),CORE(ISTART),
     &                 IADZ3,IADZ3EXP,IADT2,LENT2,
     &                 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &                 IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MLEFT)
        ENDIF
       endif
C
        CALL EXPSC2(CORE(I020),CORE(I030),IADZ3,IADZ3EXP,IRPABC)
        DIJK = EVAL(ABSOCC(I,IRPI,1),1) + EVAL(ABSOCC(J,IRPJ,1),1)
     1       + EVAL(ABSOCC(K,IRPK,2),2) + ROOT
        CALL RMD314(CORE(I020),CORE(I060),LENZ3(9),DIJK)
C
C     Evaluate overlap.
C
      IF(IUHF.EQ.0)THEN
C       CALL  ICOPY(CORE(I020),CORE(I090),LENZ3(9)*IINTFP)
       CALL  DCOPY(LENZ3(9), CORE(I020), 1, CORE(I090), 1)
       CALL LSC1SC2(CORE(I090),CORE(I020),IADZ3,IRPABC)
       CALL   VADD(CORE(I020),CORE(I020),CORE(I090),LENZ3(9),
     &             ONETHD)
      ENDIF
C
       DO 961 IABC=1,LENZ3(9)
       OVRLAP = OVRLAP + CORE(I000 - 1 + IABC) * CORE(I020 - 1 + IABC)
       AABL3R3 = AABL3R3 + CORE(I000 - 1 +IABC)*CORE(I020 - 1 + IABC)
  961  CONTINUE
C
      IF(IUHF.EQ.0)THEN
       CALL   VADD(CORE(I020),CORE(I020),CORE(I090),LENZ3(9),-ONETHD)
      ENDIF
C
C     R3 at I020, L3 at I000.
C
          DIJK = EVAL(ABSOCC(I,IRPI,1),1) + EVAL(ABSOCC(J,IRPJ,1),1)
     &         + EVAL(ABSOCC(K,IRPK,2),2)
        DO 962 IABC = 1,LENZ3(9)
        CORE(I020 - 1 + IABC) = CORE(I020 - 1 + IABC) * 
     &                          (CORE(I060 - 1 + IABC) - DIJK)
  962 CONTINUE
C
      IF(IUHF.EQ.0)THEN
C       CALL  ICOPY(CORE(I020),CORE(I090),LENZ3(9)*IINTFP)
       CALL  DCOPY(LENZ3(9), CORE(I020), 1, CORE(I090), 1)
       CALL LSC1SC2(CORE(I090),CORE(I020),IADZ3,IRPABC)
       CALL   VADD(CORE(I020),CORE(I020),CORE(I090),LENZ3(9),ONETHD)
      ENDIF
C
        DO 963 IABC=1,LENZ3(9)
        EL3R3 = EL3R3 + CORE(I000 - 1 + IABC) * CORE(I020 - 1 + IABC)
  963   CONTINUE
C
      ENDIF
C
C     Expand alpha*t3(c).
C
      IF(INT1.OR.INT2)THEN
      CALL SYMTRW2(CORE(I000),CORE(I010),CORE(ISTART),
     1             IADZ3,IADW2,IRPABC,1,2)
      ENDIF
C
C     If this is rhf, compute all alpha amplitudes and put at ISTART.
C     Then form effective amplitude at I000.
C
      IF(IUHF.EQ.0)THEN
C       CALL  ICOPY(CORE(I000),CORE(I090),LENZ3(9)*IINTFP)
       CALL  DCOPY(LENZ3(9), CORE(I000), 1, CORE(I090), 1)
       CALL LSC1SC2(CORE(I090),CORE(I000),IADZ3,IRPABC)
       CALL VADD(CORE(I000),CORE(I000),CORE(I090),LENZ3(9), 1.0D+00)
      ENDIF

      IF(INT1)THEN
      CALL  E_T1T32N(D1T1A,D1T1B,S1A,S1B,
     1             CORE(I000),CORE(I010),
     1             SCR1,SCR2,SCR3,
     1             IADZ3,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1             I,J,K,1,2,
     1             NONHF,IRREPX)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
       CALL VADD(CORE(I000),CORE(I000),CORE(I090),LENZ3(9),-1.0D+00)
      ENDIF
C
      IF(INT2)THEN
C
C     --- D2T2 = F T3 for non-Hartree-Fock cases ---
C
        IF( ((IFLAGS2(124).EQ.1.OR.IFLAGS2(124).EQ.2) .AND. NONHF) .OR.
     &        IFLAGS2(124).GE.3                                   )THEN
          IF(IFLAGS2(124).LE.2) LFOFF = 2
          IF(IFLAGS2(124).GE.3) LFOFF = 0
          CALL T2FT323(CORE(I000),CORE(I010),CORE(ISTART),
     &                 IADZ3,IADW2,
     &                 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,
     &                 IUHF,IRREPX,406, 93,LFOFF,1.0D+00,.TRUE.,.FALSE.)
        ENDIF
C
C     --- D2T2 = W T3 ---
C
      CALL E_T2T32(CORE(I000),CORE(I010),CORE(ISTART),
     1           IADZ3,IADW2,IRPI,IRPJ,IRPK,
     1           IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,IUHF,
     1           SCR1,SCR2,SCR3,IRREPX)
C
      CALL E_T2T32O(CORE(I000),CORE(I010),CORE(ISTART),IADZ3,IADW2,
     1              IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1              IUHF,SCR1,SCR2,SCR3,IRREPX)
      ENDIF
C
  880 CONTINUE
  890 CONTINUE
  900 CONTINUE
C
  980 CONTINUE
  990 CONTINUE
C
 1000 CONTINUE
C
CSSS       call checksum("NT3AAB", CORE(I000), LENZ3(9))
C
      IF(INT1)THEN
      CALL GETLST(SCR1,1,1,1,5,410)
      CALL   VADD(SCR1,SCR1,D1T1A,IRPDPD(IRREPX, 9),1.0D+00)
CSSS      call sumblk(scr1,irpdpd(irrepx,9))
      CALL PUTLST(SCR1,1,1,1,5,410)
      CALL GETLST(SCR1,1,1,1,5+IUHF,410)
      CALL   VADD(SCR1,SCR1,D1T1B,IRPDPD(IRREPX,10),1.0D+00)
      CALL PUTLST(SCR1,1,1,1,5+IUHF,410)
      ENDIF
C
      if(iuhf.eq.0)then
       el3r3 = 2.0D+00 * el3r3
       ovrlap = 2.0D+00 * ovrlap
      endif
      if(iside.eq.2)then
      write(6,9020) EL3R3
      write(6,9030) ovrlap
      endif
      if(iside.eq.2)then
       write(6,*) ' overlap contibution from aab ',aabl3r3
      endif
      RETURN
 9010 FORMAT(' @NT3AAB-I, Insufficient memory to continue. ')
 9020 FORMAT(' @NT3AAB-I, el3r3   ',F20.12)
 9030 FORMAT(' @NT3AAB-I, overlap ',F20.12)
      END
