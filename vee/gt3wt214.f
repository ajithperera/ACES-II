      SUBROUTINE GT3WT214(Z,CORE,ISPIN,IADZ,
     &                    I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &                    IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,MAXCOR)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION Z,CORE
      INTEGER ISPIN,IADZ,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &        IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,MAXCOR
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IOFFVV,IOFFOO,IOFFVO,INDEX
C-----------------------------------------------------------------------
      INTEGER LISTWO,LISTWV,DSZ,DSZEXP,IADT2,IADW,LENW
      INTEGER IRPA,IRPE,IRPEA,IRPEI,IRPEJ,IRPEK,IRPBC,
     &        IRPM,IRPMA,IRPMI,IRPMJ,IRPMK,M,MI,MJ,MK,
     &        IJ,IK,JK,I000,I010,I020,NEED,IRREP
C-----------------------------------------------------------------------
      DIMENSION Z(1),CORE(1)
      DIMENSION IADZ(8),IADT2(8),IADW(8),LENW(8)
C-----------------------------------------------------------------------
C
C     New version uses Z, IADZ, CORE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
      LISTWV = LWVOFF + ISPIN
      LISTWO = LWOOFF + ISPIN

CSSS      Write(6,*) "LISTWV,LISTVO", LISTWV,LISTWO
C
C-----------------------------------------------------------------------
C     - T(EA,IJ) * W(BC,EK) + T(EB,IJ) * W(AC,EK) + T(EC,IJ) * W(BA,EK)
C-----------------------------------------------------------------------
C
      IRPEA  = DIRPRD(IRREPT,IRPIJ)
      DSZ    = IRPDPD(IRPEA,     ISPIN)
      DSZEXP = IRPDPD(IRPEA,18 + ISPIN)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEA,18 + ISPIN)
C
      DO   10 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = I000
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     &                                VRT(DIRPRD(IRREP-1,IRPEA),ISPIN)
      ENDIF
   10 CONTINUE
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      CALL  GETLST(CORE(I000),IJ,1,1,IRPIJ,LIST2OFF + ISPIN)
CSSS      call checksum("IM50 making", core(i000), IRPDPD(IRPEA,18 + ISPIN))
CSSS      Write(6,*) "LIST2OFF+ISPIN", LIST2OFF+ISPIN
      CALL SYMEXP2(IRPEA,VRT(1,ISPIN),DSZEXP,DSZ,1,CORE(I000),
     &                                             CORE(I000))
C
      DO   20 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 20
C
      IRPEK = DIRPRD(IRPE,IRPK)
      IRPBC = DIRPRD(IRREPW,IRPEK)
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &            IOFFVO(IRPK,IRPEK,ISPIN)+(K-1)*VRT(IRPE,ISPIN)+1,
     &            VRT(IRPE,ISPIN),2,IRPEK,LISTWV)
C
      IRPA   = DIRPRD(IRPE,IRPEA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),VRT(IRPE,ISPIN),
     &           -1.0D+00,
     &           CORE(I010)       ,IRPDPD(IRPBC,ISPIN),
     &           CORE(IADT2(IRPA)),VRT(IRPE,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
   20 CONTINUE
C
C-----------------------------------------------------------------------
C     - T(EA,JK) * W(BC,EI) + T(EB,JK) * W(AC,EI) + T(EC,JK) * W(BA,EI)
C-----------------------------------------------------------------------
C
      IRPEA  = DIRPRD(IRREPT,IRPJK)
      DSZ    = IRPDPD(IRPEA,     ISPIN)
      DSZEXP = IRPDPD(IRPEA,18 + ISPIN)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEA,18 + ISPIN)
C
      DO  110 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = I000
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     &                                VRT(DIRPRD(IRREP-1,IRPEA),ISPIN)
      ENDIF
  110 CONTINUE
C
      IF(IRPJ.EQ.IRPK)THEN
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
C
      CALL  GETLST(CORE(I000),JK,1,1,IRPJK,LIST2OFF + ISPIN)
      CALL SYMEXP2(IRPEA,VRT(1,ISPIN),DSZEXP,DSZ,1,CORE(I000),
     &                                             CORE(I000))
C
      DO  120 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 120
C
      IRPEI = DIRPRD(IRPE,IRPI)
      IRPBC = DIRPRD(IRREPW,IRPEI)
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &            IOFFVO(IRPI,IRPEI,ISPIN)+(I-1)*VRT(IRPE,ISPIN)+1,
     &            VRT(IRPE,ISPIN),2,IRPEI,LISTWV)
C
      IRPA   = DIRPRD(IRPE,IRPEA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),VRT(IRPE,ISPIN),
     &           -1.0D+00,
     &           CORE(I010)       ,IRPDPD(IRPBC,ISPIN),
     &           CORE(IADT2(IRPA)),VRT(IRPE,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
  120 CONTINUE
C
C-----------------------------------------------------------------------
C       T(EA,IK) * W(BC,EJ) - T(EB,IK) * W(AC,EJ) - T(EC,IK) * W(BA,EJ)
C-----------------------------------------------------------------------
C
      IRPEA  = DIRPRD(IRREPT,IRPIK)
      DSZ    = IRPDPD(IRPEA,     ISPIN)
      DSZEXP = IRPDPD(IRPEA,18 + ISPIN)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEA,18 + ISPIN)
C
      DO  210 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = I000
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     &                                VRT(DIRPRD(IRREP-1,IRPEA),ISPIN)
      ENDIF
  210 CONTINUE
C
      IF(IRPI.EQ.IRPK)THEN
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      CALL  GETLST(CORE(I000),IK,1,1,IRPIK,LIST2OFF + ISPIN)
      CALL SYMEXP2(IRPEA,VRT(1,ISPIN),DSZEXP,DSZ,1,CORE(I000),
     &                                             CORE(I000))
C
      DO  220 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 220
C
      IRPEJ = DIRPRD(IRPE,IRPJ)
      IRPBC = DIRPRD(IRREPW,IRPEJ)
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &            IOFFVO(IRPJ,IRPEJ,ISPIN)+(J-1)*VRT(IRPE,ISPIN)+1,
     &            VRT(IRPE,ISPIN),2,IRPEJ,LISTWV)
C
      IRPA   = DIRPRD(IRPE,IRPEA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),VRT(IRPE,ISPIN),
     &            1.0D+00,
     &           CORE(I010)       ,IRPDPD(IRPBC,ISPIN),
     &           CORE(IADT2(IRPA)),VRT(IRPE,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
  220 CONTINUE
C
C-----------------------------------------------------------------------
C     T(BC,MK) * W(MA,IJ) - T(AC,MK) * W(MB,IJ) - T(BA,MK) * W(MC,IJ)
C-----------------------------------------------------------------------
C
      IRPMA = DIRPRD(IRREPW,IRPIJ)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,8 + ISPIN)
C
      DO  310 IRPA=1,NIRREP
      IRPM = DIRPRD(IRPA,IRPMA)
      LENW(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
      IF(IRPA.EQ.1)THEN
      IADW(IRPA) = I000
      ELSE
      IADW(IRPA) = IADW(IRPA-1) + LENW(IRPA-1)
      ENDIF
  310 CONTINUE
C
      IF(IRPIJ.EQ.1)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LISTWO)
C
      DO  350 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 350
C
      IRPMK = DIRPRD(IRPM,IRPK)
      IRPBC = DIRPRD(IRREPT,IRPMK)
C
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      IF(IRPM.LT.IRPK)THEN
      DO  320 M=1,POP(IRPM,ISPIN)
      MK = IOFFOO(IRPK,IRPMK,ISPIN) + (K-1)*POP(IRPM,ISPIN) + M
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MK,1,1,IRPMK,LIST2OFF + ISPIN)
  320 CONTINUE
      ENDIF
C
      IF(IRPM.GT.IRPK)THEN
      DO  330 M=1,POP(IRPM,ISPIN)
      MK = IOFFOO(IRPM,IRPMK,ISPIN) + (M-1)*POP(IRPK,ISPIN) + K
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MK,1,1,IRPMK,LIST2OFF + ISPIN)
  330 CONTINUE
      CALL VMINUS(CORE(I010),POP(IRPM,ISPIN)*IRPDPD(IRPBC,ISPIN))
      ENDIF
C
      IF(IRPM.EQ.IRPK)THEN
      DO  340 M=1,POP(IRPM,ISPIN)
         IF(M.LT.K)THEN
         MK = IOFFOO(IRPK,IRPMK,ISPIN) + INDEX(K-1) + M
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MK,1,1,IRPMK,LIST2OFF + ISPIN)
         ENDIF
C
         IF(M.GT.K)THEN
         MK = IOFFOO(IRPM,IRPMK,ISPIN) + INDEX(M-1) + K
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MK,1,1,IRPMK,LIST2OFF + ISPIN)
         CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
C
         IF(M.EQ.K)THEN
         CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
  340 CONTINUE
C
      ENDIF
C
      IRPA   = DIRPRD(IRPM,IRPMA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),POP(IRPM,ISPIN),
     &            1.0D+00,
     &           CORE(I010),IRPDPD(IRPBC,ISPIN),
     &           CORE(IADW(IRPA)),POP(IRPM,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
  350 CONTINUE
C
C-----------------------------------------------------------------------
C     T(BC,MI) * W(MA,JK) - T(AC,MI) * W(MB,JK) - T(BA,MI) * W(MC,JK)
C-----------------------------------------------------------------------
C
      IRPMA = DIRPRD(IRREPW,IRPJK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,8 + ISPIN)
C
      DO  410 IRPA=1,NIRREP
      IRPM = DIRPRD(IRPA,IRPMA)
      LENW(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
      IF(IRPA.EQ.1)THEN
      IADW(IRPA) = I000
      ELSE
      IADW(IRPA) = IADW(IRPA-1) + LENW(IRPA-1)
      ENDIF
  410 CONTINUE
C
      IF(IRPJK.EQ.1)THEN
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
C
      CALL GETLST(CORE(I000),JK,1,2,IRPJK,LISTWO)
C
      DO  450 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 450
C
      IRPMI = DIRPRD(IRPM,IRPI)
      IRPBC = DIRPRD(IRREPT,IRPMI)
C
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      IF(IRPM.LT.IRPI)THEN
      DO  420 M=1,POP(IRPM,ISPIN)
      MI = IOFFOO(IRPI,IRPMI,ISPIN) + (I-1)*POP(IRPM,ISPIN) + M
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MI,1,1,IRPMI,LIST2OFF + ISPIN)
  420 CONTINUE
      ENDIF
C
      IF(IRPM.GT.IRPI)THEN
      DO  430 M=1,POP(IRPM,ISPIN)
      MI = IOFFOO(IRPM,IRPMI,ISPIN) + (M-1)*POP(IRPI,ISPIN) + I
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MI,1,1,IRPMI,LIST2OFF + ISPIN)
  430 CONTINUE
      CALL VMINUS(CORE(I010),POP(IRPM,ISPIN)*IRPDPD(IRPBC,ISPIN))
      ENDIF
C
      IF(IRPM.EQ.IRPI)THEN
      DO  440 M=1,POP(IRPM,ISPIN)
         IF(M.LT.I)THEN
         MI = IOFFOO(IRPI,IRPMI,ISPIN) + INDEX(I-1) + M
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MI,1,1,IRPMI,LIST2OFF + ISPIN)
         ENDIF
C
         IF(M.GT.I)THEN
         MI = IOFFOO(IRPM,IRPMI,ISPIN) + INDEX(M-1) + I
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MI,1,1,IRPMI,LIST2OFF + ISPIN)
         CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
C
         IF(M.EQ.I)THEN
         CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
  440 CONTINUE
C
      ENDIF
C
      IRPA   = DIRPRD(IRPM,IRPMA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),POP(IRPM,ISPIN),
     &            1.0D+00,
     &           CORE(I010),IRPDPD(IRPBC,ISPIN),
     &           CORE(IADW(IRPA)),POP(IRPM,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
  450 CONTINUE
C
C-----------------------------------------------------------------------
C   - T(BC,MJ) * W(MA,IK) + T(AC,MJ) * W(MB,IK) + T(BA,MJ) * W(MC,IK)
C-----------------------------------------------------------------------
C
      IRPMA = DIRPRD(IRREPW,IRPIK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,8 + ISPIN)
C
      DO  510 IRPA=1,NIRREP
      IRPM = DIRPRD(IRPA,IRPMA)
      LENW(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
      IF(IRPA.EQ.1)THEN
      IADW(IRPA) = I000
      ELSE
      IADW(IRPA) = IADW(IRPA-1) + LENW(IRPA-1)
      ENDIF
  510 CONTINUE
C
      IF(IRPIK.EQ.1)THEN
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      CALL GETLST(CORE(I000),IK,1,2,IRPIK,LISTWO)
C
      DO  550 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 550
C
      IRPMJ = DIRPRD(IRPM,IRPJ)
      IRPBC = DIRPRD(IRREPT,IRPMJ)
C
      I020 = I010 + IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
      NEED = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT214',NEED,MAXCOR)
      ENDIF
C
      IF(IRPM.LT.IRPJ)THEN
      DO  520 M=1,POP(IRPM,ISPIN)
      MJ = IOFFOO(IRPJ,IRPMJ,ISPIN) + (J-1)*POP(IRPM,ISPIN) + M
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MJ,1,1,IRPMJ,LIST2OFF + ISPIN)
  520 CONTINUE
      ENDIF
C
      IF(IRPM.GT.IRPJ)THEN
      DO  530 M=1,POP(IRPM,ISPIN)
      MJ = IOFFOO(IRPM,IRPMJ,ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
      CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &            MJ,1,1,IRPMJ,LIST2OFF + ISPIN)
  530 CONTINUE
      CALL VMINUS(CORE(I010),POP(IRPM,ISPIN)*IRPDPD(IRPBC,ISPIN))
      ENDIF
C
      IF(IRPM.EQ.IRPJ)THEN
      DO  540 M=1,POP(IRPM,ISPIN)
         IF(M.LT.J)THEN
         MJ = IOFFOO(IRPJ,IRPMJ,ISPIN) + INDEX(J-1) + M
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MJ,1,1,IRPMJ,LIST2OFF + ISPIN)
         ENDIF
C
         IF(M.GT.J)THEN
         MJ = IOFFOO(IRPM,IRPMJ,ISPIN) + INDEX(M-1) + J
         CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               MJ,1,1,IRPMJ,LIST2OFF + ISPIN)
         CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
C
         IF(M.EQ.J)THEN
         CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPBC,ISPIN)),
     &               IRPDPD(IRPBC,ISPIN))
         ENDIF
  540 CONTINUE
C
      ENDIF
C
      IRPA   = DIRPRD(IRPM,IRPMA)
      CALL XGEMM('N','N',
     &           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),POP(IRPM,ISPIN),
     &           -1.0D+00,
     &           CORE(I010),IRPDPD(IRPBC,ISPIN),
     &           CORE(IADW(IRPA)),POP(IRPM,ISPIN),
     &            1.0D+00,
     &           Z(IADZ(IRPA)),IRPDPD(IRPBC,ISPIN))
  550 CONTINUE
C
      RETURN
 1010 FORMAT(' @GT3WT214-F, Insufficient memory to continue. ')
      END
