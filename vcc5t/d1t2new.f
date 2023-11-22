      SUBROUTINE D1T2NEW(T3,W,
     1                T2IJAA,T2ABJK,T2ABIK,
     1                 ICORE,IUHF,
     1                 IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1                 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1                 LNVVIJ,LNVVJK,LNVVIK,
     1                 SCR1,SCR2,SCR3,MAXCOR,T2ABA,T2ABBJ,T2ABBI,
     1                 T2AAAJ,T2AAAI)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION T2IJAA(LNVVIJ,1),
     1                 T2ABJK(LNVVJK,1),T2ABIK(LNVVIK,1)
C
      DOUBLE PRECISION T3(1),W(1)
      DOUBLE PRECISION ICORE(1)
      DOUBLE PRECISION SCR1(1),SCR2(1),SCR3(1)
      DOUBLE PRECISION T2ABA(1),T2ABBJ(1),T2ABBI(1),T2AAAJ(1),
     1                                              T2AAAI(1)
      DIMENSION IADBLK(8),LENBLK(8)
      DIMENSION IADW(8)  ,LENW(8)
      DIMENSION IADT2(8) ,LENT2(8),
     1          IADV(8)  ,LENV(8)
C
C     LOCAL ARRAYS
C
      DIMENSION IADVT(8),IADTT(8)
C
C     END OF LOCAL ARRAYS
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
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
      LISWVA = LWIC15
      LISWVB = LWIC16
      LISWVC = LWIC17
      LISWVD = LWIC18
C
      LISWOA = LWIC11
      LISWOB = LWIC12
      LISWOC = LWIC13
      LISWOD = LWIC14
C
C      AE                BE
C     T   <BC//EK>  -   T   <AC//EK>
C      IJ                IJ
C
C      AA  AB  AB        AA  AB  AB
C      AA                AA
C
      DO   10 IRREP=1,NIRREP
      LENT2(IRREP) = VRT(IRREP,1) * VRT(DIRPRD(IRREP,IRPIJ),1)
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + 
     1               VRT(IRREP-1,1) * 
     1               VRT(DIRPRD(IRREP-1,IRPIJ),1)
      ENDIF
   10 CONTINUE
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = INDEX(J-1) + I
      ELSE
      IJ = (J-1)*POP(IRPI,1) + I
      ENDIF
C
      DO   20 IRPE=1,NIRREP
      IF(VRT(IRPE,1).EQ.0) GOTO 20
      IRPBC = DIRPRD(IRPE,IRPK)
      LENGTHV = IRPDPD(IRPBC,13) * VRT(IRPE,1)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T2',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEK = DIRPRD(IRPE,IRPK)
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPK,IRPEK,4)+(K-1)*VRT(IRPE,1)+1,
     1            VRT(IRPE,1),2,IRPEK,LISWVD)
C    1            VRT(IRPE,1),2,IRPEK,30)
C
C     AT THIS POINT V MATRIX IS ORDERED (BC,E), hence use transpose
C     facility of XGEMM.
C
      IRPA   = DIRPRD(IRPE,IRPIJ)
      CALL XGEMM('N','T',VRT(IRPA,1),IRPDPD(IRPBC,13),VRT(IRPE,1),
     1 1.0D+00,T2IJAA(IADT2(IRPE),IJ),VRT(IRPA,1),
     1         ICORE(1),IRPDPD(IRPBC,13),0.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
   20 CONTINUE
C
C        AE                BE
C     - T   <BC//IE>   +  T   <AC//IE>
C        JK                JK
C
C        AB  AB  AB        AB  AB  AB
C        AB                AB
C
      DO  110 IRREP=1,NIRREP
      LENT2(IRREP) = VRT(IRREP,2) * VRT(DIRPRD(IRREP,IRPJK),1)
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + 
     1               VRT(IRREP-1,2) * 
     1               VRT(DIRPRD(IRREP-1,IRPJK),1)
      ENDIF
  110 CONTINUE
C
      DO  120 IRPE=1,NIRREP
      IF(VRT(IRPE,2).EQ.0) GOTO 120
      IRPBC = DIRPRD(IRPE,IRPI)
      LENGTHV = IRPDPD(IRPBC,13) * VRT(IRPE,2)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T2',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEI = DIRPRD(IRPE,IRPI)
C
      IF(IUHF.EQ.1)THEN
C            Anything other than RHF closed-shell (IUHF = 1)
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPI,IRPEI,3) + (I-1)*VRT(IRPE,2) + 1,
     1            VRT(IRPE,2),2,IRPEI,LISWVC)
C    1            VRT(IRPE,2),2,IRPEI,29)
      ELSE
C            RHF CASE (IUHF = 0)
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPI,IRPEI,4) + (I-1)*VRT(IRPE,2) + 1,
     1            VRT(IRPE,2),2,IRPEI,LISWVD)
C    1            VRT(IRPE,2),2,IRPEI,30)
      CALL SYMTR3(IRPEI,VRT(1,1),VRT(1,2),IRPDPD(IRPEI,13),
     1            VRT(IRPE,2),ICORE(1),SCR1,SCR2,SCR3)
      ENDIF
C
      IRPA   = DIRPRD(IRPE,IRPJK)
      CALL XGEMM('N','T',VRT(IRPA,1),IRPDPD(IRPBC,13),VRT(IRPE,2),
     1        -1.0D+00,
     1         T2ABJK(IADT2(IRPE),(K-1)*POP(IRPJ,1)+J),VRT(IRPA,1),
     1         ICORE(1),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
  120 CONTINUE
C
C      AE                 BE
C     T   <BC//JE>   -   T   <AC//JE>
C      IK                 IK
C
C      AB  AB  AB         AB  AB  AB
C      AB                 AB
C
C
      DO  210 IRREP=1,NIRREP
      LENT2(IRREP) = VRT(IRREP,2) * VRT(DIRPRD(IRREP,IRPIK),1)
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + 
     1               VRT(IRREP-1,2) * 
     1               VRT(DIRPRD(IRREP-1,IRPIK),1)
      ENDIF
  210 CONTINUE
C
      DO  220 IRPE=1,NIRREP
      IF(VRT(IRPE,2).EQ.0) GOTO 220
      IRPBC = DIRPRD(IRPE,IRPJ)
      LENGTHV = IRPDPD(IRPBC,13) * VRT(IRPE,2)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T2',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEJ = DIRPRD(IRPE,IRPJ)
C
      IF(IUHF.EQ.1)THEN
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPJ,IRPEJ,3) + (J-1)*VRT(IRPE,2) + 1,
     1            VRT(IRPE,2),2,IRPEJ,LISWVC)
C    1            VRT(IRPE,2),2,IRPEJ,29)
      ELSE
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPJ,IRPEJ,4) + (J-1)*VRT(IRPE,2) + 1,
     1            VRT(IRPE,2),2,IRPEJ,LISWVD)
C    1            VRT(IRPE,2),2,IRPEJ,30)
      CALL SYMTR3(IRPEJ,VRT(1,1),VRT(1,2),IRPDPD(IRPEJ,13),
     1            VRT(IRPE,2),ICORE(1),SCR1,SCR2,SCR3)
      ENDIF
C
      IRPA   = DIRPRD(IRPE,IRPIK)
      CALL XGEMM('N','T',VRT(IRPA,1),IRPDPD(IRPBC,13),VRT(IRPE,2),
     1         1.0D+00,
     1         T2ABIK(IADT2(IRPE),(K-1)*POP(IRPI,1)+I),VRT(IRPA,1),
     1         ICORE(1),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
  220 CONTINUE
C
C               EC
C     <AB//EJ> T
C               IK
C
C      AA  AA   AB
C               AB
C
      DO  310 IRREP=1,NIRREP
      LENT2(IRREP) = VRT(IRREP,2) * VRT(DIRPRD(IRREP,IRPIK),1)
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + 
     1               VRT(IRREP-1,2) * 
     1               VRT(DIRPRD(IRREP-1,IRPIK),1)
      ENDIF
  310 CONTINUE
C
C     SET ADDRESSES ETC FOR V
C
      DO  320 IRPE=1,NIRREP
      IF(VRT(IRPE,1).EQ.0) GOTO 320
      IRPAB = DIRPRD(IRPE,IRPJ)
        IF(IUHF.EQ.0)THEN
         LENGTHV1 = IRPDPD(IRPAB, 1) * VRT(IRPE,1)
         LENGTHV2 = IRPDPD(IRPAB,13) * VRT(IRPE,1)
         LENGTHV  = LENGTHV1 + LENGTHV2
        ELSE
         LENGTHV  = IRPDPD(IRPAB, 1) * VRT(IRPE,1)
        ENDIF
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T2',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEJ = DIRPRD(IRPE,IRPJ)
        IF(IUHF.EQ.0)THEN
        CALL GETLST(ICORE(1 + LENGTHV1),
     1              IOFFVO(IRPJ,IRPEJ,1) + (J-1)*VRT(IRPE,1) + 1,
     1              VRT(IRPE,1),2,IRPEJ,LISWVD)
        CALL ASSYM3(IRPEJ,VRT(1,1),IRPDPD(IRPAB,1),IRPDPD(IRPAB,13),
     1              VRT(IRPE,1),ICORE(1),ICORE(1+LENGTHV1),
     1              IOFFVV,1)
        ELSE
        CALL GETLST(ICORE(1),
     1              IOFFVO(IRPJ,IRPEJ,1) + (J-1)*VRT(IRPE,1) + 1,
     1              VRT(IRPE,1),2,IRPEJ,LISWVA)
C    1              VRT(IRPE,1),2,IRPEJ,27)
        ENDIF
C
      IRPC   = DIRPRD(IRPE,IRPIK)
      CALL XGEMM('N','N',IRPDPD(IRPAB,1),VRT(IRPC,2),VRT(IRPE,1),
     1            1.0D+00,
     1            ICORE(1),IRPDPD(IRPAB,1),
     1            T2ABIK(IADT2(IRPC),(K-1)*POP(IRPI,1)+I),VRT(IRPE,1),
     1            0.0D+00,
     1            T3(IADBLK(IRPC)),IRPDPD(IRPAB,1))
C  
  320 CONTINUE
C
C                 EC
C     - <AB//EI> T
C                 JK
C
C        AA  AA   AB
C                 AB
C
      DO  410 IRREP=1,NIRREP
      LENT2(IRREP) = VRT(IRREP,2) * VRT(DIRPRD(IRREP,IRPJK),1)
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + 
     1               VRT(IRREP-1,2) * 
     1               VRT(DIRPRD(IRREP-1,IRPJK),1)
      ENDIF
  410 CONTINUE
C
      DO  420 IRPE=1,NIRREP
      IF(VRT(IRPE,1).EQ.0) GOTO 420
      IRPAB = DIRPRD(IRPE,IRPI)
        IF(IUHF.EQ.0)THEN
         LENGTHV1 = IRPDPD(IRPAB, 1) * VRT(IRPE,1)
         LENGTHV2 = IRPDPD(IRPAB,13) * VRT(IRPE,1)
         LENGTHV  = LENGTHV1 + LENGTHV2
        ELSE
         LENGTHV  = IRPDPD(IRPAB, 1) * VRT(IRPE,1)
        ENDIF
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T2',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEI = DIRPRD(IRPE,IRPI)
        IF(IUHF.EQ.0)THEN
        CALL GETLST(ICORE(1 + LENGTHV1),
     1              IOFFVO(IRPI,IRPEI,1) + (I-1)*VRT(IRPE,1) + 1,
     1              VRT(IRPE,1),2,IRPEI,LISWVD)
        CALL ASSYM3(IRPEI,VRT(1,1),IRPDPD(IRPAB,1),IRPDPD(IRPAB,13),
     1              VRT(IRPE,1),ICORE(1),ICORE(1+LENGTHV1),
     1              IOFFVV,1)
        ELSE
       CALL GETLST(ICORE(1),
     1            IOFFVO(IRPI,IRPEI,1) + (I-1)*VRT(IRPE,1) + 1,
     1            VRT(IRPE,1),2,IRPEI,LISWVA)
C    1            VRT(IRPE,1),2,IRPEI,27)
        ENDIF
C
      IRPC   = DIRPRD(IRPE,IRPJK)
      CALL XGEMM('N','N',IRPDPD(IRPAB,1),VRT(IRPC,2),VRT(IRPE,1),
     1           -1.0D+00,
     1            ICORE(1),IRPDPD(IRPAB,1),
     1            T2ABJK(IADT2(IRPC),(K-1)*POP(IRPJ,1)+J),VRT(IRPE,1),
     1            1.0D+00,
     1            T3(IADBLK(IRPC)),IRPDPD(IRPAB,1))
  420 CONTINUE
C
C               BC               AC
C     <IJ//MA> T    -  <IJ//MB> T
C               MK               MK
C
C      AA  AA   AB      AA  AA   AB
C               AB               AB
C
C
      IF(IUHF.EQ.0)THEN
C
c      IF(IRPI.EQ.IRPJ)THEN
c      IJ = IOFFOO(IRPJ,IRPIJ,1) + INDEX(J-1) + I
c      JI = IOFFOO(IRPI,IRPIJ,1) + INDEX(I-1) + J
c      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,5) + (J-1) * POP(IRPI,1) + I
      JI = IOFFOO(IRPI,IRPIJ,5) + (I-1) * POP(IRPJ,1) + J
c      ENDIF
      IOFFV  = 1
      IOFFV2 = IOFFV + IRPDPD(IRPIJ,9)
      CALL GETLST(ICORE(IOFFV ),IJ,1,2,IRPIJ,LWIC14)
      CALL GETLST(ICORE(IOFFV2),JI,1,2,IRPIJ,LWIC14)
      CALL VADD(ICORE(IOFFV),ICORE(IOFFV),ICORE(IOFFV2),
     1          IRPDPD(IRPIJ,9),-1.0D+00)
      ELSE
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,1) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,1) + (J-1) * POP(IRPI,1) + I
      ENDIF
      IOFFV  = 1
      CALL GETLST(ICORE(IOFFV ),IJ,1,2,IRPIJ,LWIC11)
      ENDIF
C
      DO  510 IRPM=1,NIRREP
      IRPBC = DIRPRD(IRPM,IRPK)
      LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,1)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
  510 CONTINUE
C
      DO  520 IRPA=1,NIRREP
C
      IRPM  = DIRPRD(IRPA,IRPIJ)
      IRPBC = DIRPRD(IRPM,IRPK)
C
      IF(POP(IRPM,1).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 520
C      LENGTHT = IRPDPD(IRPBC,13) * POP(IRPM,1)
CC
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
C
C      IRPMK = DIRPRD(IRPM,IRPK)
C      CALL GETLST(ICORE(1),
C     1            IOFFOO(IRPK,IRPMK,5) + (K-1)*POP(IRPM,1) + 1,
C     1            POP(IRPM,1),1,IRPMK,46)
C
C     T matrix is ordered (Bc,M), V matrix is (M,A).
C
      CALL XGEMM('T','T',VRT(IRPA,1),IRPDPD(IRPBC,13),POP(IRPM,1),
     1         1.0D+00,
     1         ICORE(IOFFV),POP(IRPM,1),
     1         T2ABA(IADT2(IRPM)),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
C
      IOFFV = IOFFV + POP(IRPM,1) * VRT(IRPA,1)
  520 CONTINUE
C
C
C               BC               AC
C     <JK//AM> T    -  <JK//BM> T
C               IM               IM
C
C      AB  AB   AB      AB  AB   AB
C               AB               AB
C
C     For IUHF=0, integrals are from list LWIC14 and are (m,A), and the
C     record is Kj rather than Jk.
C     For IUHF=1, integrals are from list LWIC13 and are (A,m).
C
      IF(IUHF.EQ.1)THEN
C
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,1) + J
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),JK,1,2,IRPJK,LWIC13)
C
      DO  620 IRPM=1,NIRREP
C
      IRPBC = DIRPRD(IRPM,IRPI)
      IRPMI = IRPBC
      LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
C
      IRPA = DIRPRD(IRPM,IRPJK)
      IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 620
C
C      LENGTHT = IRPDPD(IRPBC,13) * POP(IRPM,2)
C
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C      DO  610 M=1,POP(IRPM,2)
C      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,13)),
C     1            IOFFOO(IRPM,IRPMI,5) + (M-1)*POP(IRPI,1) + I,
C     1            1,1,IRPMI,46)
C  610 CONTINUE
C
      CALL XGEMM('N','T',VRT(IRPA,1),IRPDPD(IRPBC,13),POP(IRPM,2),
     1         1.0D+00,
     1         ICORE(IOFFV),VRT(IRPA,1),
     1         T2ABBI(IADT2(IRPM)),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
C
      IOFFV = IOFFV + VRT(IRPA,1) * POP(IRPM,2)
  620 CONTINUE
C
      ELSE
C
C     RHF case (IUHF=0)
C
      JK = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,2) + K
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),JK,1,2,IRPJK,LWIC14)
C
      DO  630 IRPM=1,NIRREP
      IRPBC = DIRPRD(IRPM,IRPI)
      LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
  630 CONTINUE
C
      DO  650 IRPA=1,NIRREP
C
      IRPM  = DIRPRD(IRPA,IRPJK)
      IRPBC = DIRPRD(IRPM,IRPI)
      IRPMI = IRPBC
C
      IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 650
C
C      LENGTHT = IRPDPD(IRPBC,13) * POP(IRPM,2)
C
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C      DO  640 M=1,POP(IRPM,2)
C      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,13)),
C     1            IOFFOO(IRPM,IRPMI,5) + (M-1)*POP(IRPI,1) + I,
C     1            1,1,IRPMI,46)
C  640 CONTINUE
C
      CALL XGEMM('T','T',VRT(IRPA,1),IRPDPD(IRPBC,13),POP(IRPM,2),
     1         1.0D+00,
     1         ICORE(IOFFV),POP(IRPM,2),
     1         T2ABBI(IADT2(IRPM)),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
C
      IOFFV = IOFFV + VRT(IRPA,1) * POP(IRPM,2)
  650 CONTINUE
C
      ENDIF
C
C                 BC               AC
C     - <IK//AM> T    +  <IK//BM> T
C                 JM               JM
C
C        AB  AB   AB      AB  AB   AB
C                 AB               AB
C
C     For IUHF=0, integrals are from list LWIC14 and are (m,A), and the
C     record is Ki rather than Ik.
C     For IUHF=1, integrals are from list LWIC13 and are (A,m).
C
      IF(IUHF.EQ.1)THEN
C
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,1) + I
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),IK,1,2,IRPIK,LWIC13)
C
      DO  720 IRPM=1,NIRREP
C
      IRPBC = DIRPRD(IRPM,IRPJ)
      IRPMJ = IRPBC
      LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
C
      IRPA = DIRPRD(IRPM,IRPIK)
      IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 720
C
C      LENGTHT = IRPDPD(IRPBC,13) * POP(IRPM,2)
C
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C      DO  710 M=1,POP(IRPM,2)
C      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,13)),
C     1            IOFFOO(IRPM,IRPMJ,5) + (M-1)*POP(IRPJ,1) + J,
C     1            1,1,IRPMJ,46)
C  710 CONTINUE
C
      CALL XGEMM('N','T',VRT(IRPA,1),IRPDPD(IRPBC,13),POP(IRPM,2),
     1        -1.0D+00,
     1         ICORE(IOFFV),VRT(IRPA,1),
     1         T2ABBJ(IADT2(IRPM)),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
C
      IOFFV = IOFFV + VRT(IRPA,1) * POP(IRPM,2)
  720 CONTINUE
C
      ELSE
C
C     RHF case (IUHF=0)
C
      IK = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,2) + K
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),IK,1,2,IRPIK,LWIC14)
C
      DO  730 IRPM=1,NIRREP
      IRPBC = DIRPRD(IRPM,IRPJ)
      LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
  730 CONTINUE
C
      DO  750 IRPA=1,NIRREP
C
      IRPM  = DIRPRD(IRPA,IRPIK)
      IRPBC = DIRPRD(IRPM,IRPJ)
      IRPMJ = IRPBC
C
      IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 750
C
C      LENGTHT = IRPDPD(IRPBC,13) * POP(IRPM,2)
C
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C      DO  740 M=1,POP(IRPM,2)
C      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,13)),
C     1            IOFFOO(IRPM,IRPMJ,5) + (M-1)*POP(IRPJ,1) + J,
C     1            1,1,IRPMJ,46)
C  740 CONTINUE
C
      CALL XGEMM('T','T',VRT(IRPA,1),IRPDPD(IRPBC,13),POP(IRPM,2),
     1        -1.0D+00,
     1         ICORE(IOFFV),POP(IRPM,2),
     1         T2ABBJ(IADT2(IRPM)),IRPDPD(IRPBC,13),1.0D+00,
     1         W(IADW(IRPBC)),VRT(IRPA,1))
C
      IOFFV = IOFFV + VRT(IRPA,1) * POP(IRPM,2)
  750 CONTINUE
C
      ENDIF
C
C        AB
C     - T   <JK//MC>
C        IM
C
C        AA  AB  AB
C        AA
C
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,1) + J
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),JK,1,2,IRPJK,LWIC14)
C
      DO  810 IRPM=1,NIRREP
      IRPAB = DIRPRD(IRPM,IRPI)
      LENT2(IRPM) = IRPDPD(IRPAB,1) * POP(IRPM,1)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
  810 CONTINUE
C
      DO  850 IRPC=1,NIRREP
C
      IRPM = DIRPRD(IRPC,IRPJK)
C
      IRPAB = DIRPRD(IRPM,IRPI)
      IRPMI = DIRPRD(IRPM,IRPI)
C
      IF(POP(IRPM,1).EQ.0.OR.VRT(IRPC,2).EQ.0) GOTO 850
C
C      LENGTHT = IRPDPD(IRPAB,1) * POP(IRPM,1)
C
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C         IF(IRPI.GT.IRPM)THEN
C         CALL GETLST(ICORE(1),
C     1               IOFFOO(IRPI,IRPMI,1) + (I-1)*POP(IRPM,1) + 1,
C     1               POP(IRPM,1),1,IRPMI,44)
C         CALL VMINUS(ICORE(1),IRPDPD(IRPMI,1)*POP(IRPM,1))
C         ENDIF
C         IF(IRPI.LT.IRPM)THEN
C         DO  820 M=1,POP(IRPM,1)
C         CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMI,1)),
C     1               IOFFOO(IRPM,IRPMI,1) + (M-1)*POP(IRPI,1) + I,
C     1               1,1,IRPMI,44)
C  820    CONTINUE
C         ENDIF
C         IF(IRPI.EQ.IRPM)THEN
C         DO  830 M=1,POP(IRPM,1)
C            IF(I.GT.M)THEN
C            IM = INDEX(I-1) + M
C            CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMI,1)),
C     1                  IOFFOO(IRPI,IRPMI,1) + IM,
C     1                  1,1,IRPMI,44)
C            CALL VMINUS(ICORE(1 + (M-1)*IRPDPD(IRPMI,1)),
C     1                  IRPDPD(IRPMI,1))
C            ENDIF
C            IF(I.LT.M)THEN
C            IM = INDEX(M-1) + I
C            CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMI,1)),
C     1                  IOFFOO(IRPM,IRPMI,1) + IM,
C     1                  1,1,IRPMI,44)
C            ENDIF
C            IF(M.EQ.I)THEN
C            CALL   ZERO(ICORE(1 + (M-1)*IRPDPD(IRPMI,1)),
C     1            IRPDPD(IRPMI,1))
C            ENDIF
C  830    CONTINUE
C         ENDIF
C
      CALL XGEMM('N','N',IRPDPD(IRPAB,1),VRT(IRPC,2),POP(IRPM,1),
     1           -1.0D+00,
     1            T2AAAI(IADT2(IRPM)),IRPDPD(IRPAB,1),
     1            ICORE(IOFFV),POP(IRPM,1),
     1            1.0D+00,
     1            T3(IADBLK(IRPC)),IRPDPD(IRPAB,1))
C
      IOFFV = IOFFV + POP(IRPM,1) * VRT(IRPC,2)
  850 CONTINUE
C
C      AB
C     T   <IK//MC>
C      JM
C
C      AA  AB  AB
C      AA
C
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,1) + I
      IOFFV = 1
      CALL GETLST(ICORE(IOFFV),IK,1,2,IRPIK,LWIC14)
C
      DO  910 IRPM=1,NIRREP
      IRPAB = DIRPRD(IRPM,IRPJ)
      LENT2(IRPM) = IRPDPD(IRPAB,1) * POP(IRPM,1)
      IF(IRPM.EQ.1)THEN
      IADT2(IRPM) = 1
      ELSE
      IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
      ENDIF
  910 CONTINUE
C
      DO  950 IRPC=1,NIRREP
C
      IRPM  = DIRPRD(IRPC,IRPIK)
C
      IRPAB = DIRPRD(IRPM,IRPJ)
      IRPMJ = DIRPRD(IRPM,IRPJ)
C
      IF(POP(IRPM,1).EQ.0.OR.VRT(IRPC,2).EQ.0) GOTO 950
C
C      LENGTHT = IRPDPD(IRPAB,1) * POP(IRPM,1)
CC
C      IF(LENGTHT.GT.MAXCOR)THEN
C      WRITE(6,1010)
C      CALL INSMEM('D1T2',LENGTHT,MAXCOR)
C      ENDIF
CC
C         IF(IRPJ.GT.IRPM)THEN
C         CALL GETLST(ICORE(1),
C     1               IOFFOO(IRPJ,IRPMJ,1) + (J-1)*POP(IRPM,1) + 1,
C     1               POP(IRPM,1),1,IRPMJ,44)
C         CALL VMINUS(ICORE(1),IRPDPD(IRPMJ,1)*POP(IRPM,1))
C         ENDIF
C         IF(IRPJ.LT.IRPM)THEN
C         DO  920 M=1,POP(IRPM,1)
C         CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMJ,1)),
C     1               IOFFOO(IRPM,IRPMJ,1) + (M-1)*POP(IRPJ,1) + J,
C     1               1,1,IRPMJ,44)
C  920    CONTINUE
C         ENDIF
C         IF(IRPJ.EQ.IRPM)THEN
C         DO  930 M=1,POP(IRPM,1)
C            IF(J.GT.M)THEN
C            JM = INDEX(J-1) + M
C            CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMJ,1)),
C     1                  IOFFOO(IRPJ,IRPMJ,1) + JM,
C     1                  1,1,IRPMJ,44)
C            CALL VMINUS(ICORE(1 + (M-1)*IRPDPD(IRPMJ,1)),
C     1                  IRPDPD(IRPMJ,1))
C            ENDIF
C            IF(J.LT.M)THEN
C            JM = INDEX(M-1) + J
C            CALL GETLST(ICORE(1 + (M-1)*IRPDPD(IRPMJ,1)),
C     1                  IOFFOO(IRPM,IRPMJ,1) + JM,
C     1                  1,1,IRPMJ,44)
C            ENDIF
C            IF(M.EQ.J)THEN
C            CALL   ZERO(ICORE(1 + (M-1)*IRPDPD(IRPMJ,1)),
C     1            IRPDPD(IRPMJ,1))
C            ENDIF
C  930    CONTINUE
C         ENDIF
CC
      CALL XGEMM('N','N',IRPDPD(IRPAB,1),VRT(IRPC,2),POP(IRPM,1),
     1            1.0D+00,
     1            T2AAAJ(IADT2(IRPM)),IRPDPD(IRPAB,1),
     1            ICORE(IOFFV),POP(IRPM,1),
     1            1.0D+00,
     1            T3(IADBLK(IRPC)),IRPDPD(IRPAB,1))
C
      IOFFV = IOFFV + POP(IRPM,1) * VRT(IRPC,2)
  950 CONTINUE
      RETURN
 1010 FORMAT(' @D1T2-I, Insufficient memory to continue. ')
      END
