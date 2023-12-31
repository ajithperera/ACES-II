      SUBROUTINE GVT32(T3,W,ICORE,
     1                 IADT3,IADW,
     1                 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1                 I,J,K,SCR1,SCR2,SCR3,IUHF,
     1                 T2IJAA,T2JKAB,T2IKAB,
     1                 LNVVIJ,LNVVJK,LNVVIK)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TOL
      DOUBLE PRECISION T3(1),W(1),ICORE(1),SCR1(1),SCR2(1),SCR3(1)
      DOUBLE PRECISION T2IJAA(LNVVIJ,1),
     1                 T2JKAB(LNVVJK,1),T2IKAB(LNVVIK,1)
      DIMENSION IADT3(8),IADW(8)
      DIMENSION IADT2(8),LENT2(8)
      DIMENSION IADT2E(8),LENT2E(8),IADG(8),LENG(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
C
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
      ISPIN1 = 1
      ISPIN2 = 2
C
      IF(IUHF.EQ.0) GOTO 100
C
C     ROUTINE TO COMPUTE CONTRIBUTION OF T3 AAB TO G(AI,BC).
C
C                        BCD  AD
C     G(AI,BC) = - SUM  T    T
C                  J<K   IJK  JK
C                   D
C       AA AA            AAB  AB
C
C                        BCD  AD
C     G(AJ,BC) =   SUM  T    T
C                  I<K   IJK  IK
C                   D
C       AA AA            AAB  AB
C
C
      JK =                        (K-1)*POP(IRPJ,1) + J
C
      DO   10 IRPD=1,NIRREP
      IRPA = DIRPRD(IRPD,IRPJK)
      LENT2(IRPD) = VRT(IRPA,ISPIN1) * VRT(IRPD,ISPIN2)
      IF(IRPD.EQ.1)THEN
      IADT2(IRPD) = 1
      ELSE
      IADT2(IRPD) = IADT2(IRPD-1) + LENT2(IRPD-1)
      ENDIF
   10 CONTINUE
C
C     SET ADDRESSES FOR GAMMA.
C
      DO   20 IRPA=1,NIRREP
      IRPBC = DIRPRD(IRPA,IRPI)
      LENG(IRPA) = VRT(IRPA,ISPIN1) * IRPDPD(IRPBC,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
   20 CONTINUE
C
      DO   30 IRPD=1,NIRREP
      IRPA  = DIRPRD(IRPD,IRPJK)
      IF(VRT(IRPD,ISPIN2).EQ.0.OR.LENG(IRPA).EQ.0) GOTO  30
      IRPBC = DIRPRD(IRPA,IRPI)
      IRPAI = DIRPRD(IRPA,IRPI)
      AI = IOFFVO(IRPI,IRPAI,ISPIN1) + (I-1)*VRT(IRPA,ISPIN1) + 1
      CALL GETLST(ICORE(IADG(IRPA)),
     1            AI,VRT(IRPA,ISPIN1),2,IRPAI,LISGV1)
      CALL XGEMM('N','T',IRPDPD(IRPBC,ISPIN1),VRT(IRPA,ISPIN1),
     1           VRT(IRPD,ISPIN2),
     1           -1.0D+00,T3(IADT3(IRPD)),IRPDPD(IRPBC,ISPIN1),
     1            T2JKAB(IADT2(IRPD),JK),VRT(IRPA,ISPIN1),
     1            1.0D+00,ICORE(IADG(IRPA)),IRPDPD(IRPBC,ISPIN1))
      CALL PUTLST(ICORE(IADG(IRPA)),
     1            AI,VRT(IRPA,ISPIN1),2,IRPAI,LISGV1)
   30 CONTINUE
C
      IK =                        (K-1)*POP(IRPI,1) + I
C
      DO  110 IRPD=1,NIRREP
      IRPA = DIRPRD(IRPD,IRPIK)
      LENT2(IRPD) = VRT(IRPA,ISPIN1) * VRT(IRPD,ISPIN2)
      IF(IRPD.EQ.1)THEN
      IADT2(IRPD) = 1
      ELSE
      IADT2(IRPD) = IADT2(IRPD-1) + LENT2(IRPD-1)
      ENDIF
  110 CONTINUE
C
C     SET ADDRESSES FOR GAMMA.
C
      DO  120 IRPA=1,NIRREP
      IRPBC = DIRPRD(IRPA,IRPJ)
      LENG(IRPA) = VRT(IRPA,ISPIN1) * IRPDPD(IRPBC,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  120 CONTINUE
C
      DO  130 IRPD=1,NIRREP
      IRPA  = DIRPRD(IRPD,IRPIK)
      IF(VRT(IRPD,ISPIN2).EQ.0.OR.LENG(IRPA).EQ.0) GOTO 130
      IRPBC = DIRPRD(IRPA,IRPJ)
      IRPAJ = DIRPRD(IRPA,IRPJ)
      AJ = IOFFVO(IRPJ,IRPAJ,ISPIN1) + (J-1)*VRT(IRPA,ISPIN1) + 1
      CALL GETLST(ICORE(IADG(IRPA)),
     1            AJ,VRT(IRPA,ISPIN1),2,IRPAJ,LISGV1)
      CALL XGEMM('N','T',IRPDPD(IRPBC,ISPIN1),VRT(IRPA,ISPIN1),
     1           VRT(IRPD,ISPIN2),
     1            1.0D+00,T3(IADT3(IRPD)),IRPDPD(IRPBC,ISPIN1),
     1            T2IKAB(IADT2(IRPD),IK),VRT(IRPA,ISPIN1),
     1            1.0D+00,ICORE(IADG(IRPA)),IRPDPD(IRPBC,ISPIN1))
      CALL PUTLST(ICORE(IADG(IRPA)),
     1            AJ,VRT(IRPA,ISPIN1),2,IRPAJ,LISGV1)
  130 CONTINUE
C
  100 CONTINUE
C
C                        DBC  DA
C     G(AI,BC) =   SUM  T    T
C                  J<K   IJK  JK
C                   D
C       BA AB            AAB  AB
C                        AAB  AB
C
C                        DBC  DA
C     G(AJ,BC) = - SUM  T    T
C                  I<K   IJK  IK
C                   D
C       BA AB            AAB  AB
C                        AAB  AB
C
C     W IS W(B,C,D) WITH IRREP LABEL D. W IS GENERATED FROM A "PARTIAL"
C     SYMTRW1 : T3(D<B,C) ---> W(B,C,D).
C
C
C     NOTE : THERE IS A PHASE DILEMMA HERE. THUS, WHILE
C
C                     <AB//CI>  = - <AB/IC>
C                      AB  BA        AB AB
C
C            IT IS <AB/IC> WHICH IS STORED. HENCE, TO BE CONSISTENT,
C            WE INTRODUCE A MINUS SIGN WHICH DOES NOT APPEAR IN THE
C            GAMMA EQUATIONS.
C
C
      JK =                        (K-1)*POP(IRPJ,1) + J
C
C     SET T2 ADDRESSES.
C
      DO   210 IRPA=1,NIRREP
      IRPD = DIRPRD(IRPA,IRPJK)
      LENT2(IRPA) = VRT(IRPA,ISPIN2) * VRT(IRPD,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADT2(IRPA) = 1
      ELSE
      IADT2(IRPA) = IADT2(IRPA-1) + LENT2(IRPA-1)
      ENDIF
  210 CONTINUE
C
C     SET GAMMA ADDRESSES.
C
      LENGTHG = 0
      DO  220 IRPA=1,NIRREP
      IRPBC = DIRPRD(IRPA,IRPI)
      LENG(IRPA) = VRT(IRPA,ISPIN2) * IRPDPD(IRPBC,13)
      LENGTHG = LENGTHG + LENG(IRPA)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  220 CONTINUE
C
C     CONTRACT TO GET CONTRIBUTION TO GAMMA.
C
      CALL ZERO(ICORE,LENGTHG)
      DO  230 IRPD=1,NIRREP
      IRPA  = DIRPRD(IRPD,IRPJK)
      IF(VRT(IRPD,ISPIN1).EQ.0.OR.LENG(IRPA).EQ.0) GOTO 230
      IRPBC = DIRPRD(IRPA,IRPI)
      CALL XGEMM('N','N',IRPDPD(IRPBC,13),VRT(IRPA,ISPIN2),
     1                                    VRT(IRPD,ISPIN1),
     1           -1.0D+00,
     1            W(IADW(IRPD)),IRPDPD(IRPBC,13),
     1            T2JKAB(IADT2(IRPA),JK),VRT(IRPD,ISPIN1),
     1            0.0D+00,ICORE(IADG(IRPA)),IRPDPD(IRPBC,13))
  230 CONTINUE
C
C     WRITE GAMMA TO LIST.
C
      IF(IUHF.EQ.0)THEN
      DO  240 IRPA=1,NIRREP
      IRPAI = DIRPRD(IRPA,IRPI)
      IF(LENG(IRPA).EQ.0) GOTO 240
C
      CALL SYMTR3(IRPAI,VRT(1,2),VRT(1,1),IRPDPD(IRPAI,13),
     1            VRT(IRPA,ISPIN2),ICORE(IADG(IRPA)),
     1              SCR1,SCR2,SCR3)
C
      AI = IOFFVO(IRPI,IRPAI,3) + (I-1)*VRT(IRPA,ISPIN2) + 1
      CALL GETLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AI,VRT(IRPA,ISPIN2),2,IRPAI,LISGV4)
C
      CALL VADD(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(IRPA)),LENG(IRPA), 1.0D+00)
C
      CALL PUTLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AI,VRT(IRPA,ISPIN2),2,IRPAI,LISGV4)
  240 CONTINUE
      ELSE
      DO  250 IRPA=1,NIRREP
      IRPAI = DIRPRD(IRPA,IRPI)
      IF(LENG(IRPA).EQ.0) GOTO 250
      AI = IOFFVO(IRPI,IRPAI,3) + (I-1)*VRT(IRPA,ISPIN2) + 1
      CALL GETLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AI,VRT(IRPA,ISPIN2),2,IRPAI,LISGV3)
C
      CALL VADD(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(IRPA)),LENG(IRPA), 1.0D+00)
C
      CALL PUTLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AI,VRT(IRPA,ISPIN2),2,IRPAI,LISGV3)
  250 CONTINUE
      ENDIF
C
      IK =                        (K-1)*POP(IRPI,1) + I
C
C     SET T2 ADDRESSES.
C
      DO   310 IRPA=1,NIRREP
      IRPD = DIRPRD(IRPA,IRPIK)
      LENT2(IRPA) = VRT(IRPA,ISPIN2) * VRT(IRPD,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADT2(IRPA) = 1
      ELSE
      IADT2(IRPA) = IADT2(IRPA-1) + LENT2(IRPA-1)
      ENDIF
  310 CONTINUE
C
C     SET GAMMA ADDRESSES.
C
      LENGTH = 0
      DO  320 IRPA=1,NIRREP
      IRPBC = DIRPRD(IRPA,IRPJ)
      LENG(IRPA) = VRT(IRPA,ISPIN2) * IRPDPD(IRPBC,13)
      LENGTHG = LENGTHG + LENG(IRPA)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  320 CONTINUE
C
C     CONTRACT TO GET CONTRIBUTION TO GAMMA.
C
      CALL ZERO(ICORE,LENGTHG)
      DO  330 IRPD=1,NIRREP
      IRPA  = DIRPRD(IRPD,IRPIK)
      IF(VRT(IRPD,ISPIN1).EQ.0.OR.LENG(IRPA).EQ.0) GOTO 330
      IRPBC = DIRPRD(IRPA,IRPJ)
      CALL XGEMM('N','N',IRPDPD(IRPBC,13),VRT(IRPA,ISPIN2),
     1                                    VRT(IRPD,ISPIN1),
     1            1.0D+00,
     1            W(IADW(IRPD)),IRPDPD(IRPBC,13),
     1            T2IKAB(IADT2(IRPA),IK),VRT(IRPD,ISPIN1),
     1            0.0D+00,ICORE(IADG(IRPA)),IRPDPD(IRPBC,13))
  330 CONTINUE
C
C     WRITE GAMMA TO LIST.
C
      IF(IUHF.EQ.0)THEN
      DO  340 IRPA=1,NIRREP
      IRPAJ = DIRPRD(IRPA,IRPJ)
      IF(LENG(IRPA).EQ.0) GOTO 340
C
      CALL SYMTR3(IRPAJ,VRT(1,2),VRT(1,1),IRPDPD(IRPAJ,13),
     1            VRT(IRPA,ISPIN2),ICORE(IADG(IRPA)),
     1              SCR1,SCR2,SCR3)
C
      AJ = IOFFVO(IRPJ,IRPAJ,3) + (J-1)*VRT(IRPA,ISPIN2) + 1
      CALL GETLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AJ,VRT(IRPA,ISPIN2),2,IRPAJ,LISGV4)
C
      CALL VADD(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(IRPA)),LENG(IRPA), 1.0D+00)
C
      CALL PUTLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AJ,VRT(IRPA,ISPIN2),2,IRPAJ,LISGV4)
  340 CONTINUE
      ELSE
      DO  350 IRPA=1,NIRREP
      IRPAJ = DIRPRD(IRPA,IRPJ)
      IF(LENG(IRPA).EQ.0) GOTO 350
      AJ = IOFFVO(IRPJ,IRPAJ,3) + (J-1)*VRT(IRPA,ISPIN2) + 1
      CALL GETLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AJ,VRT(IRPA,ISPIN2),2,IRPAJ,LISGV3)
C
      CALL VADD(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1          ICORE(IADG(IRPA)),LENG(IRPA), 1.0D+00)
C
      CALL PUTLST(ICORE(IADG(NIRREP) + LENG(NIRREP)),
     1            AJ,VRT(IRPA,ISPIN2),2,IRPAJ,LISGV3)
  350 CONTINUE
      ENDIF
C
C                        DBC  DA
C     G(AK,BC) =   SUM  T    T
C                  I<J   IJK  IJ
C                   D
C       AB AB            AAB  AA
C                        AAB  AA
C
C     SET ADDRESSES FOR EXPANDED T2(D,A).
C
      DO  410 IRPA=1,NIRREP
      IRPD = DIRPRD(IRPA,IRPIJ)
      LENT2E(IRPA) = VRT(IRPA,ISPIN1)*VRT(IRPD,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADT2E(IRPA) = 1
      ELSE
      IADT2E(IRPA) = IADT2E(IRPA-1) + LENT2E(IRPA-1)
      ENDIF
  410 CONTINUE
C
      IF(IRPIJ.EQ.1)THEN
      IJ =                             INDEX(J-1) + I
      ELSE
      IJ =                             (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
C
C     SET GAMMA ADDRESSES.
C
      DO  420 IRPA=1,NIRREP
      IRPBC = DIRPRD(IRPA,IRPK)
      LENG(IRPA) = VRT(IRPA,ISPIN1) * IRPDPD(IRPBC,13)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  420 CONTINUE
C
C     CONTRACT TO GET CONTRIBUTION TO GAMMA.
C
      DO  430 IRPD=1,NIRREP
      IRPA  = DIRPRD(IRPD,IRPIJ)
      IF(VRT(IRPD,ISPIN1).EQ.0.OR.LENG(IRPA).EQ.0) GOTO 430
      IRPBC = DIRPRD(IRPA,IRPK)
      IRPAK = DIRPRD(IRPA,IRPK)
      AK = IOFFVO(IRPK,IRPAK,4) + (K-1)*VRT(IRPA,ISPIN1) + 1
      CALL GETLST(ICORE(IADG(IRPA)),
     1            AK,VRT(IRPA,ISPIN1),2,IRPAK,LISGV4)
      CALL XGEMM('N','N',IRPDPD(IRPBC,13),VRT(IRPA,ISPIN1),
     1           VRT(IRPD,ISPIN1),
     1            1.0D+00,W(IADW(IRPD)),IRPDPD(IRPBC,13),
     1           T2IJAA(IADT2E(IRPA),IJ),VRT(IRPD,ISPIN1),
     1            1.0D+00,
     1            ICORE(IADG(IRPA)),IRPDPD(IRPBC,13))
      CALL PUTLST(ICORE(IADG(IRPA)),
     1            AK,VRT(IRPA,ISPIN1),2,IRPAK,LISGV4)
  430 CONTINUE
      RETURN
      END
