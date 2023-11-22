      SUBROUTINE WOT323(T3,W,CORE,IADT3,IADW,
     1                  IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1                  I,J,K,ISPIN1,ISPIN2,IUHF,IMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZ,POP,VRT,DIRPRD
      DIMENSION T3(1),W(1),CORE(1)
C     DOUBLE PRECISION T2ABA(1),T2ABBJ(1),T2ABBI(1),
C    1                          T2AAAJ(1),T2AAAI(1)
      DIMENSION IADT3(8),IADW(8)
      DIMENSION IADT2(8),LENT2(8),IADG(8),LENG(8)
      DIMENSION LISTG(4),LISTT(3)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FILES/  LUOUT,MOINTS
C
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                                            LISGV3,LISGV4
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      INDEX(I) = I*(I-1)/2
C
C     Generic subroutine to calculate the contribution of AAB/BBA T3
C     amplitudes to either ijka gammas or Hbar intermediates. The spin
C     orbital equations are :
C
C     W(la,ij) =   Sum  t(abc,ijk) * <bc//kl>     (IMODE = 1)
C                  b<c
C                   k
C
C     G(la,ij) = - Sum  t(abc,ijk) * t(bc,kl)     (IMODE = 2)
C                  b<c
C                   k
C
      IF(IMODE.EQ.1)THEN
      SIGN  =  1.0D+00
      LISTG(1) =  LWIC11
      LISTG(2) =  LWIC12
      LISTG(3) =  LWIC13
      LISTG(4) =  LWIC14
      LISTT(1) =  14
      LISTT(2) =  15
      LISTT(3) =  16
      ENDIF
      IF(IMODE.EQ.2)THEN
      SIGN  = -1.0D+00
      LISTG(1) =  LISGO1
      LISTG(2) =  LISGO2
      LISTG(3) =  LISGO3
      LISTG(4) =  LISGO4
      LISTT(1) =  44
      LISTT(2) =  45
      LISTT(3) =  46
      ENDIF
      IF(IMODE.NE.1.AND.IMODE.NE.2)THEN
      WRITE(LUOUT,1010) IMODE
 1010 FORMAT(' @WOT323-F, Invalid value of IMODE. Given ',I10)
      STOP
      ENDIF
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1)*POP(IRPI,ISPIN1) + I
      ENDIF
      IK =    IOFFOO(IRPK,IRPIK,4+ISPIN1) + (K-1)*POP(IRPI,ISPIN1) + I
      JK =    IOFFOO(IRPK,IRPJK,4+ISPIN1) + (K-1)*POP(IRPJ,ISPIN1) + J
      IKRHF = IOFFOO(IRPI,IRPIK,6) + (I-1)*POP(IRPK,ISPIN2) + K
      JKRHF = IOFFOO(IRPJ,IRPJK,6) + (J-1)*POP(IRPK,ISPIN2) + K
C
C
C                        ABC  AB
C     G(IK,LC) =   SUM  T    T
C                  A<B   IJK  JL
C                   J
C       AB AB            AAB  AA
C
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,C), LABELLED BY IRPC.
C
      DO  120 IRPC=1,NIRREP
      IRPL = DIRPRD(IRPC,IRPIK)
      LENG(IRPC) = POP(IRPL,ISPIN1) * VRT(IRPC,ISPIN2)
      IF(IRPC.EQ.1)THEN
      IADG(IRPC) = 1
      ELSE
      IADG(IRPC) = IADG(IRPC-1) + LENG(IRPC-1)
      ENDIF
  120 CONTINUE
      CALL GETLST(CORE(1),IK,1,2,IRPIK,LISTG(5-ISPIN1))
C
C     SET ADDRESSES FOR T2.
C
      DO  130 IRPL=1,NIRREP
      IRPJL = DIRPRD(IRPL,IRPJ)
      IRPAB = IRPJL
      LENT2(IRPL) = IRPDPD(IRPAB,ISPIN1) * POP(IRPL,ISPIN1)
      IF(IRPL.EQ.1)THEN
C      IADT2(IRPL) = 1
      IADT2(IRPL) = IADG(NIRREP) + LENG(NIRREP)
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
  130 CONTINUE
C
C     READ T2(A<B,L).
C
      DO  170 IRPL=1,NIRREP
      IRPC = DIRPRD(IRPIK,IRPL)
      IF(POP(IRPL,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 170
      IRPJL = DIRPRD(IRPJ,IRPL)
C
      IF(IRPJ.LT.IRPL)THEN
      DO  140    L=1,POP(IRPL,ISPIN1)
      JL = IOFFOO(IRPL,IRPJL,ISPIN1) + (L-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,ISPIN1)),
     1            JL,1,1,IRPJL,LISTT(ISPIN1))
  140 CONTINUE
      ENDIF
C
      IF(IRPJ.GT.IRPL)THEN
      JL = IOFFOO(IRPJ,IRPJL,ISPIN1) + (J-1)*POP(IRPL,ISPIN1) + 1
      CALL GETLST(CORE(IADT2(IRPL)),
     1            JL,POP(IRPL,ISPIN1),1,IRPJL,LISTT(ISPIN1))
      CALL VMINUS(CORE(IADT2(IRPL)),LENT2(IRPL))
      ENDIF
C
      IF(IRPJ.EQ.IRPL)THEN
      DO  160    L=1,POP(IRPL,ISPIN1)
C
      IF(J.LT.L)THEN
      JL = IOFFOO(IRPL,IRPJL,ISPIN1) + INDEX(L-1) + J
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,ISPIN1)),
     1            JL,1,1,IRPJL,LISTT(ISPIN1))
      ENDIF
C
      IF(J.GT.L)THEN
      JL = IOFFOO(IRPJ,IRPJL,ISPIN1) + INDEX(J-1) + L
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,ISPIN1)),
     1            JL,1,1,IRPJL,LISTT(ISPIN1))
      CALL VMINUS(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,ISPIN1)),
     1            IRPDPD(IRPJL,ISPIN1))
      ENDIF
C
      IF(J.EQ.L)THEN
      CALL ZERO(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,ISPIN1)),
     1          IRPDPD(IRPJL,ISPIN1))
      ENDIF
  160 CONTINUE
C
      ENDIF
C
C     WE NOW HAVE T2 IN THE FORM T2(A<B,L) AND T3 IN THE FORM
C     T3(A<B,C). CONTRACT TOGETHER TO FORM G(L,C), WHICH CAN
C     BE SUMMED INTO GAMMA(LC | IK).
C
      DISSIZ = POP(IRPL,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = IRPDPD(IRPJL,ISPIN1)
      CALL XGEMM('T','N',DISSIZ,NDIS,NSUM, -SIGN,
     1           CORE(IADT2(IRPL)),NSUM,T3(IADT3(IRPC)),NSUM, 1.0D+00,
     1           CORE(IADG(IRPC)),DISSIZ)
  170 CONTINUE
      CALL PUTLST(CORE(1),IK,1,2,IRPIK,LISTG(5-ISPIN1))
C
C                        ABC  AB
C     G(JK,LC) = - SUM  T    T
C                  A<B   IJK  IL
C                   I
C       AB AB            AAB  AA
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,C), LABELLED BY IRPC.
C
      DO  220 IRPC=1,NIRREP
      IRPL = DIRPRD(IRPC,IRPJK)
      LENG(IRPC) = POP(IRPL,ISPIN1) * VRT(IRPC,ISPIN2)
      IF(IRPC.EQ.1)THEN
      IADG(IRPC) = 1
      ELSE
      IADG(IRPC) = IADG(IRPC-1) + LENG(IRPC-1)
      ENDIF
  220 CONTINUE
      CALL GETLST(CORE(1),JK,1,2,IRPJK,LISTG(5-ISPIN1))
C
C     SET ADDRESSES FOR T2.
C
      DO  230 IRPL=1,NIRREP
      IRPIL = DIRPRD(IRPL,IRPI)
      IRPAB = IRPIL
      LENT2(IRPL) = IRPDPD(IRPAB,ISPIN1) * POP(IRPL,ISPIN1)
      IF(IRPL.EQ.1)THEN
      IADT2(IRPL) = IADG(NIRREP) + LENG(NIRREP)
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
  230 CONTINUE
C
C     READ T2(A<B,L).
C
      DO  270 IRPL=1,NIRREP
      IRPC = DIRPRD(IRPJK,IRPL)
      IF(POP(IRPL,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 270
      IRPIL = DIRPRD(IRPI,IRPL)
C
      IF(IRPI.LT.IRPL)THEN
      DO  240    L=1,POP(IRPL,ISPIN1)
      IL = IOFFOO(IRPL,IRPIL,ISPIN1) + (L-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,ISPIN1)),
     1            IL,1,1,IRPIL,LISTT(ISPIN1))
  240 CONTINUE
      ENDIF
C
      IF(IRPI.GT.IRPL)THEN
      IL = IOFFOO(IRPI,IRPIL,ISPIN1) + (I-1)*POP(IRPL,ISPIN1) + 1
      CALL GETLST(CORE(IADT2(IRPL)),
     1            IL,POP(IRPL,ISPIN1),1,IRPIL,LISTT(ISPIN1))
      CALL VMINUS(CORE(IADT2(IRPL)),LENT2(IRPL))
      ENDIF
C
      IF(IRPI.EQ.IRPL)THEN
      DO  260    L=1,POP(IRPL,ISPIN1)
C
      IF(I.LT.L)THEN
      IL = IOFFOO(IRPL,IRPIL,ISPIN1) + INDEX(L-1) + I
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,ISPIN1)),
     1            IL,1,1,IRPIL,LISTT(ISPIN1))
      ENDIF
C
      IF(I.GT.L)THEN
      IL = IOFFOO(IRPI,IRPIL,ISPIN1) + INDEX(I-1) + L
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,ISPIN1)),
     1            IL,1,1,IRPIL,LISTT(ISPIN1))
      CALL VMINUS(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,ISPIN1)),
     1            IRPDPD(IRPIL,ISPIN1))
      ENDIF
C
      IF(I.EQ.L)THEN
      CALL ZERO(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,ISPIN1)),
     1          IRPDPD(IRPIL,ISPIN1))
      ENDIF
  260 CONTINUE
C
      ENDIF
C
C     WE NOW HAVE T2 IN THE FORM T2(A<B,L) AND T3 IN THE FORM
C     T3(A<B,C). CONTRACT TOGETHER TO FORM G(L,C), WHICH CAN
C     BE SUMMED INTO GAMMA(LC | JK).
C
      DISSIZ = POP(IRPL,ISPIN1)
      NDIS   = VRT(IRPC,ISPIN2)
      NSUM   = IRPDPD(IRPIL,ISPIN1)
      CALL XGEMM('T','N',DISSIZ,NDIS,NSUM,  SIGN,
     1           CORE(IADT2(IRPL)),NSUM,T3(IADT3(IRPC)),NSUM, 1.0D+00,
     1           CORE(IADG(IRPC)),DISSIZ)
  270 CONTINUE
      CALL PUTLST(CORE(1),JK,1,2,IRPJK,LISTG(5-ISPIN1))
C
C
CJ      IF(IUHF.GT.0) GOTO 1000
C
C
C                        ABC  BC
C     G(IK,LA) =   SUM  T    T
C                  B,C   IJK  JL
C                   J
C       AB BA            AAB  AB
C
C
C
C     NOTE : THERE IS A PHASE DILEMMA HERE. THUS, WHILE
C
C                     <IJ//KA>  = - <IJ/AK>
C                      AB  BA        AB AB
C
C            IT IS <IJ/AK> WHICH IS STORED. HENCE, TO BE CONSISTENT,
C            WE INTRODUCE A MINUS SIGN WHICH DOES NOT APPEAR IN THE
C            GAMMA EQUATIONS.
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,A), LABELLED BY IRPA.
C
CJ      DO  1320 IRPA=1,NIRREP
CJ      IRPL = DIRPRD(IRPA,IRPIK)
CJ      LENG(IRPA) = POP(IRPL,ISPIN2) * VRT(IRPA,ISPIN1)
CJ      IF(IRPA.EQ.1)THEN
CJ      IADG(IRPA) = 1
CJ      ELSE
CJ      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
CJ      ENDIF
CJ 1320 CONTINUE
C
C     SET ADDRESSES FOR T2.
C
CJ      DO  1330 IRPL=1,NIRREP
CJ      IRPJL = DIRPRD(IRPL,IRPJ)
CJ      IRPBC = IRPJL
CJ      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN2)
CJ      IF(IRPL.EQ.1)THEN
CJ      IADT2(IRPL) = 1
CJ      ELSE
CJ      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
CJ      ENDIF
CJ 1330 CONTINUE
C
C     READ T2(BC,L).
C
CJ      DO 1370 IRPL=1,NIRREP
CJ      IF(POP(IRPL,ISPIN2).EQ.0) GOTO 1370
CJ      IRPJL = DIRPRD(IRPJ,IRPL)
C
CCJ      DO  1340    L=1,POP(IRPL,ISPIN2)
CCJ      JL = IOFFOO(IRPL,IRPJL,5) + (L-1)*POP(IRPJ,ISPIN1) + J
CCJ      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,13)),
CCJ     1            JL,1,1,IRPJL,46)
CCJ 1340 CONTINUE
C
CJ      IRPA = DIRPRD(IRPIK,IRPL)
CJ      CALL XGEMM('T','N',
CJ     1           POP(IRPL,ISPIN2),VRT(IRPA,ISPIN1),IRPDPD(IRPJL,13),
CJC    1           -1.0D+00,CORE(IADT2(IRPL)),IRPDPD(IRPJL,13),
CJ     1           -1.0D+00,T2ABBJ(IADT2(IRPL)),IRPDPD(IRPJL,13),
CJ     1           W(IADW(IRPA)),IRPDPD(IRPJL,13),
CJ     1            1.00D+00,
CJ     1 GOOOV(OOOVAD(IRPIK,4)+(IKRHF-1)*IRPDPD(IRPIK,11) +
CJ     1         IADG(IRPA) - 1,IGPOS4),
CJ     1           POP(IRPL,ISPIN2))
CJ 1370 CONTINUE
C
C                        ABC  BC
C     G(JK,LA) =   SUM  T    T
C                  B,C   IJK  IL
C                   I
C       AB BA            AAB  AB
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,A), LABELLED BY IRPA.
C
CJ      DO 1420 IRPA=1,NIRREP
CJ      IRPL = DIRPRD(IRPA,IRPJK)
CJ      LENG(IRPA) = POP(IRPL,ISPIN2) * VRT(IRPA,ISPIN1)
CJ      IF(IRPA.EQ.1)THEN
CJ      IADG(IRPA) = 1
CJ      ELSE
CJ      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
CJ      ENDIF
CJ 1420 CONTINUE
C
C     SET ADDRESSES FOR T2.
C
CJ      DO  1430 IRPL=1,NIRREP
CJ      IRPIL = DIRPRD(IRPL,IRPI)
CJ      IRPBC = IRPIL
CJ      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN2)
CJ      IF(IRPL.EQ.1)THEN
CJ      IADT2(IRPL) = 1
CJ      ELSE
CJ      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
CJ      ENDIF
CJ 1430 CONTINUE
C
C     READ T2(BC,L).
C
CJ      DO 1470 IRPL=1,NIRREP
CJ      IF(POP(IRPL,ISPIN2).EQ.0) GOTO 1470
CJ      IRPIL = DIRPRD(IRPI,IRPL)
C
CCJ      DO  1440    L=1,POP(IRPL,ISPIN2)
CCJ      IL = IOFFOO(IRPL,IRPIL,5) + (L-1)*POP(IRPI,ISPIN1) + I
CCJ      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,13)),
CCJ     1            IL,1,1,IRPIL,46)
CCJ 1440 CONTINUE
CC
C     WE NOW HAVE T2 IN THE FORM T2(BC,L) AND T3 IN THE FORM
C     W(BC,A). CONTRACT TOGETHER TO FORM G(L,A), WHICH CAN
C     BE SUMMED INTO GAMMA(LA | JK).
C
CJ      IRPA = DIRPRD(IRPJK,IRPL)
CJ      CALL XGEMM('T','N',
CJ     1           POP(IRPL,ISPIN2),VRT(IRPA,ISPIN1),IRPDPD(IRPIL,13),
CJC    1            1.0D+00,CORE(IADT2(IRPL)),IRPDPD(IRPIL,13),
CJ     1            1.0D+00,T2ABBI(IADT2(IRPL)),IRPDPD(IRPIL,13),
CJ     1           W(IADW(IRPA)),IRPDPD(IRPIL,13),
CJ     1            1.00D+00,
CJ     1 GOOOV(OOOVAD(IRPJK,4)+(JKRHF-1)*IRPDPD(IRPJK,11) +
CJ     1         IADG(IRPA) - 1,IGPOS4),
CJ     1           POP(IRPL,ISPIN2))
CJ 1470 CONTINUE
C
c      IF(IUHF.EQ.0) RETURN
C
c 1000 CONTINUE
C
C                        ABC  BC
C     G(IK,LA) =   SUM  T    T
C                  B,C   IJK  JL
C                   J
C       AB BA            AAB  AB
C
C
C
C     NOTE : THERE IS A PHASE DILEMMA HERE. THUS, WHILE
C
C                     <IJ//KA>  = - <IJ/AK>
C                      AB  BA        AB AB
C
C            IT IS <IJ/AK> WHICH IS STORED. HENCE, TO BE CONSISTENT,
C            WE INTRODUCE A MINUS SIGN WHICH DOES NOT APPEAR IN THE
C            GAMMA EQUATIONS.
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,A), LABELLED BY IRPA.
C
      KI = IOFFOO(IRPI,IRPIK,7-ISPIN1) + (I-1)*POP(IRPK,ISPIN2) + K
c     CALL GETLST(CORE(1),KI,1,2,IRPIK,LISTG(2 + ISPIN1))
      CALL GETLST(CORE(1),KI,1,2,IRPIK,LISTG(2 + ISPIN1 + 1 - IUHF))
C
      DO  320 IRPA=1,NIRREP
      IRPL = DIRPRD(IRPA,IRPIK)
      LENG(IRPA) = POP(IRPL,ISPIN2) * VRT(IRPA,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  320 CONTINUE
C
C     SET ADDRESSES FOR T2.
C
      DO  330 IRPL=1,NIRREP
      IRPJL = DIRPRD(IRPL,IRPJ)
      IRPBC = IRPJL
      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN2)
      IF(IRPL.EQ.1)THEN
      IADT2(IRPL) = IADG(NIRREP) + LENG(NIRREP)
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
  330 CONTINUE
C
C     READ T2(BC,L).
C
      DO  370 IRPL=1,NIRREP
      IRPA = DIRPRD(IRPIK,IRPL)
      IF(POP(IRPL,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 370
      IRPJL = DIRPRD(IRPJ,IRPL)
C
      DO  340    L=1,POP(IRPL,ISPIN2)
      JL = IOFFOO(IRPL,IRPJL,4+ISPIN1) + (L-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPJL,13)),
     1            JL,1,1,IRPJL,LISTT(3))
  340 CONTINUE
C
C     WE NOW HAVE T2 IN THE FORM T2(BC,L) AND T3 IN THE FORM
C     W(BC,A). CONTRACT TOGETHER TO FORM G(L,A), WHICH CAN
C     BE SUMMED INTO GAMMA(LA | IK).
C
      DISSIZ = POP(IRPL,ISPIN2)
      NDIS   = VRT(IRPA,ISPIN1)
      NSUM   = IRPDPD(IRPJL,13)
      CALL XGEMM('T','N',DISSIZ,NDIS,NSUM, SIGN,
     1           CORE(IADT2(IRPL)),NSUM,W(IADW(IRPA)),NSUM,1.0D+00,
     1           CORE(IADG(IRPA)),DISSIZ)
  370 CONTINUE
c      CALL PUTLST(CORE(1),KI,1,2,IRPIK,LISTG(2 + ISPIN1))
      CALL PUTLST(CORE(1),KI,1,2,IRPIK,LISTG(2 + ISPIN1 + 1 - IUHF))
C
C                        ABC  BC
C     G(JK,LA) =   SUM  T    T
C                  B,C   IJK  IL
C                   I
C       AB BA            AAB  AB
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,A), LABELLED BY IRPA.
C
      KJ = IOFFOO(IRPJ,IRPJK,7-ISPIN1) + (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(1),KJ,1,2,IRPJK,LISTG(2 + ISPIN1 + 1 - IUHF))
C
      DO  420 IRPA=1,NIRREP
      IRPL = DIRPRD(IRPA,IRPJK)
      LENG(IRPA) = POP(IRPL,ISPIN2) * VRT(IRPA,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
  420 CONTINUE
C
C     SET ADDRESSES FOR T2.
C
      DO  430 IRPL=1,NIRREP
      IRPIL = DIRPRD(IRPL,IRPI)
      IRPBC = IRPIL
      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN2)
      IF(IRPL.EQ.1)THEN
C     IADT2(IRPL) = 1
      IADT2(IRPL) = IADG(NIRREP) + LENG(NIRREP)
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
  430 CONTINUE
C
C     READ T2(BC,L).
C
      DO  470 IRPL=1,NIRREP
      IRPA = DIRPRD(IRPJK,IRPL)
      IF(POP(IRPL,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 470
      IRPIL = DIRPRD(IRPI,IRPL)
C
      DO  440    L=1,POP(IRPL,ISPIN2)
      IL = IOFFOO(IRPL,IRPIL,4+ISPIN1) + (L-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(IADT2(IRPL) + (L-1)*IRPDPD(IRPIL,13)),
     1            IL,1,1,IRPIL,LISTT(3))
  440 CONTINUE
C
C     WE NOW HAVE T2 IN THE FORM T2(BC,L) AND T3 IN THE FORM
C     W(BC,A). CONTRACT TOGETHER TO FORM G(L,A), WHICH CAN
C     BE SUMMED INTO GAMMA(LA | JK).
C
      DISSIZ = POP(IRPL,ISPIN2)
      NDIS   = VRT(IRPA,ISPIN1)
      NSUM   = IRPDPD(IRPIL,13)
      CALL XGEMM('T','N',DISSIZ,NDIS,NSUM, -SIGN,
     1           CORE(IADT2(IRPL)),NSUM,W(IADW(IRPA)),NSUM, 1.0D+00,
     1           CORE(IADG(IRPA)),DISSIZ)
  470 CONTINUE
      CALL PUTLST(CORE(1),KJ,1,2,IRPJK,LISTG(2 + ISPIN1 + 1 - IUHF))
C
      IF(IUHF.EQ.0) RETURN
C
C                        ABC  BC
C     G(IJ,LA) =   SUM  T    T
C                  B<C   IJK  LK
C                   K
C       AA AA            AAB  AB
C
C     W IS (B,C,A) REPRESENTATION OF T3 (A<B,C).
C
C     SET ADDRESSES FOR GAMMA INCREMENTS (L,A), LABELLED BY IRPA.
C
      DO   20 IRPA=1,NIRREP
      IRPL = DIRPRD(IRPA,IRPIJ)
      LENG(IRPA) = POP(IRPL,ISPIN1) * VRT(IRPA,ISPIN1)
      IF(IRPA.EQ.1)THEN
      IADG(IRPA) = 1
      ELSE
      IADG(IRPA) = IADG(IRPA-1) + LENG(IRPA-1)
      ENDIF
   20 CONTINUE
      CALL GETLST(CORE(1),IJ,1,2,IRPIJ,LISTG(ISPIN1))
C
C     SET ADDRESSES FOR T2.
C
      DO   30 IRPL=1,NIRREP
      IRPKL = DIRPRD(IRPL,IRPK)
      IRPBC = IRPKL
      LENT2(IRPL) = IRPDPD(IRPBC,13) * POP(IRPL,ISPIN1)
      IF(IRPL.EQ.1)THEN
      IADT2(IRPL) = IADG(NIRREP) + LENG(NIRREP)
      ELSE
      IADT2(IRPL) = IADT2(IRPL-1) + LENT2(IRPL-1)
      ENDIF
   30 CONTINUE
C
C     READ T2(BC,L) FOR EACH IRPL.
C
      DO   70 IRPL=1,NIRREP
      IRPA = DIRPRD(IRPIJ,IRPL)
      IF(POP(IRPL,ISPIN1).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 70
      IRPKL = DIRPRD(IRPK,IRPL)
C
      KL = IOFFOO(IRPK,IRPKL,4+ISPIN1) + (K-1)*POP(IRPL,ISPIN1) + 1
      CALL GETLST(CORE(IADT2(IRPL)),
     1            KL,POP(IRPL,ISPIN1),1,IRPKL,LISTT(3))
C
C     WE NOW HAVE T2 IN THE FORM T2(BC,L) AND T3 IN THE FORM
C     T3(BC,A). CONTRACT TOGETHER TO FORM G(L,A), WHICH CAN
C     BE SUMMED INTO GAMMA(LA | IJ).
C
      DISSIZ = POP(IRPL,ISPIN1)
      NDIS   = VRT(IRPA,ISPIN1)
      NSUM   = IRPDPD(IRPKL,13)
      CALL XGEMM('T','N',DISSIZ,NDIS,NSUM, -SIGN,
     1           CORE(IADT2(IRPL)),NSUM,W(IADW(IRPA)),NSUM, 1.0D+00,
     1           CORE(IADG(IRPA)),DISSIZ)
   70 CONTINUE
      CALL PUTLST(CORE(1),IJ,1,2,IRPIJ,LISTG(ISPIN1))
C
      RETURN
      END
