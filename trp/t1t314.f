      SUBROUTINE T1T314(D1T1,S1,W,VJK,VIK,VIJ,
     1                 IADW,
     1                 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1                 I,J,K,IOFFOO,ISPIN,NONHF,LAMBDA)
      IMPLICIT INTEGER (A-Z)
      LOGICAL NONHF,LAMBDA
      DOUBLE PRECISION D1T1(1),S1(1)
      DOUBLE PRECISION W(1)
      DOUBLE PRECISION VJK(1),VIK(1),VIJ(1)
      DIMENSION IADW(8)
      DIMENSION IOFFOO(8,8,2)
      DIMENSION IOFFT1(8),LENT1(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
C
      INDEX(I) = I*(I-1)/2
C     List offset 
      IF (LAMBDA) THEN
         LISTOFF=100
      ELSE
         LISTOFF=0
      END IF
C
C     ROUTINE TO COMPUTE INCLUSION OF T3 AAA OR T3 BBB IN T1 A 
C     OR T1 B.
C
C     IRPI, IRPJ ARE SYMMETRY BLOCKS OF T1 A OR B.
C        I,    J ARE VALUES OF SPECIFIC ELEMENTS.
C
C     GET INTEGRALS
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      ELSE
      IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
      CALL GETLST(VIJ,IJ,1,2,IRPIJ,13+ISPIN)
C
      IF(IRPI.EQ.IRPK)THEN
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      ELSE
      IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
      CALL GETLST(VIK,IK,1,2,IRPIK,13+ISPIN)
C
      IF(IRPJ.EQ.IRPK)THEN
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + INDEX(K-1) + J
      ELSE
      JK = IOFFOO(IRPK,IRPJK,ISPIN) + (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
      CALL GETLST(VJK,JK,1,2,IRPJK,13+ISPIN)
C
C     COMPUTE T1 OFFSETS
C
      DO  30 IRREP=1,NIRREP
      LENT1(IRREP) = POP(IRREP,ISPIN) * VRT(IRREP,ISPIN)
      IF(IRREP.EQ.1)THEN
      IOFFT1(IRREP) = 0
      ELSE
      IOFFT1(IRREP) = IOFFT1(IRREP-1) + LENT1(IRREP-1)
      ENDIF
   30 CONTINUE
C
C      A  A                     ABC
C     D  T  =    SUM  <BC//JK> T        (W(BC,A))
C      I  I      B<C            IJK
C                J<K
C
      IRPA  = IRPI
      IRPBC = DIRPRD(IRPA,IRPIJK)
C      DO  110    A=1,VRT(IRPA,ISPIN)
C      DO  100   BC=1,IRPDPD(IRPBC,ISPIN)
CC
C       D1T1(IOFFT1(IRPA) + (I-1)*VRT(IRPA,ISPIN) + A) =
C     1 D1T1(IOFFT1(IRPA) + (I-1)*VRT(IRPA,ISPIN) + A) +
C     1 VJK(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
C     1                        - 1 + BC)
C  100 CONTINUE
C  110 CONTINUE
C
      CALL VECMAT(VJK,W(IADW(IRPA)),
     1            D1T1(IOFFT1(IRPA) + (I-1)*VRT(IRPA,ISPIN) + 1),
     1            IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),0,0)
C
C      A  A                     ABC
C     D  T  =  - SUM  <BC//IK> T        (W(BC,A))
C      J  J      B<C            IJK
C                I<K
C
      IRPA  = IRPJ
      IRPBC = DIRPRD(IRPA,IRPIJK)
C      DO  210    A=1,VRT(IRPA,ISPIN)
C      DO  200   BC=1,IRPDPD(IRPBC,ISPIN)
CC
C       D1T1(IOFFT1(IRPA) + (J-1)*VRT(IRPA,ISPIN) + A) =
C     1 D1T1(IOFFT1(IRPA) + (J-1)*VRT(IRPA,ISPIN) + A) -
C     1 VIK(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
C     1                        - 1 + BC)
C  200 CONTINUE
C  210 CONTINUE
      CALL VECMAT(VIK,W(IADW(IRPA)),
     1            D1T1(IOFFT1(IRPA) + (J-1)*VRT(IRPA,ISPIN) + 1),
     1            IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),0,1)
C
C      A  A                     ABC
C     D  T  =    SUM  <BC//IJ> T        (W(BC,A))
C      K  K      B<C            IJK
C                I<J
C
      IRPA  = IRPK
      IRPBC = DIRPRD(IRPA,IRPIJK)
C      DO  310    A=1,VRT(IRPA,ISPIN)
C      DO  300   BC=1,IRPDPD(IRPBC,ISPIN)
CC
C       D1T1(IOFFT1(IRPA) + (K-1)*VRT(IRPA,ISPIN) + A) =
C     1 D1T1(IOFFT1(IRPA) + (K-1)*VRT(IRPA,ISPIN) + A) +
C     1 VIJ(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
C     1                        - 1 + BC)
C  300 CONTINUE
C  310 CONTINUE
      CALL VECMAT(VIJ,W(IADW(IRPA)),
     1            D1T1(IOFFT1(IRPA) + (K-1)*VRT(IRPA,ISPIN) + 1),
     1            IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),0,0)
C
      IF(.NOT.NONHF) RETURN
C
C     TRY TO DO THE FUNNY FOURTH-ORDER PIECE ("IT WAS MY IDEA").
C
C     GET AMPLITUDES
C
      CALL GETLST(VIJ,IJ,1,2,IRPIJ,43+ISPIN+LISTOFF)
      CALL GETLST(VIK,IK,1,2,IRPIK,43+ISPIN+LISTOFF)
      CALL GETLST(VJK,JK,1,2,IRPJK,43+ISPIN+LISTOFF)
C
C        A             BC      ABC
C       S  =    SUM   T       T        (W(BC,A))
C        I      B<C    JK      IJK
C               J<K
C
      IRPA  = IRPI
      IRPBC = DIRPRD(IRPA,IRPIJK)
      DO  410    A=1,VRT(IRPA,ISPIN)
      DO  400   BC=1,IRPDPD(IRPBC,ISPIN)
C
       S1(IOFFT1(IRPA) + (I-1)*VRT(IRPA,ISPIN) + A) =
     1 S1(IOFFT1(IRPA) + (I-1)*VRT(IRPA,ISPIN) + A) +
     1 VJK(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                        - 1 + BC)
  400 CONTINUE
  410 CONTINUE
C
C        A              BC      ABC
C       S  =   - SUM   T       T        (W(BC,A))
C        J       B<C    IK      IJK
C                I<K
C
      IRPA  = IRPJ
      IRPBC = DIRPRD(IRPA,IRPIJK)
      DO  510    A=1,VRT(IRPA,ISPIN)
      DO  500   BC=1,IRPDPD(IRPBC,ISPIN)
C
       S1(IOFFT1(IRPA) + (J-1)*VRT(IRPA,ISPIN) + A) =
     1 S1(IOFFT1(IRPA) + (J-1)*VRT(IRPA,ISPIN) + A) -
     1 VIK(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                        - 1 + BC)
  500 CONTINUE
  510 CONTINUE
C
C        A              BC      ABC
C       S  =     SUM   T       T        (W(BC,A))
C        K       B<C    IJ      IJK
C                I<J
C
      IRPA  = IRPK
      IRPBC = DIRPRD(IRPA,IRPIJK)
      DO  610    A=1,VRT(IRPA,ISPIN)
      DO  600   BC=1,IRPDPD(IRPBC,ISPIN)
C
       S1(IOFFT1(IRPA) + (K-1)*VRT(IRPA,ISPIN) + A) =
     1 S1(IOFFT1(IRPA) + (K-1)*VRT(IRPA,ISPIN) + A) +
     1 VIJ(BC) * W(IADW(IRPA) + (A-1)*IRPDPD(IRPBC,ISPIN)
     1                        - 1 + BC)
  600 CONTINUE
  610 CONTINUE
      RETURN
      END
