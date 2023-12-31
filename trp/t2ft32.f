      SUBROUTINE T2FT32(T3,W,CORE,IADT3,IADW,NLEFT,
     1                 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1                 I,J,K,ICLLVL,IUHF,IMODE,LAMBDA)
      IMPLICIT INTEGER (A-Z)
      LOGICAL LAMBDA
      DOUBLE PRECISION T3(1),W(1),CORE(1)
CSSS      DOUBLE PRECISION SCR1(20000),SCR2(20000),SCR3(20000)
      DIMENSION IADT3(8),IADW(8)
      DIMENSION IADFVO(8),LENFVO(8)
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
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
C
      INDEX(I) = I*(I-1)/2
C     List offset
      IF (LAMBDA) THEN
         LISTOFF=100
         LISTOFF2 = 59
      ELSE
         LISTOFF=0
         LISTOFF2=0
      ENDIF

C
C     ROUTINE TO COMPUTE INCLUSION OF F*T3 AAB IN T2 AA OR T2 AB.
C
      ISPIN1 = 1
      ISPIN2 = 2
C
      IF(IUHF.EQ.0) GOTO 100
C
C     --- AA D2T2 INCREMENTS ---
C
C      AB  AB                      ABC
C     D   T    =    SUM  F(C,K)   T          T3(A<B,C) * F(C,K)
C      IJ  IJ       K,C            IJK
C
C      AA  AA              B B     AAB
C
C     --- LIST NUMBERS FOR D2T2 INCREMENTS. ---
C
C     --- 1. CCSDT-n Energy Calculations. ---
C
      IF((ICLLVL.GE.13.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      LIST2I = 60 + ISPIN1
      ENDIF
C
C     --- 2. CCSD(T), QCISD(T) gradients. ---
C
      IF(ICLLVL.EQ.21. OR.ICLLVL.EQ.22)THEN
      IF(IMODE.EQ.1)THEN
      LIST2I = 113 + ISPIN1 + LISTOFF2
      ENDIF
      IF(IMODE.EQ.2)THEN
      LIST2I = LIST2I1 - 1 + ISPIN1
      ENDIF
      ENDIF
C
C     --- TOTAL LENGTH OF INTERMEDIATE. ---
C
      LENF = NTBB
C
C     --- ADDRESSES AND LENGTHS OF SYMMETRY BLOCKS OF INTERMEDIATE. ---
C
      DO   10 IRPC=1,NIRREP
      LENFVO(IRPC) = VRT(IRPC,ISPIN2) * POP(IRPC,ISPIN2)
      IF(IRPC.EQ.1)THEN
      IADFVO(IRPC) = 1
      ELSE
      IADFVO(IRPC) = IADFVO(IRPC-1) + LENFVO(IRPC-1)
      ENDIF
   10 CONTINUE
C
C     --- READ FOCK MATRIX/INTERMEDIATE. ---
C
      IF(ICLLVL.EQ.13)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,2 + ISPIN2,93)
      ENDIF
      IF((ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,    ISPIN2,93)
      ENDIF
      IF(ICLLVL.EQ.21. OR.ICLLVL.EQ.22)THEN
         LIST1ORL1 = 90+LISTOFF
      IF(IMODE.EQ.1)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,    ISPIN2,LIST1ORL1)
      ENDIF
      IF(IMODE.EQ.2)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,2 + ISPIN2,93)
      CALL SSCAL(LENF, 0.5D+00,CORE(IADFVO(1)),1)
      ENDIF
      ENDIF
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPIJ,ISPIN1)
C
      IRPC  = IRPK
      IRPAB = DIRPRD(IRPIJK,IRPC)
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(T3(IADT3(IRPC)),
     1            CORE(IADFVO(IRPC) + (K-1)*VRT(IRPC,ISPIN2)),
     1            CORE(IADT2),LENT2,VRT(IRPC,ISPIN2),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      IF(IRPIJ.NE.1)THEN
      IJ = (J-1)*POP(IRPI,ISPIN1) + I
      ELSE
      IJ = INDEX(J-1) + I
      ENDIF
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN1)+IJ,
     1            1,1,IRPIJ,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            IRPDPD(IRPIJ,ISPIN1), 1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN1)+IJ,
     1            1,1,IRPIJ,LIST2I)
C
  100 CONTINUE
C
C     --- AB D2T2 INCREMENTS ---
C
C     IF(ISPIN.EQ.1)THEN
C     LIST2I = LIST2I1
C     ELSE
C     LIST2I = LIST2I2
C     ENDIF
C
C     --- LIST NUMBERS FOR D2T2 INCREMENTS. ---
C
C     --- 1. CCSDT-n Energy Calculations. ---
C
      IF((ICLLVL.GE.13.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      LIST2I = 63
      ENDIF
C
C     --- 2. CCSD(T), QCISD(T) gradients. ---
C
      IF(ICLLVL.EQ.21. OR.ICLLVL.EQ.22)THEN
      IF(IMODE.EQ.1)THEN
      LIST2I = 116 + LISTOFF2
      ENDIF
      IF(IMODE.EQ.2)THEN
      LIST2I = LIST2I3
      ENDIF
      ENDIF
C
C     --- TOTAL LENGTH OF INTERMEDIATE. ---
C
      LENF = NTAA
C
C     --- ADDRESSES AND LENGTHS OF SYMMETRY BLOCKS OF INTERMEDIATE. ---
C
C     IADFVO = 1
      DO  110 IRPC=1,NIRREP
      LENFVO(IRPC) = VRT(IRPC,ISPIN1) * POP(IRPC,ISPIN1)
      IF(IRPC.EQ.1)THEN
      IADFVO(IRPC) = 1
      ELSE
      IADFVO(IRPC) = IADFVO(IRPC-1) + LENFVO(IRPC-1)
      ENDIF
  110 CONTINUE
C
C     --- READ FOCK MATRIX/INTERMEDIATE. ---
C
C     IADFVO = 1
      IF(ICLLVL.EQ.13)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,2 + ISPIN1,93)
      ENDIF
      IF((ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,    ISPIN1,93)
      ENDIF
      IF(ICLLVL.EQ.21. OR.ICLLVL.EQ.22)THEN
      IF(IMODE.EQ.1)THEN
         LIST1ORL1 = 90+LISTOFF
      CALL GETLST(CORE(IADFVO(1)),1,1,2,    ISPIN1,LIST1ORL1)
      ENDIF
      IF(IMODE.EQ.2)THEN
      CALL GETLST(CORE(IADFVO(1)),1,1,2,2 + ISPIN1,93)
      CALL SSCAL(LENF, 0.5D+00,CORE(IADFVO(1)),1)
      ENDIF
      ENDIF
C
C      BC  BC                      ABC
C     D   T    =  - SUM  F(A,J)   T          W(A,B,C) * F(A,J)
C      IK  IK       J,A            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPIK,13)
C
      IRPA  = IRPJ
      IRPBC = DIRPRD(IRPIJK,IRPA)
      CALL ZERO(CORE(IADT2),LENT2)
C     CALL VECMAT(CORE(IADFVO(IRPA) + (J-1)*VRT(IRPA,ISPIN1)),
C    1            W(IADW(IRPA)),CORE(IADT2),
C    1            VRT(IRPA,ISPIN1),LENT2,0,0)
      CALL MATVEC(W(IADW(IRPA)),
     1            CORE(IADFVO(IRPA) + (J-1)*VRT(IRPA,ISPIN1)),
     1            CORE(IADT2),
     1            LENT2,VRT(IRPA,ISPIN1),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      IK = (K-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            LENT2,-1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I)
C
C Extract the last address of the CORE and starting from there
C allocate space for scr arrays that are being used in SYMTR3 calls.
C
      ILAST_OF_CORE  = IADT2I + LENT2 
      IADDR_4_SCR1   = ILAST_OF_CORE + 1
      IADDR_4_SCR2   = IADDR_4_SCR1  + LENT2
      IADDR_4_SCR3   = IADDR_4_SCR2  + LENT2
      ILAST_OF_SCR   = IADDR_4_SCR3  + LENT2
      IF (ILAST_OF_SCR .GE. NLEFT) CALL INSMEM("@T2FT32",
     &                             ILAST_OF_SCR,NLEFT)
C
      IF(IUHF.EQ.0)THEN
      CALL SYMTR3(IRPIK,VRT(1,2),VRT(1,1),IRPDPD(IRPIK,13),1,
     1            CORE(IADT2),CORE(IADDR_4_SCR1),
     1            CORE(IADDR_4_SCR2),CORE(IADDR_4_SCR3))
      IK = (I-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPI,IRPIK,6)+IK,
     1            1,1,IRPIK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            LENT2,-1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPI,IRPIK,6)+IK,
     1            1,1,IRPIK,LIST2I)
      ENDIF
C
C      BC  BC                      ABC
C     D   T    =    SUM  F(A,I)   T          W(A,B,C) * F(A,I)
C      JK  JK       I,A            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPJK,13)
C
      IRPA  = IRPI
      IRPBC = DIRPRD(IRPIJK,IRPA)
      CALL ZERO(CORE(IADT2),LENT2)
C     CALL VECMAT(CORE(IADFVO(IRPA) + (I-1)*VRT(IRPA,ISPIN1)),
C    1            W(IADW(IRPA)),CORE(IADT2),
C    1            VRT(IRPA,ISPIN1),LENT2,0,0)
      CALL MATVEC(W(IADW(IRPA)),
     1            CORE(IADFVO(IRPA) + (I-1)*VRT(IRPA,ISPIN1)),
     1            CORE(IADT2),
     1            LENT2,VRT(IRPA,ISPIN1),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      JK = (K-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            LENT2, 1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I)
C
C Extract the last address of the CORE and starting from there
C allocate space for scr arrays that are being used in SYMTR3 calls.
C
      ILAST_OF_CORE  = IADT2I + LENT2
      IADDR_4_SCR1   = ILAST_OF_CORE + 1
      IADDR_4_SCR2   = IADDR_4_SCR1  + LENT2
      IADDR_4_SCR3   = IADDR_4_SCR2  + LENT2
      ILAST_OF_SCR   = IADDR_4_SCR3  + LENT2
      IF (ILAST_OF_SCR .GE. NLEFT) CALL INSMEM("@T2FT32",
     &                             ILAST_OF_SCR,NLEFT)
C
      IF(IUHF.EQ.0)THEN
      CALL SYMTR3(IRPJK,VRT(1,2),VRT(1,1),IRPDPD(IRPJK,13),1,
     1            CORE(IADT2),CORE(IADDR_4_SCR1),
     1            CORE(IADDR_4_SCR2),CORE(IADDR_4_SCR3))
      JK = (J-1)*POP(IRPK,ISPIN2) + K
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPJ,IRPJK,6)+JK,
     1            1,1,IRPJK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            LENT2, 1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPJ,IRPJK,6)+JK,
     1            1,1,IRPJK,LIST2I)
      ENDIF
      RETURN
      END
