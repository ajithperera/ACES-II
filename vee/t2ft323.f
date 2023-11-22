      SUBROUTINE T2FT323(T3,W,CORE,IADT3,IADW,
     &                   IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &                   I,J,K,IUHF,IRREPX,LZOFF,LISTF,LFOFF,SCFACF,
     &                   AAB,BBA)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION T3(1),W(1),CORE(1),SCFACF
      INTEGER IADT3,IADW,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,IUHF,
     &        IRREPX,LZOFF,LISTF,LFOFF
      LOGICAL AAB,BBA
C-----------------------------------------------------------------------
C     Common blocks.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,
     &        IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IOFFVV,IOFFOO,IOFFVO
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION ZILCH,ONE,ONEM
      INTEGER IADFVO,LENFVO,IADT2,LENT2,IADT2I,LENF,ISPIN1,ISPIN2,
     &        IJ,JK,IK,IRPA,IRPC,IRPAB,IRPBC,ISCR1,ISCR2,ISCR3,IEND,
     &        INDEX
C-----------------------------------------------------------------------
C
      DIMENSION IADT3(8),IADW(8)
      DIMENSION IADFVO(8),LENFVO(8)
C-----------------------------------------------------------------------
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C-----------------------------------------------------------------------
C
      DATA ZILCH,ONE,ONEM / 0.0D+00, 1.0D+00, -1.0D+00/
C-----------------------------------------------------------------------
C
      INDEX(I) = I*(I-1)/2
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Generic routine to compute
C
C     Z(ab,ij) = Z(ab,ij) + F_{kc} * T(abc,ijk)
C
C     F is assumed to be overall symmetry 1, while Z and T have symmetry
C     IRREPX.
C
C     LZOFF  : LZOFF+1, LZOFF+2, and LZOFF+3 are lists for Z2 increment.
C     LISTF  : List for F_{kc} quantity (storage (c,k)).
C     LFOFF  : Offset on LISTF, i.e. alpha F is part LFOFF+1, beta F is
C              on part LFOFF+2.
C     SCFACF : Scale factor for F. Usually 1.
C     AAB    : If this is .TRUE. it means we are doing AAB triples.
C     BBA    : If this is .TRUE. it means we are doing BBA triples.
C-----------------------------------------------------------------------
C
      IF((AAB.AND.BBA) .OR. (.NOT.AAB .AND. .NOT.BBA))THEN
        WRITE(6,*) ' @T2FT323-F, Incompatible values of AAB, BBA. '
        CALL ERREX
      ENDIF
C
      IF(AAB)THEN
        ISPIN1 = 1
        ISPIN2 = 2
      ENDIF
C
      IF(BBA)THEN
        ISPIN1 = 2
        ISPIN2 = 1
      ENDIF
C
      IF(IUHF.NE.0)THEN
C
C       Z(AB,IJ) = SUM_{ck} F(c,k) * T(A<B,c;I<J,k)
C
C     --- TOTAL LENGTH OF F. ---
C
        LENF = NTBB
C
C     --- ADDRESSES AND LENGTHS OF SYMMETRY BLOCKS OF F. ---
C
        DO   10 IRPC=1,NIRREP
        LENFVO(IRPC) = VRT(IRPC,ISPIN2) * POP(IRPC,ISPIN2)
        IF(IRPC.EQ.1)THEN
          IADFVO(IRPC) = 1
        ELSE
          IADFVO(IRPC) = IADFVO(IRPC-1) + LENFVO(IRPC-1)
        ENDIF
   10   CONTINUE
C
C     --- Read F and scale it. ---
C
        CALL GETLST(CORE(IADFVO(1)),1,1,2,LFOFF + ISPIN2,LISTF)
        CALL SSCAL(LENF,SCFACF,CORE(IADFVO(1)),1)
C
C     --- Set address for Z. ---
C
        IRPC  = IRPK
        IRPAB = DIRPRD(IRREPX,IRPIJ)
        IADT2 = IADFVO(1) + LENF
        LENT2 = IRPDPD(IRPAB,ISPIN1)
C
        CALL ZERO(CORE(IADT2),LENT2)
        CALL MATVEC(T3(IADT3(IRPC)),
     &              CORE(IADFVO(IRPC) + (K-1)*VRT(IRPC,ISPIN2)),
     &              CORE(IADT2),LENT2,VRT(IRPC,ISPIN2),0,0)
C
        IADT2I = IADT2 + LENT2
        IF(IRPIJ.NE.1)THEN
          IJ = (J-1)*POP(IRPI,ISPIN1) + I
        ELSE
          IJ = INDEX(J-1) + I
        ENDIF
        CALL GETLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN1)+IJ,
     &              1,1,IRPIJ,LZOFF + ISPIN1)
        CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),LENT2,ONE)
        CALL PUTLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN1)+IJ,
     &              1,1,IRPIJ,LZOFF + ISPIN1)
C
      ENDIF
C
C-----------------------------------------------------------------------
C     Z(Ab,Ij).
C-----------------------------------------------------------------------
C
C     --- TOTAL LENGTH OF F. ---
C
      LENF = NTAA
C
C     --- ADDRESSES AND LENGTHS OF SYMMETRY BLOCKS OF F. ---
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
C     --- Read and scale F. ---
C
      CALL GETLST(CORE(IADFVO(1)),1,1,2,LFOFF + ISPIN1,LISTF)
      CALL SSCAL(LENF,SCFACF,CORE(IADFVO(1)),1)
C
C      BC  BC                      ABC
C     D   T    =  - SUM  F(A,J)   T          W(A,B,C) * F(A,J)
C      IK  IK       J,A            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IRPA  = IRPJ
      IRPBC = DIRPRD(IRREPX,IRPIK)
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPBC,13)
C
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(W(IADW(IRPA)),
     &            CORE(IADFVO(IRPA) + (J-1)*VRT(IRPA,ISPIN1)),
     &            CORE(IADT2),
     &            LENT2,VRT(IRPA,ISPIN1),0,0)
C
      IADT2I = IADT2  + LENT2
      ISCR1  = IADT2I + LENT2
      ISCR2  = ISCR1  + LENT2
      ISCR3  = ISCR2  + LENT2
      IEND   = ISCR3  + LENT2
C
      IF(AAB)THEN
        IK = (K-1)*POP(IRPI,ISPIN1) + I
        CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,5)+IK,
     &              1,1,IRPIK,LZOFF + 3)
        CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),LENT2,ONEM)
        CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,5)+IK,
     &              1,1,IRPIK,LZOFF + 3)
      ENDIF
C
      IF(BBA .OR. IUHF.EQ.0)THEN
        IF(IUHF.EQ.0)THEN
          CALL SYMTR3(IRPBC,VRT(1,2),VRT(1,1),IRPDPD(IRPBC,13),1,
     &                CORE(IADT2),CORE(ISCR1),CORE(ISCR2),CORE(ISCR3))
        ENDIF
        IK = (I-1)*POP(IRPK,ISPIN2) + K
        CALL GETLST(CORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     &              1,1,IRPIK,LZOFF + 3)
        CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),LENT2,ONEM)
        CALL PUTLST(CORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     &              1,1,IRPIK,LZOFF + 3)
      ENDIF
C
C      BC  BC                      ABC
C     D   T    =    SUM  F(A,I)   T          W(A,B,C) * F(A,I)
C      JK  JK       I,A            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IRPA  = IRPI
      IRPBC = DIRPRD(IRREPX,IRPJK)
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPBC,13)
C
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(W(IADW(IRPA)),
     &            CORE(IADFVO(IRPA) + (I-1)*VRT(IRPA,ISPIN1)),
     &            CORE(IADT2),
     &            LENT2,VRT(IRPA,ISPIN1),0,0)
C
      IADT2I = IADT2  + LENT2
      ISCR1  = IADT2I + LENT2
      ISCR2  = ISCR1  + LENT2
      ISCR3  = ISCR2  + LENT2
      IEND   = ISCR3  + LENT2
C
      IF(AAB)THEN
        JK = (K-1)*POP(IRPJ,ISPIN1) + J
        CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,5)+JK,
     &              1,1,IRPJK,LZOFF + 3)
        CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),LENT2,ONE)
        CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,5)+JK,
     &              1,1,IRPJK,LZOFF + 3)
      ENDIF
C
      IF(BBA .OR. IUHF.EQ.0)THEN
        IF(IUHF.EQ.0)THEN
          CALL SYMTR3(IRPBC,VRT(1,2),VRT(1,1),IRPDPD(IRPBC,13),1,
     &                CORE(IADT2),CORE(ISCR1),CORE(ISCR2),CORE(ISCR3))
        ENDIF
        JK = (J-1)*POP(IRPK,ISPIN2) + K
        CALL GETLST(CORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     &            1,1,IRPJK,LZOFF + 3)
        CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),LENT2,ONE)
        CALL PUTLST(CORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     &              1,1,IRPJK,LZOFF + 3)
      ENDIF
      RETURN
      END
