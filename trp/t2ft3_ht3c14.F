      SUBROUTINE T2FT3_HT3C14(W,CORE,IADW,IRPI,IRPJ,IRPK,IRPIJ,
     &                        IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,
     &                        ICLLVL,IMODE,LAMBDA)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W(1),CORE(1)
      LOGICAL LAMBDA
      DIMENSION IADW(8)
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
         LISTOFF2=59
      ELSE
         LISTOFF=0
         LISTOFF2=0
      ENDIF
      
C     ROUTINE TO COMPUTE INCLUSION OF F*T3 AAA OR F*T3 BBB IN T2 AA 
C     OR T2 BB.
C
C     W IS (C<D,A) REPRESENTATION OF T3 (A<C<D). NOTE THAT THIS
C     ROUTINE MUST BE CALLED BEFORE T2T314 SINCE THE LATTER TRANSPOSES
C     AND DESTROYS W.
C
C     IF(ISPIN.EQ.1)THEN
C     LIST2I = LIST2I1
C     ELSE
C     LIST2I = LIST2I2
C     ENDIF
C
C     --- TOTAL LENGTH OF INTERMEDIATE. ---
C
      IF(ISPIN.EQ.1)THEN
      LENF = NTAA
      ELSE
      LENF = NTBB
      ENDIF
C
C     --- ADDRESSES AND LENGTHS OF SYMMETRY BLOCKS OF INTERMEDIATE. ---
C
C     IADFVO = 1
      DO   10 IRPC=1,NIRREP
      LENFVO(IRPC) = VRT(IRPC,ISPIN) * POP(IRPC,ISPIN)
      IF(IRPC.EQ.1)THEN
      IADFVO(IRPC) = 1
      ELSE
      IADFVO(IRPC) = IADFVO(IRPC-1) + LENFVO(IRPC-1)
      ENDIF
   10 CONTINUE
C
C     --- Read Fock matrix/intermediate/T1 vector. ---
C
      CALL GETLST(CORE(IADFVO(1)),1,1,2,    ISPIN,93)
C
C      AB  AB                      ABC
C     D   T    =    SUM  F(C,K)   T       I<J<K   W(A<B,C) * F(C,K)
C      IJ  IJ       K,C            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPIJ,ISPIN)
C
      IRPC  = IRPK
      IRPAB = DIRPRD(IRPIJK,IRPC)
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(W(IADW(IRPC)),
     1            CORE(IADFVO(IRPC) + (K-1)*VRT(IRPC,ISPIN)),
     1            CORE(IADT2),LENT2,VRT(IRPC,ISPIN),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      IF(IRPIJ.NE.1)THEN
      IJ = (J-1)*POP(IRPI,ISPIN) + I
      ELSE
      IJ = INDEX(J-1) + I
      ENDIF
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN)+IJ,
     1            1,1,IRPIJ,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            IRPDPD(IRPIJ,ISPIN), 1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN)+IJ,
     1            1,1,IRPIJ,LIST2I)
C
C
C      AB  AB                      ABC
C     D   T    =  - SUM  F(C,J)   T       I<J<K   W(A<B,C) * F(C,J)
C      IK  IK       J,C            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPIK,ISPIN)
C
      IRPC  = IRPJ
      IRPAB = DIRPRD(IRPIJK,IRPC)
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(W(IADW(IRPC)),
     1            CORE(IADFVO(IRPC) + (J-1)*VRT(IRPC,ISPIN)),
     1            CORE(IADT2),LENT2,VRT(IRPC,ISPIN),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      IF(IRPIK.NE.1)THEN
      IK = (K-1)*POP(IRPI,ISPIN) + I
      ELSE
      IK = INDEX(K-1) + I
      ENDIF
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,ISPIN)+IK,
     1            1,1,IRPIK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            IRPDPD(IRPIK,ISPIN),-1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPIK,ISPIN)+IK,
     1            1,1,IRPIK,LIST2I)
C
C      AB  AB                      ABC
C     D   T    =    SUM  F(C,I)   T       I<J<K   W(A<B,C) * F(C,I)
C      JK  JK       I,C            IJK
C
C     --- SET ADDRESS FOR D2T2 MATRIX. ---
C
      IADT2 = IADFVO(1) + LENF
      LENT2 = IRPDPD(IRPJK,ISPIN)
C
      IRPC  = IRPI
      IRPAB = DIRPRD(IRPIJK,IRPC)
      CALL ZERO(CORE(IADT2),LENT2)
      CALL MATVEC(W(IADW(IRPC)),
     1            CORE(IADFVO(IRPC) + (I-1)*VRT(IRPC,ISPIN)),
     1            CORE(IADT2),LENT2,VRT(IRPC,ISPIN),0,0)
C
C     --- READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK ---
C     ---                      TO DISK.                            ---
C
      IADT2I = IADT2 + LENT2
      IF(IRPJK.NE.1)THEN
      JK = (K-1)*POP(IRPJ,ISPIN) + J
      ELSE
      JK = INDEX(K-1) + J
      ENDIF
      CALL GETLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,ISPIN)+JK,
     1            1,1,IRPJK,LIST2I)
      CALL   VADD(CORE(IADT2I),CORE(IADT2I),CORE(IADT2),
     1            IRPDPD(IRPJK,ISPIN), 1.0D+00)
      CALL PUTLST(CORE(IADT2I),IOFFOO(IRPK,IRPJK,ISPIN)+JK,
     1            1,1,IRPJK,LIST2I)
C
      RETURN
      END
