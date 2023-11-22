      SUBROUTINE E_T2T314(W,ICORE,
     1                    IADW,
     1                    IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1                    IRPACD,I,J,K,ISPIN,IRREPX)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W(1),ICORE(1)
      INTEGER DIRPRD,POP,VRT
      DIMENSION IADW(8)
      DIMENSION IADT2E(8),LENT2E(8),IADV(8),LENV(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     1                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
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
      IF(ISPIN.EQ.1)THEN
      LIST2I = LIST2I1
      LISTWI = LWIC25
      ELSE
      LIST2I = LIST2I2
      LISTWI = LWIC26
      ENDIF
C
C     ROUTINE TO COMPUTE INCLUSION OF T3 AAA OR T3 BBB IN T2 AA 
C     OR T2 BB.
C
C     TRANSPOSE W. W IS (C<D,A) REPRESENTATION OF T3 (A<C<D).
C
      DO   10 IRPA=1,NIRREP
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO 10
      IRPCD = DIRPRD(IRPA,IRPACD)
      CALL TRANSP(W(IADW(IRPA)),ICORE(IADW(IRPA)),
     1            VRT(IRPA,ISPIN),IRPDPD(IRPCD,ISPIN))
   10 CONTINUE
C
C     PUT TRANSPOSE IN W LOCATION (IE W IS LOST).
C
      DO   20 IRPA=1,NIRREP
      IF(VRT(IRPA,ISPIN).EQ.0) GOTO 20
      IRPCD = DIRPRD(IRPA,IRPACD)
C      CALL ICOPY(ICORE(IADW(IRPA)),W(IADW(IRPA)),
C     1           IINTFP*VRT(IRPA,ISPIN)*IRPDPD(IRPCD,ISPIN))
      CALL DCOPY(VRT(IRPA,ISPIN)*IRPDPD(IRPCD,ISPIN), ICORE(IADW(IRPA)),
     &           1, W(IADW(IRPA)), 1)

   20 CONTINUE
C
C      AB  AB                      ACD
C     D   T    =    SUM  <BK//CD> T       I<J<K   W(A,C<D) * V(C<D,B)
C      IJ  IJ       C<D            IJK
C                    K
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO   30 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPK)
      LENV(IRPB) = VRT(IRPB,ISPIN)*IRPDPD(IRPCD,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADV(IRPB) = 1
      ELSE
      IADV(IRPB) = IADV(IRPB-1) + LENV(IRPB-1)
      ENDIF
   30 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES
C
      DO   40 IRPB=1,NIRREP
      IF(VRT(IRPB,ISPIN).EQ.0) GOTO 40
      IRPBK = DIRPRD(IRPB,IRPK)
      CALL GETLST(ICORE(IADV(IRPB)),
     1     IOFFVO(IRPK,IRPBK,ISPIN)+(K-1)*VRT(IRPB,ISPIN)+1,
     1     VRT(IRPB,ISPIN),1,IRPBK,LISTWI)
   40 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      IRPAB   = DIRPRD(IRREPX,IRPIJ)
      DO   50 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPAB)
      LENT2E(IRPB) = VRT(IRPA,ISPIN) * VRT(IRPB,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
   50 CONTINUE
C
C     FORM EXPANDED T2 MATRICES. (A,B)
C
      DO   60 IRPCD=1,NIRREP
      IRPA = DIRPRD(IRPCD,IRPACD)
      IRPB = DIRPRD(IRPCD,IRPK)
      IF(VRT(IRPA,ISPIN).EQ.0.OR.VRT(IRPB,ISPIN).EQ.0) GOTO 60
      CALL XGEMM('N','N',VRT(IRPA,ISPIN),VRT(IRPB,ISPIN),
     1           IRPDPD(IRPCD,ISPIN), 1.0D+00,
     1           W(IADW(IRPA)),VRT(IRPA,ISPIN),
     1           ICORE(IADV(IRPB)),IRPDPD(IRPCD,ISPIN),
     1           0.0D+00,ICORE(IADT2E(IRPB)),VRT(IRPA,ISPIN))
   60 CONTINUE
C
C     ANTISYMMETRIZE THE EXPANDED T2 TO GIVE THE T2 INCREMENT.
C     (A<B) = (A,B) - (B,A).
C
      IADT2 = IADT2E(NIRREP) + LENT2E(NIRREP)
      CALL ZERO(ICORE(IADT2),IRPDPD(IRPAB,ISPIN))
C
      CALL EXPVV(ICORE(IADT2),ICORE(IADT2E(1)),IADT2E,
     1           IOFFVV(1,1,ISPIN),IRPAB,ISPIN)
C
C     READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK TO
C     DISK.
C
      IADT2I = IADT2 + IRPDPD(IRPAB,ISPIN)
      IF(IRPIJ.NE.1)THEN
      IJ = (J-1)*POP(IRPI,ISPIN) + I
      ELSE
      IJ = INDEX(J-1) + I
      ENDIF
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN)+IJ,
     1            1,1,IRPIJ,LIST2I)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2),
     1            IRPDPD(IRPAB,ISPIN), 1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPIJ,ISPIN)+IJ,
     1            1,1,IRPIJ,LIST2I)
C
C      AB  AB                      ACD
C     D   T    =  - SUM  <BJ//CD> T       I<J<K   W(A,C<D) * V(C<D,B)
C      IK  IK       C<D            IJK
C                    J
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  130 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPJ)
      LENV(IRPB) = VRT(IRPB,ISPIN)*IRPDPD(IRPCD,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADV(IRPB) = 1
      ELSE
      IADV(IRPB) = IADV(IRPB-1) + LENV(IRPB-1)
      ENDIF
  130 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES
C
      DO  140 IRPB=1,NIRREP
      IF(VRT(IRPB,ISPIN).EQ.0) GOTO 140
      IRPBJ = DIRPRD(IRPB,IRPJ)
      CALL GETLST(ICORE(IADV(IRPB)),
     1     IOFFVO(IRPJ,IRPBJ,ISPIN)+(J-1)*VRT(IRPB,ISPIN)+1,
     1     VRT(IRPB,ISPIN),1,IRPBJ,LISTWI)
  140 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      IRPAB   = DIRPRD(IRREPX,IRPIK)
      DO  150 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPAB)
      LENT2E(IRPB) = VRT(IRPA,ISPIN) * VRT(IRPB,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  150 CONTINUE
C
C     FORM EXPANDED T2 MATRICES. (A,B)
C
      DO  160 IRPCD=1,NIRREP
      IRPA = DIRPRD(IRPCD,IRPACD)
      IRPB = DIRPRD(IRPCD,IRPJ)
      IF(VRT(IRPA,ISPIN).EQ.0.OR.VRT(IRPB,ISPIN).EQ.0) GOTO 160
      CALL XGEMM('N','N',VRT(IRPA,ISPIN),VRT(IRPB,ISPIN),
     1           IRPDPD(IRPCD,ISPIN), 1.0D+00,
     1           W(IADW(IRPA)),VRT(IRPA,ISPIN),
     1           ICORE(IADV(IRPB)),IRPDPD(IRPCD,ISPIN),
     1           0.0D+00,ICORE(IADT2E(IRPB)),VRT(IRPA,ISPIN))
  160 CONTINUE
C
C     ANTISYMMETRIZE THE EXPANDED T2 TO GIVE THE T2 INCREMENT.
C     (A<B) = (A,B) - (B,A).
C
      IADT2 = IADT2E(NIRREP) + LENT2E(NIRREP)
      CALL ZERO(ICORE(IADT2),IRPDPD(IRPAB,ISPIN))
      CALL EXPVV(ICORE(IADT2),ICORE(IADT2E(1)),IADT2E,
     1           IOFFVV(1,1,ISPIN),IRPAB,ISPIN)
C
C     READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK TO
C     DISK.
C
      IADT2I = IADT2 + IRPDPD(IRPAB,ISPIN)
      IF(IRPIK.NE.1)THEN
      IK = (K-1)*POP(IRPI,ISPIN) + I
      ELSE
      IK = INDEX(K-1) + I
      ENDIF
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPK,IRPIK,ISPIN)+IK,
     1            1,1,IRPIK,LIST2I)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2),
     1            IRPDPD(IRPAB,ISPIN),-1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPK,IRPIK,ISPIN)+IK,
     1            1,1,IRPIK,LIST2I)
C
C      AB  AB                      ACD
C     D   T    =    SUM  <BI//CD> T       I<J<K   W(A,C<D) * V(C<D,B)
C      JK  JK       C<D            IJK
C                    I
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  230 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPI)
      LENV(IRPB) = VRT(IRPB,ISPIN)*IRPDPD(IRPCD,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADV(IRPB) = 1
      ELSE
      IADV(IRPB) = IADV(IRPB-1) + LENV(IRPB-1)
      ENDIF
  230 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES
C
      DO  240 IRPB=1,NIRREP
      IF(VRT(IRPB,ISPIN).EQ.0) GOTO 240
      IRPBI = DIRPRD(IRPB,IRPI)
      CALL GETLST(ICORE(IADV(IRPB)),
     1     IOFFVO(IRPI,IRPBI,ISPIN)+(I-1)*VRT(IRPB,ISPIN)+1,
     1     VRT(IRPB,ISPIN),1,IRPBI,LISTWI)
  240 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      IRPAB   = DIRPRD(IRREPX,IRPJK)
      DO  250 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPAB)
      LENT2E(IRPB) = VRT(IRPA,ISPIN) * VRT(IRPB,ISPIN)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  250 CONTINUE
C
C     FORM EXPANDED T2 MATRICES. (A,B)
C
      DO  260 IRPCD=1,NIRREP
      IRPA = DIRPRD(IRPCD,IRPACD)
      IRPB = DIRPRD(IRPCD,IRPI)
      IF(VRT(IRPA,ISPIN).EQ.0.OR.VRT(IRPB,ISPIN).EQ.0) GOTO 260
      CALL XGEMM('N','N',VRT(IRPA,ISPIN),VRT(IRPB,ISPIN),
     1           IRPDPD(IRPCD,ISPIN), 1.0D+00,
     1           W(IADW(IRPA)),VRT(IRPA,ISPIN),
     1           ICORE(IADV(IRPB)),IRPDPD(IRPCD,ISPIN),
     1           0.0D+00,ICORE(IADT2E(IRPB)),VRT(IRPA,ISPIN))
  260 CONTINUE
C
C     ANTISYMMETRIZE THE EXPANDED T2 TO GIVE THE T2 INCREMENT.
C     (A<B) = (A,B) - (B,A).
C
      IADT2 = IADT2E(NIRREP) + LENT2E(NIRREP)
      CALL ZERO(ICORE(IADT2),IRPDPD(IRPAB,ISPIN))
      CALL EXPVV(ICORE(IADT2),ICORE(IADT2E(1)),IADT2E,
     1           IOFFVV(1,1,ISPIN),IRPAB,ISPIN)
C
C     READ T2 INCREMENT, ADD CURRENT PIECE, AND WRITE SUM BACK TO
C     DISK.
C
      IADT2I = IADT2 + IRPDPD(IRPAB,ISPIN)
      IF(IRPJK.NE.1)THEN
      JK = (K-1)*POP(IRPJ,ISPIN) + J
      ELSE
      JK = INDEX(K-1) + J
      ENDIF
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPK,IRPJK,ISPIN)+JK,
     1            1,1,IRPJK,LIST2I)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2),
     1            IRPDPD(IRPAB,ISPIN), 1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPK,IRPJK,ISPIN)+JK,
     1            1,1,IRPJK,LIST2I)
      RETURN
      END
