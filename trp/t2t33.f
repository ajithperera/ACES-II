      SUBROUTINE T2T33(T3,W,ICORE,
     1                 IADT3,IADW,
     1                 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1                 I,J,K,LAMBDA)
      IMPLICIT INTEGER (A-Z)
      LOGICAL LAMBDA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION T3(1),W(1),ICORE(1)
      DIMENSION IADW(8),IADT3(8)
      DIMENSION IADT2E(8),LENT2E(8),IADV(8),LENV(8)
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
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
C
      INDEX(I) = I*(I-1)/2
C     List offset
      IF (LAMBDA) THEN
         LISTOFF=200
      ELSE
         LISTOFF=0
      ENDIF
C
C     ROUTINE TO COMPUTE INCLUSION OF T3 BBA IN T2 BB AND T2 AB
C
C      AB  AB                      ACD
C     D   T    =    SUM  <BK//CD> T       I<J<K   W(A,CD) * V(CD,B)
C      IJ  IJ       C<D            IJK
C                    K
C
C      BB  BB             BA  BA   BBA
C      BB  BB                      BBA
C
C
C      NOTE : W ARRAY IS W(A,CD), LABELLED BY IRPA.
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO   30 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPK)
      LENV(IRPB) = VRT(IRPB,2)*IRPDPD(IRPCD,13)
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
      IF(VRT(IRPB,2).EQ.0) GOTO 40
      IRPBK = DIRPRD(IRPB,IRPK)
      CALL GETLST(ICORE(IADV(IRPB)),
     1     IOFFVO(IRPK,IRPBK,3)+(K-1)*VRT(IRPB,2)+1,
     1     VRT(IRPB,2),2,IRPBK,LWIC27)
C    1     VRT(IRPB,2),2,IRPBK,29)
   40 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      DO   50 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPIJ)
      LENT2E(IRPB) = VRT(IRPA,2) * VRT(IRPB,2)
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
      IRPA = DIRPRD(IRPCD,IRPIJK)
      IRPB = DIRPRD(IRPCD,IRPK)
      IF(VRT(IRPA,2).EQ.0.OR.VRT(IRPB,2).EQ.0) GOTO 60
      CALL XGEMM('T','N',VRT(IRPA,2),VRT(IRPB,2),IRPDPD(IRPCD,13),
     1           1.0D+00,
     1           W(IADW(IRPA)),IRPDPD(IRPCD,13),
     1           ICORE(IADV(IRPB)),IRPDPD(IRPCD,13),
     1           0.0D+00,ICORE(IADT2E(IRPB)),VRT(IRPA,2))
   60 CONTINUE
C
C     ANTISYMMETRIZE THE EXPANDED T2 TO GIVE THE T2 INCREMENT.
C     (A<B) = (A,B) - (B,A).
C     LOOKS LIKE WE NEED A NEW ROUTINE FOR THIS. MODEL AFTER
C     EXPSC2. INDEED, WE SHOULD BE ABLE TO HACK UP EXPSC2 AND
C     IT'S DEPENDENTS !
C
      IADT2 = IADT2E(NIRREP) + LENT2E(NIRREP)
      CALL ZERO(ICORE(IADT2),IRPDPD(IRPIJ,2))
      CALL EXPVV(ICORE(IADT2),ICORE(IADT2E(1)),IADT2E,
     1           IOFFVV(1,1,2),IRPIJ,2)
C
C     READ CURRENT T2 INCREMENT, ADD THIS PIECE TO IT, AND WRITE SUM
C     BACK TO DISK.
C
      IADT2I = IADT2 + IRPDPD(IRPIJ,2)
      IF(IRPIJ.NE.1)THEN
      IJ = (J-1)*POP(IRPI,2) + I
      ELSE
      IJ = INDEX(J-1) + I
      ENDIF
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPIJ,2)+IJ,
     1            1,1,IRPIJ,LIST2I2+LISTOFF)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2),
     1            IRPDPD(IRPIJ,2), 1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPIJ,2)+IJ,
     1            1,1,IRPIJ,LIST2I2+LISTOFF)
C
C
C      AB  AB                      ACD
C     D   T    =    SUM  <JB//CD> T
C      IK  IK       C<D            IJK
C                    J
C      BA  BA             BA  BA   BBA
C      BA  BA                      BBA
C
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  130 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPJ)
      LENV(IRPB) = VRT(IRPB,1)*IRPDPD(IRPCD,13)
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
      IF(VRT(IRPB,1).EQ.0) GOTO 140
      IRPBJ = DIRPRD(IRPB,IRPJ)
      IRPCD = IRPBJ
      CALL GETLST(ICORE(IADV(IRPB)),
     1            IOFFVO(IRPJ,IRPBJ,4)+(J-1)*VRT(IRPB,1)+1,
     1            VRT(IRPB,1),2,IRPBJ,LWIC28)
C    1            VRT(IRPB,1),2,IRPBJ,30)
  140 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      DO  150 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPIK)
      LENT2E(IRPB) = VRT(IRPA,2) * VRT(IRPB,1)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  150 CONTINUE
C
C     FORM T2 MATRICES. (A,B) (BETA,ALPHA)
C
      DO  160 IRPCD=1,NIRREP
      IRPA = DIRPRD(IRPCD,IRPIJK)
      IRPB = DIRPRD(IRPCD,IRPJ)
      IF(VRT(IRPA,2).EQ.0.OR.VRT(IRPB,1).EQ.0) GOTO 160
      CALL XGEMM('T','N',VRT(IRPA,2),VRT(IRPB,1),IRPDPD(IRPCD,13),
     1          1.0D+00,
     1          W(IADW(IRPA)),IRPDPD(IRPCD,13),
     1          ICORE(IADV(IRPB)),IRPDPD(IRPCD,13),
     1          0.0D+00,ICORE(IADT2E(IRPB)),VRT(IRPA,2))
  160 CONTINUE
C
C     READ CURRENT T2 INCREMENT, ADD THIS PIECE TO IT, AND WRITE SUM
C     BACK TO DISK.
C
      IADT2I = IADT2E(NIRREP) + LENT2E(NIRREP)
      IK = (I-1)*POP(IRPK,1) + K
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I3+LISTOFF)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2E(1)),
     1            IRPDPD(IRPIK,13), 1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I3+LISTOFF)
C
C
C
C      AB  AB                      ACD
C     D   T    =  - SUM  <IB//CD> T
C      JK  JK       C<D            IJK
C                    I
C      BA  BA             BA  BA   BBA
C      BA  BA                      BBA
C
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  230 IRPB=1,NIRREP
      IRPCD = DIRPRD(IRPB,IRPI)
      LENV(IRPB) = VRT(IRPB,1)*IRPDPD(IRPCD,13)
      IF(IRPB.EQ.1)THEN
      IADV(IRPB) = 1
      ELSE
      IADV(IRPB) = IADV(IRPB-1) + LENV(IRPB-1)
      ENDIF
  230 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES.
C
      DO  240 IRPB=1,NIRREP
      IF(VRT(IRPB,1).EQ.0) GOTO 240
      IRPBI = DIRPRD(IRPB,IRPI)
      IRPCD = IRPBI
      CALL GETLST(ICORE(IADV(IRPB)),
     1     IOFFVO(IRPI,IRPBI,4)+(I-1)*VRT(IRPB,1)+1,
     1     VRT(IRPB,1),2,IRPBI,LWIC28)
C    1     VRT(IRPB,1),2,IRPBI,30)
  240 CONTINUE
C
C     SET ADDRESSES FOR T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      DO  250 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPJK)
      LENT2E(IRPB) = VRT(IRPA,2) * VRT(IRPB,1)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  250 CONTINUE
C
C     FORM T2 MATRICES. (A,B)
C
      DO  260 IRPCD=1,NIRREP
      IRPA = DIRPRD(IRPCD,IRPIJK)
      IRPB = DIRPRD(IRPCD,IRPI)
      IF(VRT(IRPA,2).EQ.0.OR.VRT(IRPB,1).EQ.0) GOTO 260
      CALL XGEMM('T','N',VRT(IRPA,2),VRT(IRPB,1),IRPDPD(IRPCD,13),
     1           1.0D+00,
     1           W(IADW(IRPA)),IRPDPD(IRPCD,13),
     1           ICORE(IADV(IRPB)),IRPDPD(IRPCD,13),
     1           0.0D+00,
     1           ICORE(IADT2E(IRPB)),VRT(IRPA,2))
  260 CONTINUE
C
C     READ CURRENT T2 INCREMENT, ADD THIS PIECE TO IT, AND WRITE SUM
C     BACK TO DISK.
C
      IADT2I = IADT2E(NIRREP) + LENT2E(NIRREP)
      JK = (J-1)*POP(IRPK,1) + K
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I3+LISTOFF)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2E(1)),
     1            IRPDPD(IRPJK,13),-1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I3+LISTOFF)
C
C
C      AB  AB                      CDB
C     D   T    =    SUM  <AJ//CD> T
C      IK  IK       C<D            IJK
C                    J
C      BA  BA             BB  BB   BBA
C      BA  BA                      BBA
C
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  330 IRPA=1,NIRREP
      IRPCD = DIRPRD(IRPA,IRPJ)
      LENV(IRPA) = VRT(IRPA,2)*IRPDPD(IRPCD,2)
      IF(IRPA.EQ.1)THEN
      IADV(IRPA) = 1
      ELSE
      IADV(IRPA) = IADV(IRPA-1) + LENV(IRPA-1)
      ENDIF
  330 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES
C
      IADVT = IADV(NIRREP) + LENV(NIRREP)
C
      DO  340 IRPA=1,NIRREP
      IF(VRT(IRPA,2).EQ.0) GOTO 340
      IRPAJ = DIRPRD(IRPA,IRPJ)
      IRPCD = IRPAJ
      CALL GETLST(ICORE(IADV(IRPA)),
     1     IOFFVO(IRPJ,IRPAJ,2) + (J-1)*VRT(IRPA,2) + 1,
     1     VRT(IRPA,2),2,IRPAJ,LWIC26)
C    1     VRT(IRPA,2),2,IRPAJ,28)
C
C     TRANSPOSE BLOCKS OF INTEGRALS OUT-OF-PLACE AND COPY BACK
C     TO ORIGINAL POSITION.
C
CJ      CALL TRANSP(ICORE(IADV(IRPA)),ICORE(IADVT),
CJ     1            VRT(IRPA,2),IRPDPD(IRPCD,2))
CJ      CALL  ICOPY(IINTFP*VRT(IRPA,2)*IRPDPD(IRPCD,2),
CJ     1            ICORE(IADVT),1,ICORE(IADV(IRPA)),1)
  340 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      DO  350 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPIK)
      LENT2E(IRPB) = VRT(IRPA,2) * VRT(IRPB,1)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  350 CONTINUE
C
C     FORM T2 MATRICES. (A,B)
C
      DO  360 IRPCD=1,NIRREP
      IRPB = DIRPRD(IRPCD,IRPIJK)
      IRPA = DIRPRD(IRPCD,IRPJ)
      IF(VRT(IRPA,2).EQ.0.OR.VRT(IRPB,1).EQ.0) GOTO 360
      CALL XGEMM('T','N',VRT(IRPA,2),VRT(IRPB,1),IRPDPD(IRPCD,2),
     1           1.0D+00,
     1           ICORE(IADV(IRPA)),IRPDPD(IRPCD,2),
     1           T3(IADT3(IRPB)),IRPDPD(IRPCD,2),
     1           0.0D+00,
     1           ICORE(IADT2E(IRPB)),VRT(IRPA,2))
  360 CONTINUE
C
C     READ CURRENT T2 INCREMENT, ADD THIS PIECE TO IT, AND WRITE SUM
C     BACK TO DISK.
C
      IADT2I = IADT2E(NIRREP) + LENT2E(NIRREP)
      IK = (I-1)*POP(IRPK,1) + K
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I3+LISTOFF)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2E(1)),
     1            IRPDPD(IRPIK,13), 1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPI,IRPIK,5)+IK,
     1            1,1,IRPIK,LIST2I3+LISTOFF)
C
C      AB  AB                      CDB
C     D   T    =  - SUM  <AI//CD> T
C      JK  JK       C<D            IJK
C                    I
C      BA  BA             BB  BB   BBA
C      BA  BA                      BBA
C
C     SET ADDRESSES FOR INTEGRALS/INTERMEDIATES.
C
      DO  430 IRPA=1,NIRREP
      IRPCD = DIRPRD(IRPA,IRPI)
      LENV(IRPA) = VRT(IRPA,2)*IRPDPD(IRPCD,2)
      IF(IRPA.EQ.1)THEN
      IADV(IRPA) = 1
      ELSE
      IADV(IRPA) = IADV(IRPA-1) + LENV(IRPA-1)
      ENDIF
  430 CONTINUE
C
C     READ INTEGRALS/INTERMEDIATES
C
      IADVT = IADV(NIRREP) + LENV(NIRREP)
C
      DO  440 IRPA=1,NIRREP
      IF(VRT(IRPA,1).EQ.0) GOTO 440
      IRPAI = DIRPRD(IRPA,IRPI)
      IRPCD = IRPAI
      CALL GETLST(ICORE(IADV(IRPA)),
     1     IOFFVO(IRPI,IRPAI,2) + (I-1)*VRT(IRPA,2) + 1,
     1     VRT(IRPA,2),2,IRPAI,LWIC26)
C    1     VRT(IRPA,2),2,IRPAI,28)
C
C     TRANSPOSE BLOCKS OF INTEGRALS OUT-OF-PLACE AND COPY BACK
C     TO ORIGINAL POSITION.
C
CJ      CALL TRANSP(ICORE(IADV(IRPA)),ICORE(IADVT),
CJ     1            VRT(IRPA,2),IRPDPD(IRPCD,2))
CJ      CALL  ICOPY(IINTFP*VRT(IRPA,2)*IRPDPD(IRPCD,2),
CJ     1            ICORE(IADVT),1,ICORE(IADV(IRPA)),1)
  440 CONTINUE
C
C     SET ADDRESSES FOR EXPANDED T2 MATRICES. (A,B)
C
      IOFFT2E = IADV(NIRREP) + LENV(NIRREP)
      DO  450 IRPB=1,NIRREP
      IRPA = DIRPRD(IRPB,IRPJK)
      LENT2E(IRPB) = VRT(IRPA,2) * VRT(IRPB,1)
      IF(IRPB.EQ.1)THEN
      IADT2E(IRPB) = IOFFT2E
      ELSE
      IADT2E(IRPB) = IADT2E(IRPB-1) + LENT2E(IRPB-1)
      ENDIF
  450 CONTINUE
C
C     FORM T2 MATRICES. (A,B)
C
      DO  460 IRPCD=1,NIRREP
      IRPB = DIRPRD(IRPCD,IRPIJK)
      IRPA = DIRPRD(IRPCD,IRPI)
      IF(VRT(IRPA,2).EQ.0.OR.VRT(IRPB,1).EQ.0) GOTO 460
      CALL XGEMM('T','N',VRT(IRPA,2),VRT(IRPB,1),IRPDPD(IRPCD,2),
     1           1.0D+00,
     1           ICORE(IADV(IRPA)),IRPDPD(IRPCD,2),
     1           T3(IADT3(IRPB)),IRPDPD(IRPCD,2),
     1           0.0D+00,
     1           ICORE(IADT2E(IRPB)),VRT(IRPA,2))
  460 CONTINUE
C
C     READ CURRENT T2 INCREMENT, ADD THIS PIECE TO IT, AND WRITE SUM
C     BACK TO DISK.
C
      IADT2I = IADT2E(NIRREP) + LENT2E(NIRREP)
      JK = (J-1)*POP(IRPK,1) + K
      CALL GETLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I3+LISTOFF)
      CALL   VADD(ICORE(IADT2I),ICORE(IADT2I),ICORE(IADT2E(1)),
     1            IRPDPD(IRPJK,13),-1.0D+00)
      CALL PUTLST(ICORE(IADT2I),IOFFOO(IRPJ,IRPJK,5)+JK,
     1            1,1,IRPJK,LIST2I3+LISTOFF)
C
      RETURN
      END