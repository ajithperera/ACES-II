
#define _TJK
#define _TIK
#define _TIJ
#define _HHJK
#define _HHIK
#define _HHIJ

      SUBROUTINE D1T14(W,T2IJ,T2IK,T2JK,
     1                 VIJ,VIK,VJK,
     1                 ICORE,ISPIN,
     1                 IADBLK,LENBLK,IADT2,LENT2,IADV,LENV,
     1                 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1                 DSZEIJ,DSZEIK,DSZEJK,
     1                 LNOVIJ,LNOVIK,LNOVJK,MAXCOR,LAMBDA)
      IMPLICIT INTEGER (A-Z)
      LOGICAL LAMBDA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION W(1),
     1                 T2IJ(DSZEIJ,1),
     1                 T2IK(DSZEIK,1),
     1                 T2JK(DSZEJK,1)
      DOUBLE PRECISION  VIJ(LNOVIJ,1),
     1                  VIK(LNOVIK,1),
     1                  VJK(LNOVJK,1)
      DOUBLE PRECISION ICORE(1)
C
      DIMENSION IADBLK(8),LENBLK(8)
      DIMENSION IADT2(8) ,LENT2(8),
     1          IADV(8)  ,LENV(8)
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
C     List offset
      IF (LAMBDA) THEN
         LISTOFF=100
      ELSE
         LISTOFF=0
      ENDIF
      TOL = 1.0D-12
C
      IF(ISPIN.EQ.1)THEN
      LISTWV = LWIC15
      LISTWO = LWIC11
      ELSE
      LISTWV = LWIC16
      LISTWO = LWIC12
      ENDIF
C
C               EA              EC              EB
C   - <BC//EK> T    - <AB//EK> T    + <AC//EK> T
C               IJ              IJ              IJ
C
      DO   10 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     1                  VRT(DIRPRD(IRREP-1,IRPIJ),ISPIN)
      ENDIF
   10 CONTINUE
C
      IF(IRPI.EQ.IRPJ)THEN
      IJ = INDEX(J-1) + I
      ELSE
      IJ = (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
C     Loop over symmetry blocks of E. For each non-empty block, read
C     integral/intermediate, determine which blocks of A and BC will be
C     computed and compute them.
C
      DO   20 IRPE=1,NIRREP
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 20
      IRPBC = DIRPRD(IRPE,IRPK)
      LENGTHV = IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T14',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEK = DIRPRD(IRPE,IRPK)
      EKEFF = IOFFVO(IRPK,IRPEK,ISPIN)+(K-1)*VRT(IRPE,ISPIN)+1
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPK,IRPEK,ISPIN)+(K-1)*VRT(IRPE,ISPIN)+1,
     1            VRT(IRPE,ISPIN),2,IRPEK,LISTWV)
C    1            VRT(IRPE,ISPIN),2,IRPEK,26+ISPIN)
C
      IRPA   = DIRPRD(IRPE,IRPIJ)
#ifdef _TIJ
      CALL XGEMM('N','N',IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),
     1          VRT(IRPE,ISPIN),-1.0D+00,
     1          ICORE(1),IRPDPD(IRPBC,ISPIN),
     1          T2IJ(IADT2(IRPA),IJ),VRT(IRPE,ISPIN),
     1          1.0D+00,W(IADBLK(IRPA)),IRPDPD(IRPBC,ISPIN))
#endif
   20 CONTINUE
C
C               EA              EC              EB
C   - <BC//EI> T    - <AB//EI> T    + <AC//EI> T
C               JK              JK              JK
C
      DO  110 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     1                  VRT(DIRPRD(IRREP-1,IRPJK),ISPIN)
      ENDIF
  110 CONTINUE
C
      IF(IRPJ.EQ.IRPK)THEN
      JK = INDEX(K-1) + J
      ELSE
      JK = (K-1)*POP(IRPJ,ISPIN) + J
      ENDIF
C
      DO  120 IRPE=1,NIRREP
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 120
      IRPBC = DIRPRD(IRPE,IRPI)
      LENGTHV = IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T14',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEI = DIRPRD(IRPE,IRPI)
      EIEFF = IOFFVO(IRPI,IRPEI,ISPIN)+(I-1)*VRT(IRPE,ISPIN)+1
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPI,IRPEI,ISPIN)+(I-1)*VRT(IRPE,ISPIN)+1,
     1            VRT(IRPE,ISPIN),2,IRPEI,LISTWV)
C    1            VRT(IRPE,ISPIN),2,IRPEI,26+ISPIN)
C
      IRPA   = DIRPRD(IRPE,IRPJK)
#ifdef _TJK
      CALL XGEMM('N','N',IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),
     1          VRT(IRPE,ISPIN),-1.0D+00,
     1          ICORE(1),IRPDPD(IRPBC,ISPIN),
     1          T2JK(IADT2(IRPA),JK),VRT(IRPE,ISPIN),
     1          1.0D+00,W(IADBLK(IRPA)),IRPDPD(IRPBC,ISPIN))
#endif
  120 CONTINUE
C
C                 EA              EC              EB
C     - <BC//EJ> T    - <AB//EJ> T    + <AC//EJ> T
C                 IK              IK              IK
C
      DO  210 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADT2(IRREP) = 1
      ELSE
      IADT2(IRREP) = IADT2(IRREP-1) + VRT(IRREP-1,ISPIN) * 
     1                  VRT(DIRPRD(IRREP-1,IRPIK),ISPIN)
      ENDIF
  210 CONTINUE
C
      IF(IRPI.EQ.IRPK)THEN
      IK = INDEX(K-1) + I
      ELSE
      IK = (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      DO  220 IRPE=1,NIRREP
      IF(VRT(IRPE,ISPIN).EQ.0) GOTO 220
      IRPBC = DIRPRD(IRPE,IRPJ)
      LENGTHV = IRPDPD(IRPBC,ISPIN) * VRT(IRPE,ISPIN)
C
      IF(LENGTHV.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T14',LENGTHV,MAXCOR)
      ENDIF
C
      IRPEJ = DIRPRD(IRPE,IRPJ)
      EJEFF = IOFFVO(IRPJ,IRPEJ,ISPIN)+(J-1)*VRT(IRPE,ISPIN)+1
      CALL GETLST(ICORE(1),
     1            IOFFVO(IRPJ,IRPEJ,ISPIN)+(J-1)*VRT(IRPE,ISPIN)+1,
     1            VRT(IRPE,ISPIN),2,IRPEJ,LISTWV)
C    1            VRT(IRPE,ISPIN),2,IRPEJ,26+ISPIN)
C
      IRPA   = DIRPRD(IRPE,IRPIK)
#ifdef _TIK
      CALL XGEMM('N','N',IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),
     1          VRT(IRPE,ISPIN), 1.0D+00,
     1          ICORE(1),IRPDPD(IRPBC,ISPIN),
     1          T2IK(IADT2(IRPA),IK),VRT(IRPE,ISPIN),
     1          1.0D+00,W(IADBLK(IRPA)),IRPDPD(IRPBC,ISPIN))
#endif
  220 CONTINUE
C
C
C     NOW FOR THE M SUM.
C
C      BC                  AB                   AC
C     T    <IJ//MA>   +   T    <IJ//MC>    -   T    <IJ//MB>
C      MK                  MK                   MK
C
C     COMPUTE LENGTHS AND HENCE ADDRESSES OF V SUBBLOCKS. THEY
C     WILL BE (M,A). LABEL BY A.
C
      DO  310 IRPA=1,NIRREP
      IRPM = DIRPRD(IRPA,IRPIJ)
      LENV(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
      IF(IRPA.EQ.1)THEN
      IADV(IRPA) = 1
      ELSE
      IADV(IRPA) = IADV(IRPA-1) + LENV(IRPA-1)
      ENDIF
  310 CONTINUE
C
      IF(IRPIJ.EQ.1)THEN
C     IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + INDEX(J-1) + I
      IJ =                            INDEX(J-1) + I
      ELSE
C     IJ = IOFFOO(IRPJ,IRPIJ,ISPIN) + (J-1)*POP(IRPI,ISPIN) + I
      IJ =                            (J-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      DO  350 IRPM=1,NIRREP
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 350
      IRPBC = DIRPRD(IRPM,IRPK)
      LENGTHT = IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
C
      IF(LENGTHT.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T14',LENGTHT,MAXCOR)
      ENDIF
C
      IRPMK = DIRPRD(IRPM,IRPK)
C
      IF(IRPM.LT.IRPK)THEN
      DO  320 M=1,POP(IRPM,ISPIN)
      MK = IOFFOO(IRPK,IRPMK,ISPIN) + (K-1)*POP(IRPM,ISPIN) + M
      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1            MK,1,1,IRPMK,43+ISPIN+LISTOFF)
  320 CONTINUE
      ENDIF
C
      IF(IRPM.GT.IRPK)THEN
      DO  330 M=1,POP(IRPM,ISPIN)
      MK = IOFFOO(IRPM,IRPMK,ISPIN) + (M-1)*POP(IRPK,ISPIN) + K
      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1            MK,1,1,IRPMK,43+ISPIN+LISTOFF)
  330 CONTINUE
      CALL VMINUS(ICORE(1),POP(IRPM,ISPIN)*IRPDPD(IRPMK,ISPIN))
      ENDIF
C
      IF(IRPM.EQ.IRPK)THEN
      DO  340 M=1,POP(IRPM,ISPIN)
         IF(M.LT.K)THEN
         MK = IOFFOO(IRPK,IRPMK,ISPIN) + INDEX(K-1) + M
         CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1               MK,1,1,IRPMK,43+ISPIN+LISTOFF)
         ENDIF
C
         IF(M.GT.K)THEN
         MK = IOFFOO(IRPM,IRPMK,ISPIN) + INDEX(M-1) + K
         CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1               MK,1,1,IRPMK,43+ISPIN+LISTOFF)
         CALL VMINUS(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1               IRPDPD(IRPMK,ISPIN))
         ENDIF
C
         IF(M.EQ.K)THEN
         CALL    ZERO(ICORE(1+(M-1)*IRPDPD(IRPMK,ISPIN)),
     1                IRPDPD(IRPMK,ISPIN))
         ENDIF
  340 CONTINUE
C
      ENDIF
C
      IRPA   = DIRPRD(IRPM,IRPIJ)
#ifdef _HHIJ
      CALL XGEMM('N','N',IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),
     1           POP(IRPM,ISPIN), 1.0D+00,
     1           ICORE(1),IRPDPD(IRPBC,ISPIN),
     1           VIJ(IADV(IRPA),IJ),POP(IRPM,ISPIN),
     1            1.0D+00,W(IADBLK(IRPA)),IRPDPD(IRPBC,ISPIN))
#endif
  350 CONTINUE
C
C      BC                  AB                   AC
C     T    <JK//MA>   +   T    <JK//MC>    -   T    <JK//MB>
C      MI                  MI                   MI
C
C     COMPUTE LENGTHS AND HENCE ADDRESSES OF V SUBBLOCKS. THEY
C     WILL BE (M,A). LABEL BY A.
C
      DO IRPA = 1, NIRREP
         IRPM = DIRPRD(IRPA,IRPJK)
         LENV(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
         IF (IRPA.EQ.1) THEN
            IADV(IRPA) = 1
         ELSE
            IADV(IRPA) = IADV(IRPA-1) + LENV(IRPA-1)
         END IF
      END DO
      IF (IRPJK.EQ.1) THEN
         JK = INDEX(K-1) + J
      ELSE
         JK = (K-1)*POP(IRPJ,ISPIN) + J
      END IF
      DO IRPM = 1, NIRREP
      IF (POP(IRPM,ISPIN).NE.0) THEN
      IRPBC = DIRPRD(IRPM,IRPI)
      LENGTHT = IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
      IF (LENGTHT.GT.MAXCOR) THEN
         WRITE(6,1010)
         CALL INSMEM('D1T14',LENGTHT,MAXCOR)
      END IF
      IRPMI = DIRPRD(IRPM,IRPI)
      IF (IRPM.LT.IRPI) THEN
         DO M = 1, POP(IRPM,ISPIN)
            MI = IOFFOO(IRPI,IRPMI,ISPIN) + (I-1)*POP(IRPM,ISPIN) + M
            CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                  MI,1,1,IRPMI,43+ISPIN+LISTOFF)
         END DO
      END IF
      IF (IRPM.GT.IRPI) THEN
         DO M = 1, POP(IRPM,ISPIN)
            MI = IOFFOO(IRPM,IRPMI,ISPIN) + (M-1)*POP(IRPI,ISPIN) + I
            CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                  MI,1,1,IRPMI,43+ISPIN+LISTOFF)
         END DO
         CALL VMINUS(ICORE(1),POP(IRPM,ISPIN)*IRPDPD(IRPMI,ISPIN))
      END IF
      IF (IRPM.EQ.IRPI) THEN
         DO M = 1, POP(IRPM,ISPIN)
            IF (M.LT.I) THEN
               MI = IOFFOO(IRPI,IRPMI,ISPIN) + INDEX(I-1) + M
               CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                     MI,1,1,IRPMI,43+ISPIN+LISTOFF)
            END IF
            IF (M.GT.I) THEN
               MI = IOFFOO(IRPM,IRPMI,ISPIN) + INDEX(M-1) + I
               CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                     MI,1,1,IRPMI,43+ISPIN+LISTOFF)
               CALL VMINUS(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                     IRPDPD(IRPMI,ISPIN))
            END IF
            IF (M.EQ.I) THEN
               CALL ZERO(ICORE(1+(M-1)*IRPDPD(IRPMI,ISPIN)),
     1                   IRPDPD(IRPMI,ISPIN))
            END IF
         END DO
      END IF
      IRPA = DIRPRD(IRPM,IRPJK)
#ifdef _HHJK
      CALL XGEMM('N','N',
     1           IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),POP(IRPM,ISPIN),
     1           1.0D0,ICORE(1),          IRPDPD(IRPBC,ISPIN),
     1                 VJK(IADV(IRPA),JK),POP(IRPM,ISPIN),
     1           1.0D0,W(IADBLK(IRPA)),   IRPDPD(IRPBC,ISPIN))
#endif
      END IF
      END DO
C
C
C      BC                  AB                   AC
C   - T    <IK//MA>   -   T    <IK//MC>    +   T    <IK//MB>
C      MJ                  MJ                   MJ
C
C     COMPUTE LENGTHS AND HENCE ADDRESSES OF V SUBBLOCKS. THEY
C     WILL BE (M,A). LABEL BY A.
C
      DO  510 IRPA=1,NIRREP
      IRPM = DIRPRD(IRPA,IRPIK)
      LENV(IRPA) = POP(IRPM,ISPIN) * VRT(IRPA,ISPIN)
      IF(IRPA.EQ.1)THEN
      IADV(IRPA) = 1
      ELSE
      IADV(IRPA) = IADV(IRPA-1) + LENV(IRPA-1)
      ENDIF
  510 CONTINUE
C
      IF(IRPIK.EQ.1)THEN
C     IK = IOFFOO(IRPK,IRPIK,ISPIN) + INDEX(K-1) + I
      IK =                            INDEX(K-1) + I
      ELSE
C     IK = IOFFOO(IRPK,IRPIK,ISPIN) + (K-1)*POP(IRPI,ISPIN) + I
      IK =                            (K-1)*POP(IRPI,ISPIN) + I
      ENDIF
C
      DO  550 IRPM=1,NIRREP
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 550
      IRPBC = DIRPRD(IRPM,IRPJ)
      LENGTHT = IRPDPD(IRPBC,ISPIN) * POP(IRPM,ISPIN)
C
      IF(LENGTHT.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('D1T14',LENGTHT,MAXCOR)
      ENDIF
C
      IRPMJ = DIRPRD(IRPM,IRPJ)
C
      IF(IRPM.LT.IRPJ)THEN
      DO  520    M=1,POP(IRPM,ISPIN)
      MJ = IOFFOO(IRPJ,IRPMJ,ISPIN) + (J-1)*POP(IRPM,ISPIN) + M
      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1            MJ,1,1,IRPMJ,43+ISPIN+LISTOFF)
  520 CONTINUE
      ENDIF
C
      IF(IRPM.GT.IRPJ)THEN
      DO  530    M=1,POP(IRPM,ISPIN)
      MJ = IOFFOO(IRPM,IRPMJ,ISPIN) + (M-1)*POP(IRPJ,ISPIN) + J
      CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1            MJ,1,1,IRPMJ,43+ISPIN+LISTOFF)
  530 CONTINUE
      CALL VMINUS(ICORE(1),POP(IRPM,ISPIN)*IRPDPD(IRPMJ,ISPIN))
      ENDIF
C
      IF(IRPM.EQ.IRPJ)THEN
      DO  540    M=1,POP(IRPM,ISPIN)
         IF(M.LT.J)THEN
         MJ = IOFFOO(IRPJ,IRPMJ,ISPIN) + INDEX(J-1) + M
         CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1               MJ,1,1,IRPMJ,43+ISPIN+LISTOFF)
         ENDIF
C
         IF(M.GT.J)THEN
         MJ = IOFFOO(IRPM,IRPMJ,ISPIN) + INDEX(M-1) + J
         CALL GETLST(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1               MJ,1,1,IRPMJ,43+ISPIN+LISTOFF)
         CALL VMINUS(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1               IRPDPD(IRPMJ,ISPIN))
         ENDIF
C
         IF(M.EQ.J)THEN
         CALL    ZERO(ICORE(1+(M-1)*IRPDPD(IRPMJ,ISPIN)),
     1                IRPDPD(IRPMJ,ISPIN))
         ENDIF
  540 CONTINUE
C
      ENDIF
C
      IRPA   = DIRPRD(IRPM,IRPIK)
#ifdef _HHIK
      CALL XGEMM('N','N',IRPDPD(IRPBC,ISPIN),VRT(IRPA,ISPIN),
     1           POP(IRPM,ISPIN),-1.0D+00,
     1           ICORE(1),IRPDPD(IRPBC,ISPIN),
     1           VIK(IADV(IRPA),IK),POP(IRPM,ISPIN),
     1            1.0D+00,W(IADBLK(IRPA)),IRPDPD(IRPBC,ISPIN))
#endif
  550 CONTINUE
      RETURN
 1010 FORMAT(' @D1T14-I, Insufficient memory to continue. ')
      END
