      SUBROUTINE EL2R1(CORE,MAXCOR,IUHF,IRREPW,IRREPT,
     &                 LWVOFF,LWOOFF,LT1,LT1OFF,LT2OFF,LZOFF,E2W1)
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C     Argument declarations.
C-----------------------------------------------------------------------
      DOUBLE PRECISION CORE(1),E2W1
      INTEGER MAXCOR,IUHF,IRREPW,IRREPT,LWVOFF,LWOOFF,LT1,LT1OFF,
     &        LT2OFF,LZOFF
C
C-----------------------------------------------------------------------
C     Common block declarations.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER NSTART,NIRREP,IRREPA,IRREPB,DIRPRD
      INTEGER IRPDPD,ISYTYP,ID
      INTEGER POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB
      INTEGER ISYMOFF
C
C-----------------------------------------------------------------------
C     Local variable declarations.
C-----------------------------------------------------------------------
      INTEGER I0T1(2),I000,I010,I020,I030,I040,I050
      INTEGER IRREPZ,DSZZ,NDSZ,NDSZEXP,I,J,IJ,IJEXP,F,FI,FJ
      INTEGER ISPIN,IRPIJ,IRPAB,IRPI,IRPJ,IRPF,IRPFI,IRPFJ
C
      INTEGER IOFFZ,IOFFW,IOFFT,IRPA,IRPB,IRPM,IRPMA,IRPMB,DSZZEXP,
     &        DSZW,NDSW
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DOUBLE PRECISION E2W1AA,E2W1BB,E2W1AB,ESPAD
      LOGICAL SKIPVVVO,SKIPOOOV
C
C-----------------------------------------------------------------------
C     In-line function.
C-----------------------------------------------------------------------
      INTEGER INDEX
C-----------------------------------------------------------------------
C     External function.
C-----------------------------------------------------------------------
      DOUBLE PRECISION SDOT
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C-----------------------------------------------------------------------
      DATA ONE,ONEM,ZILCH / 1.0D+00,-1.0D+00,0.0D+00/
C-----------------------------------------------------------------------
      INDEX(I) = I*(I-1)/2
C-----------------------------------------------------------------------
C
C     - \sum_f W(abfi) * t(f,j) + \sum_f W(abfj) * t(f,i)
C
C     IRREPT --- Irrep of T_1/R_1 vector.
C     IRREPW --- Irrep of W(abfi).
C     IRREPZ --- Irrep of W * T. This must also be the L2 symmetry for
C                this routine to given an energy.
C                On the whole we expect IRREPT=IRREPX, IRREPW=1, but we
C                try to code for more general case.
C
      IRREPZ = DIRPRD(IRREPW,IRREPT)
C
      E2W1   = ZILCH
      E2W1AA = ZILCH
      E2W1BB = ZILCH
      E2W1AB = ZILCH
      ESPAD  = ZILCH
C
C-----------------------------------------------------------------------
C     Get T1 vectors.
C-----------------------------------------------------------------------
      I0T1(1) = 1
      I0T1(2) = I0T1(1) + IRPDPD(IRREPT, 9)
      I000    = I0T1(2) + IRPDPD(IRREPT,10)
      CALL GETLST(CORE(I0T1(1)),1,1,1,LT1OFF+1     ,LT1)
      CALL GETLST(CORE(I0T1(2)),1,1,1,LT1OFF+1+IUHF,LT1)
C-----------------------------------------------------------------------
C
      SKIPVVVO = .FALSE.
      SKIPOOOV = .FALSE.
      IF(SKIPVVVO) GOTO 300
C
      IF(IUHF.NE.0)THEN
C
      DO 100 ISPIN=1,2
      DO  90 IRPIJ=1,NIRREP
      IRPAB = DIRPRD(IRREPZ,IRPIJ)
      DSZZ = IRPDPD(IRPAB,ISPIN)
      NDSZ = IRPDPD(IRPIJ,ISPIN + 2)
      IF(DSZZ.EQ.0.OR.NDSZ.EQ.0) GOTO 90
C
      NDSZEXP = IRPDPD(IRPIJ,20+ISPIN)
      IRPFI   = DIRPRD(IRREPW,IRPAB)
C
      I010 = I000 + DSZZ * NDSZEXP
      I020 = I010 + DSZZ
      CALL ZERO(CORE(I000),DSZZ*NDSZEXP)
C
      DO  80 IRPI=1,NIRREP
      IRPF    = DIRPRD(IRPFI,IRPI)
C
      DO  70    I=1,POP(IRPI,ISPIN)
      DO  60    F=1,VRT(IRPF,ISPIN)
C
      FI = ISYMOFF(IRPI,IRPFI,8+ISPIN) + (I-1)*VRT(IRPF,ISPIN) + F - 1
      CALL GETLST(CORE(I010),FI,1,2,IRPFI,LWVOFF+ISPIN)
C
      IRPJ    = DIRPRD(IRPI,IRPIJ)
      DO  50    J=1,POP(IRPJ,ISPIN)
      IJEXP = ISYMOFF(IRPJ,IRPIJ,20+ISPIN) + (J-1)*POP(IRPI,ISPIN) +I-1
C
      FJ  = ISYMOFF(IRPJ,IRREPT,8+ISPIN) + (J-1)*VRT(IRPF,ISPIN) +F-1
C
      CALL SAXPY(DSZZ,-CORE(I0T1(ISPIN) - 1 + FJ),CORE(I010),1,
     &           CORE(I000 + (IJEXP-1)*DSZZ),1)
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C
C-----------------------------------------------------------------------
C     We have Z(A<B;I,J) at I000. Try to generate 
C     Z(A<B;I<J) = Z(A<B;I,J) - Z(A<B;J,I) at I010. We can dump it to
C     disk to be read later or we can immediately contract with L2 to
C     get energy.
C-----------------------------------------------------------------------
C
      I020 = I010 + DSZZ * NDSZ
      CALL ASSYM(IRPIJ,POP(1,ISPIN),DSZZ,DSZZ,CORE(I010),CORE(I000))
C
      CALL GETLST(CORE(I000),1,NDSZ,2,IRPIJ,LT2OFF+ISPIN)
      IF(ISPIN.EQ.1)THEN
        E2W1AA = E2W1AA + SDOT(DSZZ*NDSZ,CORE(I000),1,CORE(I010),1)
      ELSE
        E2W1BB = E2W1BB + SDOT(DSZZ*NDSZ,CORE(I000),1,CORE(I010),1)
      ENDIF
C
   90 CONTINUE
  100 CONTINUE
      write(6,*) ' @EL2R1-I, AA energy after VVVO is ',E2W1AA
      write(6,*) ' @EL2R1-I, BB energy after VVVO is ',E2W1BB
      ENDIF
C
      DO 290 IRPIJ=1,NIRREP
      IRPAB = DIRPRD(IRREPZ,IRPIJ)
      DSZZ = IRPDPD(IRPAB,13)
      NDSZ = IRPDPD(IRPIJ,14)
      IF(DSZZ.EQ.0.OR.NDSZ.EQ.0) GOTO 290
C
C     -W(AbfI) * T(f,j).
C     UHF :   + W(Ab;If) * T(f,j) ; W from LWVOFF+3.
C     RHF :   + W(bA;fI) = W(Ba;Fi); W from LWVOFF+4; T(f,j)=T(F,J).
C             => reorder W. Or reorder Z ?
C
      I010 = I000 + DSZZ * NDSZ
      CALL ZERO(CORE(I000),DSZZ*NDSZ)
C
      I020 = I010 + DSZZ
      IRPFI   = DIRPRD(IRREPW,IRPAB)
C
      DO 180 IRPI=1,NIRREP
      IRPF    = DIRPRD(IRPFI,IRPI)
C
      DO 170    I=1,POP(IRPI,1)
      DO 160    F=1,VRT(IRPF,2)
C
      IF(IUHF.EQ.0)THEN
      FI = ISYMOFF(IRPI,IRPFI,11) + (I-1)*VRT(IRPF,2) + F - 1
      CALL GETLST(CORE(I010),FI,1,2,IRPFI,LWVOFF+4)
      ELSE
      FI = ISYMOFF(IRPF,IRPFI,18) + (F-1)*POP(IRPI,1) + I - 1
      CALL GETLST(CORE(I010),FI,1,2,IRPFI,LWVOFF+3)
      ENDIF
C
      IRPJ    = DIRPRD(IRPI,IRPIJ)
      DO 150    J=1,POP(IRPJ,2)
      IJ = ISYMOFF(IRPJ,IRPIJ,14) + (J-1)*POP(IRPI,1) + I - 1
C
      FJ = ISYMOFF(IRPJ,IRREPT,10) + (J-1)*VRT(IRPF,2) + F - 1
C
      CALL SAXPY(DSZZ,CORE(I0T1(2) - 1 + FJ),CORE(I010),1,
     &           CORE(I000 + (IJ-1)*DSZZ),1)
  150 CONTINUE
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
C
C-----------------------------------------------------------------------
C     RHF : we have Z(bA;I,j) at I000. Try to generate Z(Ab;Ij).
C     UHF : we already have what we want.
C-----------------------------------------------------------------------
C
      IF(IUHF.EQ.0)THEN
        I020 = I010 + MAX(DSZZ,NDSZ)
        I030 = I020 + MAX(DSZZ,NDSZ)
        I040 = I030 + MAX(DSZZ,NDSZ)
        CALL SYMTR3(IRPAB,VRT(1,2),VRT(1,1),DSZZ,NDSZ,CORE(I000),
     &              CORE(I010),CORE(I020),CORE(I030))
      ENDIF
C
C      W(AbFj) * T(F,I).
C
      I020 = I010 + DSZZ
      IRPFJ = DIRPRD(IRREPW,IRPAB)
C
      DO 280 IRPJ=1,NIRREP
      IRPF    = DIRPRD(IRPFJ,IRPJ)
C
      DO 270    J=1,POP(IRPJ,2)
      DO 260    F=1,VRT(IRPF,1)
C
      FJ = ISYMOFF(IRPJ,IRPFJ,11) + (J-1)*VRT(IRPF,1) + F - 1
      CALL GETLST(CORE(I010),FJ,1,2,IRPFJ,LWVOFF+4)
C
      IRPI    = DIRPRD(IRPJ,IRPIJ)
      DO 250    I=1,POP(IRPI,1)
      IJ = ISYMOFF(IRPJ,IRPIJ,14) + (J-1)*POP(IRPI,1) + I - 1
C
      FI = ISYMOFF(IRPI,IRREPT, 9) + (I-1)*VRT(IRPF,1) + F - 1
C
      CALL SAXPY(DSZZ,CORE(I0T1(1) - 1 + FI),CORE(I010),1,
     &           CORE(I000 + (IJ-1)*DSZZ),1)
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(Ab;Ij) is at I000. Dump to disk and/or contract with L(Ab;Ij).
C     Spin-adapt to include alpha-alpha energy in RHF case.
C-----------------------------------------------------------------------
C
      CALL GETLST(CORE(I010),1,NDSZ,2,IRPIJ,LT2OFF+3)
      E2W1AB = E2W1AB + SDOT(DSZZ*NDSZ,CORE(I000),1,CORE(I010),1)
C
      IF(IUHF.EQ.0)THEN
        I020 = I010 + DSZZ * NDSZ
        I030 = I020 + MAX(DSZZ,NDSZ)
        I040 = I030 + MAX(DSZZ,NDSZ)
        CALL SPINAD3(IRPAB,VRT,DSZZ,NDSZ,CORE(I000),
     &               CORE(I020),CORE(I030))
        ESPAD = ESPAD + SDOT(DSZZ*NDSZ,CORE(I000),1,CORE(I010),1)
      ENDIF
C
  290 CONTINUE
      write(6,*) ' @EL2R1-I, AB energy after VVVO is ',E2W1AB
      write(6,*) ' @EL2R1-I, Spin-adapted energy after VVVO is ',ESPAD
C
  300 CONTINUE
      IF(SKIPOOOV) RETURN
C
C-----------------------------------------------------------------------
C     Z(ab,ij) = \sum_m W(maij) T(b,m) - W(mbij) T(a,m)
C-----------------------------------------------------------------------
C
      IF(IUHF.NE.0)THEN
C
      DO 400 ISPIN=1,2
C
      DO 390 IRPIJ=1,NIRREP
C
      IRPAB = DIRPRD(IRPIJ,IRREPZ)
      DSZZ  = IRPDPD(IRPAB,ISPIN)
      NDSZ  = IRPDPD(IRPIJ,ISPIN + 2)
C
      IF(DSZZ.EQ.0.OR.NDSZ.EQ.0) GOTO 390
C
      IRPMA = DIRPRD(IRPIJ,IRREPW)
C
C     Dimensions on disk --- (ij,ma) list. Read and transpose symmetry
C     block of W. End up with W(ma,ij) at I020.
C
      DSZW  = NDSZ
      NDSW  = IRPDPD(IRPMA,15+ISPIN)
      DSZZEXP = IRPDPD(IRPAB,18+ISPIN)
C
      I010 = I000 + DSZZ
      I020 = I010 + DSZZEXP
      I030 = I020 + NDSW * DSZW
      I040 = I030 + DSZW * NDSW
C
      CALL GETLST(CORE(I030),1,NDSW,2,IRPMA,LWOOFF+ISPIN)
      CALL TRANSP(CORE(I030),CORE(I020),NDSW,DSZW)
C
C     I000 --- Zsqueezed(a<b) = Z(a,b) - Z(b,a).
C     I010 --- Z(a,b)
C     I020 --- W(ma,ij) (one symmetry block)
C     I030 --- W(ij,ma) (one symmetry block)
C
      DO 380 IJ=1,NDSZ
C
      CALL ZERO(CORE(I000),DSZZ)
      CALL ZERO(CORE(I010),DSZZEXP)
C
      DO 370 IRPM=1,NIRREP
C
      IF(POP(IRPM,ISPIN).EQ.0) GOTO 370
C
      IRPA = DIRPRD(IRPM,IRPMA)
      IRPB = DIRPRD(IRPM,IRREPT)
C
      IF(VRT(IRPA,ISPIN).EQ.0.OR.VRT(IRPB,ISPIN).EQ.0) GOTO 370
C
      IOFFZ = I010 + ISYMOFF(IRPB,IRPAB,18+ISPIN) - 1
      IOFFW = I020 + (IJ-1)*NDSW + ISYMOFF(IRPA,IRPMA,15+ISPIN) - 1
      IOFFT = I0T1(ISPIN) + ISYMOFF(IRPM,IRREPT,8+ISPIN) - 1
C
      CALL XGEMM('T','T',
     &           VRT(IRPA,ISPIN),VRT(IRPB,ISPIN),POP(IRPM,ISPIN),
     &           ONE,
     &           CORE(IOFFW),POP(IRPM,ISPIN),
     &           CORE(IOFFT),VRT(IRPB,ISPIN),ZILCH,
     &           CORE(IOFFZ),VRT(IRPA,ISPIN))
  370 CONTINUE
C
C-----------------------------------------------------------------------
C     For each IJ we have Z(a,b) at I010. Compute Z(a,b)-Z(b,a) and
C     squeeze it to Z(a<b) at I000. Then dump to disk and/or contract
C      with L2.
C-----------------------------------------------------------------------
C
      CALL ASSYM4(IRPAB,VRT(1,ISPIN),DSZZ,DSZZEXP,1,CORE(I000),
     &            CORE(I010),ISPIN,18+ISPIN)
C
      CALL GETLST(CORE(I010),IJ,1,2,IRPIJ,LT2OFF+ISPIN)
      IF(ISPIN.EQ.1)THEN
        E2W1AA = E2W1AA + SDOT(DSZZ,CORE(I000),1,CORE(I010),1)
      ELSE
        E2W1BB = E2W1BB + SDOT(DSZZ,CORE(I000),1,CORE(I010),1)
      ENDIF
C
  380 CONTINUE
  390 CONTINUE
  400 CONTINUE
      write(6,*) ' @EL2R1-I, AA energy after OOOV is ',E2W1AA
      write(6,*) ' @EL2R1-I, BB energy after OOOV is ',E2W1BB
C
      ENDIF
C
      DO 590 IRPIJ=1,NIRREP
C
      IRPAB = DIRPRD(IRPIJ,IRREPZ)
      DSZZ  = IRPDPD(IRPAB,13)
      NDSZ  = IRPDPD(IRPIJ,14)
C
      IF(DSZZ.EQ.0.OR.NDSZ.EQ.0) GOTO 590
C
      IRPMA = DIRPRD(IRPIJ,IRREPW)
C
C     UHF : W(mAIj) = -W(AmIj). List LWOOFF+3.
C     RHF : W(mAIj) = -W(mAjI) = -W(MaJi). Read Ji, list LWOOFF+4. 
C     Dimensions on disk --- (ij,ma) list. Read and transpose symmetry
C     block of W. End up with W(ma,ij) at I010.
C     UHF order : A,m. RHF order : M,a.
C
      DSZW  = NDSZ
      NDSW  = IRPDPD(IRPMA,11)
C
C     I000 --- Z(a,b)
C     I010 --- W(ma,ij) (one symmetry block); A,m (UHF); m,A (RHF).
C     I020 --- W(ij,ma) (one symmetry block)
C
      I010 = I000 + DSZZ
      I020 = I010 + NDSW * DSZW
      I030 = I020 + DSZW * NDSW
C
      CALL GETLST(CORE(I020),1,NDSW,2,IRPMA,LWOOFF+4-IUHF)
      CALL TRANSP(CORE(I020),CORE(I010),NDSW,DSZW)
C
C     In RHF we flip i and j. In UHF we flip a and m. In this way we
C     always end up with W(mA,Ij) order.
C
      I030 = I020 + MAX(NDSW,DSZW)
      I040 = I030 + MAX(NDSW,DSZW)
      I050 = I040 + MAX(NDSW,DSZW)
      IF(IUHF.EQ.0)THEN
        CALL SYMTR1(IRPIJ,POP(1,2),POP(1,1),NDSW,CORE(I010),
     &              CORE(I020),CORE(I030),CORE(I040))
      ELSE
        CALL SYMTR3(IRPMA,VRT(1,1),POP(1,2),NDSW,DSZW,CORE(I010),
     &              CORE(I020),CORE(I030),CORE(I040))
      ENDIF
C
      DO 480 IJ=1,NDSZ
C
      CALL ZERO(CORE(I000),DSZZ)
C
      DO 470 IRPM=1,NIRREP
C
      IF(POP(IRPM,2).EQ.0) GOTO 470
C
      IRPA = DIRPRD(IRPM,IRPMA)
      IRPB = DIRPRD(IRPM,IRREPT)
C
      IF(VRT(IRPA,1).EQ.0.OR.VRT(IRPB,2).EQ.0) GOTO 470
C
      IOFFZ = I000 + ISYMOFF(IRPB,IRPAB,13) - 1
      IOFFW = I010 + (IJ-1)*NDSW + ISYMOFF(IRPA,IRPMA,25) - 1
      IOFFT = I0T1(2) + ISYMOFF(IRPM,IRREPT,10) - 1
C
      CALL XGEMM('T','T',
     &           VRT(IRPA,1),VRT(IRPB,2),POP(IRPM,2),
     &           ONEM,
     &           CORE(IOFFW),POP(IRPM,2),
     &           CORE(IOFFT),VRT(IRPB,2),ZILCH,
     &           CORE(IOFFZ),VRT(IRPA,1))
  470 CONTINUE
C
C-----------------------------------------------------------------------
C     For each Ij we have Z(A,b) at I000. Dump to disk and/or contract
C     with L2.
C-----------------------------------------------------------------------
C
      CALL GETLST(CORE(I020),IJ,1,2,IRPIJ,LT2OFF+3)
      E2W1AB = E2W1AB + SDOT(DSZZ,CORE(I000),1,CORE(I020),1)
C
  480 CONTINUE
C
C-----------------------------------------------------------------------
C     Z(Ab,Ij) = - W(Mb,Ij) T(A,M)
C-----------------------------------------------------------------------

      IRPMB = DIRPRD(IRPIJ,IRREPW)
C
C     UHF, RHF : W(MbIj). List LWOOFF+4.
C     Dimensions on disk --- (ij,mb) list. Read and transpose symmetry
C     block of W. End up with W(mb,ij) at I010.
C
      DSZW  = NDSZ
      NDSW  = IRPDPD(IRPMB,18)
C
C     I000 --- Z(A,b)
C     I010 --- W(Mb,Ij) (one symmetry block).
C     I020 --- W(Ij,Mb) (one symmetry block)
C
      I010 = I000 + DSZZ
      I020 = I010 + NDSW * DSZW
      I030 = I020 + DSZW * NDSW
C
      CALL GETLST(CORE(I020),1,NDSW,2,IRPMB,LWOOFF+4)
      CALL TRANSP(CORE(I020),CORE(I010),NDSW,DSZW)
C
      DO 580 IJ=1,NDSZ
C
      CALL ZERO(CORE(I000),DSZZ)
C
      DO 570 IRPM=1,NIRREP
C
      IF(POP(IRPM,1).EQ.0) GOTO 570
C
      IRPB = DIRPRD(IRPM,IRPMB)
      IRPA = DIRPRD(IRPM,IRREPT)
C
      IF(VRT(IRPA,1).EQ.0.OR.VRT(IRPB,2).EQ.0) GOTO 570
C
      IOFFZ = I000 + ISYMOFF(IRPB,IRPAB,13) - 1
      IOFFW = I010 + (IJ-1)*NDSW + ISYMOFF(IRPB,IRPMB,18) - 1
      IOFFT = I0T1(1) + ISYMOFF(IRPM,IRREPT, 9) - 1
C
      CALL XGEMM('N','N',
     &           VRT(IRPA,1),VRT(IRPB,2),POP(IRPM,1),
     &           ONEM,
     &           CORE(IOFFT),VRT(IRPA,2),
     &           CORE(IOFFW),POP(IRPM,1),ZILCH,
     &           CORE(IOFFZ),VRT(IRPA,1))
  570 CONTINUE
C
C-----------------------------------------------------------------------
C     For each Ij we have Z(A,b) at I000. Dump to disk and/or contract
C     with L2. Spin-adapt to include alpha-alpha energy in RHF case.
C-----------------------------------------------------------------------
C
      CALL GETLST(CORE(I020),IJ,1,2,IRPIJ,LT2OFF+3)
      E2W1AB = E2W1AB + SDOT(DSZZ,CORE(I000),1,CORE(I020),1)
C
      IF(IUHF.EQ.0)THEN
        I030 = I020 + DSZZ
        I040 = I030 + DSZZ
        I050 = I040 + DSZZ
        CALL SPINAD3(IRPAB,VRT,DSZZ,1,CORE(I000),
     &               CORE(I030),CORE(I040))
        ESPAD = ESPAD + SDOT(DSZZ,CORE(I000),1,CORE(I020),1)
      ENDIF
C
  580 CONTINUE
C
  590 CONTINUE
      write(6,*) ' @EL2R1-I, AB energy after OOOV is ',E2W1AB
      write(6,*) ' @EL2R1-I, Spin-adapted energy after OOOV is ',ESPAD
      IF(IUHF.NE.0)THEN
        E2W1 = E2W1AA + E2W1BB + E2W1AB
      ELSE
        E2W1 =                   ESPAD
      ENDIF
C
      RETURN
      END
