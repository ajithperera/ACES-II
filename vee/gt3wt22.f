      SUBROUTINE GT3WT22(T3,W,CORE,
     &                 IADT3,IADW,IADT2,LENT2,
     &                 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,
     &                 IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,IUHF,MAXCOR,
     &                 total,irrepx)
C-----------------------------------------------------------------------
C
C     This is moving toward being a generic subroutine for computing
C     contractions of the form Z3 = W*T2 for the AAB spin case of Z3.
C
C     T3    Z3(A<B,c). Labelled by IRPC. Must be initialized by caller.
C     W     Z3(A,Bc).  Labelled by IRPBC. Must be initialized by caller.
C           [N.B. W is not the "W" which is contracted with T2]
C     CORE  Scratch space.
C     IADT3 Addresses for symmetry blocks of T3. Must be supplied by
C           caller.
C     IADW  Addresses for symmetry blocks of W. Must be supplied by
C           caller.
C     IRREPT Overall symmetry of T2.
C     IRREPW Overall symmetry of W in W*T2.
C     
C     Notes : list expectations (Watch for diversions from standards)
C
C     1-NIRREP  LIST2OFF + 1     T2(AB,IJ)      A<B ; I<J     ***
C     1-NIRREP  LIST2OFF + 2     T2(ab,ij)      a<b ; i<j     ***
C     1-NIRREP  LIST2OFF + 3     T2(Ab,Ij)      A,b ; I;j
C
C     1-NIRREP  LWOOFF   + 1      W(KA,IJ)      K,A ; I<J     ***
C     1-NIRREP  LWOOFF   + 2      W(ka,ij)      k,a ; i<j     ***
C     1-NIRREP  LWOOFF   + 3      W(Ak,Ij)      A,k ; I,j     ***
C     1-NIRREP  LWOOFF   + 4      W(Ka,Ij)      K,a ; I,j
C
C     1-NIRREP  LWVOFF   + 1      W(AB,CI)      A<B ; C,I     ***
C     1-NIRREP  LWVOFF   + 2      W(ab,ci)      a<b ; c,i     ***
C     1-NIRREP  LWVOFF   + 3      W(Ab,Ic)      A,b ; c,I     ***
C     1-NIRREP  LWVOFF   + 4      W(Ab,Ci)      A,b ; C,i
C
C    *** means that the list is not used by this routine in RHF calcul-
C        ations, regardless of whether it is present or not. This helps
C        generality --- e.g. all alpha T2 vectors are present on list
C        44, but all alpha R2 vectors are not present on list 444 in
C        RHF EOM-CC calculations. Similarly, all alpha Hbar elements
C        are not usually available in RHF calculations.
C
C-----------------------------------------------------------------------
C
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION T3,W,CORE, SUM, total
      INTEGER IADT3,IADW,IADT2,LENT2,I,J,K,IRPI,IRPJ,IRPK,
     &        IRPIJ,IRPIK,IRPJK,IRREPT,LIST2OFF,IRREPW,LWOOFF,LWVOFF,
     &        IUHF,MAXCOR,LENGTH,kk
C-----------------------------------------------------------------------
C     Common blocks.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        IOFFVV,IOFFOO,IOFFVO
C-----------------------------------------------------------------------
C     Functions.
C-----------------------------------------------------------------------
      INTEGER INDEX, IDSYMSZ
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION ONE,ONEM
      INTEGER LISWOA,LISWOB,LISWOC,LISWOD,LISWVA,LISWVB,LISWVC,LISWVD,
     &        IRPA,IRPAB,IRPAE,IRPBC,IRPC,
     &        IRPE,IRPEI,IRPEJ,IRPEK,IRPEC,
     &        IRPM,IRPIM,IRPMI,IRPJM,IRPMJ,IRPMK,IRPMA,IRPMC,
     &        IJ,JI,IJAB,JIAB,IK,JK,M,IM,IMAB,MIAB,JM,JMAB,KM,MJAB,
     &        IMOFF,JMOFF,
     &        ISPIN1,ISPIN2,
     &        I000,I010,I020,I030,I040,I050,NEED,
     &        DSZ,DSZEXP,LENV,IOFFV,irrepx
C
C-----------------------------------------------------------------------
      DIMENSION T3(1),W(1),CORE(1),IADT3(8),IADW(8),IADT2(8),LENT2(8)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &                NF1BB,NF2BB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C-----------------------------------------------------------------------
      DATA ONE /1.0D+00/, ONEM /-1.0D+00/
C-----------------------------------------------------------------------
      INDEX(I) = I*(I-1)/2
C-----------------------------------------------------------------------
C
      LISWVA = LWVOFF + 1
      LISWVB = LWVOFF + 2
      LISWVC = LWVOFF + 3
      LISWVD = LWVOFF + 4
C
      LISWOA = LWOOFF + 1
      LISWOB = LWOOFF + 2
      LISWOC = LWOOFF + 3
      LISWOD = LWOOFF + 4
C
      ISPIN1 = 1
      ISPIN2 = 2
C
C-----------------------------------------------------------------------
C      AE                BE
C     T   <BC//EK>  -   T   <AC//EK>
C      IJ                IJ
C
C      AA  AB  AB        AA  AB  AB
C      AA                AA
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPIJ)
C     
      DSZ    = IRPDPD(IRPAE,     ISPIN1)
      DSZEXP = IRPDPD(IRPAE,18 + ISPIN1)
C   
CSSS      call zero(CORE(I000), IDSYMSZ(irrepx,ISYTYP(1,46),ISYTYP(2,46)))
CSSS      call getall(CORE(I000),1,irrepx, 403)
CSSS      call checksum("403 in GT3w", CORE(I000),
CSSS     &               IDSYMSZ(irrepx,ISYTYP(1,46),ISYTYP(2,46)),sum)

      IF(IUHF.EQ.0)THEN
        I000 = 1
        I010 = I000 + DSZEXP
        I020 = I010 + DSZEXP
        I030 = I020 + DSZEXP
C
        IJAB = IOFFOO(IRPJ,IRPIJ,5) + (J-1)*POP(IRPI,ISPIN1) + I
        JIAB = IOFFOO(IRPI,IRPIJ,5) + (I-1)*POP(IRPJ,ISPIN1) + J
        CALL GETLST(CORE(I010),IJAB,1,1,IRPIJ,LIST2OFF + 3)

CSSS        Write(6,*) IJAB, IRPIJ
CSSS        Write(6, "(5(F10.6))") (CORE(I010+kk), kk=0, DSZEXP-1)

        CALL GETLST(CORE(I020),JIAB,1,1,IRPIJ,LIST2OFF + 3)

CSSS        Write(6, "(5(F10.6))") (CORE(I020+kk), kk=0, DSZEXP-1)
CSSS        call checksum("DR3=WT2:403", CORE(I010),DSZEXP,sum)

        CALL   VADD(CORE(I010),CORE(I010),CORE(I020),DSZEXP,ONEM)
        CALL  DCOPY(DSZEXP, CORE(I010), 1, CORE(I000), 1)
     
CSSS        Write(6,*) 
CSSS        Write(6, "(5(F10.6))") (CORE(I000+kk), kk=0, DSZEXP-1)
CSSS        call checksum("DR3=WT2:403", CORE(I000),DSZEXP,sum)

      ELSE
        I000 = 1
        I010 = I000 + DSZEXP
        I020 = I010 + DSZ

        IF(IRPI.EQ.IRPJ)THEN
        IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
        ELSE
        IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1)*POP(IRPI,ISPIN1) + I
        ENDIF
        CALL GETLST(CORE(I010),IJ,1,1,IRPIJ,LIST2OFF + ISPIN1)
        CALL SYMEXP2(IRPAE,VRT(1,ISPIN1),DSZEXP,DSZ,1,CORE(I000),
     &                                                CORE(I010))

CSSS        call checksum("DR3=WT2:403", CORE(I000),DSZEXP,sum)
CSSS        total = total + sum
CSSS        Write(6, "(F15.10)") total

      ENDIF
C
      DO   10 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN1) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
   10 CONTINUE
C
      DO   20 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 20

      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEK = DIRPRD(IRPE,IRPK)
      IRPBC = DIRPRD(IRREPW,IRPEK)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN1)
      I020  = I010 + LENV
      NEED  = IINTFP * I020
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     1            IOFFVO(IRPK,IRPEK,4)+(K-1)*VRT(IRPE,ISPIN1)+1,
     1            VRT(IRPE,ISPIN1),2,IRPEK,LISWVD)
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN1),
     &           ONE,CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &               CORE(I010),IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
   20 CONTINUE

CSSS       length = VRT(IRPA,ISPIN1)*IRPDPD(IRPBC,13)
CSSS       call checksum("DR3=WT2", W(IADW(IRPBC)), 
CSSS     &  LENGTH,sum)

CSSS       total = total + sum
CSSS       Write(6,*) total
C
C-----------------------------------------------------------------------
C        AE                BE
C     - T   <BC//IE>   +  T   <AC//IE>
C        JK                JK
C
C        AB  AB  AB        AB  AB  AB
C        AB                AB
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPJK)
C
      DO  110 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN2) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
  110 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPAE,13)
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(I000),JK,1,1,IRPJK,LIST2OFF + 3)
C
      DO  120 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 120
C
      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEI = DIRPRD(IRPE,IRPI)
      IRPBC = DIRPRD(IRREPW,IRPEI)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
      I040 = I030 + IRPDPD(IRPBC,13)
      I050 = I040 + LENV
      NEED = IINTFP * I050
C
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1010)
      CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
        CALL GETLST(CORE(I040),
     &              IOFFVO(IRPI,IRPEI,4) + (I-1)*VRT(IRPE,ISPIN2) + 1,
     &              VRT(IRPE,ISPIN2),2,IRPEI,LISWVD)
        CALL SYMTR3(IRPBC,VRT(1,ISPIN1),VRT(1,ISPIN2),IRPDPD(IRPBC,13),
     &              VRT(IRPE,ISPIN2),
     &              CORE(I040),CORE(I010),CORE(I020),CORE(I030))
      ELSE
        CALL GETLST(CORE(I040),
     &              IOFFVO(IRPI,IRPEI,3) + (I-1)*VRT(IRPE,ISPIN2) + 1,
     &              VRT(IRPE,ISPIN2),2,IRPEI,LISWVC)
      ENDIF
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN2),
     &           ONEM,
     &           CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &           CORE(I040),IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
  120 CONTINUE
C
C-----------------------------------------------------------------------
C      AE                 BE
C     T   <BC//JE>   -   T   <AC//JE>
C      IK                 IK
C
C      AB  AB  AB         AB  AB  AB
C      AB                 AB
C-----------------------------------------------------------------------
C
      IRPAE = DIRPRD(IRREPT,IRPIK)
C
      DO  210 IRPE=1,NIRREP
        IRPA = DIRPRD(IRPE,IRPAE)
        LENT2(IRPE) = VRT(IRPE,ISPIN2) * VRT(IRPA,ISPIN1)
          IF(IRPE.EQ.1)THEN
            IADT2(IRPE) = 1
          ELSE
            IADT2(IRPE) = IADT2(IRPE-1) + LENT2(IRPE-1)
          ENDIF
  210 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPAE,13)
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(I000),IK,1,1,IRPIK,LIST2OFF + 3)
C
      DO  220 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN2).EQ.0) GOTO 220
C
      IRPA  = DIRPRD(IRPE,IRPAE)
      IRPEJ = DIRPRD(IRPE,IRPJ)
      IRPBC = DIRPRD(IRREPW,IRPEJ)
      LENV  = IRPDPD(IRPBC,13) * VRT(IRPE,ISPIN2)
C
      I020 = I010 + IRPDPD(IRPBC,13)
      I030 = I020 + IRPDPD(IRPBC,13)
      I040 = I030 + IRPDPD(IRPBC,13)
      I050 = I040 + LENV
      NEED = IINTFP * I050
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
        CALL GETLST(CORE(I040),
     &              IOFFVO(IRPJ,IRPEJ,4) + (J-1)*VRT(IRPE,ISPIN2) + 1,
     &              VRT(IRPE,ISPIN2),2,IRPEJ,LISWVD)
        CALL SYMTR3(IRPBC,VRT(1,ISPIN1),VRT(1,ISPIN2),IRPDPD(IRPBC,13),
     &              VRT(IRPE,ISPIN2),
     &              CORE(I040),CORE(I010),CORE(I020),CORE(I030))
      ELSE
        CALL GETLST(CORE(I040),
     &              IOFFVO(IRPJ,IRPEJ,3) + (J-1)*VRT(IRPE,ISPIN2) + 1,
     &              VRT(IRPE,ISPIN2),2,IRPEJ,LISWVC)
      ENDIF
C
      CALL XGEMM('N','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),VRT(IRPE,ISPIN2),
     &           ONE,
     &           CORE(I000 + IADT2(IRPE) - 1),VRT(IRPA,ISPIN1),
     &           CORE(I040),IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
  220 CONTINUE
C
C-----------------------------------------------------------------------
C               EC
C     <AB//EJ> T
C               IK
C
C      AA  AA   AB
C               AB
C-----------------------------------------------------------------------
C
      IRPEC = DIRPRD(IRREPT,IRPIK)
C
      DO  310 IRPC=1,NIRREP
        IRPE = DIRPRD(IRPEC,IRPC)
        LENT2(IRPC) = VRT(IRPC,ISPIN2) * VRT(IRPE,ISPIN1)
          IF(IRPC.EQ.1)THEN
            IADT2(IRPC) = 1
          ELSE
            IADT2(IRPC) = IADT2(IRPC-1) + LENT2(IRPC-1)
          ENDIF
  310 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEC,13)
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(I000),IK,1,1,IRPIK,LIST2OFF + 3)
C
      DO  320 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 320
C
      IRPC   = DIRPRD(IRPE,IRPEC)
      IRPEJ = DIRPRD(IRPE,IRPJ)
      IRPAB = DIRPRD(IRREPW,IRPEJ)
C
        IF(IUHF.EQ.0)THEN
         I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
         I030 = I020 + IRPDPD(IRPAB,13)     * VRT(IRPE,ISPIN1)
         NEED = IINTFP * I030
        ELSE
         I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
         NEED = IINTFP * I020
        ENDIF
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
C
        CALL GETLST(CORE(I020),
     &              IOFFVO(IRPJ,IRPEJ,1) + (J-1)*VRT(IRPE,ISPIN1) + 1,
     &              VRT(IRPE,ISPIN1),2,IRPEJ,LISWVD)
        CALL ASSYM3(IRPAB,VRT(1,ISPIN1),IRPDPD(IRPAB,ISPIN1),
     &              IRPDPD(IRPAB,13),
     &              VRT(IRPE,ISPIN1),CORE(I010),CORE(I020),IOFFVV,1)
      ELSE
C
        CALL GETLST(CORE(I010),
     &         IOFFVO(IRPJ,IRPEJ,ISPIN1) + (J-1)*VRT(IRPE,ISPIN1) + 1,
     &              VRT(IRPE,ISPIN1),2,IRPEJ,LISWVA)
      ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),VRT(IRPE,ISPIN1),
     &           ONE,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IADT2(IRPC) - 1),VRT(IRPE,ISPIN1),
     &            ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C  
  320 CONTINUE
C
C-----------------------------------------------------------------------
C                 EC
C     - <AB//EI> T
C                 JK
C
C        AA  AA   AB
C                 AB
C-----------------------------------------------------------------------
C
      IRPEC = DIRPRD(IRREPT,IRPJK)
C
      DO  410 IRPC=1,NIRREP
        IRPE = DIRPRD(IRPEC,IRPC)
        LENT2(IRPC) = VRT(IRPC,ISPIN2) * VRT(IRPE,ISPIN1)
          IF(IRPC.EQ.1)THEN
            IADT2(IRPC) = 1
          ELSE
            IADT2(IRPC) = IADT2(IRPC-1) + LENT2(IRPC-1)
          ENDIF
  410 CONTINUE
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPEC,13)
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(I000),JK,1,1,IRPJK,LIST2OFF + 3)
C
      DO  420 IRPE=1,NIRREP
C
      IF(VRT(IRPE,ISPIN1).EQ.0) GOTO 420
C
      IRPC   = DIRPRD(IRPE,IRPEC)
      IRPEI = DIRPRD(IRPE,IRPI)
      IRPAB = DIRPRD(IRREPW,IRPEI)
C
        IF(IUHF.EQ.0)THEN
         I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
         I030 = I020 + IRPDPD(IRPAB,13)     * VRT(IRPE,ISPIN1)
         NEED = IINTFP * I030
        ELSE
         I020 = I010 + IRPDPD(IRPAB,ISPIN1) * VRT(IRPE,ISPIN1)
         NEED = IINTFP * I020
        ENDIF
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      IF(IUHF.EQ.0)THEN
C
        CALL GETLST(CORE(I020),
     &              IOFFVO(IRPI,IRPEI,1) + (I-1)*VRT(IRPE,ISPIN1) + 1,
     &              VRT(IRPE,ISPIN1),2,IRPEI,LISWVD)
        CALL ASSYM3(IRPAB,VRT(1,ISPIN1),IRPDPD(IRPAB,ISPIN1),
     &              IRPDPD(IRPAB,13),
     &              VRT(IRPE,ISPIN1),CORE(I010),CORE(I020),IOFFVV,1)
      ELSE
C
        CALL GETLST(CORE(I010),
     &         IOFFVO(IRPI,IRPEI,ISPIN1) + (I-1)*VRT(IRPE,ISPIN1) + 1,
     &              VRT(IRPE,ISPIN1),2,IRPEI,LISWVA)
      ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),VRT(IRPE,ISPIN1),
     &           ONEM,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IADT2(IRPC) - 1),VRT(IRPE,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C  
  420 CONTINUE
C
C
C-----------------------------------------------------------------------
C               BC               AC
C     <IJ//MA> T    -  <IJ//MB> T
C               MK               MK
C
C      AA  AA   AB      AA  AA   AB
C               AB               AB
C-----------------------------------------------------------------------
C
      IRPMA = DIRPRD(IRREPW,IRPIJ)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,9)
C
      IF(IUHF.EQ.0)THEN
C
        IJ = IOFFOO(IRPJ,IRPIJ,5) + (J-1) * POP(IRPI,1) + I
        JI = IOFFOO(IRPI,IRPIJ,5) + (I-1) * POP(IRPJ,1) + J
        CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LISWOD)
        CALL GETLST(CORE(I010),JI,1,2,IRPIJ,LISWOD)
        CALL   VADD(CORE(I000),CORE(I000),CORE(I010),IRPDPD(IRPMA,9),
     &              ONEM)
C
      ELSE
C
        IF(IRPI.EQ.IRPJ)THEN
          IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + INDEX(J-1) + I
        ELSE
          IJ = IOFFOO(IRPJ,IRPIJ,ISPIN1) + (J-1) * POP(IRPI,1) + I
        ENDIF
          CALL GETLST(CORE(I000),IJ,1,2,IRPIJ,LISWOA)
C
      ENDIF
C
      DO  510 IRPM=1,NIRREP
        IRPMK = DIRPRD(IRPM,IRPK)
        IRPBC = DIRPRD(IRREPT,IRPMK)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,1)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
  510 CONTINUE
C
      IOFFV = 0
      DO  520 IRPA=1,NIRREP
C
      IRPM  = DIRPRD(IRPA,IRPMA)
      IRPMK = DIRPRD(IRPM,IRPK)
      IRPBC = DIRPRD(IRREPT,IRPMK)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 520
C
      I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN1)
      NEED = I020 * IINTFP 
C
      IF(NEED.GT.MAXCOR)THEN
        WRITE(6,1010)
        CALL INSMEM('GT3WT22',NEED,MAXCOR)
      ENDIF
C
      CALL GETLST(CORE(I010),
     &            IOFFOO(IRPK,IRPMK,5) + (K-1)*POP(IRPM,ISPIN1) + 1,
     &            POP(IRPM,ISPIN1),1,IRPMK,LIST2OFF + 3)
C
C     T matrix is ordered (Bc,M), V matrix is (M,A).
C
      CALL XGEMM('T','T',
     &           VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN1),
     &           ONE,
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           CORE(I010)        ,IRPDPD(IRPBC,13),ONE,
     &           W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPA,ISPIN1)
  520 CONTINUE
C
C
C-----------------------------------------------------------------------
C               BC               AC
C     <JK//AM> T    -  <JK//BM> T
C               IM               IM
C
C      AB  AB   AB      AB  AB   AB
C               AB               AB
C-----------------------------------------------------------------------
C
C     For IUHF=0, integrals are from list LWIC14 and are (m,A), and the
C     record is Kj rather than Jk.
C     For IUHF=1, integrals are from list LWIC13 and are (A,m).
C
      IRPMA = DIRPRD(IRREPW,IRPJK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,11)
C
      IF(IUHF.EQ.1)THEN
C
        JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,1) + J
        CALL GETLST(CORE(I000),JK,1,2,IRPJK,LISWOC)
C
        IOFFV = 0
        DO  620 IRPM=1,NIRREP
C
        IRPMI = DIRPRD(IRPM,IRPI)
        IRPBC = DIRPRD(IRREPT,IRPMI)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
C
        IRPA = DIRPRD(IRPM,IRPMA)
        IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 620
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
          IF(NEED.GT.MAXCOR)THEN
            WRITE(6,1010)
            CALL INSMEM('GT3WT22',NEED,MAXCOR)
          ENDIF
C
        DO  610 M=1,POP(IRPM,ISPIN2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPM,IRPMI,5) + (M-1)*POP(IRPI,ISPIN1) + I,
     &              1,1,IRPMI,LIST2OFF + 3)
  610   CONTINUE
C
        CALL XGEMM('N','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONE,
     &             CORE(I000 + IOFFV),VRT(IRPA,ISPIN1),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  620   CONTINUE
C
      ELSE
C
C     RHF case (IUHF=0)
C
        JK = IOFFOO(IRPJ,IRPJK,5) + (J-1)*POP(IRPK,2) + K
        CALL GETLST(CORE(I000),JK,1,2,IRPJK,LISWOD)
C
        DO  630 IRPM=1,NIRREP
          IRPMI = DIRPRD(IRPM,IRPI)
          IRPBC = DIRPRD(IRREPT,IRPMI)
          LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
            IF(IRPM.EQ.1)THEN
              IADT2(IRPM) = 1
            ELSE
              IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
            ENDIF
  630   CONTINUE
C
        IOFFV = 0
        DO  650 IRPA=1,NIRREP
C
        IRPM  = DIRPRD(IRPA,IRPMA)
        IRPMI = DIRPRD(IRPM,IRPI)
        IRPBC = DIRPRD(IRREPT,IRPMI)
C
        IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 650
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        DO  640 M=1,POP(IRPM,2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPM,IRPMI,5) + (M-1)*POP(IRPI,ISPIN1) + I,
     &            1,1,IRPMI,LIST2OFF + 3)
  640   CONTINUE
C
        CALL XGEMM('T','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONE,
     &             CORE(I000 + IOFFV),POP(IRPM,ISPIN2),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  650   CONTINUE
C
      ENDIF
C
C
C-----------------------------------------------------------------------
C                 BC               AC
C     - <IK//AM> T    +  <IK//BM> T
C                 JM               JM
C
C        AB  AB   AB      AB  AB   AB
C                 AB               AB
C-----------------------------------------------------------------------
C
C     For IUHF=0, integrals are from list LWIC14 and are (m,A), and the
C     record is Ki rather than Ik.
C     For IUHF=1, integrals are from list LWIC13 and are (A,m).
C
      IRPMA = DIRPRD(IRREPW,IRPIK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMA,11)
C
      IF(IUHF.EQ.1)THEN
C
        IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,1) + I
        CALL GETLST(CORE(I000),IK,1,2,IRPIK,LISWOC)
C
        IOFFV = 0
        DO  720 IRPM=1,NIRREP
C
        IRPMJ = DIRPRD(IRPM,IRPJ)
        IRPBC = DIRPRD(IRREPT,IRPMJ)
        LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,2)
          IF(IRPM.EQ.1)THEN
            IADT2(IRPM) = 1
          ELSE
            IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
          ENDIF
C
        IRPA = DIRPRD(IRPM,IRPMA)
        IF(POP(IRPM,2).EQ.0.OR.VRT(IRPA,1).EQ.0) GOTO 720
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
          IF(NEED.GT.MAXCOR)THEN
            WRITE(6,1010)
            CALL INSMEM('GT3WT22',NEED,MAXCOR)
          ENDIF
C
        DO  710 M=1,POP(IRPM,ISPIN2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPM,IRPMJ,5) + (M-1)*POP(IRPJ,ISPIN1) + J,
     &              1,1,IRPMJ,LIST2OFF + 3)
  710   CONTINUE
C
        CALL XGEMM('N','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONEM,
     &             CORE(I000 + IOFFV),VRT(IRPA,ISPIN1),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  720   CONTINUE
C
      ELSE
C
C     RHF case (IUHF=0)
C
        IK = IOFFOO(IRPI,IRPIK,5) + (I-1)*POP(IRPK,2) + K
        CALL GETLST(CORE(I000),IK,1,2,IRPIK,LISWOD)
C
        DO  730 IRPM=1,NIRREP
          IRPMJ = DIRPRD(IRPM,IRPJ)
          IRPBC = DIRPRD(IRREPT,IRPMJ)
          LENT2(IRPM) = IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
            IF(IRPM.EQ.1)THEN
              IADT2(IRPM) = 1
            ELSE
              IADT2(IRPM) = IADT2(IRPM-1) + LENT2(IRPM-1)
            ENDIF
  730   CONTINUE
C
        IOFFV = 0
        DO  750 IRPA=1,NIRREP
C
        IRPM  = DIRPRD(IRPA,IRPMA)
        IRPMJ = DIRPRD(IRPM,IRPJ)
        IRPBC = DIRPRD(IRREPT,IRPMJ)
C
        IF(POP(IRPM,ISPIN2).EQ.0.OR.VRT(IRPA,ISPIN1).EQ.0) GOTO 750
C
        I020 = I010 + IRPDPD(IRPBC,13) * POP(IRPM,ISPIN2)
        NEED = IINTFP * I020
C
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        DO  740 M=1,POP(IRPM,2)
        CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPBC,13)),
     &              IOFFOO(IRPM,IRPMJ,5) + (M-1)*POP(IRPJ,ISPIN1) + J,
     &            1,1,IRPMJ,LIST2OFF + 3)
  740   CONTINUE
C
        CALL XGEMM('T','T',
     &             VRT(IRPA,ISPIN1),IRPDPD(IRPBC,13),POP(IRPM,ISPIN2),
     &             ONEM,
     &             CORE(I000 + IOFFV),POP(IRPM,ISPIN2),
     &             CORE(I010        ),IRPDPD(IRPBC,13),ONE,
     &             W(IADW(IRPBC)),VRT(IRPA,ISPIN1))
C
        IOFFV = IOFFV + VRT(IRPA,ISPIN1) * POP(IRPM,ISPIN2)
  750   CONTINUE
C
      ENDIF
C
C
C-----------------------------------------------------------------------
C        AB
C     - T   <JK//MC>
C        IM
C
C        AA  AB  AB
C        AA
C-----------------------------------------------------------------------
C
      IRPMC = DIRPRD(IRREPW,IRPJK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMC,12)
      JK = IOFFOO(IRPK,IRPJK,5) + (K-1)*POP(IRPJ,ISPIN1) + J
      CALL GETLST(CORE(I000),JK,1,2,IRPJK,LISWOD)
C
      IOFFV = 0
      DO  850 IRPC=1,NIRREP
C
      IRPM = DIRPRD(IRPC,IRPMC)
C
      IRPIM = DIRPRD(IRPM,IRPI)
      IRPAB = DIRPRD(IRREPT,IRPIM)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 850
C
      IF(IUHF.EQ.0)THEN
C
        I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
        I030 = I020 + IRPDPD(IRPAB,13)
        I040 = I030 + IRPDPD(IRPAB,13)
        NEED = IINTFP * I040
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        DO  810    M=1,POP(IRPM,ISPIN1)
C
        IMAB = IOFFOO(IRPM,IRPIM,5) + (M-1)*POP(IRPI,ISPIN1) + I
        MIAB = IOFFOO(IRPI,IRPIM,5) + (I-1)*POP(IRPM,ISPIN1) + M
        CALL GETLST(CORE(I020),IMAB,1,1,IRPIM,LIST2OFF + 3)
        CALL GETLST(CORE(I030),MIAB,1,1,IRPIM,LIST2OFF + 3)
        CALL   VADD(CORE(I020),CORE(I020),CORE(I030),
     &              IRPDPD(IRPAB,13),ONEM)
        CALL  SQSYM(IRPAB,VRT(1,ISPIN1),IRPDPD(IRPAB,ISPIN1),
     &              IRPDPD(IRPAB,13),1,
     &              CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),CORE(I020))
C
        IF(IRPI.EQ.IRPM.AND.I.EQ.M)THEN
          CALL ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &              IRPDPD(IRPAB,ISPIN1))
        ENDIF
C
  810   CONTINUE
C
      ELSE
C
        I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
        NEED = IINTFP * I020
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        IF(IRPI.GT.IRPM)THEN
          IMOFF = IOFFOO(IRPI,IRPIM,ISPIN1) + (I-1)*POP(IRPM,ISPIN1) + 1
          CALL GETLST(CORE(I010),IMOFF,POP(IRPM,ISPIN1),1,IRPIM,
     &                LIST2OFF + ISPIN1)
          CALL VMINUS(CORE(I010),IRPDPD(IRPAB,ISPIN1)*POP(IRPM,ISPIN1))
        ENDIF
C
        IF(IRPI.LT.IRPM)THEN
          DO  820 M=1,POP(IRPM,ISPIN1)
          IM = IOFFOO(IRPM,IRPIM,ISPIN1) + (M-1)*POP(IRPI,ISPIN1) + I
          CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                IM,1,1,IRPIM,LIST2OFF + ISPIN1)
  820     CONTINUE
        ENDIF
C
        IF(IRPI.EQ.IRPM)THEN
          DO  830 M=1,POP(IRPM,ISPIN1)
C
            IF(I.GT.M)THEN
              IM = IOFFOO(IRPI,IRPIM,ISPIN1) + INDEX(I-1) + M
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),IM,
     &                    1,1,IRPIM,LIST2OFF + ISPIN1)
              CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
C
            IF(I.LT.M)THEN
              IM = IOFFOO(IRPM,IRPIM,ISPIN1) + INDEX(M-1) + I
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),IM,
     &                    1,1,IRPIM,LIST2OFF + ISPIN1)
             ENDIF
C
            IF(I.EQ.M)THEN
              CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
  830     CONTINUE
        ENDIF
C
      ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),POP(IRPM,ISPIN1),
     &           ONEM,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPC,ISPIN2)
  850 CONTINUE
C
C
C-----------------------------------------------------------------------
C      AB
C     T   <IK//MC>
C      JM
C
C      AA  AB  AB
C      AA
C-----------------------------------------------------------------------
C
      IRPMC = DIRPRD(IRREPW,IRPIK)
C
      I000 = 1
      I010 = I000 + IRPDPD(IRPMC,12)
      IK = IOFFOO(IRPK,IRPIK,5) + (K-1)*POP(IRPI,ISPIN1) + I
      CALL GETLST(CORE(I000),IK,1,2,IRPIK,LISWOD)
C
      IOFFV = 0
      DO  950 IRPC=1,NIRREP
C
      IRPM = DIRPRD(IRPC,IRPMC)
C
      IRPJM = DIRPRD(IRPM,IRPJ)
      IRPAB = DIRPRD(IRREPT,IRPJM)
C
      IF(POP(IRPM,ISPIN1).EQ.0.OR.VRT(IRPC,ISPIN2).EQ.0) GOTO 950
C
      IF(IUHF.EQ.0)THEN
C
        I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
        I030 = I020 + IRPDPD(IRPAB,13)
        I040 = I030 + IRPDPD(IRPAB,13)
        NEED = IINTFP * I040
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        DO  910    M=1,POP(IRPM,ISPIN1)
C
        JMAB = IOFFOO(IRPM,IRPJM,5) + (M-1)*POP(IRPJ,ISPIN1) + J
        MJAB = IOFFOO(IRPJ,IRPJM,5) + (J-1)*POP(IRPM,ISPIN1) + M
        CALL GETLST(CORE(I020),JMAB,1,1,IRPJM,LIST2OFF + 3)
        CALL GETLST(CORE(I030),MJAB,1,1,IRPJM,LIST2OFF + 3)
        CALL   VADD(CORE(I020),CORE(I020),CORE(I030),
     &              IRPDPD(IRPAB,13),ONEM)
        CALL  SQSYM(IRPAB,VRT(1,ISPIN1),IRPDPD(IRPAB,ISPIN1),
     &              IRPDPD(IRPAB,13),1,
     &              CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),CORE(I020))
C
        IF(IRPJ.EQ.IRPM.AND.J.EQ.M)THEN
          CALL ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &              IRPDPD(IRPAB,ISPIN1))
        ENDIF
C
  910   CONTINUE
C
      ELSE
C
        I020 = I010 + IRPDPD(IRPAB,ISPIN1) * POP(IRPM,ISPIN1)
        NEED = IINTFP * I020
        IF(NEED.GT.MAXCOR)THEN
          WRITE(6,1010)
          CALL INSMEM('GT3WT22',NEED,MAXCOR)
        ENDIF
C
        IF(IRPJ.GT.IRPM)THEN
          JMOFF = IOFFOO(IRPJ,IRPJM,ISPIN1) + (J-1)*POP(IRPM,ISPIN1) + 1
          CALL GETLST(CORE(I010),JMOFF,POP(IRPM,ISPIN1),1,IRPJM,
     &                LIST2OFF + ISPIN1)
          CALL VMINUS(CORE(I010),IRPDPD(IRPAB,ISPIN1)*POP(IRPM,ISPIN1))
        ENDIF
C
        IF(IRPJ.LT.IRPM)THEN
          DO  920 M=1,POP(IRPM,ISPIN1)
          JM = IOFFOO(IRPM,IRPJM,ISPIN1) + (M-1)*POP(IRPJ,ISPIN1) + J
          CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                JM,1,1,IRPJM,LIST2OFF + ISPIN1)
  920     CONTINUE
        ENDIF
C
        IF(IRPJ.EQ.IRPM)THEN
          DO  930 M=1,POP(IRPM,ISPIN1)
C
            IF(J.GT.M)THEN
              JM = IOFFOO(IRPJ,IRPJM,ISPIN1) + INDEX(J-1) + M
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),JM,
     &                    1,1,IRPJM,LIST2OFF + ISPIN1)
              CALL VMINUS(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
C
            IF(J.LT.M)THEN
              JM = IOFFOO(IRPM,IRPJM,ISPIN1) + INDEX(M-1) + J
              CALL GETLST(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),JM,
     &                    1,1,IRPJM,LIST2OFF + ISPIN1)
             ENDIF
C
            IF(J.EQ.M)THEN
              CALL   ZERO(CORE(I010 + (M-1)*IRPDPD(IRPAB,ISPIN1)),
     &                    IRPDPD(IRPAB,ISPIN1))
            ENDIF
  930     CONTINUE
        ENDIF
C
      ENDIF
C
      CALL XGEMM('N','N',
     &           IRPDPD(IRPAB,ISPIN1),VRT(IRPC,ISPIN2),POP(IRPM,ISPIN1),
     &           ONE,
     &           CORE(I010),IRPDPD(IRPAB,ISPIN1),
     &           CORE(I000 + IOFFV),POP(IRPM,ISPIN1),
     &           ONE,
     &           T3(IADT3(IRPC)),IRPDPD(IRPAB,ISPIN1))
C
      IOFFV = IOFFV + POP(IRPM,ISPIN1) * VRT(IRPC,ISPIN2)
  950 CONTINUE
C
      RETURN
 1010 FORMAT(' @GT3WT22-F, Insufficient memory to continue. ')
      END
