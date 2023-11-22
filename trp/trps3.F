      sUBROUTINE TRPS3(T1A,T1B,FOVA,FOVB,D1T1A,D1T1B,S1A,S1B,
     1                 EVAL,SCR1,SCR2,SCR3,
     1                 OOOV,OOOV2,GOOOV,BUF,
     1                 T2BAB,T2BAAI,T2BAAJ,T2BBBI,T2BBBJ,
     1                 IADBLK,LENBLK,IADT2,LENT2,
     1                 IADV,LENV,OOOVAD,LENINT,
     1                 ICORE,
     1                 LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,
     1                 INT1,INT2,NONHF,IUHF,LCCSD_T,FINAL,EXTRNL_CCSD)
      IMPLICIT INTEGER (A-Z)
      LOGICAL INT1,INT2,NONHF
      LOGICAL IJEQL,NONEQL
      LOGICAL CCSDT4,CCSDT,LCCSD_T,FINAL
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      LOGICAL EXTRNL_CCSD
      DOUBLE PRECISION T1A(1),T1B(1),FOVA(1),FOVB(1),
     1                 D1T1A(1),D1T1B(1),S1A(1),S1B(1), DDOT,
     1                 SCR1(1),SCR2(1),SCR3(1),EVAL(NOCA+NVRTA,2)
C     OOOV  ijka integrals/intermediates for D3T3 = WT2
C     OOOV2 ijka integrals/intermediates for D2T2 = WT3 (includes 1b)
C     GOOOV ijka gammas
      DOUBLE PRECISION  OOOV(LNOOOV,4),OOOV2(LNOOOV,4),
     1                 GOOOV(LNOOOV,4),BUF(1)
      DOUBLE PRECISION T2BAB(1),T2BAAI(1),T2BAAJ(1),
     1                 T2BBBI(1),T2BBBJ(1)
      DOUBLE PRECISION ICORE(1)
C
      DOUBLE PRECISION E4TAAA,E4TBBB,E4TBAA,E4TBBA,E5STA,E5STB,
     1                 E5ST2A,E5ST2B,E4STA,E4STB,E4ST2A,E4ST2B
      DOUBLE PRECISION E5TEST,SDOT,E6TEST
      DOUBLE PRECISION DOOA,DOOB,DVVA,DVVB
      DOUBLE PRECISION DIJK
C
#include "maxbasfn.par"
#include "flags.h"
      DIMENSION OOOVAD(8,4),LENINT(8,4)
      DIMENSION LEN(8,8),IADBLK(8),LENBLK(8)
      DIMENSION                        IADT2(8) ,LENT2(8),
     1                                 IADV(8)  ,LENV(8)
      DIMENSION IADW(8),LENW(8)
      DIMENSION IADW2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      EQUIVALENCE(ICLLVL,IFLAGS( 2))
      EQUIVALENCE(IDRLVL,IFLAGS( 3))
      EQUIVALENCE(IREFNC,IFLAGS(11))
      EQUIVALENCE(IQRHFP,IFLAGS(32))
      EQUIVALENCE(IQRHFM,IFLAGS(33))
      EQUIVALENCE(INONHF,IFLAGS(38))
      EQUIVALENCE(IORBTP,IFLAGS(39))
C
C     T3 ARRAYS
C
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXBASFN,8,2),ABSOCC(MAXBASFN,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TBAA,E4TBBA,E5STA,E5STB,
     1                E4STA ,E4STB
      COMMON /DENST3/ DOOA(MAXBASFN),DOOB(MAXBASFN),
     1                DVVA(MAXBASFN),DVVB(MAXBASFN)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
      COMMON /TESTEN/ E5TEST,E6TEST
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      COMMON /T3METH/ CCSDT4,CCSDT
C
      INDEX(I) = I*(I-1)/2
C
C
      WRITE(6,1015)
 1015 FORMAT(' @TRPS3-I, Spin case BBA ')
C
C     Initialize D1T1 and "S1" increments.
C
      CALL ZERO(D1T1A,NTAA)
      CALL ZERO(D1T1B,NTBB)
      CALL ZERO(S1A  ,NTAA)
      CALL ZERO(S1B  ,NTBB)
#ifdef _DEBUG_LVL0
      iprint=1
#endif
C
C     Get eigenvalues (more precisely, diagonal Fock matrix elements)
C
      CALL GETREC(20,'JOBARC','SCFEVALA',(NOCA+NVRTA)*IINTFP,
     1            EVAL(1,1))
      CALL GETREC(20,'JOBARC','SCFEVALB',(NOCB+NVRTB)*IINTFP,
     1            EVAL(1,2))
C
C     Determine the lengths and dimensions of the NIRREP**2 T2 BBA arrays.
C
C      WRITE(6,1070)
 1070 FORMAT(' @TRPS3-I, TRIPLES SYMMETRY BLOCK INFORMATION ')
      DO  100 IRPIJK=1,NIRREP
      DO   90 IRPC  =1,NIRREP
      IRPAB = DIRPRD(IRPIJK,IRPC)
      LEN(IRPC,IRPIJK) = IRPDPD(IRPAB,2) * VRT(IRPC,1)
C      WRITE(6,1080) IRPIJK,IRPC,IRPAB,LEN(IRPC,IRPIJK)
 1080 FORMAT('   FOR IJK BLOCK ',I4,/,
     1       '       C   BLOCK ',I4,' AB BLOCK ',I4,' LENGTH IS ',I6)
   90 CONTINUE
  100 CONTINUE
C
C     Compute addresses and lengths of OOOV quantities.
C
      DO  410 ISPIN=1,4
      DO  400 IRREP=1,NIRREP
      IF(ISPIN.EQ.1)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP, 3) * IRPDPD(IRREP,16)
      ENDIF
      IF(ISPIN.EQ.2)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP, 4) * IRPDPD(IRREP,17)
      ENDIF
      IF(ISPIN.EQ.3)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP,14) * IRPDPD(IRREP,11)
      ENDIF
      IF(ISPIN.EQ.4)THEN
      LENINT(IRREP,ISPIN) = IRPDPD(IRREP,14) * IRPDPD(IRREP,18)
      ENDIF
  400 CONTINUE
  410 CONTINUE
C
      OOOVAD(1,1) = 1
      OOOVAD(1,2) = 1
      OOOVAD(1,3) = 1
      OOOVAD(1,4) = 1
      DO  430 ISPIN=1,4
      IF(NIRREP.GT.1)THEN
      DO  420 IRREP = 2,NIRREP
      OOOVAD(IRREP,ISPIN) = OOOVAD(IRREP-1,ISPIN) +
     1                      LENINT(IRREP-1,ISPIN)
  420 CONTINUE
      ENDIF
  430 CONTINUE
C      WRITE(6,1200)
 1200 FORMAT(' @TRPS3-I, OOOV ARRAY LENGTHS AND ADDRESSES ')
      DO  450 ISPIN=1,4
      DO  440 IRREP=1,NIRREP
C      WRITE(6,1210) ISPIN,IRREP,LENINT(IRREP,ISPIN),OOOVAD(IRREP,ISPIN)
 1210 FORMAT(' SPIN, IRREP, LENGTH, ADDRESS ',2I4,I6,I6)
  440 CONTINUE
  450 CONTINUE
C
      DO  460 IRREP=1,NIRREP
      CALL GETTRN(OOOV(OOOVAD(IRREP,4),4),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LWIC14)
C    1            2,IRREP,10)
      CALL GETTRN(OOOV(OOOVAD(IRREP,2),2),BUF,
     1            IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     1            2,IRREP,LWIC12)
C    1            2,IRREP,8)
      CALL GETTRN(OOOV(OOOVAD(IRREP,3),3),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LWIC13)
C    1            2,IRREP,9)
  460 CONTINUE
C
      DO  471 IRPIJ=1,NIRREP
      CALL SYMTR3(IRPIJ,POP(1,2),VRT(1,2),IRPDPD(IRPIJ,17),
     1            IRPDPD(IRPIJ,4),OOOV(OOOVAD(IRPIJ,2),2),
     1            SCR1,SCR2,SCR3)
      CALL SYMTR3(IRPIJ,VRT(1,1),POP(1,2),IRPDPD(IRPIJ,11),
     1            IRPDPD(IRPIJ,14),OOOV(OOOVAD(IRPIJ,3),3),
     1            SCR1,SCR2,SCR3)
      CALL SYMTR3(IRPIJ,POP(1,1),VRT(1,2),IRPDPD(IRPIJ,12),
     1            IRPDPD(IRPIJ,14),OOOV(OOOVAD(IRPIJ,4),4),
     1            SCR1,SCR2,SCR3)
  471 CONTINUE
C
      IF(INT2)THEN
C
      DO  9460 IRREP=1,NIRREP
      CALL GETTRN(OOOV2(OOOVAD(IRREP,4),4),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LWIC24)
      CALL GETTRN(OOOV2(OOOVAD(IRREP,2),2),BUF,
     1            IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     1            2,IRREP,LWIC22)
      CALL GETTRN(OOOV2(OOOVAD(IRREP,3),3),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LWIC23)
 9460 CONTINUE
C
      DO  9471 IRPIJ=1,NIRREP
      CALL SYMTR3(IRPIJ,POP(1,2),VRT(1,2),IRPDPD(IRPIJ,17),
     1            IRPDPD(IRPIJ,4),OOOV2(OOOVAD(IRPIJ,2),2),
     1            SCR1,SCR2,SCR3)
      CALL SYMTR3(IRPIJ,POP(1,1),VRT(1,2),IRPDPD(IRPIJ,12),
     1            IRPDPD(IRPIJ,14),OOOV2(OOOVAD(IRPIJ,4),4),
     1            SCR1,SCR2,SCR3)
 9471 CONTINUE
      ENDIF
C
      IF(IDRLVL.GT.0)THEN
C
      DO  9461 IRREP=1,NIRREP
      CALL GETTRN(GOOOV(OOOVAD(IRREP,4),4),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LISGO4)
      CALL GETTRN(GOOOV(OOOVAD(IRREP,2),2),BUF,
     1            IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     1            2,IRREP,LISGO2)
      CALL GETTRN(GOOOV(OOOVAD(IRREP,3),3),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LISGO3)
      CALL SYMTR3(IRREP,VRT(1,1),POP(1,2),IRPDPD(IRREP,11),
     1            IRPDPD(IRREP,14),GOOOV(OOOVAD(IRREP,3),3),
     1            SCR1,SCR2,SCR3)
 9461 CONTINUE
      ENDIF
      IF (EXTRNL_CCSD) THEN
         
         DO IRREP=1,NIRREP
            CALL GETTRN(GOOOV(OOOVAD(IRREP,4),4),BUF,
     &                  IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     &                  2,IRREP,110)
            CALL GETTRN(GOOOV(OOOVAD(IRREP,2),2),BUF,
     &                  IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     &                  2,IRREP,108)
            CALL GETTRN(GOOOV(OOOVAD(IRREP,3),3),BUF,
     &                  IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     &                  2,IRREP,109)
            CALL SYMTR3(IRREP,VRT(1,1),POP(1,2),IRPDPD(IRREP,11),
     &                  IRPDPD(IRREP,14),GOOOV(OOOVAD(IRREP,3),3),
     &                  SCR1,SCR2,SCR3)
         ENDDO
      ENDIF 
C
C     Compute number of free double words.
      MCORE = NFREE / IINTFP
cYAU      WRITE(6,1230) NFREE,MCORE,IINTFP
 1230 FORMAT(' @TRPS3-I, Number of free integer words ',I10,/,
     1       '           Number of free double  words ',I10,/,
     1       '           Ratio                        ',I10)
C
      ISPIN1 = 2
      ISPIN2 = 1
C
      DO 1000 IRPIJK=1,NIRREP
C
      IJKVAL = 0
C
      DO  990 IRPK  =1,NIRREP
      IF(POP(IRPK,ISPIN2).EQ.0) GOTO 990
C
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN2)
C
      DO  980 IRPJ=1,NIRREP
      IF(POP(IRPJ,ISPIN1).EQ.0) GOTO 980
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IF(IRPI.GT.IRPJ) GOTO 980
C
      IF(POP(IRPI,ISPIN1).EQ.0) GOTO 980
C
      IJEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ)THEN
      IJEQL  = .TRUE.
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      NIJ = (POP(IRPJ,ISPIN1) * (POP(IRPJ,ISPIN1)-1))/2
      ELSE
      NONEQL = .TRUE.
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN1)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN1)
      NIJ =  POP(IRPI,ISPIN1) * POP(IRPJ,ISPIN1)
      ENDIF
C
      NIK = POP(IRPI,ISPIN1) * POP(IRPK,ISPIN2)
      NJK = POP(IRPJ,ISPIN1) * POP(IRPK,ISPIN2)
C
      IF(IJEQL.AND.POP(IRPJ,ISPIN1).LT.2) GOTO 980
C
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
C     COMPUTE ADDRESSES AND LENGTHS OF ABC SYMMETRY BLOCKS FOR
C     THIS IJK. THIS LOOP IS EFFECTIVELY OVER IRREPS OF C (IE
C     COLUMNS OF THE EVENTUAL TARGET).
C
C     MAIN CORE VECTOR ADDRESSING : (DOUBLE PRECISION WORDS)
C
C           I000   SYMMETRY PACKED T3 FOR THIS IJK BLOCK
C           I010                   W
C           I020   UNEXPANDED T2(  ,IJ) BBBB FOR THIS IJ BLOCK
C           I030   EXPANDED   T2(  ,IJ) BBBB
C           I040   T2(  ,JK) ABAB FOR THIS JK BLOCK
C           I050   T2(  ,IK) ABAB FOR THIS IK BLOCK
C           ISTART CURRENT STARTING ADDRESS OF REUSABLE CORE
C
      I000 = 1
      DO  600 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADBLK(IRREP) = 1
      ELSE
      IADBLK(IRREP) = IADBLK(IRREP-1) + LEN(IRREP-1,IRPIJK)
      ENDIF
      LENBLK(IRREP) = LEN(IRREP,IRPIJK)
  600 CONTINUE
C
      LENABC=0
      DO  610 IRREP=1,NIRREP
      LENABC = LENABC + LENBLK(IRREP)
  610 CONTINUE
C
      I010 = I000 + LENABC
      DO  630 IRREP=1,NIRREP
      IRPA = DIRPRD(IRPIJK,IRREP)
      LENW(IRREP) = VRT(IRPA,2) * IRPDPD(IRREP,13)
  630 CONTINUE
C
      DO  640 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADW(IRREP) = 1
      ELSE
      IADW(IRREP) = IADW(IRREP-1) + LENW(IRREP-1)
      ENDIF
  640 CONTINUE
C
      IOFFT2 = I010
      DO  641 IRREP=1,NIRREP
      IOFFT2 = IOFFT2 + LENW(IRREP)
  641 CONTINUE
C
      I020 = IOFFT2
C     GET VARIOUS "SMALL" BLOCKS OF T2 AND <IJ//KA> IN CORE, ORGANIZED
C     APPROPRIATELY.
C
C     BLOCK OF BBBB T(  ,IJ) WHERE I AND J ARE DIFFERENT SYMMETRY
C     BLOCKS. LENGTH IS IRPDPD(IRPIJ,2) * NIJ.
C
      CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,IRPIJ,45)
C
C     EXPAND INTO A DIFFERENT REGION OF CORE SO WE KEEP UNTOUCHED T2
C     AS WELL AS EXPANDED ONE. LENGTH IS DSZEXP * NIJ.
C
      DSZ   = IRPDPD(IRPIJ,2)
      DSZEXP= IRPDPD(IRPIJ,18+ISPIN1)
      I030 = I020 + DSZ    * NIJ
      I040 = I030 + DSZEXP * NIJ
      CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1             ICORE(I020))
C
C     BLOCK OF ABAB T(  ,JK). LENGTH IS IRPDPD(IRPJK,13) * NJK
C
      I050 = I040 + IRPDPD(IRPJK,13) * NJK
      CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,IRPJK,46)
C
C     BLOCK OF ABAB T(  ,IK). LENGTH IS IRPDPD(IRPIK,13) * NIK
C
      CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1,IRPIK,46)
C
      I060 = I050 + IRPDPD(IRPIK,13) * NIK
C
      LNVOIJ  = IRPDPD(IRPIJ,10)
      LNVOIK  = IRPDPD(IRPIK,11)
      LNVOJK  = IRPDPD(IRPJK,11)
      LNVOIK2 = IRPDPD(IRPIK,12)
      LNVOJK2 = IRPDPD(IRPJK,12)
C
C     SET A STARTING ADDRESS FOR REST OF CORE, THAT PART WHICH CAN BE
C     REUSED AFTER D1T3.
C
      I070 = I060 + LENABC
C
C     If this is a CCSD(T) or QCISD(T) gradient calculation, allocate
C     memory for "disconnected" triple amplitudes.
C
CSSS      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
      I080 = I070 + LENABC
      I090 = I080 + I020 - I010
      I100 = I090 + LENABC
      ISTART = I100
CSSS      ELSE
CSSS      ISTART = I070
CSSS      ENDIF
      IF (LCCSD_T) THEN
C
C IL00 is L3 (symmetry packed)
C
          IL00 = ISTART
          IL10 = IL00 + LENABC
          ISTART=IL10+(I020-I010)
      ENDIF
      NLEFT = MCORE - ISTART
C
      IF(ISTART+LENABC.GE.MCORE)THEN
         WRITE(6,9010)
         CALL INSMEM('TRPS3',ISTART+LENABC,MCORE)
      ENDIF
C
C     COMPUTE DENOMINATOR ARRAY

      CALL MKD33(ICORE(I060),EVAL,IADBLK,IRPIJK,NOCA,NOCB,NVRTA,NVRTB)
C
      DO  649    K=KLOW,KHIGH
      DO  648    J=JLOW,JHIGH
      IF(IJEQL) IHIGH = J-1
      DO  647    I=ILOW,IHIGH
C
      CALL ZERO(ICORE(I000),LENABC)
C
      LENGTHW = 0
      DO 9647 IRREP=1,NIRREP
      LENGTHW = LENGTHW + LENW(IRREP)
 9647 CONTINUE
      CALL ZERO(ICORE(I010),LENGTHW)
C
      IF (.NOT. FINAL) THEN
      CALL D1T3( ICORE(I000),ICORE(I010),
     1           ICORE(I030),ICORE(I040),ICORE(I050),
     1 OOOV(OOOVAD(IRPIJ,2)+IOFFOO(IRPJ,IRPIJ,2)*IRPDPD(IRPIJ,17),2),
     1 OOOV(OOOVAD(IRPJK,4)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,12),4),
     1 OOOV(OOOVAD(IRPIK,4)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,12),4),
     1 OOOV(OOOVAD(IRPJK,3)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,11),3),
     1 OOOV(OOOVAD(IRPIK,3)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,11),3),
     1 ICORE(ISTART),IUHF,
     1 IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1 DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1 IRPDPD(IRPIJ,17),IRPDPD(IRPJK,12),IRPDPD(IRPIK,12),
     1 IRPDPD(IRPJK,11),IRPDPD(IRPIK,11),
     1 SCR1,SCR2,SCR3,.FALSE.)
C
c      IJK = IJKPOS(IRPI,IRPJ,IRPK,2)
c      IOFF = IJKOFF(IJK,IRPIJK,ISPIN1+1)
c      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ + INDEX(J-1) + I
c      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN1) + I
      IJKVAL = IJKVAL + 1
C
C     CCSDT or non-HF (non-canonical) if semicanonical orbitals not
C     used. (ii) Occupied orbital sum.
C
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT323V(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              2,1,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
      ENDIF
      IF(CCSDT4.OR.CCSDT)THEN
C     CCSDT abef intermediate contribution.
      CALL T3WT323P(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              2,1,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
CC     CCSDT maei intermediate contribution.
      CALL T3WT323R(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,2,1,
     1              IADBLK,IADW,IRPI,IRPJ,IRPK,I,J,K,IRPIJK,IUHF)
      ENDIF
C
C     CCSDT or non-HF (non-canonical) if semicanonical orbitals not
C     used. (ii) Occupied orbital sum.
C
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT323O(ICORE(I000),ICORE(ISTART),NFREE,IUHF,2,1,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
      IF(CCSDT4.OR.CCSDT)THEN
      CALL T3WT323HH(ICORE(I000),ICORE(ISTART),NFREE,IUHF,2,1,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
C
      CALL EXPSC3(ICORE(I000),ICORE(I010),IADBLK,IADW,IRPIJK)
C
C     REMOVE DENOMINATORS
C
      DIJK = EVAL(ABSOCC(I,IRPI,2),2) + EVAL(ABSOCC(J,IRPJ,2),2)
     1     + EVAL(ABSOCC(K,IRPK,1),1)
      CALL RMD314(ICORE(I000),ICORE(I060),LENABC,DIJK)
#ifdef _DEBUG_LVLM`
      write(*,*)
      call checksum("ITER-TRPS3-0", ICORE(I000), LENABC)
#endif
      ELSE
C
C If this iterative NONHF-CCSD(T), converged amplitudes are already
C on the disk.
C
         IJKVAL = IJKVAL + 1
#ifdef _NOSKIP
         IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
         CALL ZERO(ICORE(I000),LENABC)
         CALL T3FT323V(ICORE(I000),ICORE(I010),ICORE(ISTART),NFREE,IUHF,
     1              2,1,LENABC,LENGTHW,IRPIJK,IJKVAL,IADBLK,IADW)
         CALL T3FT323O(ICORE(I000),ICORE(ISTART),NFREE,IUHF,2,1,LENABC,
     1              IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
         ENDIF
         IF (ISTART+LENABC .GT. NLEFT) THEN
            CALL INSMEM("@-TRPS3",ISTART+LENABC,NLEFT)
         ENDIF

         CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,2 + 1 + 4)
         E4TBBA = E4TBBA + Ddot(LENABC, ICORE(ISTART), 1, 
     &                            ICORE(I000),1) 
#endif
      CALL ZERO(ICORE(I000),LENABC)
      CALL ZERO(ICORE(I010),LENGTHW)
      CALL D1T3( ICORE(I000),ICORE(I010),
     1           ICORE(I030),ICORE(I040),ICORE(I050),
     1 OOOV(OOOVAD(IRPIJ,2)+IOFFOO(IRPJ,IRPIJ,2)*IRPDPD(IRPIJ,17),2),
     1 OOOV(OOOVAD(IRPJK,4)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,12),4),
     1 OOOV(OOOVAD(IRPIK,4)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,12),4),
     1 OOOV(OOOVAD(IRPJK,3)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,11),3),
     1 OOOV(OOOVAD(IRPIK,3)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,11),3),
     1 ICORE(ISTART),IUHF,
     1 IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1 I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1 DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1 IRPDPD(IRPIJ,17),IRPDPD(IRPJK,12),IRPDPD(IRPIK,12),
     1 IRPDPD(IRPJK,11),IRPDPD(IRPIK,11),
     1 SCR1,SCR2,SCR3,.FALSE.)

      CALL EXPSC3(ICORE(I000),ICORE(I010),IADBLK,IADW,IRPIJK)
      CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,2 + 1 + 4)

CSSS      Write(*,*)
CSSS      Call checksum("T3(C)+T3D-(D)", ICORE(IStart), Lenabc)

      E4TBBA = E4TBBA + Ddot(LENABC, ICORE(ISTART), 1, 
     &                       ICORE(I000),1) 
      E6TEST = E6TEST +  Ddot(LENABC, ICORE(ISTART), 1,
     &                       ICORE(I000),1)

CSSS      WRITE(6,*) ' @TRPS2-I, E6TEST from disconnected triples ',E6TEST
CSSS        CALL E4TIJK(ICORE(ISTART),ICORE(I060),LENABC,DIJK,E4TBBA,
CSSS     &        __FILE__)

      DIJK = EVAL(ABSOCC(I,IRPI,2),2) + EVAL(ABSOCC(J,IRPJ,2),2)
     1     + EVAL(ABSOCC(K,IRPK,1),1)
      CALL RMD314(ICORE(I000),ICORE(I060),LENABC,DIJK)
CSSS      CALL DAXPY(LENABC, -1.0D0, ICORE(I000), 1, ICORE(ISTART), 1)
      CALL DCOPY(LENABC, ICORE(ISTART), 1, ICORE(I070), 1)

      ENDIF 
C
C     COMPUTE T(CCSD) ENERGY
C
C
C For L based CCSD(T), we also need to compute L3(ijk,abc)= L2(ij,{le}a)<bc||{le}k>
C and have both T3 and L3 contributions together. The current implementation
C requires two arrays of o^3v^3 vectors in memory.
C There is the possibility for an ilternative implementation
C T3(ijk,abc)o<ek||bc> = I(ij,ae) and E4 = I(ij,ae)oL(ij,ae). The memory
C requirments for the second schemem is slightly less but it is not as stright
C forward. The alternative implementaion can be considered later if the
C memory requirements become excessive. Ajith Perera, 10/05.
C     
      IF (LCCSD_T) THEN

         CALL ZERO(ICORE(IL00), LENABC)
         CALL ZERO(ICORE(IL10),LENGTHW)
C
C Get various "small" blocks of L2 And <IJ||KA> in core, organized
C appropriately.
C           
C Block of bbbb L2(  ,ij) where I and J are different symmetrY
C blocks. LENGTH is IRPDPD(IRPIJ,2) * NIJ.
C           
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,IRPIJ,
     &               145)
C
C Expand into a different region of core so we keep untouched L2
C as well as expanded one. LENGTH is DSZEXP * NIJ.
C
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
C
C Block of ABAB T(  ,JK). LENGTH is IRPDPD(IRPJK,13) * NJK
C
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,IRPJK,
     &               146)
C
C Block of ABAB T(  ,IK). LENGTH is IRPDPD(IRPIK,13) * NIK
C
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1,IRPIK,
     &               146)

         CALL D1T3( ICORE(IL00),ICORE(IL10),
     1    ICORE(I030),ICORE(I040),ICORE(I050),
     1    OOOV(OOOVAD(IRPIJ,2)+IOFFOO(IRPJ,IRPIJ,2)*IRPDPD(IRPIJ,17),2),
     1    OOOV(OOOVAD(IRPJK,4)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,12),4),
     1    OOOV(OOOVAD(IRPIK,4)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,12),4),
     1    OOOV(OOOVAD(IRPJK,3)+IOFFOO(IRPJ,IRPJK,5)*IRPDPD(IRPJK,11),3),
     1    OOOV(OOOVAD(IRPIK,3)+IOFFOO(IRPI,IRPIK,5)*IRPDPD(IRPIK,11),3),
     1    ICORE(ISTART),IUHF,
     1    IADBLK,LENBLK,IADW,LENW,IADT2,LENT2,IADV,LENV,
     1    I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPJK,IRPIK,IRPIJK,
     1    DSZEXP,IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1    IRPDPD(IRPIJ,17),IRPDPD(IRPJK,12),IRPDPD(IRPIK,12),
     1    IRPDPD(IRPJK,11),IRPDPD(IRPIK,11),
     1    SCR1,SCR2,SCR3,.TRUE.)
         
         CALL EXPSC3(ICORE(IL00),ICORE(IL10),IADBLK,IADW,IRPIJK)

         CALL RMD314(ICORE(IL00),ICORE(I060),LENABC,DIJK)
         
C
         CALL E4TIJKL(ICORE(I000),ICORE(IL00),ICORE(I060),LENABC,
     &                DIJK, E4TBBA,__FILE__)
C
C Reload the T2 vectors back in to the appropriate locations, so any methods
C other than CCSD(T) can  smoothly proceed.
C
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,IRPIJ,
     &               45)
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,IRPJK,
     &               46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1,IRPIK,
     &               46)
C
      ELSE
C
         IF (.NOT. FINAL) 
     &   CALL E4TIJK(ICORE(I000),ICORE(I060),LENABC,DIJK,E4TBBA,
     &        __FILE__)
         IF (EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVL0
         If (iprint .eq.1) write(6,"(2a)") " In trps3.F,Entered",
     &                                     " EXTRNL_CC block"
         iprint=-1
#endif 
         ENDIF 
C This endif is for LCCSD_T
C
      ENDIF 
C
C     If this is a CCSD(T) or QCISD(T) gradient calculation, compute
C     "disconnected" triple amplitudes.  For LCCSD(T) compute
C     disconnected lambda triple amplitude.
C
CAP      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN

      IF((IDRLVL.GT.0.AND.
     1  (ICLLVL.EQ.21.OR.ICLLVL.EQ.22).OR.EXTRNL_CCSD))THEN

      IF (.NOT. FINAL) THEN
      CALL ZERO(ICORE(I070),LENABC)
      CALL ZERO(ICORE(I010),IOFFT2 - I010)
      CALL S1S223(ICORE(ISTART),ICORE(I070),ICORE(I010),IADBLK,IADW,
     1            2,1,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1            IUHF,1,LCCSD_T)

#ifdef _DEBUG_LVLM
      write(*,*)
      call checksum("ITER-TRPS3-1", ICORE(I070), LENABC)
#endif
      IF(NONHF)THEN
      CALL ZERO(ICORE(I090),LENABC)
      CALL ZERO(ICORE(I010),IOFFT2 - I010)
      CALL S1S223(ICORE(ISTART),ICORE(I090),ICORE(I010),IADBLK,IADW,
     1            2,1,I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1            IUHF,2,LCCSD_T)
      CALL VADD(ICORE(I070),ICORE(I070),ICORE(I090),LENABC,1.0D+00)
      ENDIF
C
#ifdef _DEBUG_NONHFTRP
      CALL RMD314(ICORE(I070),ICORE(I060),LENABC,DIJK)
      CALL DAXPY(LENABC, 1.0D0, ICORE(I000), 1, ICORE(I070),1)
      CALL E4TCTDIJK(ICORE(I000), ICORE(I070), ICORE(I060),
     &               LENABC,DIJK, E5TEST,__FILE__)
      WRITE(6,*) ' @TRPS2-I, E5ST from disconnected triples ',E5TEST
      Write(*,*)
      Call checksum("T3(C)+T3D-(I)", ICORE(I070), Lenabc)
#endif
      ENDIF
      IF (.NOT. FINAL) 
     &  CALL RMD314(ICORE(I070),ICORE(I060),LENABC,DIJK)
      ENDIF
C
C
C     ICORE(I000)     A<B,C T3(C)
C     ICORE(IL00)     A<B,C L3(C)
C     ICORE(I070)     A<B,C T3(D)/L3(D)
C
C     If this is a QCISD(T) calculation, multiply T3(C) by 2.
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC,2.0D+00,ICORE(I000),1)
      ENDIF
C
C     Expand alpha*t3(c).
C
CAP      IF(INT1.OR.INT2.OR.IDRLVL.GT.0)THEN

      IF(INT1.OR.INT2.OR.IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL SYMTRW2(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1             IADBLK,IADW2,IRPIJK,2,1)
      IF (LCCSD_T)
C     Expand alpha*l3(c)
     &   CALL SYMTRW2(ICORE(IL00),ICORE(IL10),ICORE(ISTART),
     &                IADBLK,IADW2,IRPIJK,2,1)
      ENDIF
C
      IF(INT1)THEN
      CALL  T1T32N(D1T1B,D1T1A,S1B,S1A,
     1            ICORE(I000),ICORE(I010),
     1            ICORE(I020),ICORE(I040),ICORE(I050),
     1            SCR1,SCR2,SCR3,
     1            IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     1            IRPIJK,I,J,K,2,1,
     1            IRPDPD(IRPIJ,2),IRPDPD(IRPJK,13),IRPDPD(IRPIK,13),
     1            NONHF,LCCSD_T)
#ifdef _DEBUG_LVLM
      call checksum("D1T1A:   ",D1T1A,NTAA)
      call checksum("D1T1B:   ",D1T1B,NTBB)
#endif
      ENDIF
C
C     Form T3(C) contribution to GAMMA(IJ,AB). For external CCSD
C     this contribution is included in FORM_HT3C3 (see below).
C
      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
      CALL T2FT33(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            SCR1,SCR2,SCR3,ICLLVL,1,LCCSD_T)
      IF(NONHF)THEN
      CALL T2FT33(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            SCR1,SCR2,SCR3,ICLLVL,2,LCCSD_T)
      ENDIF
      ENDIF
C
C     If this is a QCISD(T) calculation, divide 2*T3(C) by 2.
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC   ,0.5D+00,ICORE(I000),1)
      CALL SSCAL(I020-I010,0.5D+00,ICORE(I010),1)
      ENDIF
C
C     Form alpha*t3(d)/l3(d) and expand it.
C
CAP      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN

      IF((IDRLVL.GT.0.AND.
     1  (ICLLVL.EQ.21.OR.ICLLVL.EQ.22).OR.EXTRNL_CCSD))THEN

      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC,2.0D+00,ICORE(I070),1)
      ENDIF
      CALL SYMTRW2(ICORE(I070),ICORE(I080),ICORE(ISTART),
     1             IADBLK,IADW2,IRPIJK,2,1)
      ENDIF
C
C     We now have :
C                     I000 t3(c)
C                     I010 expanded t3(c)
C                     IL00 l3(c)
C                     IL10 expanded l3(c)
C                     I070 alpha*t3(d)/l3(d)
C                     I080 expanded alpha*t3(d)/l3(d)
C
CAP      IF(IDRLVL.GT.0)THEN

      IF (IDRLVL.GT.0 .OR. EXTRNL_CCSD) THEN
C
C     GRADIENTS : DENSITY EVALUATION (DIAGONAL ONLY).
C
C     Put t3(c)+alpha*t3(d) or l3(c)+alpha*l3(d) at I070 and I080
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
         IF (LCCSD_T) THEN
            CALL VADD(ICORE(I070),ICORE(IL00),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I080),ICORE(IL10),ICORE(I080),I020-I010
     &         ,1.0D+00)
         ELSE
            IF (.NOT. FINAL) THEN
            CALL VADD(ICORE(I070),ICORE(I000),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I080),ICORE(I010),ICORE(I080),I020-I010
     &         ,1.0D+00)
            ENDIF
      ENDIF
      ENDIF
C
      IF(ICLLVL.EQ.4.OR.ICLLVL.EQ.9.OR.ICLLVL.EQ.11)THEN
      CALL DOT33(DOOA,DOOB,ICORE(I000),ICORE(I000),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2))
C
      CALL DVT33(DVVA,DVVB,ICORE(I000),ICORE(I010),
     1                     ICORE(I000),ICORE(I010),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2))
      ENDIF
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
      CALL DOT33(DOOA,DOOB,ICORE(I000),ICORE(I070),
     1           LENABC,IRPI,IRPJ,IRPK,
     1           I,J,K,NOCCO(1),NOCCO(2))
C
#ifdef _DEBUG_LVLM
      call checksum("DOOA", DOOA, NOCCO(1)*NOCCO(1))
#endif
      CALL DVT33(DVVA,DVVB,ICORE(I000),ICORE(I010),
     1                     ICORE(I070),ICORE(I080),IADBLK,IADW2,
     1           IRPIJK,NVRTO(1),NVRTO(2))
#ifdef _DEBUG_LVLM
      call checksum("DVVA", DVVA, NVRTO(1)*NVRTO(1))
#endif
      ENDIF
C
C     GRADIENTS : GAMMA EVALUATION.
C
C     Put t3(c)+(alpha/2)*t3(d) at I000 and I010.
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
         IF (LCCSD_T) THEN
C     In the LCCSD(T) case, put l3(c)+alpha*l3(d) at IL00 and IL10
             call dcopy(lenabc,icore(i070),1,icore(il00),1)
             call dcopy(i020-i010,icore(i080),1,icore(il10),1)
         ELSE
            IF (.NOT. FINAL) THEN
            CALL VADD(ICORE(I000),ICORE(I000),ICORE(I070),LENABC
     &         ,1.0D+00)
            CALL VADD(ICORE(I010),ICORE(I010),ICORE(I080),I020-I010
     &         ,1.0D+00)
      CALL SSCAL(LENABC   ,0.5D+00,ICORE(I000),1)
      CALL SSCAL(I020-I010,0.5D+00,ICORE(I010),1)
      ENDIF

      ENDIF
      ENDIF
C
      IF (LCCSD_T) THEN
C     Add contributions from l3(c)+alpha*l3(d)
C Scale Lambda3 and T3(c) because there is an implicit factor of two
C for gamma
         call dscal(lenabc,0.5d0,icore(i000),1)
         call dscal(lenabc,0.5d0,icore(il00),1)
         call dscal(i020-i010,0.5d0,icore(i010),1)
         call dscal(i020-i010,0.5d0,icore(il10),1)
         CALL GVT33(ICORE(IL00),ICORE(IL10),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ICORE(I030),
     &      ICORE(I040),ICORE(I050),DSZEXP,IRPDPD(IRPJK,13),
     &      IRPDPD(IRPIK,13))
         CALL GOT33(GOOOV,LNOOOV,OOOVAD,ICORE(IL00),ICORE(IL10),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &      I,J,K,.FALSE.)
C     Load Lambda2 elements into memory
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,IRPIJ,
     &               145)
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,IRPJK,
     &               146)
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1,IRPIK,
     &               146)
C     Add contributions with t3(c) to the gamma matrices
         CALL GVT33(ICORE(I000),ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ICORE(I030),
     &      ICORE(I040),ICORE(I050),DSZEXP,IRPDPD(IRPJK,13),
     &      IRPDPD(IRPIK,13))
C Reload the T2 vectors back in to the appropriate locations, so any methods
C other than LCCSD(T) can  smoothly proceed.
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,IRPIJ,
     &               45)
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,ICORE(I030),
     1                ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,IRPJK,
     &               46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1,IRPIK,
     &               46)
         CALL GOT33(GOOOV,LNOOOV,OOOVAD,ICORE(I000),ICORE(I010),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &      I,J,K,.TRUE.)
C Scale Lambda3 and T3(c) because there is an implicit factor of two
C for gamma
         call dscal(lenabc,2.0d0,icore(i000),1)
         call dscal(lenabc,2.0d0,icore(il00),1)
         call dscal(i020-i010,2.0d0,icore(i010),1)
         call dscal(i020-i010,2.0d0,icore(il10),1)

      ELSE
         IF (.NOT. EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVLM
         call checksum("T33(C)---",icore(i000),lenabc)
         call checksum("T33(d)---",icore(i010),iofft2-i010)
#endif
         CALL GVT33(ICORE(I000),ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ICORE(I030),
     &      ICORE(I040),ICORE(I050),DSZEXP,IRPDPD(IRPJK,13),
     &      IRPDPD(IRPIK,13))
         CALL GOT33(GOOOV,LNOOOV,OOOVAD,ICORE(I000),ICORE(I010),
     &      ICORE(ISTART),IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &      I,J,K,.FALSE.)
         ENDIF 
      ENDIF
C
      ENDIF

      IF (EXTRNL_CCSD) THEN
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,
     &               IRPIJ,15)
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,
     &                ICORE(I030),ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,
     &               IRPJK,16)
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1, 
     &               IRPIK,16)
         CALL SYMTRW2(ICORE(I000),ICORE(I010),ICORE(ISTART),
     &                IADBLK,IADW2,IRPIJK,2,1)
#ifdef _DEBUG_LVLM
         call checksum("T33(C)---",icore(i000),lenabc)
         call checksum("T33(d)---",icore(i010),iofft2-i010)
#endif
         CALL T32ABCI3(ICORE(I000),ICORE(I010),ICORE(ISTART),
     &                 IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,
     &                 IRPJK,IRPIJK,I,J,K,ICORE(I030),
     &                 ICORE(I040),ICORE(I050),DSZEXP,
     &                 IRPDPD(IRPJK,13),IRPDPD(IRPIK,13))
         CALL GETLST(ICORE(I020),IOFFOO(IRPJ,IRPIJ,2)+1,NIJ,1,
     &               IRPIJ,45)
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN1),DSZEXP,DSZ,NIJ,
     &                ICORE(I030),ICORE(I020))
         CALL GETLST(ICORE(I040),IOFFOO(IRPJ,IRPJK,5)+1,NJK,1,
     &               IRPJK,46)
         CALL GETLST(ICORE(I050),IOFFOO(IRPI,IRPIK,5)+1,NIK,1, 
     &               IRPIK,46)

         CALL T32IJKA3(GOOOV,LNOOOV,OOOVAD,ICORE(I000),
     &                ICORE(I010),ICORE(ISTART),IADBLK,IADW2,
     &                IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,
     &               .FALSE.)
      ENDIF 
C
      IF(INT2)THEN
C
C     --- D2T2 = F T3 for non-Hartree-Fock cases ---
C
      IF((ICLLVL.EQ.13.AND.NONHF).OR. 
     1   (ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL T2FT33(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            SCR1,SCR2,SCR3,ICLLVL,0,.FALSE.)
      ENDIF
      IF (EXTRNL_CCSD .AND. NONHF) THEN
      CALL T2FT33(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1            IADBLK,IADW2,NLEFT,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1            ICLLVL,IUHF,0,.FALSE.)
      ENDIF 
C
C     --- D2T2 = W T3 ---
      CALL T2T33(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1           IADBLK,IADW2,IRPI,IRPJ,IRPK,
     1           IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,.FALSE.)
C
      CALL T2T33O(ICORE(I000),ICORE(I010),ICORE(ISTART),
     1 OOOV2(OOOVAD(IRPIJ,2)+IOFFOO(IRPJ,IRPIJ,2)*IRPDPD(IRPIJ,17),2),
     1 OOOV2(OOOVAD(IRPIK,3) + IOFFOO(IRPI,IRPIK,5)*LNVOIK,3),
     1 OOOV2(OOOVAD(IRPJK,3) + IOFFOO(IRPJ,IRPJK,5)*LNVOJK,3),
     1 OOOV2(OOOVAD(IRPIK,4) + IOFFOO(IRPI,IRPIK,5)*LNVOIK2,4),
     1 OOOV2(OOOVAD(IRPJK,4) + IOFFOO(IRPJ,IRPJK,5)*LNVOJK2,4),
     1 IADBLK,IADW2,
     1 IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1 LNVOIJ,LNVOIK,LNVOJK,LNVOIK2,LNVOJK2,.FALSE.)
C
      IF (LCCSD_T) THEN
C     --- D2L2 = L3 W ---
         CALL T2T33(ICORE(IL00),ICORE(IL10),ICORE(ISTART),IADBLK,IADW2,
     &      IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,.TRUE.)
C
         CALL T2T33O(ICORE(IL00),ICORE(IL10),ICORE(ISTART),OOOV2(
     &      OOOVAD(IRPIJ,2)+IOFFOO(IRPJ,IRPIJ,2)*IRPDPD(IRPIJ,17),2),
     &      OOOV2(OOOVAD(IRPIK,3) + IOFFOO(IRPI,IRPIK,5)*LNVOIK,3),
     &      OOOV2(OOOVAD(IRPJK,3) + IOFFOO(IRPJ,IRPJK,5)*LNVOJK,3),
     &      OOOV2(OOOVAD(IRPIK,4) + IOFFOO(IRPI,IRPIK,5)*LNVOIK2,4),
     &      OOOV2(OOOVAD(IRPJK,4) + IOFFOO(IRPJ,IRPJK,5)*LNVOJK2,4),
     &      IADBLK,IADW2,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     &      LNVOIJ,LNVOIK,LNVOJK,LNVOIK2,LNVOJK2,.TRUE.)
      ENDIF

      ENDIF
C This is no longer needed since T3 -> T2 are computed by default since
C INT1 and INT2 is true for extrnl_ccsd
CSSS      IF (EXTRNL_CCSD) THEN
CSSS         CALL FORM_HT3C3(ICORE(I000),ICORE(I010),SCR1,SCR2,SCR3,
CSSS     &                   ICORE(ISTART),OOOV,OOOV2,OOOVAD,LNOOOV,
CSSS     &                   IADBLK,IADW,IADW2,IRPI,IRPJ,IRPK,IRPIJ,
CSSS     &                   IRPIK,IRPJK,IRPIJK,I,J,K,LNVOIJ,LNVOIK,
CSSS     &                   LNVOJK,LNVOIK2,LNVOJK2,ICLLVL,.FALSE.,
CSSS     &                   NONHF)
CSSS      ENDIF 
C
c      IJK = IJKPOS(IRPI,IRPJ,IRPK,2)
c      IOFF = IJKOFF(IJK,IRPIJK,ISPIN1+1)
c      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ + INDEX(J-1) + I
c      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN1) + I
C     write(6,*) ' @TRPS3-I, sym, record ',IRPIJK,IJKVAL
      IF(T3STOR)THEN
      CALL PUTLIST(ICORE(I000),IJKVAL,1,1,IRPIJK,ISPIN1 + 1)
      ENDIF
#ifdef _DEBUG_LVLM
      write(*,*)
      call checksum("TRPS3", ICORE(I000), LENABC)
#endif
  647 CONTINUE
  648 CONTINUE
  649 CONTINUE
C
  980 CONTINUE
  990 CONTINUE
C
 1000 CONTINUE
      IF (.NOT. FINAL) THEN
      WRITE(6,3000) E4TBBA
 3000 FORMAT(' @TRPS3-I, E4TBBA ',F20.12)

      ELSE 

      WRITE(6,3001) E4TBBA
 3001 FORMAT(' @TRPS3-I, E4TBBA+E5ST+E4ST ',F20.12)

      ENDIF 

      IF (EXTRNL_CCSD) THEN
         CALL GETLST(SCR1,1,1,1,3,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1A:   ",SCR1,NTAA)
#endif 
         CALL VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
         CALL PUTLST(SCR1,1,1,1,3,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1A:   ",D1T1A,NTAA)
         call checksum("D1T1A:   ",SCR1,NTAA)
#endif
         CALL GETLST(SCR1,1,1,1,4,90)
#ifdef _DEBUG_LVLM
         call checksum("D1T1B:   ",SCR1,NTBB)
#endif 
         CALL VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
#ifdef _DEBUG_LVLM
         call checksum("D1T1B:   ",D1T1B,NTBB)
         call checksum("D1T1B:   ",SCR1,NTBB)
#endif 
         CALL PUTLST(SCR1,1,1,1,4,90)
         DO IRREP=1,NIRREP
            CALL PUTTRN(GOOOV(OOOVAD(IRREP,2),2),BUF,
     &                  IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     &                  2,IRREP,108)
            CALL SYMTR3(IRREP,POP(1,2),VRT(1,1),IRPDPD(IRREP,11),
     &                  IRPDPD(IRREP,14),GOOOV(OOOVAD(IRREP,3),3),
     &                  SCR1,SCR2,SCR3)
            CALL PUTTRN(GOOOV(OOOVAD(IRREP,3),3),BUF,
     &                  IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     &                  2,IRREP,109)
            CALL PUTTRN(GOOOV(OOOVAD(IRREP,4),4),BUF,
     &                  IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     &                  2,IRREP,110)
        ENDDO
#ifdef _DEBUG_LVL0
        Write(6,"(a)") " ---After trps3,Lists:61-63,90(3,4)---"
CSSS        call check_ints(icore(istart),nleft,iuhf,.false.)
        call check_t32w(icore(istart),nleft,iuhf)
        call check_t2(icore(istart),nleft,iuhf)
C        irrepx=1
C        Lenab=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
C        iend = istart+lenab
C        nleft = nleft - iend
C        call Form_t2aaaa(icore(istart),icore(iend),nleft,lenab)
#endif
      ENDIF

      IF(IDRLVL.GT.0)THEN
C
      DO  9463 IRREP=1,NIRREP
      CALL PUTTRN(GOOOV(OOOVAD(IRREP,2),2),BUF,
     1            IRPDPD(IRREP, 4),IRPDPD(IRREP,17),
     1            2,IRREP,LISGO2)
      CALL SYMTR3(IRREP,POP(1,2),VRT(1,1),IRPDPD(IRREP,11),
     1            IRPDPD(IRREP,14),GOOOV(OOOVAD(IRREP,3),3),
     1            SCR1,SCR2,SCR3)
      CALL PUTTRN(GOOOV(OOOVAD(IRREP,3),3),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,11),
     1            2,IRREP,LISGO3)
      CALL PUTTRN(GOOOV(OOOVAD(IRREP,4),4),BUF,
     1            IRPDPD(IRREP,14),IRPDPD(IRREP,18),
     1            2,IRREP,LISGO4)
 9463 CONTINUE
C
      ENDIF
C
      IF((ICLLVL.GE.13.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
      CALL GETLST(SCR1,1,1,1,3,90)
#ifdef _DEBUG_LVL0
      call checksum("D1T1A:   ",SCR1,NTAA)
#endif 
      CALL   VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
#ifdef _DEBUG_LVL0
      call checksum("D1T1A:   ",D1T1A,NTAA)
      call checksum("D1T1A:   ",SCR1,NTAA)
#endif 
      CALL PUTLST(SCR1,1,1,1,3,90)
      CALL GETLST(SCR1,1,1,1,4,90)
#ifdef _DEBUG_LVL0
      call checksum("D1T1B:   ",SCR1,NTBB)
#endif 
      CALL   VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
#ifdef _DEBUG_LVL0
      call checksum("D1T1B:   ",D1T1B,NTBB)
      call checksum("D1T1B:   ",SCR1,NTBB)
#endif 
      CALL PUTLST(SCR1,1,1,1,4,90)
      ENDIF
C
C     If this is CCSD(T) or QCISD(T) gradient, write D1T1 triples contribution
C     to disk.
C
      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
C List 104 has been used eleswhere.
CSSS        CALL GETLST(SCR1,1,1,1,1,104)
        CALL GETLST(SCR1,1,1,1,1,4)
        CALL   VADD(SCR1,SCR1,D1T1A,NTAA,1.0D+00)
        CALL PUTLST(SCR1,1,1,1,1,4)
        CALL GETLST(SCR1,1,1,1,2,4)
        CALL   VADD(SCR1,SCR1,D1T1B,NTBB,1.0D+00)
        CALL PUTLST(SCR1,1,1,1,2,4)
        IF(NONHF)THEN
          CALL GETLST(SCR1,1,1,1,5,160)
          CALL   VADD(SCR1,SCR1,S1A,NTAA,0.5D+00)
          CALL PUTLST(SCR1,1,1,1,5,160)
          CALL GETLST(SCR1,1,1,1,6,160)
          CALL   VADD(SCR1,SCR1,S1B,NTBB,0.5D+00)
          CALL PUTLST(SCR1,1,1,1,6,160)
        ENDIF
      ENDIF
C
      IF(ICLLVL.EQ.9.OR.ICLLVL.EQ.13.OR.ICLLVL.EQ.14.OR.ICLLVL.EQ.15.OR.
     1   ICLLVL.EQ.16.OR.ICLLVL.EQ.17.OR.ICLLVL.EQ.18.OR.
     1   ICLLVL.EQ.33.OR.ICLLVL.EQ.34) RETURN
      IF(ICLLVL.EQ.4.AND.IREFNC.LE.1.AND..NOT.NONHF) RETURN
C
C     COMPUTE CCSD(T) OR QCISD(T) FIFTH-ORDER PIECE OF ENERGY
C
C For L-based CCSD instead of T1, we use L1 to compute the
C 5th order triples contribution. Ajith Perera, 10/05.
C
      IF (.NOT. FINAL) THEN

      IF (LCCSD_T) THEN
          CALL GETLST(T1A,1,1,1,1,190)
          CALL GETLST(T1B,1,1,1,2,190)
      ELSE
          CALL GETLST(T1A,1,1,1,1,90)
          CALL GETLST(T1B,1,1,1,2,90)
      ENDIF 
C
      E5ST2A = 0.0D+00
      E5ST2B = 0.0D+00
      DO  3040 AI=1,NTAA
      E5ST2A = E5ST2A + T1A(AI) * D1T1A(AI)
 3040 CONTINUE
C
      DO  3050 AI=1,NTBB
      E5ST2B = E5ST2B + T1B(AI) * D1T1B(AI)
 3050 CONTINUE
C
      E5STA = E5STA + E5ST2A
      E5STB = E5STB + E5ST2B
      WRITE(6,6020) E5STA,E5STB
 6020 FORMAT(' @TRPS3-I, E5ST A ',F20.15,/,
     1       ' @TRPS3-I, E5ST B ',F20.15)
C
C     Return if reference is Hartree-Fock.
C

#ifdef _DEBUG_LVL0
      if (.not. extrnl_ccsd) then
      Write(6,"(a)") " ---After trps3,Lists:61-63,90(3,4)---"
CSSS      call check_ints(icore(istart),nleft,iuhf,.false.) 
CSSS      call check_t32w_grad(icore(istart),nleft,iuhf)
CSSS      call check_t32w(icore(istart),nleft,iuhf)
      call check_t2(icore(istart),nleft,iuhf)
      call check_t2grd(icore(istart),nleft,iuhf)
C      irrepx=1
C      Lenab=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
C      iend = istart+lenab
C      nleft = nleft - iend
C      call Form_t2aaaa(icore(istart),icore(iend),nleft,lenab)
      endiF
#endif
      IF(.NOT.NONHF) RETURN
C
C     COMPUTE NON-HF CCSD(T) OR QCISD(T) FOURTH-ORDER PIECE.
C

      CALL GETLST(FOVA,1,1,1,3,93)
      CALL GETLST(FOVB,1,1,1,4,93)
C
      E4ST2A = 0.0D+00
      E4ST2B = 0.0D+00
      DO  4040 AI=1,NTAA
      E4ST2A = E4ST2A + FOVA(AI) * S1A(AI)
 4040 CONTINUE
C
      DO  4050 AI=1,NTBB
      E4ST2B = E4ST2B + FOVB(AI) * S1B(AI)
 4050 CONTINUE
C
      E4STA = E4STA + E4ST2A
      E4STB = E4STB + E4ST2B
      WRITE(6,7020) E4STA,E4STB
 7020 FORMAT(' @TRPS3-I, E4ST A ',F20.15,/,
     1       ' @TRPS3-I, E4ST B ',F20.15)
      ENDIF 

      RETURN
 9010 FORMAT(' TRPS3-I, Insufficient memory to continue. ')
      
      END
