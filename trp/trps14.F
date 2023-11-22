      SUBROUTINE TRPS14(OOOV,OOOV2,GOOOV,BUF,EVAL,
     1                  IADBLK,LENBLK,IADT2,LENT2,IADV,LENV,
     1                  OOOVAD,LENINT,LNOOOV,
     1                  ICORE,NFREE,IUHF,ISPIN,NO,NV,INT1,INT2,NONHF,
     1                  LCCSD_T,FINAL,EXTRNL_CCSD)
      IMPLICIT INTEGER (A-Z)
#include "maxbasfn.par"
C
C The three most relevant papers for all of the standard triple
C excitation methods: Watts et al., JCP, 98, 8718, 1993 (JCP2)
C for L-CCSD(T), Crawford and Stanton, IJQC,70,601,Eqns. A16 (IJQC)
C
      LOGICAL INT1,INT2,NONHF,FINAL
      LOGICAL IJKEQL,IJEQL,JKEQL,NONEQL
      LOGICAL CCSDT4,CCSDT,LCCSD_T
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      LOGICAL EXTRNL_CCSD
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION DIJK,DOOA,DOOB,DVVA,DVVB
      DOUBLE PRECISION EVAL(NO+NV), DDOT
      DOUBLE PRECISION E4TAAA,E4TBBB,E4TBAA,E4TABB,E5STA,E5STB,
     1                 E4STA ,E4STB ,E4T14,E5ST14,E4ST14
      DOUBLE PRECISION E5TEST,SDOT,E6TEST
C     ijka integrals/intermediates used in D3T3 = WT2
      DOUBLE PRECISION OOOV(LNOOOV),BUF(NO*(NO-1)/2)
C     ijka integrals/intermediates used in D2T2 = WT3 (includes 1b)
      DOUBLE PRECISION OOOV2(LNOOOV)
      DIMENSION OOOVAD(8),LENINT(8)
C
      DOUBLE PRECISION GOOOV(LNOOOV)
      DOUBLE PRECISION ICORE(1)
C
      DIMENSION LEN(8,8),IADBLK(8),LENBLK(8)
      DIMENSION                        IADT2(8) ,LENT2(8),
     1                                 IADV(8)  ,LENV(8)
      DIMENSION POPT3(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
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
C     TRIPLES COMMON BLOCKS
C
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXBASFN,8,2),ABSOCC(MAXBASFN,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TBAA,E4TABB,E5STA,E5STB,
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
      IF(ISPIN.EQ.1)THEN
      NT1 = NTAA
      ELSE
      NT1 = NTBB
      ENDIF
#ifdef _DEBUG_LVL0
      iprint=1
#endif 
C
C     MAIN CORE VECTOR ADDRESSING
C
      IM90 = 1
      IM80 = IM90 + NT1
      IM70 = IM80 + NT1
      IM60 = IM70 + NT1
      IM50 = IM60 + NT1
C
      IF(ISPIN.EQ.1)THEN
      WRITE(6,1015)
 1015 FORMAT(' @TRPS14-I, Spin case AAA ')
      ELSE
      WRITE(6,1020)
 1020 FORMAT(' @TRPS14-I, Spin case BBB ')
      ENDIF
C
C     INITIALIZE T1 INCREMENTS
C
      CALL ZERO(ICORE(IM90),NT1)
C
C     INITIALIZE "F" INCREMENTS
C
      CALL ZERO(ICORE(IM70),NT1)
C
C     GET EIGENVALUES
C
      IF(ISPIN.EQ.1)THEN
      CALL GETREC(20,'JOBARC','SCFEVALA',(NO+NV)*IINTFP,EVAL)
      ELSE
      CALL GETREC(20,'JOBARC','SCFEVALB',(NO+NV)*IINTFP,EVAL)
      ENDIF
C
C     Compute (B<C,A) lengths for given IRPIJK (=IRPABC)
C
      DO   90 IRPIJK=1,NIRREP
      DO   80 IRPA  =1,NIRREP
      IRPBC = DIRPRD(IRPIJK,IRPA)
      LEN(IRPA,IRPIJK) = IRPDPD(IRPBC,ISPIN) * VRT(IRPA,ISPIN)
   80 CONTINUE
   90 CONTINUE
C
C     COMPUTE LENGTHS OF SYMMETRY PACKED ABC TRIPLES FOR EACH
C     IRREP
C
      DO  220 IRREP=1,NIRREP
      POPT3(IRREP) = 0
  220 CONTINUE
C
      DO 233 IRPC=1,NIRREP
      DO 232 IRPB=1,IRPC
      DO 231 IRPA=1,IRPB
      IRPBC  = DIRPRD(IRPB,IRPC)
      IRPABC = DIRPRD(IRPA,IRPBC)
C
C     Compute length of the block. There are four possible cases.
C
      IF(IRPA.EQ.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) * (VRT(IRPC,ISPIN) - 1) *
     1                            (VRT(IRPC,ISPIN) - 2)) / 6
      ENDIF
C
      IF(IRPA.EQ.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) *  VRT(IRPB,ISPIN)      *
     1                            (VRT(IRPB,ISPIN) - 1)) / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.EQ.IRPC)THEN
      LENGTH = (VRT(IRPC,ISPIN) * (VRT(IRPC,ISPIN) - 1) *
     1                             VRT(IRPA,ISPIN)    )  / 2
      ENDIF
C
      IF(IRPA.NE.IRPB.AND.IRPB.NE.IRPC)THEN
      LENGTH =  VRT(IRPC,ISPIN) *  VRT(IRPB,ISPIN)      *
     1                             VRT(IRPA,ISPIN)
      ENDIF
C
      POPT3(IRPABC) = POPT3(IRPABC) + LENGTH
  231 CONTINUE
  232 CONTINUE
  233 CONTINUE
C

c      DO  230 C=3,NV
c      DO  230 B=2,C-1
c      DO  230 A=1,B-1
c      IRPA = IRREPS(NO + A,ISPIN)
c      IRPB = IRREPS(NO + B,ISPIN)
c      IRPC = IRREPS(NO + C,ISPIN)
c      IRPBC = DIRPRD(IRPB,IRPC)
c      IRPABC = DIRPRD(IRPBC,IRPA)
c      POPT3(IRPABC) = POPT3(IRPABC) + 1
c  230 CONTINUE
C
      DO  240 IRREP=1,NIRREP
C     WRITE(6,1200) IRREP,POPT3(IRREP)
  240 CONTINUE
 1200 FORMAT(' @TRPS14-I, IRREP, A<B<C POPULATION ',I4,I10)
C
C     COMPUTE ADDRESSES FOR BLOCKS OF OOOV INTEGRALS. PRESENTLY
C     THESE ARE HELD IN CORE.
C
      DO  250 IRREP=1,NIRREP
      LENINT(IRREP) = IRPDPD(IRREP,2+ISPIN) * IRPDPD(IRREP,15+ISPIN)
  250 CONTINUE
C
      OOOVAD(1) = 1
      IF(NIRREP.GT.1)THEN
      DO  260 IRREP = 2,NIRREP
      OOOVAD(IRREP) = OOOVAD(IRREP-1) + LENINT(IRREP-1)
  260 CONTINUE
      ENDIF
C      WRITE(6,1210)
 1210 FORMAT(' @TRPS14-I, OOOV ARRAY LENGTHS AND ADDRESSES ')
      DO  270 IRREP=1,NIRREP
C      WRITE(6,1220) IRREP,LENINT(IRREP),OOOVAD(IRREP)
 1220 FORMAT(' IRREP, LENGTH, ADDRESS ',I4,I6,I6)
  270 CONTINUE
C
      IF(ISPIN.EQ.1)THEN
      LISTWO = LWIC11
      ELSE
      LISTWO = LWIC12
      ENDIF
C
      DO  280 IRREP=1,NIRREP
      CALL GETTRN(OOOV(OOOVAD(IRREP)),BUF,IRPDPD(IRREP,2+ISPIN),
     1            IRPDPD(IRREP,15+ISPIN),2,IRREP,LISTWO)
  280 CONTINUE
C
C     Read integrals/intermediates for CCSDT-1B piece
C
      IF(INT2)THEN
C
      IF(ISPIN.EQ.1)THEN
      DO  281 IRREP=1,NIRREP
      CALL GETTRN(OOOV2(OOOVAD(IRREP)),BUF,IRPDPD(IRREP,2+ISPIN),
     1            IRPDPD(IRREP,15+ISPIN),2,IRREP,LWIC21)
  281 CONTINUE
      ELSE
      DO  282 IRREP=1,NIRREP
      CALL GETTRN(OOOV2(OOOVAD(IRREP)),BUF,IRPDPD(IRREP,2+ISPIN),
     1            IRPDPD(IRREP,15+ISPIN),2,IRREP,LWIC22)
  282 CONTINUE
      ENDIF
C
      ENDIF
C
C     Initialize OOOV gamma if this is a gradient calculation.
C
      IF((IDRLVL.GT.0 .OR. EXTRNL_CCSD) .AND.IUHF.GT.0)THEN
      CALL ZERO(GOOOV,LNOOOV)
      ENDIF
C
C     COMPUTE THE NUMBER OF FREE DOUBLE WORDS (YUK !)
      MCORE = NFREE / IINTFP
cYAU      WRITE(6,1230) NFREE,MCORE,IINTFP
 1230 FORMAT(' @TRPS14-I, Number of free integer words ',I10,/,
     1       '            Number of free double words  ',I10,/,
     1       '            Ratio                        ',I10)
C
      E4T14  = 0.0D+00
      E5ST14 = 0.0D+00
C
      DO  1000 IRPIJK=1,NIRREP
C
C The orginal IJKVAL =  0 is wrong. This variable never get updated
C before it is used. So, set it to 1 instead. Ajith Perera, 09/2006.
C 
C This is incorrect.  Need to be updating ijkval as it cycles through
C IJKVAL should start at 0 and then be added to after D1T14
C Andrew Taube 12/2009

C     IJKVAL = 1
      IJKVAL = 0
C
      DO   990   IRPK=1,NIRREP
      IF(POP(IRPK,ISPIN).EQ.0) GOTO 990
      DO   980   IRPJ=1,IRPK
      IRPJK = DIRPRD(IRPJ,IRPK)
      IRPI  = DIRPRD(IRPJK,IRPIJK)
      IRPIJ =  DIRPRD(IRPI,IRPJ)
      IRPIK =  DIRPRD(IRPI,IRPK)
C
      IF(IRPI.GT.IRPJ.OR.POP(IRPI,ISPIN).EQ.0.OR.
     1                   POP(IRPJ,ISPIN).EQ.0) GOTO 980
C
      IF(IRPI.EQ.IRPJ.AND.
     1   IRPI.EQ.IRPK.AND.POP(IRPK,ISPIN).LT.3) GOTO 980
      IF(IRPJ.EQ.IRPK.AND.POP(IRPK,ISPIN).LT.2) GOTO 980
      IF(IRPI.EQ.IRPJ.AND.POP(IRPJ,ISPIN).LT.2) GOTO 980
C
      IJKEQL = .FALSE.
      IJEQL  = .FALSE.
      JKEQL  = .FALSE.
      NONEQL = .FALSE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.EQ.IRPK) IJKEQL = .TRUE.
      IF(IRPJ.EQ.IRPK.AND.IRPI.NE.IRPK)  JKEQL = .TRUE.
      IF(IRPI.EQ.IRPJ.AND.IRPI.NE.IRPK)  IJEQL = .TRUE.
      IF(IRPI.NE.IRPJ.AND.IRPI.NE.IRPK.AND.IRPJ.NE.IRPK) NONEQL = .TRUE.
C
      IF(IJKEQL)THEN
      NIJ = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      NIK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      NJK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      ENDIF
C
      IF(JKEQL)THEN
      NIJ = POP(IRPI,ISPIN) * POP(IRPJ,ISPIN)
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = (POP(IRPK,ISPIN) * (POP(IRPK,ISPIN)-1))/2
      ENDIF
C
      IF(IJEQL)THEN
      NIJ = (POP(IRPJ,ISPIN) * (POP(IRPJ,ISPIN)-1))/2
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = POP(IRPJ,ISPIN) * POP(IRPK,ISPIN)
      ENDIF
C
      IF(NONEQL)THEN
      NIJ = POP(IRPI,ISPIN) * POP(IRPJ,ISPIN)
      NIK = POP(IRPI,ISPIN) * POP(IRPK,ISPIN)
      NJK = POP(IRPJ,ISPIN) * POP(IRPK,ISPIN)
      ENDIF
C
C     COMPUTE ADDRESSES AND LENGTHS OF ABC SYMMETRY BLOCKS FOR
C     THIS IJK. THIS LOOP IS EFFECTIVELY OVER IRREPS OF A (IE
C     COLUMNS OF THE EVENTUAL TARGET).
C
      DO  300 IRREP=1,NIRREP
      IF(IRREP.EQ.1)THEN
      IADBLK(IRREP) = 1
      ELSE
      IADBLK(IRREP) = IADBLK(IRREP-1) + LEN(IRREP-1,IRPIJK)
      ENDIF
      LENBLK(IRREP) = LEN(IRREP,IRPIJK)
  300 CONTINUE
C
      LENABC=0
      DO  310 IRREP=1,NIRREP
      LENABC = LENABC + LENBLK(IRREP)
  310 CONTINUE
C

C     IM50 is t3(a<b<c)
      I000 = IM50 + POPT3(IRPIJK)
      I010 = I000 + LENABC
C
C     SET ADDRESSES FOR T2 VECTORS
C
      I020 = I010 + IRPDPD(IRPIJ,ISPIN) * NIJ
      I030 = I020 + IRPDPD(IRPIK,ISPIN) * NIK
      I040 = I030 + IRPDPD(IRPJK,ISPIN) * NJK
C
C     DETERMINE DISTRIBUTION SIZES FOR EXPANDED T2 VECTORS
C
      DSZEIJ = IRPDPD(IRPIJ,18+ISPIN)
      DSZEIK = IRPDPD(IRPIK,18+ISPIN)
      DSZEJK = IRPDPD(IRPJK,18+ISPIN)
C
      I050 = I040 + DSZEIJ * NIJ
      I060 = I050 + DSZEIK * NIK
      I070 = I060 + DSZEJK * NJK
C
C     READ T2 VECTORS (NEED OFFSETS IN EACH IRREP_IOFFOO FROM OTHER
C     ROUTINES SHOULD SUFFICE).
C
      CALL GETLST(ICORE(I010),IOFFOO(IRPJ,IRPIJ,ISPIN)+1,
     1            NIJ,1,IRPIJ,43+ISPIN)
      CALL GETLST(ICORE(I020),IOFFOO(IRPK,IRPIK,ISPIN)+1,
     1            NIK,1,IRPIK,43+ISPIN)
      CALL GETLST(ICORE(I030),IOFFOO(IRPK,IRPJK,ISPIN)+1,
     1            NJK,1,IRPJK,43+ISPIN)
C
C     EXPAND T2 VECTORS
C
      CALL SYMEXP2(IRPIJ,VRT(1,ISPIN),DSZEIJ,IRPDPD(IRPIJ,ISPIN),
     1             NIJ,ICORE(I040),ICORE(I010))
      CALL SYMEXP2(IRPIK,VRT(1,ISPIN),DSZEIK,IRPDPD(IRPIK,ISPIN),
     1             NIK,ICORE(I050),ICORE(I020))
      CALL SYMEXP2(IRPJK,VRT(1,ISPIN),DSZEJK,IRPDPD(IRPJK,ISPIN),
     1             NJK,ICORE(I060),ICORE(I030))
C
C     SET ADDRESSES FOR TRANSPOSED BLOCKS OF OOOV INTEGRALS (USED FOR
C     INCLUSION OF T3 IN T2 (T2T314O).
C
      LNVOIJ = IRPDPD(IRPIJ,8+ISPIN)
      LNVOIK = IRPDPD(IRPIK,8+ISPIN)
      LNVOJK = IRPDPD(IRPJK,8+ISPIN)
      I080 = I070 + LNVOIJ * NIJ
      I090 = I080 + LNVOIK * NIK
      I100 = I090 + LNVOJK * NJK
C
C     DO THE OUT-OF-PLACE TRANSPOSITIONS.
C
      IF(INT2)THEN
C
      CALL SYMTRA2(IRPIJ,POP(1,ISPIN),VRT(1,ISPIN),LNVOIJ,NIJ,
     1             OOOV2(OOOVAD(IRPIJ) + 
     1                  IOFFOO(IRPJ,IRPIJ,ISPIN) * LNVOIJ),
     1             ICORE(I070))
      CALL SYMTRA2(IRPIK,POP(1,ISPIN),VRT(1,ISPIN),LNVOIK,NIK,
     1             OOOV2(OOOVAD(IRPIK) + 
     1                  IOFFOO(IRPK,IRPIK,ISPIN) * LNVOIK),
     1             ICORE(I080))
      CALL SYMTRA2(IRPJK,POP(1,ISPIN),VRT(1,ISPIN),LNVOJK,NJK,
     1             OOOV2(OOOVAD(IRPJK) + 
     1                  IOFFOO(IRPK,IRPJK,ISPIN) * LNVOJK),
     1             ICORE(I090))
      ENDIF
C
      I110 = I100 + POPT3(IRPIJK)
      I120 = I110 + IRPDPD(IRPJK,ISPIN)
      I130 = I120 + IRPDPD(IRPIK,ISPIN)
      I140 = I130 + IRPDPD(IRPIJ,ISPIN)
C
C     If this is a CCSD(T) or QCISD(T) derivative calculation, allocate
C     space for "disconnected" triples.
C
CSSS      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN

      I150 = I140 + POPT3(IRPIJK)
      ISTART = I150

CSSS      ELSE

CSSS      ISTART = I140

CSSS      ENDIF

      IF (LCCSD_T) THEN
C
C IL50 is L3(a<b<c)
C
          IL50 = ISTART
          IL00 = IL50 + POPT3(IRPIJK)
          ISTART=IL00+LENABC
      ENDIF
      NLEFT = MCORE - ISTART
C
C The following memory test is not correct. We use ICORE(ISTART) 
C to store quantities of length (a<b,c) for gradients and iterative
C triples. This might be the cause of mysterious crashes of large
C triples gradient calculations for low symmetry molecules.
C Ajith Perera 11/2005.
CSSS      IF(ISTART.GE.MCORE)THEN
CSSS      WRITE(6,9010)
CSSS      CALL INSMEM('TRPS14',ISTART,MCORE)
CSSS      ENDIF
      IF (ISTART+LENABC.GE.MCORE) THEN
         WRITE(6,9010)
         CALL INSMEM('TRPS14',ISTART+LENABC,MCORE)
      END IF
C
C     FORM ABC DENOMINATOR FOR THIS IJK SYMMETRY BLOCK.
C
      CALL MKD314(ICORE(I100),EVAL,ISPIN,NO,NV,IRPIJK,POPT3)
C
      IF(NONEQL)THEN
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 1
      JHIGH = POP(IRPJ,ISPIN)
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN)
      ENDIF
C
      IF(IJKEQL)THEN
      KLOW  = 3
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 2
      ILOW  = 1
      ENDIF
C
      IF(IJEQL)THEN
      KLOW  = 1
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 2
      JHIGH = POP(IRPJ,ISPIN)
      ILOW  = 1
      ENDIF
C
      IF(JKEQL)THEN
      KLOW  = 2
      KHIGH = POP(IRPK,ISPIN)
      JLOW  = 1
      ILOW  = 1
      IHIGH = POP(IRPI,ISPIN)
      ENDIF
C
      DO   430 K=KLOW,KHIGH
C
      IF(IJKEQL.OR.JKEQL) JHIGH = K-1
      DO   420 J=JLOW,JHIGH
C
      IF(IJKEQL.OR.IJEQL) IHIGH = J-1
      DO   410 I=ILOW,IHIGH
C
      CALL ZERO(ICORE(I000),LENABC)
      CALL ZERO(ICORE(IM50),POPT3(IRPIJK))
C
C Compute the T3(ijk,abc)= T2(ij,{le}a)<bc||{le}k> 
c
      IF (.NOT. FINAL) THEN

      CALL D1T14(ICORE(I000),ICORE(I040),ICORE(I050),ICORE(I060),
     1 OOOV(OOOVAD(IRPIJ) + 
     1      IOFFOO(IRPJ,IRPIJ,ISPIN)*IRPDPD(IRPIJ,15+ISPIN)),
     1 OOOV(OOOVAD(IRPIK) + 
     1      IOFFOO(IRPK,IRPIK,ISPIN)*IRPDPD(IRPIK,15+ISPIN)),
     1 OOOV(OOOVAD(IRPJK) + 
     1      IOFFOO(IRPK,IRPJK,ISPIN)*IRPDPD(IRPJK,15+ISPIN)),
     1           ICORE(ISTART),ISPIN,
     1           IADBLK,LENBLK,IADT2,LENT2,IADV,LENV,
     1           I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1           DSZEIJ,DSZEIK,DSZEJK,
     1 IRPDPD(IRPIJ,15+ISPIN),IRPDPD(IRPIK,15+ISPIN),
     1 IRPDPD(IRPJK,15+ISPIN),NLEFT,.FALSE.)
C Update IJKVAL (see comment at start of loop 1000)
      IJKVAL = IJKVAL + 1
C
C     CCSDT or non-HF calculations without semicanonical orbitals.
C     (i) Virtual orbital sum.
C
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT314V(ICORE(I000),ICORE(ISTART),NFREE,ISPIN,POPT3(IRPIJK),
     1              LENABC,IRPIJK,IJKVAL,IADBLK)
      ENDIF
C
      IF(CCSDT4.OR.CCSDT)THEN
C     CCSDT abef intermediate contribution.
      CALL T3WT314P(ICORE(I000),ICORE(ISTART),NFREE,ISPIN,POPT3(IRPIJK),
     1              LENABC,IRPIJK,IJKVAL,IADBLK)
CC     CCSDT maei intermediate contribution.
      CALL T3WT314R(ICORE(I000),ICORE(ISTART),NFREE,ISPIN,LENABC,
     1              IADBLK,IRPI,IRPJ,IRPK,I,J,K,IRPIJK)
      ENDIF
C
C     Convert the a<b;c T3 to a<b<c.
C
      CALL SYMCONTW(ICORE(IM50),ICORE(I000),IADBLK,ISPIN,IRPIJK)
C
C     CCSDT or non-HF calculations without semicanonical orbitals.
C     (ii) Occupied orbital sum.
C     Note : T3FT314O computes contributions to an a<b<c T3 rather than
C     to a<b;c. Therefore, it MUST be called after SYMCONTW.
C
      CALL ZERO(ICORE(ISTART),LENABC)
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL T3FT314O(ICORE(IM50),ICORE(ISTART),NFREE,ISPIN,
     1               POPT3(IRPIJK),IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
      IF(CCSDT4.OR.CCSDT)THEN
      CALL T3WT314HH(ICORE(IM50),ICORE(ISTART),NFREE,ISPIN,
     1               POPT3(IRPIJK),IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
C
      DIJK = EVAL(ABSOCC(I,IRPI,ISPIN)) + EVAL(ABSOCC(J,IRPJ,ISPIN))
     1                                  + EVAL(ABSOCC(K,IRPK,ISPIN))
      CALL RMD314(ICORE(IM50),ICORE(I100),POPT3(IRPIJK),DIJK)

      ELSE 
C
C If this iterative NONHF-CCSD(T), converged amplitudes are already
C on the disk. 

      IJKVAL = IJKVAL + 1
#ifdef _NOSKIP
      IF(ICLLVL.EQ.18.OR.(INONHF.EQ.1.AND.IORBTP.EQ.0))THEN
      CALL ZERO(ICORE(IM50), POPT3(IRPIJK))
      CALL T3FT314V(ICORE(IM50),ICORE(ISTART),NFREE,ISPIN,POPT3(IRPIJK),
     1              LENABC,IRPIJK,IJKVAL,IADBLK)
      CALL T3FT314O(ICORE(IM50),ICORE(ISTART),NFREE,ISPIN,
     1               POPT3(IRPIJK),IRPI,IRPJ,IRPK,IRPIJK,I,J,K)
      ENDIF
C     
      IF (ISTART+POPT3(IRPIJK) .GT. NLEFT) THEN
          CALL INSMEM("@-TRPS14", ISTART+POPT3(IRPIJK), NLEFT) 
      ENDIF
      CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)
      E4T14 = E4T14 + DDOT(POPT3(IRPIJK), ICORE(IM50), 1, 
     &                   ICORE(ISTART), 1)
#endif 
      CALL D1T14(ICORE(I000),ICORE(I040),ICORE(I050),ICORE(I060),
     1 OOOV(OOOVAD(IRPIJ) +
     1      IOFFOO(IRPJ,IRPIJ,ISPIN)*IRPDPD(IRPIJ,15+ISPIN)),
     1 OOOV(OOOVAD(IRPIK) +
     1      IOFFOO(IRPK,IRPIK,ISPIN)*IRPDPD(IRPIK,15+ISPIN)),
     1 OOOV(OOOVAD(IRPJK) +
     1      IOFFOO(IRPK,IRPJK,ISPIN)*IRPDPD(IRPJK,15+ISPIN)),
     1           ICORE(ISTART),ISPIN,
     1           IADBLK,LENBLK,IADT2,LENT2,IADV,LENV,
     1           I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1           DSZEIJ,DSZEIK,DSZEJK,
     1 IRPDPD(IRPIJ,15+ISPIN),IRPDPD(IRPIK,15+ISPIN),
     1 IRPDPD(IRPJK,15+ISPIN),NLEFT,.FALSE.)

      CALL SYMCONTW(ICORE(IM50),ICORE(I000),IADBLK,ISPIN,IRPIJK)
      CALL GETLIST(ICORE(ISTART),IJKVAL,1,1,IRPIJK,1 + 3*(ISPIN-1) + 4)

CSSS      Write(*,*)
cSSS      Call checksum("T3(C)+T3D-(D)", ICORE(ISTART), POPT3(IRREPIJK))

      E4T14 = E4T14 + DDOT(POPT3(IRPIJK), ICORE(IM50), 1, 
     &                   ICORE(ISTART), 1)

      E6TEST = E6TEST + DDOT(POPT3(IRPIJK), ICORE(IM50), 1,
     &                   ICORE(ISTART), 1)
CSSS      WRITE(6,*) ' @TRPS2-I, E6TEST from disconnected triples ',E6TEST

CSSS         CALL E4TIJK(ICORE(ISTART),ICORE(I100),POPT3(IRPIJK),DIJK,
CSSS     &               E4T14,__FILE__)
C
      DIJK = EVAL(ABSOCC(I,IRPI,ISPIN)) + EVAL(ABSOCC(J,IRPJ,ISPIN))
     1                                  + EVAL(ABSOCC(K,IRPK,ISPIN))
      CALL RMD314(ICORE(IM50),ICORE(I100),POPT3(IRPIJK),DIJK)
CSSS      CALL DAXPY(POPT3(IRPIJK), -1.0D0, ICORE(IM50), 1, 
CSSS     &           ICORE(ISTART), 1)
CSSS      CALL DCOPY(POPT3(IRPIJK), ICORE(ISTART), 1, ICORE(IM50), 1)
      CALL DCOPY(POPT3(IRPIJK), ICORE(ISTART), 1, ICORE(I140), 1)
      ENDIF 
C
C For L-based CCSD(T), we also need to compute
C    L3(ijk,abc)= L2(ij,{le}a)<bc||{le}k>
C and have both T3 and L3 contributions together. The current
C implementation requires two arrays of o^3v^3 vectors in memory.
C There is the possibility for an alternative implementation
C    T3(ijk,abc)o<ek||bc>=I(ij,ae) and E4=I(ij,ae)oL(ij,ae). The memory
C requirments for the second schemem is slightly less but it is not as straight
C forward. The alternative implementaion can be considered later if the
C memory requirements become excessive. (Eqns. 9 in IJQC) Ajith Perera, 10/05.
C
C
      IF (LCCSD_T) THEN
C Note that we are going to use the same space allocted for T2
C to store and manipulate L2 ans then reload T2 back for the
C rest of the code (this may not be necessary).
c Read L2 vectors (need offsets in each irrep_ioffoo from other
c routines should suffice).
C
         CALL GETLST(ICORE(I010),IOFFOO(IRPJ,IRPIJ,ISPIN)+1,
     1               NIJ,1,IRPIJ,143+ISPIN)
         CALL GETLST(ICORE(I020),IOFFOO(IRPK,IRPIK,ISPIN)+1,
     1               NIK,1,IRPIK,143+ISPIN)
         CALL GETLST(ICORE(I030),IOFFOO(IRPK,IRPJK,ISPIN)+1,
     1               NJK,1,IRPJK,143+ISPIN)
C
C Expand L2 vectors
c
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN),DSZEIJ,IRPDPD(IRPIJ,ISPIN),
     1                NIJ,ICORE(I040),ICORE(I010))
         CALL SYMEXP2(IRPIK,VRT(1,ISPIN),DSZEIK,IRPDPD(IRPIK,ISPIN),
     1                NIK,ICORE(I050),ICORE(I020))
         CALL SYMEXP2(IRPJK,VRT(1,ISPIN),DSZEJK,IRPDPD(IRPJK,ISPIN),
     1                NJK,ICORE(I060),ICORE(I030))
C
C L3(ijk,abc)= L2(ij,{le}a)<bc||{le}k>
C
         CALL ZERO(ICORE(IL00),LENABC)
         CALL ZERO(ICORE(IL50),POPT3(IRPIJK))
C
         CALL D1T14(ICORE(IL00),ICORE(I040),ICORE(I050),ICORE(I060),
     1              OOOV(OOOVAD(IRPIJ) +
     1              IOFFOO(IRPJ,IRPIJ,ISPIN)*IRPDPD(IRPIJ,15+ISPIN)),
     1              OOOV(OOOVAD(IRPIK) +
     1              IOFFOO(IRPK,IRPIK,ISPIN)*IRPDPD(IRPIK,15+ISPIN)),
     1              OOOV(OOOVAD(IRPJK) +
     1              IOFFOO(IRPK,IRPJK,ISPIN)*IRPDPD(IRPJK,15+ISPIN)),
     1              ICORE(ISTART),ISPIN,
     1              IADBLK,LENBLK,IADT2,LENT2,IADV,LENV,
     1              I,J,K,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1              DSZEIJ,DSZEIK,DSZEJK,
     1              IRPDPD(IRPIJ,15+ISPIN),IRPDPD(IRPIK,15+ISPIN),
     1              IRPDPD(IRPJK,15+ISPIN),NLEFT,.TRUE.)
C
C Convert the a<b;c L3 to a<b<c.
C 
         CALL SYMCONTW(ICORE(IL50),ICORE(IL00),IADBLK,ISPIN,IRPIJK)
C
         CALL RMD314(ICORE(IL50),ICORE(I100),POPT3(IRPIJK),DIJK)

         CALL E4TIJKL(ICORE(IM50),ICORE(IL50),ICORE(I100),
     &                POPT3(IRPIJK),DIJK,E4T14, 
     &                __FILE__)
C
C Reload the T2 vectors, so any methods other than CCSD(T) can
C smoothly proceed. 
C
         CALL GETLST(ICORE(I010),IOFFOO(IRPJ,IRPIJ,ISPIN)+1,
     1               NIJ,1,IRPIJ,43+ISPIN)
         CALL GETLST(ICORE(I020),IOFFOO(IRPK,IRPIK,ISPIN)+1,
     1               NIK,1,IRPIK,43+ISPIN)
         CALL GETLST(ICORE(I030),IOFFOO(IRPK,IRPJK,ISPIN)+1,
     1               NJK,1,IRPJK,43+ISPIN)
C
C Expand T2 vectorS
C
         CALL SYMEXP2(IRPIJ,VRT(1,ISPIN),DSZEIJ,IRPDPD(IRPIJ,ISPIN),
     1                NIJ,ICORE(I040),ICORE(I010))
         CALL SYMEXP2(IRPIK,VRT(1,ISPIN),DSZEIK,IRPDPD(IRPIK,ISPIN),
     1                NIK,ICORE(I050),ICORE(I020))
         CALL SYMEXP2(IRPJK,VRT(1,ISPIN),DSZEJK,IRPDPD(IRPJK,ISPIN),
     1                NJK,ICORE(I060),ICORE(I030))

      ELSE
         IF (.NOT. FINAL)
     &   CALL E4TIJK(ICORE(IM50),ICORE(I100),POPT3(IRPIJK),DIJK,
     &               E4T14,__FILE__)
C
C Endif for the L-CCSDT
C
      END IF
C
C     If this is a CCSD(T) or QCISD(T) gradient calculation, compute
C     "disconnected" triples at I140.  For LCCSD(T) only need the
C     disconnected lambda triples piece.  Put that at I140.
C
CAP      IF(((IDRLVL.GT.0) .AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))) THEN

      IF((IDRLVL.GT.0.AND.
     1  (ICLLVL.EQ.21.OR.ICLLVL.EQ.22).OR.EXTRNL_CCSD))THEN

      IF (.NOT. FINAL) THEN
      CALL ZERO(ICORE(I000),LENABC)
      CALL ZERO(ICORE(I140),POPT3(IRPIJK))
C
C This is T3(disconnected) and above is T3(connected), The LCCSD_T
C needs L3(disconneted) and the L3(connected) is already made (see
C above.
C
CSSS      Write(*,*) "Entering S1S214"
      CALL S1S214(ICORE(ISTART),ICORE(I000),IADBLK,ISPIN,I,J,K,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,1,LCCSD_T)
      IF(NONHF)THEN
CSSS      Write(*,*) "Entering S1S214 for NONHF"
      CALL S1S214(ICORE(ISTART),ICORE(I000),IADBLK,ISPIN,I,J,K,
     1            IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,2,LCCSD_T)
      ENDIF

      CALL SYMCONTW(ICORE(I140),ICORE(I000),IADBLK,ISPIN,IRPIJK)

#ifdef _DEBUG_NONHFTRP
C Notice that debug block can be used only derivative flags are
C set. Notice that I140 and ISTART also need to be correcttly 
C set. 
C
      CALL RMD314(ICORE(I140),ICORE(I100),POPT3(IRPIJK),DIJK)
      CALL DAXPY(POPT3(IRPIJK), 1.0D0, ICORE(IM50), 1, ICORE(I140),
     &           1)
      CALL E4TCTDIJK(ICORE(IM50), ICORE(I140), ICORE(I100),
     &               POPT3(IRPIJK),DIJK,
     &               E5TEST,__FILE__)
        WRITE(6,*) ' @TRPS14-I, E5ST from disconnected triples ',
     &                E5TEST
      Write(*,*)
      Call checksum("T3(C)+T3D-(I)", ICORE(I140), POPT3(IRREPIJK))
#endif 
      ENDIF 

      IF (.NOT. FINAL) 
     & CALL RMD314(ICORE(I140),ICORE(I100),POPT3(IRPIJK),DIJK)

C Endif for the (IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22)
      ENDIF

      IF (EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVL0
       if(iprint .eq. 1) write(6,"(2a)") " In trps14.F,Entered",
     &                                   " EXTRNL_CC block"
         iprint=-1
#endif
        ENDIF
#ifdef _DEBUG_LVLM
          call checksum("T31(C)---",icore(i000),lenabc)
#endif
C
C     ICORE(IM50) A<B<C ordered T3(C)
C     ICORE(IL50) A<B<C ordered L3(C)
C     ICORE(I000) A<B,C ordered T3(C) (scaled for QCISD(T)).
C ICORE(IL00) A<B,C ordered L3(C)
C ICORE(I140) A<B<C ordered T3(C)/L3(C) denominator removed
COLD ICORE(IL40) A<B<C ordered L3(D) denominator removed
C
C     Expand T3(C) and L3(C)
C
      IF(INT1.OR.INT2.OR.IDRLVL.GT.0.OR.EXTRNL_CCSD)THEN
      CALL ZERO(ICORE(I000),LENABC)
      CALL SYMEXPT3(ICORE(IM50),ICORE(I000),IADBLK,ISPIN,IRPIJK)

      IF (LCCSD_T) THEN
         CALL ZERO(ICORE(IL00),LENABC)
         CALL SYMEXPT3(ICORE(IL50),ICORE(IL00),IADBLK,ISPIN,IRPIJK) 
      END IF
      ENDIF
C
C     If this is a QCISD(T) calculation multiply expanded T3(C) by 2.
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC,2.0D+00,ICORE(I000),1)
      ENDIF
C
C     If this is a QCISD(T) or CCSD(T) gradient calculation, calculate
C     the contribution of T3(C) to GAMMA(IJ,AB).
C     Gamma(IJ,AB)= <T1||T3> in the case of CCSD(T) and
C     Gamma(IJ,AB)= <L1||T3> for L-CCSD(T). 
C     (Eqn. 40 in JCP2 and part of the first terms in Eqn. A22 in IJQC)
C     For external CCSD this contribution is included in FORM_HT3C3
C     (see below).

      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
CSSS      Write(*,*) "Entering T2FT314 and T2FT314 for GAMMA"

      IF(IUHF.GT.0)THEN
           CALL T2FT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,IRPJ,
     1                  IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,
     &                  ICLLVL,1,LCCSD_T)
      IF(NONHF)THEN
             CALL T2FT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,
     1                    IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,
     &                    K,ISPIN,ICLLVL,2,LCCSD_T)
      ENDIF
      ENDIF
      ENDIF
C
C     Include T3 in T1
C
C In the case of CCSD(T) gradients this WT3 (S(i,a)) can be added to
C Lambda (i.e. these terms remain constant) but if we are doing
C L-CCSD(T) gradients, these terms are the right hand side of a linear
C equation that determines a new set of amplitudes (designated as
C epsilon Crawford and Stanton, IJQC,70,601,Eqns. A16 and Eqn. 36 in
C JCP2). See below for the formation S(ij,ab).
C
      IF(INT1)THEN
CSSS      Write(*,*) "ENtering T1T314"
      CALL T1T314(ICORE(IM90),ICORE(IM70),ICORE(I000),
     1            ICORE(I110),ICORE(I120),ICORE(I130),
     1            IADBLK,IRPI,IRPJ,IRPK,
     1            IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,IOFFOO,ISPIN,
     &            NONHF,LCCSD_T)
      ENDIF
C
CAP      IF(IDRLVL.GT.0)THEN

      IF (IDRLVL.GT.0 .OR. EXTRNL_CCSD) THEN
C
C     GRADIENTS : DENSITY EVALUATION (DIAGONAL ELEMENTS).
C
C     If this is a QCISD(T) or CCSD(T) gradient calculation, form
C     expanded T3(D) in ICORE(ISTART).  For LCCSD(T) form expanded
C     L3(D) in ICORE(ISTART)
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
      CALL ZERO(ICORE(ISTART),LENABC)
      CALL SYMEXPT3(ICORE(I140),ICORE(ISTART),IADBLK,ISPIN,IRPIJK)
      END IF
C
C     If this is a QCISD(T) gradient calculation, divide expanded
C     2*T3(C) by 2 to recover expanded T3(C).
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(LENABC,0.5D+00,ICORE(I000),1)
      ENDIF
C
C     If this is a QCISD(T) or CCSD(T) gradient calculation, put expanded
C     T3(C) + alpha*T3(D) in ICORE(ISTART). If this is an MBPT(4) or
C     CCSD+T(CCSD) gradient simply copy T3(C) to ICORE(ISTART).
C
      IF(ICLLVL.EQ.21)THEN
      ALPHA = 2.0D+00
      ELSE
      IF(ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
      ALPHA = 1.0D+00
      ELSE
      IF(ICLLVL.EQ.4.OR.ICLLVL.EQ.9.OR.ICLLVL.EQ.11)THEN
      ALPHA = 0.0D+00
      ENDIF
      ENDIF
      ENDIF

      IF (LCCSD_T) THEN
C     For L-CCSD(T) gradients from  L3(C) + alpha*L3(D) in ICORE(ISTART)
         CALL VADD(ICORE(ISTART),ICORE(IL00),ICORE(ISTART),LENABC,
     &             ALPHA)
      ELSE
      IF (.NOT. FINAL) THEN
         CALL VADD(ICORE(ISTART),ICORE(I000),ICORE(ISTART),LENABC,
     &             ALPHA)
      ENDIF 

      ENDIF
C
C     For L-CCSD(T) gradients from  L3(C) + alpha*L3(D) in ICORE(ISTART)
C
CSSS      IF (LCCSD_T) CALL VADD(ICORE(ILSTART),ICORE(IL00),
CSSS     &                       ICORE(ILSTART),LENABC,ALPHA) 
C
C     Pass expanded T3(C) and T3(C) + alpha*T3(D) to virtual-virtual
C     density routine. The term computed here is <T3(C)|T3(C) + alpha*T3(D)>
C     for CCSD(T) and <L3(C)+alpha*L3(D)|T3(C)> for L-CCSD(T). This
C     is a contribution to the one particle density matrix D(a,b). (Eqn.
C     34 in JCP2 and 3rd term in Eqn. A22 in IJQC). 
C

      IF(ISPIN.EQ.1)THEN
      CALL DVT314(DVVA,ICORE(I000),ICORE(ISTART),
     1            IADBLK,IRPIJK,ISPIN,NVRTO(ISPIN))
#ifdef _DEBUG_LVLM
      Call checksum("DVVA", DVVA,NVRTO(1)*NVRTO(1))
#endif
      ELSE
      CALL DVT314(DVVB,ICORE(I000),ICORE(ISTART),
     1            IADBLK,IRPIJK,ISPIN,NVRTO(ISPIN))
      ENDIF
C
C     Repeat the procedure with the unexpanded T3's/L3's
C
      IF (LCCSD_T) THEN
         CALL VADD(ICORE(I140),ICORE(IL50),ICORE(I140),POPT3(IRPIJK),
     &      ALPHA)
      ELSE
            IF (.NOT. FINAL) THEN
            CALL VADD(ICORE(I140),ICORE(IM50),ICORE(I140),
     &                POPT3(IRPIJK), ALPHA)
            ENDIF 
      END IF

CSSS      IF (LCCSD_T) CALL VADD(ICORE(IL40),ICORE(IL50),ICORE(IL40),
CSSS     &                       POPT3(IRPIJK),ALPHA)
C
C
C     This is the contribution to the D(i,j) from the same terms
C     as above. Eqn. 35 in JCP2 and 3rd term in Eqn. A22 in IJQC.
C
      IF(ISPIN.EQ.1)THEN
      CALL DOT314(DOOA,ICORE(IM50),ICORE(I140),POPT3(IRPIJK),
     1            IRPI,IRPJ,IRPK,I,J,K,ISPIN,NOCCO(ISPIN))
#ifdef _DEBUG_LVLM
      Call checksum("DOOA", DOOA, NOCCO(1)*NOCCO(1))
#endif
      ELSE
      CALL DOT314(DOOB,ICORE(IM50),ICORE(I140),POPT3(IRPIJK),
     1            IRPI,IRPJ,IRPK,I,J,K,ISPIN,NOCCO(ISPIN))
      ENDIF
C
C     Clean-up step : rearrange so that we have the following
C
C     ICORE(IM50)   A<B<C ordered T3(C)
C     ICORE(IL50)   A<B<C ordered L3(C)
C     ICORE(I000)   A<B,C ordered T3(C)
C     ICORE(IL00)   A<B,C ordered L3(C)
C     ICORE(I140)   A<B<C ordered T3(D)/L3(D)
C     ICORE(ISTART) A<B,C ordered T3(D)/L3(D)
C
C     This step is unnecessary in MBPT(4) and CCSD+T(CCSD) gradient
C     calculations.
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
         IF (LCCSD_T) THEN
            CALL VADD(ICORE(I140),ICORE(I140),ICORE(IL50),POPT3(IRPIJK)
     &                ,-1.0D+00)
            CALL VADD(ICORE(ISTART),ICORE(ISTART),ICORE(IL00),LENABC
     &                ,-1.0D+00)
         ELSE
            IF (.NOT. FINAL) THEN
            CALL VADD(ICORE(I140),ICORE(I140),ICORE(IM50),POPT3(IRPIJK)
     &                ,-1.0D+00)
            CALL VADD(ICORE(ISTART),ICORE(ISTART),ICORE(I000),LENABC
     &                ,-1.0D+00)
            ENDIF 

         ENDIF
      ENDIF
C
CSSS      IF (LCCSD_T) THEN
CSSS         CALL VADD(ICORE(IL40),ICORE(IL40),ICORE(IL50),POPT3(IRPIJK),
CSSS     &            -1.00D0)
CSSS         CALL VADD(ICORE(ILSTART),ICORE(ILSTART), ICORE(IL00), 
CSSS     &             LENABC, -1.00D0)
CSSS      ENDIF
CSSS      ENDIF
C
C     For QCISD(T) remove factor of 2.
C
      IF(ICLLVL.EQ.21)THEN
      CALL SSCAL(POPT3(IRPIJK),0.5D+00,ICORE(I140),1)
      CALL SSCAL(LENABC,0.5D+00,ICORE(ISTART),1)
      ENDIF
C
C     Now form expanded and unexpanded T3(C) + (alpha/2)*T3(D) in T3
C     and ICORE(I000). This is only necessary for CCSD(T) or QCISD(T)
C     (for MBPT(4) and CCSD+T(CCSD) T3 and ICORE(I000) already have what
C     is needed).
C
      IF(ICLLVL.EQ.21)THEN
      ALPHA = 1.0D+00
      ENDIF
      IF(ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
      ALPHA = 0.5D+00
      ENDIF
      IF (LCCSD_T) THEN
        ALPHA = 1.0d0
      ENDIF
#ifdef _DEBUG_LVLM
          call checksum("T32(C)---",icore(i000),lenabc)
#endif

      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.EXTRNL_CCSD)THEN
C     CALL VADD(T3,T3,ICORE(I140),POPT3(IRPIJK),ALPHA)
         IF (LCCSD_T) THEN
C     In LCCSD(T) case, put l3(c)+alpha*l3(d) at IL00
            CALL VADD(ICORE(IL00),ICORE(IL00),ICORE(ISTART),LENABC
     &         ,ALPHA)
         ELSE
            IF (.NOT. FINAL) THEN
            CALL VADD(ICORE(I000),ICORE(I000),ICORE(ISTART),LENABC,
     &         ALPHA)
            ENDIF
         ENDIF
CSSS      IF (LCCSD_T) CALL VADD(ICORE(IL00),ICORE(IL00),ICORE(ILSTART),
CSSS                             LENABC,ALPHA) 
      ENDIF
#ifdef _DEBUG_LVLM
          call checksum("T33(C)---",icore(i000),lenabc)
#endif
C
C     GRADIENTS : GAMMA EVALUATION. The term computed here are <T3(C)|2T3(C)
C     + T3(D)> for CCSD(T) and <L3(C)+L3(D)|T3(C)> for L-CCSD(T) contributions
C     to the Gamma(ij,ka) and Gamma(ab,cd) (Eqns. 38 and 39 in JCP2 and second
C     term in Eqn. A22 in IJQC). 
      IF(IUHF.GT.0)THEN
         IF (LCCSD_T) THEN
C Scale Lambda3 and T3(c) because there is an implicit factor of two
C for gamma
            call dscal(lenabc,0.5d0,icore(il00),1)
            call dscal(lenabc,0.5d0,icore(i000),1)
C     Add contributions from l3(c)+alpha*l3(d)
            CALL GVT314(ICORE(IL00),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     &         IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,.FALSE.)
            CALL GOT314(GOOOV,OOOVAD,LNOOOV,ICORE(IL00),ICORE(ISTART),
     &         IADBLK,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,ISPIN
     &         ,.FALSE.)

C     Add contributions from t3(c)
            CALL GVT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     &         IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,.TRUE.)
            CALL GOT314(GOOOV,OOOVAD,LNOOOV,ICORE(I000),ICORE(ISTART),
     &         IADBLK,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,ISPIN
     &         ,.TRUE.)
C Scale Lambda3 and T3(c) because there is not an implicit factor of two
C for everything below in the calculation
            call dscal(lenabc,2.0d0,icore(il00),1)
            call dscal(lenabc,2.0d0,icore(i000),1)
         ELSE
         IF (.NOT. EXTRNL_CCSD) THEN
CSSS      Write(*,*) "Entering GVT314 and GOT314"
           CALL GVT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     &         IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,.FALSE.)
           CALL GOT314(GOOOV,OOOVAD,LNOOOV,ICORE(I000),ICORE(ISTART),
     &         IADBLK,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,I,J,K,ISPIN
     &       ,.FALSE.)
         ENDIF
      ENDIF 
      ENDIF
      ENDIF 

      IF (EXTRNL_CCSD) THEN
           CALL T32ABCI14(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,
     &                    IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,
     &                    K,ISPIN,.FALSE.)
           CALL T32IJKA14(GOOOV,OOOVAD,LNOOOV,ICORE(I000),
     &                    ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,IRPIJ,
     &                    IRPIK,IRPJK,I,J,K,ISPIN,.FALSE.)
      ENDIF
C
CSSS      Write(*,*) "INT2-- Before T3-T2 and T3-T1", INT2
      IF(INT2.AND.IUHF.GT.0)THEN
C
C     INCLUDE T3 IN T2. In the case of CCSD(T) gradients this
C     WT3 (S(ij,ab)) can be added to Lambda (ie. these terms are
C     remains constant) but if we are doing L-CCSD(T) gradients, these
C     terms are the right hand side of a linear equations 
C     that detemine a new set of amplitudes (designated as epsilon
C     Crawford and Stanton, IJQC,70,601,Eqns. A16 and Eqn. 37 in JCP2).
C

CSSS      Write(*,*) "Entering T2T314O"
      CALL T2T314O(ICORE(I000),ICORE(I070),ICORE(I080),ICORE(I090),
     1             ICORE(ISTART),
     1             IADBLK,
     1             IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,
     1             I,J,K,ISPIN,LNVOIJ,LNVOIK,LNVOJK,.FALSE.)

      IF (LCCSD_T) THEN
C     D2L2 = L3 W
         CALL T2T314O(ICORE(IL00),ICORE(I070),ICORE(I080),ICORE(I090),
     &      ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,
     &      IRPIJK,I,J,K,ISPIN,LNVOIJ,LNVOIK,LNVOJK,.TRUE.)
      ENDIF
C
C     --- D2T2 = F*T3 in non-Hartree-Fock cases ---
C
      IF((ICLLVL.EQ.13.AND.NONHF).OR. 
     1   (ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
CSSS      Write(*,*) "Entering T2FT314"
      CALL T2FT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     1             IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,
     1             ICLLVL,0,.FALSE.)
      ENDIF

      IF (EXTRNL_CCSD .AND. NONHF) THEN
         CALL T2FT314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,
     &                IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,
     &                K,ISPIN,ICLLVL,0,.FALSE.)
      ENDIF 
CSSS      Write(*,*) "Entering T2T314"
      CALL  T2T314(ICORE(I000),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     1             IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,
     1             ISPIN,.FALSE.)
      IF (LCCSD_T) THEN
C     D2L2 = L3 W
         CALL  T2T314(ICORE(IL00),ICORE(ISTART),IADBLK,IRPI,IRPJ,IRPK,
     &      IRPIJ,IRPIK,IRPJK,IRPIJK,I,J,K,ISPIN,.TRUE.)
      ENDIF
      ENDIF
C This is no longer needed since T3 -> T2 are computed by default since
C INT1 and INT2 is true for extrnl_ccsd. Note that if we need to call
C we need to store T3 to a diffrent location before T2T314 call since
C it destroys it. 
CSSS      IF (EXTRNL_CCSD) THEN
CSSS      CALL FORM_HT3C14(ICORE(I000),ICORE(I070),ICORE(I080),
CSSS     &                      ICORE(I090),ICORE(ISTART),IADBLK,IRPI,
CSSS     &                      IRPJ,IRPK,IRPIJ,IRPIK,IRPJK,IRPIJK,I,
CSSS     &                      J,K,ISPIN,LNVOIJ,LNVOIK,LNVOJK,.FALSE.,
CSSS     &                      NONHF,IUHF)
CSSS      ENDIF
C
C     Write out the amplitudes. Now the fun begins !
C
c      IJK  = IJKPOS(IRPI,IRPJ,IRPK,1)
c      IOFF = IJKOFF(IJK,IRPIJK,1+3*(ISPIN-1))
c      IF(IJKEQL) IJKVAL = IOFF + ((K-1)*(K-2)*(K-3))/6 + INDEX(J-1) + I
c      IF(IJEQL ) IJKVAL = IOFF + (K-1)*NIJ             + INDEX(J-1) + I
c      IF(JKEQL ) IJKVAL = IOFF + (INDEX(K-1)+ J-1)*POP(IRPI,ISPIN)  + I
c      IF(NONEQL) IJKVAL = IOFF + (K-1)*NIJ + (J-1)*POP(IRPI,ISPIN)  + I
C     write(6,*) ' @TRPS14-I, sym, record ',IRPIJK,IJKVAL

      IF(T3STOR)THEN
      CALL PUTLIST(ICORE(IM50),IJKVAL,1,1,IRPIJK,1+3*(ISPIN-1))
      ENDIF
#ifdef _DEBUG_LVLM
      Write(*,*)
      Call checksum("TRPS14", ICORE(IM50), POPT3(IRREPIJK))
#endif
  410 CONTINUE
  420 CONTINUE
  430 CONTINUE
C
  980 CONTINUE
  990 CONTINUE
 1000 CONTINUE
C
      IF (EXTRNL_CCSD) THEN
          CALL GETLST(ICORE(IM80),1,1,1,2+ISPIN,90)
          CALL VADD(ICORE(IM80),ICORE(IM80),ICORE(IM90),NT1,1.0D+00)
          CALL PUTLST(ICORE(IM80),1,1,1,2+ISPIN,90)

          IF (IUHF.GT.0)THEN
             IF (ISPIN.EQ.1)THEN
                DO IRREP=1,NIRREP
                   LEN = IRPDPD(IRREP,2+ISPIN)*IRPDPD(IRREP,15+ISPIN)
CSSS                   CALL DAXPY(LEN,1.0D0,OOOV2(OOOVAD(IRREP)),1,
CSSS     &                        GOOOV(OOOVAD(IRREP)),1)
                   CALL PUTTRN(GOOOV(OOOVAD(IRREP)),BUF,
     &                         IRPDPD(IRREP,2+ISPIN),
     &                         IRPDPD(IRREP,15+ISPIN),2,IRREP,107)
                ENDDO
             ELSE
                DO IRREP=1,NIRREP
                   LEN = IRPDPD(IRREP,2+ISPIN)*IRPDPD(IRREP,15+ISPIN)
CSSS                   CALL DAXPY(LEN,1.0D0,OOOV2(OOOVAD(IRREP)),1,
CSSS     &                        GOOOV(OOOVAD(IRREP)),1)
                   CALL PUTTRN(GOOOV(OOOVAD(IRREP)),BUF,
     &                         IRPDPD(IRREP,2+ISPIN),
     &                         IRPDPD(IRREP,15+ISPIN),2,IRREP,108)
                ENDDO
             ENDIF
          ENDIF
#ifdef _DEBUG_LVL0
          Write(6,"(a)") " ---After trps14,Lists:61-63,90(3,4)---"
CSSS          call check_ints(icore(istart),nleft,iuhf,.false.)
          call check_t32w(icore(istart),nleft,iuhf)
          call check_t2(icore(istart),nleft,iuhf)
#endif
      ENDIF 
C
C     If this is a gradient calculation, write ooov gamma to disk.
C
      IF(IDRLVL.GT.0.AND.IUHF.GT.0)THEN
      IF(ISPIN.EQ.1)THEN
      DO  285 IRREP=1,NIRREP
      CALL PUTTRN(GOOOV(OOOVAD(IRREP)),BUF,IRPDPD(IRREP,2+ISPIN),
     1            IRPDPD(IRREP,15+ISPIN),2,IRREP,LISGO1)
  285 CONTINUE
      ELSE
      DO  286 IRREP=1,NIRREP
      CALL PUTTRN(GOOOV(OOOVAD(IRREP)),BUF,IRPDPD(IRREP,2+ISPIN),
     1            IRPDPD(IRREP,15+ISPIN),2,IRREP,LISGO2)
  286 CONTINUE
      ENDIF
C
      ENDIF
C
      IF (.NOT. FINAL) THEN
      IF(ISPIN.EQ.1)THEN
      E4TAAA = E4T14
      WRITE(6,6000) E4TAAA
 6000 FORMAT(' @TRPS14-I, E4TAAA ',F20.12)
      ELSE
      E4TBBB = E4T14
      WRITE(6,6010) E4TBBB
 6010 FORMAT(' @TRPS14-I, E4TBBB ',F20.12)
      ENDIF
      
      ELSE

      IF(ISPIN.EQ.1)THEN
      E4TAAA = E4T14
      WRITE(6,6001) E4TAAA 
 6001 FORMAT(' @TRPS14-I, E4TAAA+E5STA+E4STA ',F20.12)
      ELSE
      E4TBBB = E4T14
      WRITE(6,6012) E4TBBB
 6012 FORMAT(' @TRPS14-I, E4TBBB+E5STA+E4STA ',F20.12)
      ENDIF

      ENDIF 
#ifdef _DEBUG_LVLM
      if (.not. extrnl_ccsd) then
      Write(6,"(a)") " ---After trps14,Lists:61-63,90(3,4)---"
CSSS      call check_ints(icore(istart),nleft,iuhf,.false.)
CSSS      call check_t32w_grad(icore(istart),nleft,iuhf)
      call check_t2(icore(istart),nleft,iuhf)
      call check_t2grd(icore(istart),nleft,iuhf)
      endif
#endif
C
C      IF THIS IS ITERATIVE TRIPLES WRITE T1 INCREMENTS TO DISK
C
      IF((ICLLVL.GE.13.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
        CALL GETLST(ICORE(IM80),1,1,1,2+ISPIN,90)
        CALL VADD(ICORE(IM80),ICORE(IM80),ICORE(IM90),NT1,1.0D+00)
        CALL PUTLST(ICORE(IM80),1,1,1,2+ISPIN,90)
      ENDIF
C
C     If this is CCSD(T) or QCISD(T) gradient, write D1T1 triples contribution
C     to disk.
C
      IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))THEN
        CALL GETLST(ICORE(IM80),1,1,1,ISPIN,4)
        CALL VADD(ICORE(IM80),ICORE(IM80),ICORE(IM90),NT1,1.0D+00)
        CALL PUTLST(ICORE(IM80),1,1,1,ISPIN,4)
        IF(NONHF)THEN
          CALL GETLST(ICORE(IM80),1,1,1,4+ISPIN,160)
          CALL VADD(ICORE(IM80),ICORE(IM80),ICORE(IM70),NT1,0.5D+00)
          CALL PUTLST(ICORE(IM80),1,1,1,4+ISPIN,160)
        ENDIF
      ENDIF
C
C     IF THIS IS MBPT(4) OR ITERATIVE TRIPLES RETURN
C
      IF(ICLLVL.EQ.9.OR.ICLLVL.EQ.13.OR.ICLLVL.EQ.14.OR.ICLLVL.EQ.15.OR.
     1   ICLLVL.EQ.16.OR.ICLLVL.EQ.17.OR.ICLLVL.EQ.18.OR.
     1   ICLLVL.EQ.33.OR.ICLLVL.EQ.34) RETURN
      IF(ICLLVL.EQ.4.AND.IREFNC.LE.1.AND..NOT.NONHF) RETURN

      IF (.NOT. FINAL) THEN

C For iterative NON-HF triples, the total contribution is
C already computed)
C
C     COMPUTE CCSD(T) OR QCISD(T) FIFTH-ORDER PIECE OF ENERGY
C
C For L-based CCSD instead of T1, we use L1 to compute the 
C 5th order triples contribution. Ajith Perera, 10/05.

      IF (LCCSD_T) THEN
         CALL GETLST(ICORE(IM80),1,1,1,ISPIN,190)
      ELSE
         CALL GETLST(ICORE(IM80),1,1,1,ISPIN,90)
      ENDIF
C
      E5ST14 = 0.0D+00
      DO  7010 AI = 1,NT1
      E5ST14 = E5ST14 + ICORE(IM90 + AI - 1) * ICORE(IM80 + AI - 1)
 7010 CONTINUE
C
      IF(ISPIN.EQ.1)THEN
      E5STA = E5ST14
      WRITE(6,6020) E5STA
 6020 FORMAT(' @TRPS14-I, E5ST A ',F20.15)
      ELSE
      E5STB = E5ST14
      WRITE(6,6030) E5STB
 6030 FORMAT(' @TRPS14-I, E5ST B ',F20.15)
      ENDIF
C
C     IF THIS IS A HARTREE-FOCK REFERENCE RETURN
C
      IF(.NOT.NONHF) RETURN
C
C     COMPUTE NON-HF FOURTH-ORDER PIECE FOR CCSD(T) OR QCISD(T).
C
      CALL GETLST(ICORE(IM60),1,1,1,2+ISPIN,93)
C
      E4ST14 = 0.0D+00
      DO  8010 AI = 1,NT1
      E4ST14 = E4ST14 + ICORE(IM70 + AI - 1) * ICORE(IM60 + AI - 1)
 8010 CONTINUE
C
      IF(ISPIN.EQ.1)THEN
      E4STA = E4ST14
      WRITE(6,8020) E4STA
 8020 FORMAT(' @TRPS14-I, E4ST A ',F20.15)
      ELSE
      E4STB = E4ST14
      WRITE(6,8030) E4STB
 8030 FORMAT(' @TRPS14-I, E4ST B ',F20.15)
      ENDIF
      
      ENDIF 

      RETURN
 9010 FORMAT(' @TRPS14-I, Insufficient memory to continue. ')
      END
