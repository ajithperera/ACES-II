
#include "flags.h"

      SUBROUTINE TRPS_PCCD(ICORE,MAXCOR,IUHF,E4T,FINAL,EXTRNL_CCSD,
     +                     PCCD)
      IMPLICIT INTEGER (A-Z)

      DOUBLE PRECISION E4T,E4TAAA,E4TBBB,E4TBAA,E4TABB,E5ST,
     1                 E5STA,E5STB,E4ST,E4STA,E4STB,ECCTOT,
     1                 ETCCSD,EPAREN
      DOUBLE PRECISION E5TEST,E6TEST
      DOUBLE PRECISION DOOA,DOOB,DVVA,DVVB
      DOUBLE PRECISION ALLT3
      DOUBLE PRECISION E45T(20)
      LOGICAL INT1,INT2,NONHF,HF,FINAL
      LOGICAL ISKP14,ISKP2,ISKP3
      LOGICAL CCSDT4,CCSDT
      LOGICAL INCORE
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      LOGICAL T3GRAD
      LOGICAL EXTRNL_CCSD
      LOGICAL PCCD
#include "maxbasfn.par"
#include "trps.par"
      DATA MAXNBF/MAXBASFN/
      DIMENSION ICORE(MAXCOR)

      integer in_year,  in_mon,  in_mday,  in_hour,  in_min,  in_sec
      integer out_year, out_mon, out_mday, out_hour, out_min, out_sec

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
      EQUIVALENCE(IPTORB,IFLAGS(64))
      COMMON /FLAGS2/  IFLAGS2(500)
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXBASFN,8,2),ABSOCC(MAXBASFN,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TBAA,E4TABB,E5STA,E5STB,
     1                E4STA ,E4STB,E4TBAA2, E4TBAA3
      COMMON /DENST3/ DOOA(MAXBASFN),DOOB(MAXBASFN),
     1                DVVA(MAXBASFN),DVVB(MAXBASFN)
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     1                LISGV3,LISGV4
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     &                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /TESTEN/ E5TEST, E6TEST
      COMMON /T3METH/ CCSDT4,CCSDT
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3CORE/ INCORE
C-----------------------------------------------------------------------
C
      WRITE(6,1005)
 1005 FORMAT(/,72('-'),/,T17,
     &       ' VCC Triple Excitation Contributions ',/,72('-'),/)
C
      CALL ZERO(E45T,20)
      E4T  = 0.0D+00
      E5ST = 0.0D+00
      E4ST = 0.0D+00
C
      CCSDT4 = ICLLVL.EQ.17
      CCSDT  = ICLLVL.EQ.18.OR.ICLLVL.EQ.29.OR.ICLLVL.EQ.30
C
      CALL DIMT3(DISTSZ(1,1),NDISTS(1,1))
      CALL DIMT3(DISTSZ(1,5),NDISTS(1,5))
C
C-----------------------------------------------------------------------
C     If triples are stored and there is enough core, put in ALLT3.
C-----------------------------------------------------------------------
      IF(T3STOR)THEN
        NT3 = NRECS * LNPHYR
         Write(6,"(a,a,a,2(1x,I10))") " The number of T3",
     &                                    " amps and the number that",
     &                                    " can be kept in memory:",
     &                                      NT3, NUMT3
        Write(6,*) 
        IF(NT3.LT.NUMT3)THEN
          INCORE = .TRUE.
          WRITE(6,9010)
        ELSE
          INCORE = .FALSE.
          WRITE(6,9020)
        ENDIF
        IF(INCORE) CALL GETALLT3
      ENDIF
 9010 FORMAT(' @TRPS-I, There is sufficient core for the T3 file. ')
 9020 FORMAT(' @TRPS-I, There is insufficient core for the T3 file. ')
C-----------------------------------------------------------------------
C
      E5TEST = 0.0D+00
      E6TEST = 0.0D+00
C
C-----------------------------------------------------------------------
C     Debug aids to skip one or more spin cases.
C-----------------------------------------------------------------------
      ISKP14 = .FALSE.
      ISKP2  = .FALSE.
      ISKP3  = .FALSE.
C
C-----------------------------------------------------------------------
C     For special noniterative HF triples we pretend we are doing
C     a simpler calculation until after xvcc, although we of course
C     store the T3 amplitudes.
C
C     CCSD+TQ*, CCSD(TQ), CCSD+TQ, CC5SD[T] ---> CCSD(T)
C     QCISD(TQ)                             ---> QCISD(T)
C     CCSDT+Q*, CCSDT+Q                     ---> CCSDT
C-----------------------------------------------------------------------
C
      IF(ICLLVL.EQ.12.OR.ICLLVL.EQ.27.OR.ICLLVL.EQ.28.OR.ICLLVL.EQ.31)
     &                                                              THEN
        METACT = ICLLVL
        ICLLVL = 22
      ELSEIF(ICLLVL.EQ.26)THEN
        METACT = ICLLVL
        ICLLVL = 21
      ELSEIF(ICLLVL.EQ.29.OR.ICLLVL.EQ.30)THEN
        METACT = ICLLVL
        ICLLVL = 18
      ENDIF
C-----------------------------------------------------------------------
C     If IFLAGS(3)=1 but this is not an analytical derivative calculat-
C     ion, reset IFLAGS(3). Also, check we have perturbed canonical orb-
C     itals (should be unnecessary to check this these days).
C-----------------------------------------------------------------------
C
      T3GRAD = .FALSE.
      IF(IDRLVL.GT.0)THEN
C
        T3GRAD = .TRUE.
C
C     If this is a numerical gradient and IDRLVL > 0, turn the
C     derivative flag off, but turn on at end.
C
CSSS        Write(*,*) "IFLAGS2(h_IFLAGS2_grad_calc)",
CSSS     &              IFLAGS2(h_IFLAGS2_grad_calc)
        IF (IFLAGS2(h_IFLAGS2_grad_calc).EQ.2) THEN
          if (iFlags2(h_IFLAGS2_geom_opt).ne.0) WRITE(6,1053)
          if (iFlags(h_IFLAGS_vib)       .eq.3) WRITE(6,1054)
          if (iFlags(h_IFLAGS_vib)       .eq.3) IDROLD = IDRLVL
          IDRLVL = 0
          T3GRAD = .FALSE.
        ENDIF
C
C     If this is a real gradient calculation, make sure we have
C     perturbed canonical orbitals.
C
        IF(IFLAGS(87).NE.3.AND.IFLAGS(64).NE.1.AND.T3GRAD)THEN
          WRITE(6,1052)
          CALL ERREX
        ENDIF
      ENDIF
 1052 FORMAT(' @TRPS-I, Gradients require perturbed canonical ',
     &       ' orbitals ! ',/,
     &       '          Please resubmit with PERT_ORB=CANONICAL ')
 1053 FORMAT(' @TRPS-I, Numerical optimization. IDRLVL set to 0. ')
 1054 FORMAT(' @TRPS-I, Energy only frequencies. IDRLVL set to 0. ')
C
C-----------------------------------------------------------------------
C     Determine kind of reference function.
C-----------------------------------------------------------------------
      NONHF = INONHF.NE.0.OR.
     1        IREFNC.EQ.1.AND.(IQRHFP.NE.0.OR.IQRHFM.NE.0).OR.
     1        IREFNC.EQ.2
      HF = .NOT.NONHF
C
      WRITE(6,1061)
      IF(NONHF)THEN
        WRITE(6,1062)
      ELSE
        WRITE(6,1063)
      ENDIF
C
      IF(IREFNC.EQ.0.AND.HF)    WRITE(6,1064)
      IF(IREFNC.EQ.0.AND.NONHF) WRITE(6,1065)
      IF(IREFNC.EQ.1.AND.HF)    WRITE(6,1066)
      IF(IREFNC.EQ.1.AND.NONHF) WRITE(6,1067)
      IF(IREFNC.EQ.2)           WRITE(6,1068)
      IF(IREFNC.EQ.1.AND.(IQRHFP.NE.0.OR.IQRHFM.NE.0)) WRITE(6,1069)
C
C     Do not attempt noniterative triple excitation calculations
C     for non-HF orbitals unless we have semicanonical orbitals.
C
      IF(TRIPNI.AND.NONHF.AND.IORBTP.EQ.0)THEN
        WRITE(6,1070)
        CALL ERREX
      ENDIF
 1070 FORMAT(' @TRPS-I, Caution : This calculation',
     1 ' is about to evaluate noniterative',/,
     1 '          triples corrections with',
     1 ' standard non-HF orbitals. ACES2 devel-',/,
     1 '          opers do not believe the',
     1 ' results are meaningful. We suggest',/,
     1       '          Please resubmit with ORBITAL=SEMICANONICAL ')
C
C-----------------------------------------------------------------------
C     Set flags for including T3 in T2 and T1 equations. Remember that
C     UCC(4) does not put T3 in the T1 equation.
C
C     Note : for CCSD(T) and QCISD(T), INT1 means "compute D1T1 increm-
C     ents but only increment lists if this is a gradient calculation".
C-----------------------------------------------------------------------
      INT1 = .FALSE.
      INT2 = .FALSE.
      INT1 = (TRIPIT.AND.ICLLVL.NE. 9)    .OR.
     &        ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.
     &       (ICLLVL.EQ.4.AND.(IREFNC.EQ.2.OR.IQRHFP.NE.0.OR.
     &                                        IQRHFM.NE.0.OR.
     &                                        NONHF))
      INT2 = TRIPIT .OR.
     &      (ICLLVL.EQ. 4.AND.IDRLVL.GT.0).OR.
     &      (ICLLVL.EQ.11.AND.IDRLVL.GT.0).OR.
     &      (ICLLVL.EQ.21.AND.IDRLVL.GT.0).OR.
     &      (ICLLVL.EQ.22.AND.IDRLVL.GT.0)

      IF (EXTRNL_CCSD .AND. .NOT. (INT1 .AND. INT2)) THEN
         INT1 = .TRUE.
         INT2 = .TRUE.
C
C At the moment external CCSD is restricted to non-iterative
C triples. Lets use 61-63 increment lists to store <ab||ij> 
C contributions from T3*(W1,W2)->W2 contributions. Similarly
C use 90(3,4) to store T3*w2-> W1 contributions. Lets zero
C them out here.

          IF (IUHF .NE. 0) THEN
             CALL ZEROLIST(ICORE,MAXCOR,61)
             CALL ZEROLIST(ICORE,MAXCOR,62)
             CALL ZEROLIST(ICORE,MAXCOR,63)
             CALL ZERLST(ICORE,NTAA,1,1,3,90)
             CALL ZERLST(ICORE,NTBB,1,1,4,90)
          ELSE 
             CALL ZEROLIST(ICORE,MAXCOR,63)
             CALL ZERLST(ICORE,NTAA,1,1,3,90)
          ENDIF 
          LIST2I1 = 61
          LIST2I2 = 62
          LIST2I3 = 63
          LISGV1  = 127
          LISGV2  = 128
          LISGV3  = 129
          LISGV4  = 130
          LISGO1  = 107
          LISGO2  = 108
          LISGO3  = 109
          LISGO4  = 110
          CALL SET_EXTRNL_CC_LISTS(ICORE,MAXCOR,IUHF,0)
          IF (IUHF .NE. 0) THEN
             CALL ZEROLIST(ICORE,MAXCOR,107)
             CALL ZEROLIST(ICORE,MAXCOR,108)
             CALL ZEROLIST(ICORE,MAXCOR,109)
             CALL ZEROLIST(ICORE,MAXCOR,110)
             CALL ZEROLIST(ICORE,MAXCOR,127)
             CALL ZEROLIST(ICORE,MAXCOR,128)
             CALL ZEROLIST(ICORE,MAXCOR,129)
             CALL ZEROLIST(ICORE,MAXCOR,130)
         ELSE
             CALL ZEROLIST(ICORE,MAXCOR,107)
             CALL ZEROLIST(ICORE,MAXCOR,110)
             CALL ZEROLIST(ICORE,MAXCOR,130)
         ENDIF 
      ENDIF 
C
cYAU      IF(     INT1) WRITE(6,1072)
cYAU      IF(.NOT.INT1) WRITE(6,1073)
cYAU      IF(     INT2) WRITE(6,1074)
cYAU      IF(.NOT.INT2) WRITE(6,1075)
      IF(     INT1) WRITE(6,1072)
      IF(.NOT.INT1) WRITE(6,1073)
      IF(     INT2) WRITE(6,1074)
      IF(.NOT.INT2) WRITE(6,1075)
C-----------------------------------------------------------------------
C
      MAXPOP = 0
      MAXVRT = 0
      MAX1   = 0
      MAX2   = 0
      MAX3   = 0
      MAX4   = 0
      MAX9   = 0
      MAX10  = 0
      MAX11  = 0
      MAX12  = 0
      MAX13  = 0
      MAX14  = 0
      DO   20 IRREP=1,NIRREP
      IF(POP(IRREP,1).GT.MAXPOP) MAXPOP = POP(IRREP,1)
      IF(POP(IRREP,2).GT.MAXPOP) MAXPOP = POP(IRREP,2)
      IF(VRT(IRREP,1).GT.MAXVRT) MAXVRT = VRT(IRREP,1)
      IF(VRT(IRREP,2).GT.MAXVRT) MAXVRT = VRT(IRREP,2)
C
      IF(IRPDPD(IRREP, 1).GT.MAX1 ) MAX1   = IRPDPD(IRREP, 1)
      IF(IRPDPD(IRREP, 2).GT.MAX2 ) MAX2   = IRPDPD(IRREP, 2)
      IF(IRPDPD(IRREP, 3).GT.MAX3 ) MAX3   = IRPDPD(IRREP, 3)
      IF(IRPDPD(IRREP, 4).GT.MAX4 ) MAX4   = IRPDPD(IRREP, 4)
      IF(IRPDPD(IRREP, 9).GT.MAX9 ) MAX9   = IRPDPD(IRREP, 9)
      IF(IRPDPD(IRREP,10).GT.MAX10) MAX10  = IRPDPD(IRREP,10)
      IF(IRPDPD(IRREP,11).GT.MAX11) MAX11  = IRPDPD(IRREP,11)
      IF(IRPDPD(IRREP,12).GT.MAX12) MAX12  = IRPDPD(IRREP,12)
      IF(IRPDPD(IRREP,13).GT.MAX13) MAX13  = IRPDPD(IRREP,13)
      IF(IRPDPD(IRREP,14).GT.MAX14) MAX14  = IRPDPD(IRREP,14)
   20 CONTINUE
      MAXOO  = MAX(MAX3,MAX4 ,MAX14)
      MAXVV  = MAX(MAX1,MAX2 ,MAX13)
      MAXOV  = MAX(MAX9,MAX10,MAX11,MAX12)
      SCRDIM = MAX(MAXOV,MAXVV,MAXOO)
C
      NOCA  = NOCCO(1)
      NOCB  = NOCCO(2)
      NVRTA = NVRTO(1)
      NVRTB = NVRTO(2)
C
C     Evaluate symmetry packed lengths of <ij//ka> integrals/intermediates.
C
      LNOOOV1 = 0
      LNOOOV2 = 0
      LNOOOV3 = 0
      LNOOOV4 = 0
      DO   30 IRPIJ=1,NIRREP
      LNOOOV1 = LNOOOV1 + IRPDPD(IRPIJ, 3) * IRPDPD(IRPIJ, 9)
      LNOOOV2 = LNOOOV2 + IRPDPD(IRPIJ, 4) * IRPDPD(IRPIJ,10)
      LNOOOV3 = LNOOOV3 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,11)
      LNOOOV4 = LNOOOV4 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,12)
   30 CONTINUE
C
      LNOOOV = MAX0(LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4)
cYAU      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
C
C     Set up orbital index tables
C
      CALL ABSORB(ABSVRT,ABSOCC,MAXNBF)
C
C     Calculate VV, OO, VO offsets
C
      CALL MKOFVV
      CALL MKOFOO
      CALL MKOFVO
C
C     Set W intermediate list numbers
C
      CALL INITRP
C
C     Form intermediates for ccsdt-1b and higher calculations.
C
      IF((ICLLVL.GE.14.AND.ICLLVL.LE.18).OR.ICLLVL.EQ.33.OR.
     &                                      ICLLVL.EQ.34)THEN
        CALL T1BINT(ICORE,MAXCOR,IUHF)
      ENDIF
C
C     Transpose records of LWIC17 and LWIC27 (the latter only if it is
C     not the same list).
C
      IF(IUHF.EQ.1)THEN
        CALL TSPABCI(ICORE,MAXCOR,1,EXTRNL_CCSD)
      ENDIF

c     Calculate triples contributions to W(abci) and W(ijka) (CCSDT).
C
      IF(CCSDT) CALL W45T3DRV(ICORE,MAXCOR,IUHF,EXTRNL_CCSD)
C
C     If this is a gradient calculation, initialize diagonal elements
C     of density matrix
C
      IF(IDRLVL.GT.0)THEN
        CALL ZERO(DOOA,MAXNBF)
        CALL ZERO(DOOB,MAXNBF)
        CALL ZERO(DVVA,MAXNBF)
        CALL ZERO(DVVB,MAXNBF)
      ENDIF
C
C     If this is a gradient calculation, set gamma list numbers
C
      IF(IDRLVL.GT.0)THEN
C
C     ---- MBPT(4), QCISD(T), UCC(4) ----
C
      IF(ICLLVL.EQ.4.OR.ICLLVL.EQ.21.OR.ICLLVL.EQ.9)THEN
      LISGO1 = 107
      LISGO2 = 108
      LISGO3 = 109
      LISGO4 = 110
      LISGV1 = 127
      LISGV2 = 128
      LISGV3 = 129
      LISGV4 = 130
      ENDIF
C
C     ---- CCSD + T(CCSD), CCSD(T) ----
C
      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.22)THEN
      LISGO1 = 162
      LISGO2 = 163
      LISGO3 = 164
      LISGO4 = 165

      LISGV1 = 166
      LISGV2 = 167
      LISGV3 = 168
      LISGV4 = 169
      ENDIF
C
      ENDIF
C
C     If this is a gradient calculation, set list numbers for T2 or
C     D2T2 increments.
C
      IF(IDRLVL.GT.0)THEN
C
C     ---- MBPT(4), UCC(4) ----
C
      IF(ICLLVL.EQ.4.OR.ICLLVL.EQ.9)THEN
      LIST2I1 = 61
      LIST2I2 = 62
      LIST2I3 = 63
      ENDIF
C
C     ---- CCSD + T(CCSD), CCSD(T), QCISD(T) ----
C
      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.22.OR.ICLLVL.EQ.21)THEN
C      
C These list  have been used elsewhere.
CSSS      LIST2I1 = 101
CSSS      LIST2I2 = 102
CSSS      LIST2I3 = 103

      LIST2I1 = 1
      LIST2I2 = 2
      LIST2I3 = 3
      ENDIF
C
      ENDIF
C
C     If this is a fourth-order or UCC(4) GRADIENT calculation, convert
C     the T2 quantities on lists LIST2I1, LIST2I2, LIST2I3 to D2T2
C     quantities.
C
      IF((ICLLVL.EQ.4.OR.ICLLVL.EQ.9).AND.IDRLVL.GT.0)THEN
        CALL DEWQ(ICORE,MAXCOR,IUHF,0)
      ENDIF
C
C     Initialize energies.
C
      E4T    = 0.0D+00
      E4TAAA = 0.0D+00
      E4TBBB = 0.0D+00
      E4TBAA = 0.0D+00
      E4TBBA = 0.0D+00
      E4TABB = 0.0D+00
      E5ST   = 0.0D+00
      E5STA  = 0.0D+00
      E5STB  = 0.0D+00
C
C     If this is an iterative triples calculation, set list numbers
C     for T2 increments. Also, if this is supposedly an ROHF-MBPT(4)
C     calculation set d2t2 increment list numbers.
C
      IF(IDRLVL.EQ.0.AND.(TRIPIT.OR.ICLLVL.EQ.4))THEN
        LIST2I1 = 61
        LIST2I2 = 62
        LIST2I3 = 63
      ENDIF
C
C     If this is CCSDT or CCSDT-4 reorder MBEJ lists.
C
      IF(CCSDT4.OR.CCSDT) CALL SETVOVO(ICORE,MAXCOR,IUHF,1)
C
C     For CCSDT-4 replace intermediates with integrals on lists 54-59.
C
      IF(CCSDT4)          CALL   CSDT4(ICORE,MAXCOR,IUHF)
C
C     For CCSDT negate certain MBEJ lists.
C
      IF(CCSDT)           CALL DMPMBEJ(ICORE,MAXCOR,IUHF)
C
C     If this is an RHF iterative triples energy calculation skip the
C     all alpha case.
C
      IF(IUHF.EQ.0 .OR. ISKP14) GOTO 1000
C
      DO  110 ISPIN=1,2
      NO = NOCCO(ISPIN)
      NV = NVRTO(ISPIN)
C
      IF(NO.LT.3.OR.NV.LT.3) GOTO 110
C
C     I000   T3   AMPLITUDES
C     I010   OOOV INTEGRALS
C     I020   OOOV intermediates
C     I030   BUF  BUFFER FOR OOOV INTEGRAL LIST READ
C     I040   EVAL EIGENVALUES
C     I050   IADBLK
C     I060   LENBLK
C     I070   IADT2
C     I080   LENT2
C     I090   IADV
C     I100   LENV
C     I110   OOOVAD
C     I120   LENINT
C     I130   ICORE (THE REST)
C
      I000 = 1
      I010 = I000
      I020 = I010 + IINTFP * LNOOOV
      I025 = I020 + IINTFP * LNOOOV
      I030 = I025 + IINTFP * LNOOOV
      I040 = I030 + IINTFP * NO * (NO-1)/2
      I050 = I040 + IINTFP * (NO + NV)
      I060 = I050 + IINTFP * 8
      I070 = I060 + IINTFP * 8
      I080 = I070 + IINTFP * 8
      I090 = I080 + IINTFP * 8
      I100 = I090 + IINTFP * 8
      I110 = I100 + IINTFP * 8
      I120 = I110 + IINTFP * 8
      I130 = I120 + IINTFP * 8
C     COMPUTE NUMBER OF FREE INTEGER WORDS
      NFREE = MAXCOR - (I130 - 1)
      NALLOC = I130 - 1
cYAU      WRITE(6,2001) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
        WRITE(6,2002)
        CALL INSMEM('TRPS',I130,MAXCOR)
      ENDIF
      call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
      CALL TRPS14_PCCD(ICORE(I010),ICORE(I020),ICORE(I025),ICORE(I030),
     1            ICORE(I040),
     1            ICORE(I050),ICORE(I060),
     1            ICORE(I070),ICORE(I080),
     1            ICORE(I090),ICORE(I100),
     1            ICORE(I110),ICORE(I120),LNOOOV,
     1            ICORE(I130),
     1            NFREE,IUHF,ISPIN,NO,NV,INT1,
     1            INT2,NONHF,.FALSE.,FINAL,EXTRNL_CCSD,PCCD)
      call wallclock(out_year,out_mon,out_mday,out_hour,out_min,out_sec)
CSSS      write(*,*) 'trps14 in  : ',
CSSS     &           in_year,in_mon,in_mday,in_hour,in_min,in_sec
CSSS      write(*,*) 'trps14 out : ',
CSSS     &           out_year,out_mon,out_mday,out_hour,out_min,out_sec
  110 CONTINUE
 1000 CONTINUE
C
      IF(NOCA .LT.2.OR.NOCB .LT.1.OR.
     1   NVRTA.LT.2.OR.NVRTB.LT.1.OR.ISKP2) GOTO 2000
C
C     Transpose certain OOOV lists.
C
CSSS      Write(*,*) "The derivative level", IDRLVL 
      IF(IDRLVL.EQ.0)THEN
       IF(IUHF.EQ.0)THEN
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
        IF(INT2.AND.LWIC24.NE.LWIC14)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
        ENDIF
       ELSE
        CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
        IF(INT2.AND.LWIC24.NE.LWIC14)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
        ENDIF
       ENDIF
      ELSE
CSSS       Write(*,*) "Entering setooov1"
       CALL SETOOOV(ICORE,MAXCOR,IUHF,1)
      ENDIF
C
C     I000   T1   alpha
C     I010   T1   beta
C     I020   FOV  alpha
C     I030   FOV  beta
C     I040   D1T1 alpha
C     I050   D1T1 beta
C     I060   S1   alpha
C     I070   S1   beta
C     I080   EVAL
C     I090   SCR1
C     I100   SCR2
C     I110   SCR3
C     I120   OOOV integrals
c     I130   OOOV intermediates
C     I140   OOOV gammas
C     I150   BUF  BUFFER FOR OOOV INTEGRAL LIST READ
C     I160   IADT3
C     I170   LENT3
C     I180   IADT2
C     I190   LENT2
C     I200   IADV
C     I210   LENV
C     I220   OOOVAD
C     I230   LENINT
C     I240   ICORE (THE REST)
C
      I000 = 1
      I010 = I000 + IINTFP * NTAA
      I020 = I010 + IINTFP * NTBB
      I030 = I020 + IINTFP * NTAA
      I040 = I030 + IINTFP * NTBB
      I050 = I040 + IINTFP * NTAA
      I060 = I050 + IINTFP * NTBB
      I070 = I060 + IINTFP * NTAA
      I080 = I070 + IINTFP * NTBB
      I090 = I080 + IINTFP *(NOCA + NVRTA) * 2
      I100 = I090 + IINTFP * SCRDIM
      I110 = I100 + IINTFP * SCRDIM
      I120 = I110 + IINTFP * SCRDIM
C     I130 = I120 + IINTFP * LNOOOV * 4
      I130 = I120
C     I140 = I130 + IINTFP * LNOOOV * 4
      I140 = I130
C
      IF(IUHF.EQ.0)THEN
      NGAMMA = 1
      ELSE
      NGAMMA = 4
      ENDIF
      IF(IDRLVL.GT.0 .OR. EXTRNL_CCSD)THEN
      I150 = I140 + IINTFP * LNOOOV * NGAMMA
      ELSE
      I150 = I140 + IINTFP * 1      * NGAMMA
      ENDIF
      I160 = I150 + IINTFP * MAXOO
      I161 = I160 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I162 = I161 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I163 = I162 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I164 = I163 + IINTFP * MAX1  * MAX(NOCA,NOCB)
      I165 = I164 + IINTFP * MAX1  * MAX(NOCA,NOCB)
      I170 = I165 + IINTFP * 8
      I180 = I170 + IINTFP * 8
      I190 = I180 + IINTFP * 8
      I200 = I190 + IINTFP * 8
      I210 = I200 + IINTFP * 8
      I220 = I210 + IINTFP * 8
      I230 = I220 + IINTFP * 8             * 4
      I240 = I230 + IINTFP * 8             * 4
C
C     Compute number of free integer words.
C
      NFREE = MAXCOR - (I240 - 1)
      NALLOC = I240 - 1
cYAU      WRITE(6,2011) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
        WRITE(6,2012)
        CALL INSMEM('TRPS',I240,MAXCOR)
      ENDIF
      call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
      CALL TRPS2_PCCD(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     1           ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     1           ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     1                                   ICORE(I140),ICORE(I150),
     1           ICORE(I160),ICORE(I161),ICORE(I162),
     1           ICORE(I163),ICORE(I164),
     1           ICORE(I165),ICORE(I170),ICORE(I180),ICORE(I190),
     1           ICORE(I200),ICORE(I210),ICORE(I220),ICORE(I230),
     1           ICORE(I240),
     1           LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,INT1,INT2,NONHF,
     1           IUHF,.FALSE.,FINAL,EXTRNL_CCSD,PCCD)
      call wallclock(out_year,out_mon,out_mday,out_hour,out_min,out_sec)

      IF(IDRLVL.EQ.0)THEN
       IF(IUHF.EQ.0)THEN
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
        IF(INT2.AND.LWIC24.NE.LWIC14)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
        ENDIF
       ELSE
        CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
        IF(INT2.AND.LWIC24.NE.LWIC14)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
        ENDIF
       ENDIF
      ELSE
       Print*, "Entering setoov2"
       CALL SETOOOV(ICORE,MAXCOR,IUHF,2)
      ENDIF
 2000 CONTINUE
C
      IF(IUHF.EQ.0.OR.
     1   NOCA .LT.1.OR.NOCB .LT.2.OR.
     1   NVRTA.LT.1.OR.NVRTB.LT.2.OR.ISKP3) GOTO 3000
C
C     Reorganize abci/ijka intermediate lists, t2 ab lists and,
C     if appropriate, d2t2 ab increment lists.
C
      IF(IUHF.EQ.1)THEN
        CALL TSPABCI2(ICORE,MAXCOR,1,EXTRNL_CCSD)
C
cold    CALL AT2IAB(ICORE,0,INT2)
        CALL AT2IAB(ICORE,46,1,0)
        IF(INT1) CALL AT2IAB(ICORE,16,1,0)
        IF(INT2) CALL AT2IAB(ICORE,LIST2I3,1,0)
        IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))
     &           CALL AT2IAB(ICORE,116,1,0)
        IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE,MAXCOR,IUHF,1)
      ENDIF
C
      I000 = 1
      I010 = I000 + IINTFP * NTAA
      I020 = I010 + IINTFP * NTBB
      I030 = I020 + IINTFP * NTAA
      I040 = I030 + IINTFP * NTBB
      I050 = I040 + IINTFP * NTAA
      I060 = I050 + IINTFP * NTBB
      I070 = I060 + IINTFP * NTAA
      I080 = I070 + IINTFP * NTBB
      I090 = I080 + IINTFP *(NOCA + NVRTA) * 2
      I100 = I090 + IINTFP * SCRDIM
      I110 = I100 + IINTFP * SCRDIM
      I120 = I110 + IINTFP * SCRDIM
      I130 = I120 + IINTFP * LNOOOV * 4
      I140 = I130 + IINTFP * LNOOOV * 4
      I150 = I140 + IINTFP * LNOOOV * 4
      I160 = I150 + IINTFP * MAXOO
      I161 = I160 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I162 = I161 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I163 = I162 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I164 = I163 + IINTFP * MAX2  * MAX(NOCA,NOCB)
      I165 = I164 + IINTFP * MAX2  * MAX(NOCA,NOCB)
      I170 = I165 + IINTFP * 8
      I180 = I170 + IINTFP * 8
      I190 = I180 + IINTFP * 8
      I200 = I190 + IINTFP * 8
      I210 = I200 + IINTFP * 8
      I220 = I210 + IINTFP * 8
      I230 = I220 + IINTFP * 8             * 4
      I240 = I230 + IINTFP * 8             * 4
C
C     Compute number of free integer words.
C
      NFREE = MAXCOR - (I240 - 1)
      NALLOC = I240 - 1
cYAU      WRITE(6,2021) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
        WRITE(6,2022)
        CALL INSMEM('TRPS',I240,MAXCOR)
      ENDIF
      call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
      CALL TRPS3_PCCD(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     1           ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     1           ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     1           ICORE(I120),ICORE(I130),ICORE(I140),ICORE(I150),
     1           ICORE(I160),ICORE(I161),ICORE(I162),
     1           ICORE(I163),ICORE(I164),
     1           ICORE(I165),ICORE(I170),ICORE(I180),ICORE(I190),
     1           ICORE(I200),ICORE(I210),ICORE(I220),ICORE(I230),
     1           ICORE(I240),
     1           LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,INT1,INT2,NONHF,
     1           IUHF,.FALSE.,FINAL,EXTRNL_CCSD,PCCD)
      call wallclock(out_year,out_mon,out_mday,out_hour,out_min,out_sec)
CSSS      write(*,*) 'trps3 in  : ',
CSSS     &           in_year,in_mon,in_mday,in_hour,in_min,in_sec
CSSS      write(*,*) 'trps3 out : ',
CSSS     &           out_year,out_mon,out_mday,out_hour,out_min,out_sec
C
C     (Beta,Aplha) --> (Alpha,Beta).
C
      IF(IUHF.EQ.1)THEN
        CALL TSPABCI2(ICORE,MAXCOR,2,EXTRNL_CCSD)
C
cold    CALL AT2IAB(ICORE,1,INT2)
        CALL AT2IAB(ICORE,46,1,1)
        IF(INT1) CALL AT2IAB(ICORE,16,1,1)
        IF(INT2) CALL AT2IAB(ICORE,LIST2I3,1,1)
        IF(IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))
     &           CALL AT2IAB(ICORE,116,1,1)
        IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE,MAXCOR,IUHF,2)
      ENDIF
C
 3000 CONTINUE
C
C     <Ab||cI> --> <Ab||Ic>.
C
      IF(IUHF.EQ.1)THEN
         CALL TSPABCI(ICORE,MAXCOR,2,EXTRNL_CCSD)
      ENDIF
C
C     For CCSDT negate certain MBEJ lists.
C
      IF(CCSDT)           CALL DMPMBEJ(ICORE,MAXCOR,IUHF)
C
C     If this is CCSDT or CCSDT-4 reorder MBEJ lists.
C
      IF(CCSDT4.OR.CCSDT) CALL SETVOVO(ICORE,MAXCOR,IUHF,2)
C
C     Attempt to make sure plain integrals end up on 231-233.
C
      IF(CCSDT)THEN
        CALL HBARCRAP(ICORE,MAXCOR,IUHF,.FALSE.)
        CALL   FORMW1(ICORE,MAXCOR,IUHF,.FALSE.)
      ENDIF
C
C     Attempt to compute the RMS error in T3.
C
c     CALL RMST3(ICORE,MAXCOR,IUHF)
C
C     Move the triples amplitudes so that the new (ie current) ones will
C     be the old ones in the next iteration.
C
      IF(T3STOR)THEN
        IF(TRIPIT)THEN
          CALL MOVLIST(ICORE,1,5)
          CALL MOVLIST(ICORE,2,6)
          IF(IUHF.GT.0)THEN
          CALL MOVLIST(ICORE,3,7)
          CALL MOVLIST(ICORE,4,8)
          ENDIF
        ENDIF
        IF(INCORE) CALL PUTALLT3
      ENDIF
C
C     Reset derivative level flag.
C
      IF(IFLAGS2(h_IFLAGS2_grad_calc).EQ.2)THEN
        if (iFlags2(h_IFLAGS2_geom_opt).ne.0) then
           IDRLVL = 1
           WRITE(6,2052)
        end if
        if (iFlags(h_IFLAGS_vib)       .eq.3) then
           IDRLVL = IDROLD
           WRITE(6,2053) IDRLVL
        end if
      ENDIF
 2052 FORMAT(' @TRPS-I, Numerical optimization. IDRLVL reset to 1. ')
 2053 FORMAT(' @TRPS-I, Energy only frequency. IDRLVL reset to ',I1)
C
C     Return in iterativetriples energy calculations.
C
      IF(TRIPIT.AND.ICLLVL.NE.9)        RETURN
      IF(ICLLVL.EQ. 9.AND.IDRLVL.EQ. 0) RETURN
C
C     If this is a gradient calculation, check sum final t3 gammas
C     and evaluate fourth-order energy if this is an MBPT(4) calcul-
C     ation.
C
      IF(IDRLVL.GT.0.AND.T3GRAD)THEN
        CALL GT3E4T(ICORE,MAXCOR,IUHF)
      ENDIF
C
C    If this is an MBPT(4) or UCC(4) GRADIENT calculation, convert the D2T2
C    quantities on lists LIST2I1, LIST2I2, LIST2I3 to T2 quantities.
C
      IF((ICLLVL.EQ.4.OR.ICLLVL.EQ.9).AND.IDRLVL.GT.0.AND.T3GRAD)THEN
        CALL DEWQ(ICORE,MAXCOR,IUHF,1)
      ENDIF
C
C     If this is a gradient calculation, dump diagonal parts of density
C     to lists.
C
      IF(IDRLVL.GT.0.AND.T3GRAD)THEN
        CALL DMPDEN(DOOA,DOOB,DVVA,DVVB,ICORE,
     &              NOCCO(1),NOCCO(2),NVRTO(1),NVRTO(2),IUHF,IINTFP)
      ENDIF
C
C     If method is CCSD+T(CCSD), CCSD(T), QCISD(T) then evaluate
C     CCSD/QCISD energy. Note F*T1 contribution must be included
C     if reference is non-Hartree-Fock.
C
      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.21.OR.ICLLVL.EQ.22.OR.
     &   ICLLVL.EQ.50.OR.ICLLVL.EQ.51)THEN
        CALL CCENRG(ICORE,IUHF,ECCTOT,NONHF)
      ENDIF
C
C     Summarize results.
C
      IF(IUHF.EQ.0)THEN
        E4T  = 2.0D+00 * (E4TAAA + E4TBAA)
        E5ST = 2.0D+00 * (E5STA  + E5STB )
        E4ST = 2.0D+00 * (E4STA  + E4STB )
      ELSE
        E4T  =            E4TAAA + E4TBBB + E4TBAA + E4TABB
        E5ST =            E5STA  + E5STB
        E4ST =            E4STA  + E4STB
      ENDIF
C
C     --- Hartree-Fock MBPT(4) ---
C
      IF(ICLLVL.EQ. 4.AND.HF)THEN
      WRITE(6,2031) E4T
 2031 FORMAT(' @TRPS-I, Hartree-Fock MBPT(4) triples energy = ',
     1       F15.12)
      ENDIF
C
C     --- ROHF, QRHF, or general non-Hartree-Fock MBPT(4) ---
C
      IF(ICLLVL.EQ. 4.AND.NONHF)THEN
      WRITE(6,2032) E4T,E5ST,E4ST,E4T,E5ST,E4ST
 2032 FORMAT(' @TRPS-I, Non-Hartree-Fock MBPT(4) triples energies',/,
     1       '          E4T           = ',F15.12,/,
     1       '          E5ST          = ',F15.12,/,
     1       '          E4ST          = ',F15.12,/,
     1       '          In other words  ',/,
     1       '          <T3[2](Dagger) D3  T3[2]> = ',F15.12,/,
     1       '          <T1   (Dagger) Q1 WT3[2]> = ',F15.12,/,
     1       '          <T2   (Dagger) Q2 FT3[2]> = ',F15.12)
      E4T = E4T + E5ST + E4ST
      ENDIF
C
      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.22.OR.ICLLVL.EQ.50.OR.ICLLVL.EQ.51)
     1THEN
      WRITE(6,2033) E4T,E5ST,E4ST,E4T,E5ST,E4ST
 2033 FORMAT(' @TRPS-I, E4T (T(CCSD)) = ',F15.12,/,
     1       '          E5ST          = ',F15.12,/,
     1       '          E4ST          = ',F15.12,/,
     1       '          In other words  ',/,
     1       '          <T3[2](Dagger) D3  T3[2]> = ',F15.12,/,
     1       '          <T1   (Dagger) Q1 WT3[2]> = ',F15.12,/,
     1       '          <T2   (Dagger) Q2 FT3[2]> = ',F15.12)
      ENDIF
C
C     --- CCSD+T(CCSD) with HF reference (not defined with Non-HF) ---
C
      IF(ICLLVL.EQ.11)THEN
      ETCCSD = ECCTOT + E4T
      WRITE(6,2034) ECCTOT,ETCCSD
 2034 FORMAT(' @TRPS-I, CCSD           = ',F20.12,/,
     1       '          CCSD + T(CCSD) = ',F20.12)
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,ETCCSD)
      ENDIF
C
C     --- CCSD(T) with HF reference ---
C
      IF((ICLLVL.EQ.22.OR.ICLLVL.EQ.50.OR.ICLLVL.EQ.51).AND.
     1    .NOT.NONHF)THEN
      ETCCSD = ECCTOT + E4T
      EPAREN = ECCTOT + E4T + E5ST
      WRITE(6,2035) ECCTOT,ETCCSD,EPAREN
 2035 FORMAT(' @TRPS-I, CCSD           = ',F20.12,/,
     1       '          CCSD + T(CCSD) = ',F20.12,/,
     1       '          CCSD(T)        = ',F20.12)
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
      CALL PUTREC(1,'JOBARC','CCSDENRG',IINTFP,ECCTOT)
      ENDIF
C
C     --- CCSD(T) with Non-HF reference ---
C
      IF((ICLLVL.EQ.22.OR.ICLLVL.EQ.50.OR.ICLLVL.EQ.51).AND.
     1    NONHF)THEN
      EPAREN = ECCTOT + E4T + E5ST + E4ST
      WRITE(6,2036) ECCTOT,EPAREN
 2036 FORMAT(' @TRPS-I, CCSD           = ',F20.12,/,
     1       '          CCSD(T)        = ',F20.12)
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
      ENDIF
C
C     --- QCISD(T)  with HF reference (not defined with Non-HF) ---
C
      IF(ICLLVL.EQ.21.OR.ICLLVL.EQ.50.OR.ICLLVL.EQ.51)THEN
      ETCCSD = ECCTOT + E4T
      EPAREN = ECCTOT + E4T + E5ST
      WRITE(6,2037) ECCTOT,EPAREN
 2037 FORMAT(' @TRPS-I, QCISD          = ',F20.12,/,
     1       '          QCISD(T)       = ',F20.12)
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
      ENDIF
C
      IF(NONHF.AND.IDRLVL.GT.0.AND.T3GRAD) WRITE(6,2039) E5TEST
C
C     Reset ICLLVL for special noniterative HF triples.
C
      IF(METACT.EQ.12.OR.(METACT.GE.26.AND.METACT.LE.31))THEN
        ICLLVL = METACT
      ENDIF
C
      E45T(1) = E4T
      E45T(2) = E5ST
      E45T(3) = E4ST
      CALL PUTREC(20,'JOBARC','NONITT3',IINTFP*20,E45T)

      IF (EXTRNL_CCSD) THEN
#ifdef _DEBUG_LVL0
          Write(6,*) 
          Write(6,"(a)") " ------Adding external CC flags---"
          Write(6,"(a)") " Iflags(h_IOPPAR_calc)  is set to 10->CCSD"
          Write(6,"(a)") " Iflags(h_IOPPAR_nonhf) is set to True" 
          Write(6,"(a)") " EXTRNLCC JOBARC record is added"
CSSS          call check_ints(icore,maxcor,iuhf,.false.)
          call check_t32w(icore(istart),nleft,iuhf)
          call check_t2(icore,maxcor,iuhf)
#endif 
          CALL PUTREC(20,"JOBARC","EXTRNLCC",1,1) 
          IF (IFLAGS(h_IOPPAR_calc) .EQ. 50) IFLAGS(h_IOPPAR_calc) = 51
CSSS          IFLAGS(h_IOPPAR_calc) = 10
CSSS          IF (IFLAGS(h_IOPPAR_nonhf) .EQ. 0) IFLAGS(h_IOPPAR_nonhf) = 1
          CALL PUTREC(20,"JOBARC","IFLAGS  ",100,IFLAGS)
      ENDIF 
C
      RETURN
 1061 FORMAT(' @TRPS-I, Reference function information : ')
 1062 FORMAT('          Reference function is general non-HF ')
 1063 FORMAT('          Reference function is Hartree-Fock ')
 1064 FORMAT('          Reference function is RHF Hartree-Fock ')
 1065 FORMAT('          Reference function is RHF non-HF ')
 1066 FORMAT('          Reference function is UHF Hartree-Fock ')
 1067 FORMAT(' @TRPS-I, Reference function is UHF non-HF ')
 1068 FORMAT(' @TRPS-I, Reference function is ROHF ')
 1069 FORMAT(' @TRPS-I, Reference function is QRHF ')
 1072 FORMAT(' @TRPS-I, INT1 is true ')
 1073 FORMAT(' @TRPS-I, INT1 is false ')
 1074 FORMAT(' @TRPS-I, INT2 is true ')
 1075 FORMAT(' @TRPS-I, INT2 is false ')
 1081 FORMAT(' @TRPS-I, Lengths of ijka integrals. #1 ',I10,/,
     1       '                                     #2 ',I10,/,
     1       '                                     #3 ',I10,/,
     1       '                                     #4 ',I10,/,
     1       '                                    Max ',I10)
 2001 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS14 ',I12,/,
     1 '          Number of free words for TRPS14 ',I12)
 2002 FORMAT(' TRPS-I, Insufficient memory for TRPS14 pre-allocation. ')
 2011 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS2 ',I12,/,
     1 '          Number of free words for TRPS2 ',I12)
 2012 FORMAT(' TRPS-I, Insufficient memory for TRPS2  pre-allocation. ')
 2021 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS3 ',I12,/,
     1 '          Number of free words for TRPS3 ',I12)
 2022 FORMAT(' TRPS-I, Insufficient memory for TRPS3  pre-allocation. ')
 2039 FORMAT(' @TRPS-I, E5TEST is ',F20.12)
      END
