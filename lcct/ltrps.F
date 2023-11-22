#include "flags.h"

      SUBROUTINE LTRPS(ICORE, MAXCOR, IUHF, E4T, FAKE_CALC_LEVEL)

      IMPLICIT INTEGER (A-Z)

      DOUBLE PRECISION E4T,E4TAAA,E4TBBB,E4TBAA,E4TABB,E5ST,
     1                 E5STA,E5STB,E4ST,E4STA,E4STB,ECCTOT,
     1                 ETCCSD,EPAREN
C
      DOUBLE PRECISION E5TEST
      DOUBLE PRECISION DOOA,DOOB,DVVA,DVVB
      DOUBLE PRECISION ALLT3
      DOUBLE PRECISION E45T(20)
C
      LOGICAL INT1,INT2,NONHF,HF
      LOGICAL ISKP14,ISKP2,ISKP3
      LOGICAL CCSDT4,CCSDT
      LOGICAL INCORE
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      LOGICAL T3GRAD
C
      CHARACTER*7  FAKE_CALC_LEVEL
C
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
     1                E4STA ,E4STB
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
      COMMON /TESTEN/ E5TEST
      COMMON /T3METH/ CCSDT4,CCSDT
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3CORE/ INCORE
C-----------------------------------------------------------------------
C Let's fake this to think that it is doing a CCSD(T) calculation. Lambda
C based CCCSD(T) is a special kind of CCSD(T) calculation. So it is fair
C and easy implement L-CCSD(T) if we fake the code to think that it is
C doing CCSD(T) and add the proper changes.
CC
C If this is a gradient calculation, set gamma list numbers
C
      IF (IDRLVL.GT.0)THEN
C
         LISGO1  = 162
         LISGO2  = 163
         LISGO3  = 164
         LISGO4  = 165
         LISGV1  = 166
         LISGV2  = 167
         LISGV3  = 168
         LISGV4  = 169
         LIST2I1 = 1
         LIST2I2 = 2
         LIST2I3 = 3
C
      ENDIF
C     List offset for Lambda
      LISTOFF=200
      IF (FAKE_CALC_LEVEL .EQ. "CCSD(T)") THEN
         ICLLVL = 22
         iflags(h_IOPPAR_calc)=22
         TRIPNI=.TRUE.
         TRIPNI1=.TRUE.
         TRIPIT=.FALSE.
         T3STOR=.FALSE.

         if (idrlvl .gt. 0) call init2l(iuhf)
      ENDIF
C
      WRITE(6,1005)
 1005 FORMAT(/,72('-'),/,T17,
     &       ' CC Triple Excitation Contributions ',/,72('-'),/)
C
      CALL ZERO(E45T,20)
      E4T  = 0.0D+00
      E5ST = 0.0D+00
      E4ST = 0.0D+00
C
      CCSDT4 = .FALSE.
      CCSDT  = .FALSE.
C
      CALL DIMT3(DISTSZ(1,1),NDISTS(1,1))
      CALL DIMT3(DISTSZ(1,5),NDISTS(1,5))
C
C-----------------------------------------------------------------------
C If triples are stored and there is enough core, put in ALLT3.
C-----------------------------------------------------------------------
      IF(T3STOR)THEN
         NT3 = NRECS * LNPHYR
         IF(NT3.LT.NUMT3)THEN
            INCORE = .TRUE.
            WRITE(6,9010)
         ELSE
            INCORE = .FALSE.
            WRITE(6,9020)
         ENDIF
         IF(INCORE) CALL GETALLT3
      ENDIF
 9010 FORMAT(' @LTRPS-I, There is sufficient core for the T3 file. ')
 9020 FORMAT(' @LTRPS-I, There is insufficient core for the T3 file. ')
C-----------------------------------------------------------------------
C
      E5TEST = 0.0D+00
C
C-----------------------------------------------------------------------
C     Debug aids to skip one or more spin cases.
C-----------------------------------------------------------------------
      ISKP14 = .FALSE.
      ISKP2  = .FALSE.
      ISKP3  = .FALSE.

      T3GRAD = .FALSE.
      IF(IDRLVL.GT.0)THEN
C
         T3GRAD = .TRUE.
C
C     If this is a numerical gradient and IDRLVL > 0, turn the
C     derivative flag off, but turn on at end.
C
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
     &   ' orbitals ! ',/,
     &   '          Please resubmit with PERT_ORB=CANONICAL ')
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
C Do not attempt noniterative triple excitation calculations
C for non-HF orbitals unless we have semicanonical orbitals.
C
      IF(TRIPNI.AND.NONHF.AND.IORBTP.EQ.0)THEN
         WRITE(6,1070)
         CALL ERREX
      ENDIF
 1070 FORMAT(' @LTRPS-I, Caution : This calculation',
     1 ' is about to evaluate noniterative',/,
     1 '          triples corrections with',
     1 ' standard non-HF orbitals. ACES2 devel-',/,
     1 '          opers do not believe the',
     1 ' results are meaningful. We suggest',/,
     1       '          Please resubmit with ORBITAL=SEMICANONICAL ')
C
C-----------------------------------------------------------------------
C Set flags for including T3 in T2 and T1 equations. Remember that
C UCC(4) does not put T3 in the T1 equation.
C
C Note : for CCSD(T) and QCISD(T), INT1 means "compute D1T1 increm-
C ents but only increment lists if this is a gradient calculation".
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

      DO 20 IRREP=1,NIRREP
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
 20   CONTINUE
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
C Evaluate symmetry packed lengths of <ij//ka> integrals/intermediates.
C
      LNOOOV1 = 0
      LNOOOV2 = 0
      LNOOOV3 = 0
      LNOOOV4 = 0
      DO 30 IRPIJ=1,NIRREP
         LNOOOV1 = LNOOOV1 + IRPDPD(IRPIJ, 3) * IRPDPD(IRPIJ, 9)
         LNOOOV2 = LNOOOV2 + IRPDPD(IRPIJ, 4) * IRPDPD(IRPIJ,10)
         LNOOOV3 = LNOOOV3 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,11)
         LNOOOV4 = LNOOOV4 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,12)
 30   CONTINUE
C
      LNOOOV = MAX0(LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4)
cYAU      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
C
C Set up orbital index tables
C
      CALL ABSORB(ABSVRT,ABSOCC,MAXNBF)
C
C Calculate VV, OO, VO offsets
C
      CALL MKOFVV
      CALL MKOFOO
      CALL MKOFVO
C
C Set W intermediate list numbers
C
      Print*, "Entering initrp"
      CALL INITRP
C
C     Transpose records of LWIC17 and LWIC27 (the latter only if it is
C     not the same list).
C
      IF(IUHF.EQ.1)THEN
        CALL TSPABCI(ICORE,MAXCOR,1)
      ENDIF

C
C If this is a gradient calculation, initialize diagonal elements
C of density matrix
C
      IF(IDRLVL.GT.0)THEN
         CALL ZERO(DOOA,MAXNBF)
         CALL ZERO(DOOB,MAXNBF)
         CALL ZERO(DVVA,MAXNBF)
         CALL ZERO(DVVB,MAXNBF)
      ENDIF


C
C Initialize energies.
C
      E4T    = 0.0D+00
      E4TAAA = 0.0D+00
      E4TBBB = 0.0D+00
      E4TBAA = 0.0D+00
      E4TABB = 0.0D+00
      E5ST   = 0.0D+00
      E5STA  = 0.0D+00
      E5STB  = 0.0D+00
C
      IF(IUHF.EQ.0 .OR. ISKP14) GOTO 1000
C
      DO 110 ISPIN=1, 2
         NO = NOCCO(ISPIN)
         NV = NVRTO(ISPIN)
C
         IF (NO.GE.3.OR.NV.GE.3) THEN
C
C I000   T3   AMPLITUDES
C I010   OOOV INTEGRALS
C I020   OOOV intermediates
C I030   BUF  BUFFER FOR OOOV INTEGRAL LIST READ
C I040   EVAL EIGENVALUES
C I050   IADBLK
C I060   LENBLK
C I070   IADT2
C I080   LENT2
C I090   IADV
C I100   LENV
C I110   OOOVAD
C I120   LENINT
C I130   ICORE (THE REST)
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
C COMPUTE NUMBER OF FREE INTEGER WORDS
            NFREE  = MAXCOR - (I130 - 1)
            NALLOC = I130 - 1
cYAU      WRITE(6,2001) MAXCOR,NALLOC,NFREE
            IF (NFREE.LE.0)THEN
               WRITE(6,2002)
               CALL INSMEM('LTRPS',I130,MAXCOR)
            ENDIF
      call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
      CALL TRPS14(ICORE(I010),ICORE(I020),ICORE(I025),ICORE(I030),
     1            ICORE(I040),
     1         ICORE(I050),ICORE(I060),
     1         ICORE(I070),ICORE(I080),
     1         ICORE(I090),ICORE(I100),
     1         ICORE(I110),ICORE(I120),LNOOOV,
     1         ICORE(I130),
     1         NFREE,IUHF,ISPIN,NO,NV,INT1,INT2,NONHF,.TRUE.,.FALSE.,
     1         .FALSE.)
C
            call wallclock(out_year,out_mon,out_mday,out_hour,
     &         out_min,out_sec)
            write(*,*) 'trps14 in  : ',
     &         in_year,in_mon,in_mday,in_hour,in_min,in_sec
            write(*,*) 'trps14 out : ',
     &           out_year,out_mon,out_mday,out_hour,out_min,out_sec
        ENDIF
 110  CONTINUE
C
C Entry point RHF calcualtions
C
 1000 CONTINUE

      IF (.NOT. (NOCA .LT.2.OR.NOCB .LT.1.OR.
     1   NVRTA.LT.2.OR.NVRTB.LT.1.OR.ISKP2)) THEN
C
C Transpose certain OOOV lists.
C
         IF (IDRLVL.EQ.0)THEN
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
               IF (INT2.AND.LWIC24.NE.LWIC14)THEN
                  CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
               ENDIF
            ENDIF
         ELSE
            CALL SETOOOV(ICORE,MAXCOR,IUHF,1)
         ENDIF
C
C I000   T1   alpha
C I010   T1   beta
C I020   FOV  alpha
C I030   FOV  beta
C I040   D1T1 alpha
C I050   D1T1 beta
C I060   S1   alpha
C I070   S1   beta
C I080   EVAL
C I090   SCR1
C I100   SCR2
C I110   SCR3
C I120   OOOV integrals
C I130   OOOV intermediates
C I140   OOOV gammas
C I150   BUF  BUFFER FOR OOOV INTEGRAL LIST READ
C I160   IADT3
C I170   LENT3
C I180   IADT2
C I190   LENT2
C I200   IADV
C I210   LENV
C I220   OOOVAD
C I230   LENINT
C I240   ICORE (THE REST)
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
         IF (IUHF.EQ.0)THEN
            NGAMMA = 1
         ELSE
            NGAMMA = 4
         ENDIF
C
         IF (IDRLVL.GT.0)THEN
            I150 = I140 + IINTFP * LNOOOV * NGAMMA
         ELSE
            I150 = I140 + IINTFP * 1      * NGAMMA
         ENDIF
C
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
C
cYAU      WRITE(6,2011) MAXCOR,NALLOC,NFREE
C
         IF(NFREE.LE.0)THEN
            WRITE(6,2012)
            CALL INSMEM('LTRPS',I240,MAXCOR)
         ENDIF
C
         call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
C
         CALL TRPS2(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     1      ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     1      ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     1      ICORE(I140),ICORE(I150),
     1      ICORE(I160),ICORE(I161),ICORE(I162),
     1      ICORE(I163),ICORE(I164),
     1      ICORE(I165),ICORE(I170),ICORE(I180),ICORE(I190),
     1      ICORE(I200),ICORE(I210),ICORE(I220),ICORE(I230),
     1      ICORE(I240),
     1      LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,INT1,INT2,NONHF,
     1      IUHF,.TRUE.,.FALSE.,.FALSE.)
C
         call wallclock(out_year,out_mon,out_mday,out_hour,out_min,
     &      out_sec)
C
         write(*,*) 'trps2 in  : ',
     &      in_year,in_mon,in_mday,in_hour,in_min,in_sec
         write(*,*) 'trps2 out : ',
     &      out_year,out_mon,out_mday,out_hour,out_min,out_sec
C
         IF (IDRLVL.EQ.0)THEN
            IF (IUHF.EQ.0)THEN
               CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
               IF (INT2.AND.LWIC24.NE.LWIC14)THEN
                  CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
               ENDIF
            ELSE
               CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
               CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
               CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
               CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
               IF (INT2.AND.LWIC24.NE.LWIC14)THEN
                  CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
                  CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
               ENDIF
            ENDIF
C
         ELSE
            CALL SETOOOV(ICORE,MAXCOR,IUHF,2)
         ENDIF
C
c 2000    CONTINUE
      endif
C
         IF (IUHF.EQ.0.OR.
     1      NOCA .LT.1.OR.NOCB .LT.2.OR.
     1      NVRTA.LT.1.OR.NVRTB.LT.2.OR.ISKP3) GOTO 3000
C
C Reorganize abci/ijka intermediate lists, t2 ab lists and,
C if appropriate, d2t2 ab increment lists.
C
         IF (IUHF.EQ.1)THEN
            CALL TSPABCI2(ICORE,MAXCOR,1)
C
cold    CALL AT2IAB(ICORE,0,INT2)
            CALL AT2IAB(ICORE,46,1,0)
            CALL AT2IAB(ICORE,146,1,0)
            IF (INT1) CALL AT2IAB(ICORE,16,1,0)
            IF (INT2) THEN
              CALL AT2IAB(ICORE,LIST2I3,1,0)
              CALL AT2IAB(ICORE,LIST2I3+LISTOFF,1,0)
            ENDIF
            IF (IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))
     &         CALL AT2IAB(ICORE,116+59,1,0)
C
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
C
cYAU      WRITE(6,2021) MAXCOR,NALLOC,NFREE
C
         IF(NFREE.LE.0)THEN
            WRITE(6,2022)
            CALL INSMEM('LTRPS',I240,MAXCOR)
         ENDIF
C
         call wallclock(in_year,in_mon,in_mday,in_hour,in_min,in_sec)
C
         CALL TRPS3(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     1      ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     1      ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     1      ICORE(I120),ICORE(I130),ICORE(I140),ICORE(I150),
     1      ICORE(I160),ICORE(I161),ICORE(I162),
     1      ICORE(I163),ICORE(I164),
     1      ICORE(I165),ICORE(I170),ICORE(I180),ICORE(I190),
     1      ICORE(I200),ICORE(I210),ICORE(I220),ICORE(I230),
     1      ICORE(I240),
     1      LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,INT1,INT2,NONHF,
     1      IUHF,.TRUE.,.FALSE.,.FALSE.)
C
         call wallclock(out_year,out_mon,out_mday,out_hour,out_min,
     &      out_sec)
C
         write(*,*) 'trps3 in  : ',
     &      in_year,in_mon,in_mday,in_hour,in_min,in_sec
         write(*,*) 'trps3 out : ',
     &      out_year,out_mon,out_mday,out_hour,out_min,out_sec
C
C (Beta,Aplha) --> (Alpha,Beta).
C
         IF (IUHF.EQ.1)THEN
            CALL TSPABCI2(ICORE,MAXCOR,2)
C
cold    CALL AT2IAB(ICORE,1,INT2)
            CALL AT2IAB(ICORE,46,1,1)
            CALL AT2IAB(ICORE,146,1,1)
            IF(INT1) CALL AT2IAB(ICORE,16,1,1)
C     Something here for Lambda?
            IF (INT2) THEN
               CALL AT2IAB(ICORE,LIST2I3,1,1)
               CALL AT2IAB(ICORE,LIST2I3+LISTOFF,1,1)
            ENDIF
            IF (IDRLVL.GT.0.AND.(ICLLVL.EQ.21.OR.ICLLVL.EQ.22))
     &         CALL AT2IAB(ICORE,116+59,1,1)
C
            IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE,MAXCOR,IUHF,2)
         ENDIF
C
 3000    CONTINUE
C
C <Ab||cI> --> <Ab||Ic>.
C
         IF(IUHF.EQ.1)THEN
            CALL TSPABCI(ICORE,MAXCOR,2)
         ENDIF
C
C Attempt to compute the RMS error in T3.
C
c     CALL RMST3(ICORE,MAXCOR,IUHF)
C
C Move the triples amplitudes so that the new (ie current) ones will
C be the old ones in the next iteration.
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
 2052    FORMAT(' @LTRPS-I, Numerical optimization. IDRLVL reset to 1.')
 2053    FORMAT(' @LTRPS-I, Energy only frequency. IDRLVL reset to ',I1)
C
C If this is a gradient calculation, check sum final t3 gammas
C and evaluate fourth-order energy if this is an MBPT(4) calculation.
C
         IF(IDRLVL.GT.0.AND.T3GRAD)THEN
            CALL GT3E4T(ICORE,MAXCOR,IUHF)
         ENDIF
C
C If this is a gradient calculation, dump diagonal parts of density
C to lists.
C
         IF (IDRLVL.GT.0.AND.T3GRAD)THEN
            CALL DMPDEN(DOOA,DOOB,DVVA,DVVB,ICORE,
     &         NOCCO(1),NOCCO(2),NVRTO(1),NVRTO(2),IUHF,IINTFP)
         ENDIF
C
C If method is CCSD+T(CCSD), CCSD(T), QCISD(T) then evaluate
C CCSD/QCISD energy. Note F*T1 contribution must be included
C if reference is non-Hartree-Fock.
C
         IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
            CALL CCENRG(ICORE,IUHF,ECCTOT,NONHF)
         ENDIF
C
C Summarize results.
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
         IF (ICLLVL.EQ.11.OR.ICLLVL.EQ.22)THEN
            WRITE(6,2033) E4T,E5ST,E4ST,E4T,E5ST,E4ST
 2033       FORMAT(' @LTRPS-I, E4T (T(CCSD)) = ',F15.12,/,
     1             '           E5ST          = ',F15.12,/,
     1             '           E4ST          = ',F15.12,/,
     1             '           In other words  ',/,
     1             '           <L3[2] D3  T3[2]> = ',F15.12,/,
     1             '           <L1    Q1 WT3[2]> = ',F15.12,/,
     1             '           <L2    Q2 FT3[2]> = ',F15.12)
         ENDIF
C
C --- CCSD(T) with HF reference ---
C
         IF (ICLLVL.EQ.22.AND..NOT.NONHF)THEN
            ETCCSD = ECCTOT + E4T
            EPAREN = ECCTOT + E4T + E5ST
            WRITE(6,2035) ECCTOT,ETCCSD,EPAREN
C
 2035       FORMAT(' @LTRPS-I,  CCSD          = ',F20.12,/,
     1             '           aCCSD[T]       = ',F20.12,/,
     1             '           aCCSD(T)       = ',F20.12)
            CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
            CALL PUTREC(1,'JOBARC','CCSDENRG',IINTFP,ECCTOT)
         ENDIF
C
C --- CCSD(T) with Non-HF reference ---
C
         IF (ICLLVL.EQ.22.AND.NONHF)THEN
            EPAREN = ECCTOT + E4T + E5ST + E4ST
            WRITE(6,2036) ECCTOT,EPAREN
 2036       FORMAT(' @LTRPS-I,  CCSD          = ',F20.12,/,
     1             '           aCCSD(T)       = ',F20.12)
            CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
         ENDIF
C
         IF(NONHF.AND.IDRLVL.GT.0.AND.T3GRAD) WRITE(6,2039) E5TEST
C
C Reset ICLLVL for special noniterative HF triples.
C
         IF(METACT.EQ.12.OR.(METACT.GE.26.AND.METACT.LE.31))THEN
            ICLLVL = METACT
         ENDIF
C
         E45T(1) = E4T
         E45T(2) = E5ST
         E45T(3) = E4ST
         CALL PUTREC(20,'JOBARC','NONITT3',IINTFP*20,E45T)
C
         RETURN
 1061    FORMAT(' @LTRPS-I, Reference function information : ')
 1062 FORMAT('          Reference function is general non-HF ')
 1063 FORMAT('          Reference function is Hartree-Fock ')
 1064 FORMAT('          Reference function is RHF Hartree-Fock ')
 1065 FORMAT('          Reference function is RHF non-HF ')
 1066 FORMAT('          Reference function is UHF Hartree-Fock ')
 1067 FORMAT(' @LTRPS-I, Reference function is UHF non-HF ')
 1068 FORMAT(' @LTRPS-I, Reference function is ROHF ')
 1069 FORMAT(' @LTRPS-I, Reference function is QRHF ')
 1072 FORMAT(' @LTRPS-I, INT1 is true ')
 1073 FORMAT(' @LTRPS-I, INT1 is false ')
 1074 FORMAT(' @LTRPS-I, INT2 is true ')
 1075 FORMAT(' @LTRPS-I, INT2 is false ')
 1081 FORMAT(' @LTRPS-I, Lengths of ijka integrals. #1 ',I10,/,
     1       '                                      #2 ',I10,/,
     1       '                                      #3 ',I10,/,
     1       '                                      #4 ',I10,/,
     1       '                                     Max ',I10)
 2001 FORMAT(' @LTRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS14 ',I12,/,
     1 '          Number of free words for TRPS14 ',I12)
 2002 FORMAT(' LTRPS-I, Insufficient memory for TRPS14 pre-allocation.')
 2011 FORMAT(' @LTRPS-I, Available memory in integer words ',I12,/,
     1 '           Number of pre-allocated words for TRPS2 ',I12,/,
     1 '           Number of free words for TRPS2 ',I12)
 2012 FORMAT(' LTRPS-I, Insufficient memory for TRPS2 pre-allocation. ')
 2021 FORMAT(' @LTRPS-I, Available memory in integer words ',I12,/,
     1 '           Number of pre-allocated words for TRPS3 ',I12,/,
     1 '           Number of free words for TRPS3 ',I12)
 2022 FORMAT(' LTRPS-I, Insufficient memory for TRPS3 pre-allocation. ')
 2039 FORMAT(' @LTRPS-I, E5TEST is ',F20.12)
      END
