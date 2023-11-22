#include "flags.h"
      Subroutine T3_guess(Icore, Maxcor, Iuhf, INIT, T3STOR, INCORE,
     +                    Extrnl_ccsd)

      Implicit Double Precision (A-H, O-Z)

C Solve D3T3 = W*T2 + FOO*T3 + FVV T3 iteratively for NON-HF (ROHF and
C QRHF or anyother) CCSD(T) when semicanonical orbitals can nor be used.
C Natuallrly built on top of what is available in triples codes. When
C INIT is true only the first term is included. 
C  
      Logical Int1, Int2, NONHF, T3STOR, PRESENT, INIT, INCORE
      Logical T3GRAD 
      Logical Extrnl_ccsd 
      Dimension Icore(Maxcor)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "info.com"
#include "maxbasfn.par"
#include "trps.par"

      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     &                NCOMB(4)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     &                LWIC15,LWIC16,LWIC17,LWIC18,
     &                LWIC21,LWIC22,LWIC23,LWIC24,
     &                LWIC25,LWIC26,LWIC27,LWIC28,
     &                LWIC31,LWIC32,LWIC33,
     &                LWIC34,LWIC35,LWIC36,
     &                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /GAMLIS/ LISGO1,LISGO2,LISGO3,LISGO4,LISGV1,LISGV2,
     &                LISGV3,LISGV4
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     &                NRECS,LUAUX 
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /EDEBUG/ E14, E2AAB, E2BBA
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/  IFLAGS2(500)

      INT1   = .FALSE.
      INT2   = .FALSE.
      NONHF  = .TRUE.

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

      ICLLVL = IFLAGS(2)
      IDRLVL = IFLAGS(3)
      IREFNC = IFLAGS(11)
      IQRHFP = IFLAGS(32)
      IQRHFM = IFLAGS(33)
      INONHF = IFLAGS(38)
      NONHF  = .TRUE.
C
      IF(IDRLVL.GT.0)THEN
C
C     If this is a numerical gradient and IDRLVL > 0, turn the
C     derivative flag off, but turn on at end.
C
        Write(*,*) IFLAGS2(h_IFLAGS2_grad_calc) 

        IF (IFLAGS2(h_IFLAGS2_grad_calc).EQ.2) THEN
          if (iFlags2(h_IFLAGS2_geom_opt).ne.0) WRITE(6,1053)
          if (iFlags(h_IFLAGS_vib)       .eq.3) WRITE(6,1054)
          if (iFlags(h_IFLAGS_vib)       .eq.3) IDROLD = IDRLVL
          IDRLVL    = 0
          IFLAGS(3) = 0
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

      Write(*,*) "@T3-GUESS", IDRLVL
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
      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
 1081 FORMAT(' @TRPS-I, Lengths of ijka integrals. #1 ',I10,/,
     &       '                                     #2 ',I10,/,
     &       '                                     #3 ',I10,/,
     &       '                                     #4 ',I10,/,
     &       '                                    Max ',I10)
C
      IF(IUHF.EQ.0) GOTO 1000

      DO  110 ISPIN=1,2

          E14 = 0.0D0

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
         NFREE = MAXCOR - (I130 - 1)
         NALLOC = I130 - 1
         IF(NFREE.LE.0)THEN
            WRITE(6,2002)
 2002 FORMAT(' T3_iter, Insufficient memory for T3it_trps14 
     & pre-allocation. ')
            CALL INSMEM('TRPS',I130,MAXCOR)
         ENDIF

         Write(*,*) "Entering T3IT_TRPS14"
         CALL T3IT_TRPS14(ICORE(I010),ICORE(I020),ICORE(I025),
     &                    ICORE(I030), ICORE(I040),
     &                    ICORE(I050),ICORE(I060),
     &                    ICORE(I070),ICORE(I080),
     &                    ICORE(I090),ICORE(I100),
     &                    ICORE(I110),ICORE(I120),LNOOOV,
     &                    ICORE(I130),
     &                    NFREE,IUHF,ISPIN,NO,NV,INT1,INT2,
     &                    NONHF,T3STOR,INIT)

  110 CONTINUE
 1000 CONTINUE

      IF(NOCA .LT.2.OR.NOCB .LT.1.OR.
     &   NVRTA.LT.2.OR.NVRTB.LT.1) GOTO 2000
C
C     Transpose certain OOOV lists.
C
       IF (IDRLVL .EQ. 0) THEN

       IF(IUHF.EQ.0)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
       ELSE
        CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)  
       ENDIF

       ELSE
         CALL SETOOOV(ICORE,MAXCOR,IUHF,1)
       ENDIF 

       E2AAB = 0.0D0

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
C     I130   OOOV intermediates
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

      I150 = I140 + IINTFP * 1      * NGAMMA
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

      IF(NFREE.LE.0)THEN
        WRITE(6,2012)
 2012 FORMAT(' T3_iter, Insufficient memory for T3IT_TRPS2
     &  pre-allocation. ')
        CALL INSMEM('TRPS',I240,MAXCOR)
      ENDIF

      CALL T3IT_TRPS2(ICORE(I000),ICORE(I010),
     &                ICORE(I020),ICORE(I030),
     &                ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),
     &                ICORE(I080),ICORE(I090),
     &                ICORE(I100),ICORE(I110),
     &                ICORE(I140),ICORE(I150),
     &                ICORE(I160),ICORE(I161),
     &                ICORE(I162),ICORE(I163),
     &                ICORE(I164),ICORE(I165),
     &                ICORE(I170),ICORE(I180),
     &                ICORE(I190),ICORE(I200),
     &                ICORE(I210),ICORE(I220),
     &                ICORE(I230),ICORE(I240),
     &                LNOOOV,NOCA,NOCB,NVRTA,NVRTB,
     &                NFREE,INT1,INT2,NONHF,IUHF,
     &                T3STOR,INIT)

      IF (IDRLVL .EQ. 0) THEN

      IF(IUHF.EQ.0)THEN
         CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
      ELSE
         CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
         CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
      ENDIF 

      ELSE
         CALL SETOOOV(ICORE,MAXCOR,IUHF,2)
      ENDIF 

 2000 CONTINUE
C
      IF(IUHF.EQ.0.OR.
     1   NOCA .LT.1.OR.NOCB .LT.2.OR.
     1   NVRTA.LT.1.OR.NVRTB.LT.2) GOTO 3000
C
C     Reorganize abci/ijka intermediate lists, t2 ab lists and,
C     if appropriate, d2t2 ab increment lists.
C
      Write(*,*) "Before TSPABCI2", IDRLVL
      IF(IUHF.EQ.1)THEN
        CALL TSPABCI2(ICORE,MAXCOR,1,Extrnl_ccsd)
        CALL AT2IAB(ICORE,46,1,0)
        CALL AT2IAB(ICORE,16,1,0)
      ENDIF
      Write(*,*) "After TSPABCI2"
C
      E2BBA = 0.0D0

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
      i140 = I130 + IINTFP * LNOOOV * 4
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
      IF(NFREE.LE.0)THEN
        WRITE(6,2022)
 2022 FORMAT(' T3_iter, Insufficient memory for T3IT_TRPS3 
     & pre-allocation. ')
        CALL INSMEM('TRPS',I240,MAXCOR)
      ENDIF

      CALL T3IT_TRPS3(ICORE(I000),ICORE(I010),
     &                ICORE(I020),ICORE(I030),
     &                ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),
     &                ICORE(I080),ICORE(I090),
     &                ICORE(I100),ICORE(I110),
     &                ICORE(I120),ICORE(I130),
     &                ICORE(I140),ICORE(I150),
     &                ICORE(I160),ICORE(I161),
     &                ICORE(I162),ICORE(I163),
     &                ICORE(I164),ICORE(I165),
     &                ICORE(I170),ICORE(I180),
     &                ICORE(I190),ICORE(I200),
     &                ICORE(I210),ICORE(I220),
     &                ICORE(I230),ICORE(I240),
     &                LNOOOV,NOCA,NOCB,NVRTA,
     &                NVRTB,NFREE,INT1,INT2,NONHF,
     &                IUHF,T3STOR,INIT)

      IF(IUHF.EQ.1)THEN
        CALL TSPABCI2(ICORE,MAXCOR,2,Extrnl_ccsd)
        CALL AT2IAB(ICORE,46,1,1)
        CALL AT2IAB(ICORE,16,1,1)
      ENDIF

 3000 CONTINUE
C 
      IF(T3STOR)THEN
        CALL MOVLIST(ICORE,1,5)
        CALL MOVLIST(ICORE,2,6)
        IF(IUHF.GT.0)THEN
           CALL MOVLIST(ICORE,3,7)
           CALL MOVLIST(ICORE,4,8)
        ENDIF
        IF(INCORE) CALL PUTALLT3
      ENDIF
C 
      Return
      End
