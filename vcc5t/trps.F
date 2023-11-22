      SUBROUTINE TRPS(ICORE,MAXCOR,IUHF,E4T,S1A,S1B)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION E4T,E4TAAA,E4TBBB,E4TAAB,E4TBBA,E5ST,
     1                 E5STA,E5STB,E4ST,E4STA,E4STB,
     1                 E5TD1,E5TD4,E5TD2,E5TD3,
     1                 E5TT1,E5TT4,E5TT2,E5TT3,
     1                 E5QT1,E5QT4,E5QT2,E5QT3,
     1                 E6TT1,E6TT4,E6TT2,E6TT3,
     1                 ECCTOT,
     1                 ETCCSD,EPAREN,PARENS
      DOUBLE PRECISION ALLT3
      DOUBLE PRECISION S1A(1),S1B(1)
      LOGICAL NONHF,HF
      LOGICAL ISKP14,ISKP2,ISKP3
      LOGICAL CCSDT4,CCSDT
      LOGICAL INCORE
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      INTEGER SCRDIM,MAXOO,MAXVV,MAXOV
#include "maxbasfn.par"
#include "trps.par"
      DATA MAXNBF/MAXBASFN/
      DIMENSION ICORE(MAXCOR)
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
      EQUIVALENCE(INONHF,IFLAGS(38))
      EQUIVALENCE(IORBTP,IFLAGS(39))
      EQUIVALENCE(IPTORB,IFLAGS(64))
C
C     TRIPLES COMMON BLOCKS
C
      COMMON /TRIPLES/ TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ACTORB/ ABSVRT(MAXBASFN,8,2),ABSOCC(MAXBASFN,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T3ENRG/ E4TAAA,E4TBBB,E4TAAB,E4TBBA,
     1                E5STA ,E5STB ,
     1                E4STA ,E4STB ,
     1                E5TD1 ,E5TD4 ,E5TD2 ,E5TD3,
     1                E5TT1 ,E5TT4 ,E5TT2 ,E5TT3,
     1                E5QT1 ,E5QT4 ,E5QT2 ,E5QT3,
     1                E6TT1 ,E6TT4 ,E6TT2 ,E6TT3
C
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /T3METH/ CCSDT4,CCSDT
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
      COMMON /T3CORE/ INCORE
C
      IF(T3STOR)THEN
C
C     Read the T3 file.
C
      NT3 = NRECS * LNPHYR
      IF(NT3.LT.NUMT3)THEN
      INCORE = .TRUE.
      WRITE(6,9010)
 9010 FORMAT(' @TRPS-I, There is sufficient core for the T3 file. ')
      ELSE
      INCORE = .FALSE.
      WRITE(6,9020)
 9020 FORMAT(' @TRPS-I, There is insufficient core for the T3 file. ')
      ENDIF
C
      IF(INCORE) CALL GETALLT3
      ENDIF
C
      ISKP14 = .FALSE.
      ISKP2  = .FALSE.
      ISKP3  = .FALSE.
C
      WRITE(6,1010)
 1010 FORMAT(' @TRPS-I, Welcome to TRPS. Watch your RUs ! ')
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
C     Transpose records of LWIC17 and LWIC27 (the latter only if it is
C     not the same list).
C
      IF(IUHF.EQ.1)THEN
      CALL TSPABCI(ICORE(1),MAXCOR,1)
      ENDIF
C
      E4T    = 0.0D+00
      E4TAAA = 0.0D+00
      E4TBBB = 0.0D+00
      E4TAAB = 0.0D+00
      E4TBBA = 0.0D+00
      E5ST   = 0.0D+00
      E5STA  = 0.0D+00
      E5STB  = 0.0D+00
      E5TD1  = 0.0D+00
      E5TD4  = 0.0D+00
      E5TD2  = 0.0D+00
      E5TD3  = 0.0D+00
      E5TT1  = 0.0D+00
      E5TT4  = 0.0D+00
      E5TT2  = 0.0D+00
      E5TT3  = 0.0D+00
      E5QT1  = 0.0D+00
      E5QT4  = 0.0D+00
      E5QT2  = 0.0D+00
      E5QT3  = 0.0D+00
      E6TT1  = 0.0D+00
      E6TT4  = 0.0D+00
      E6TT2  = 0.0D+00
      E6TT3  = 0.0D+00
C
C     If this is CCSDT or CCSDT-4 reorder MBEJ lists.
C
      IF(CCSDT4.OR.CCSDT) CALL SETVOVO(ICORE(1),MAXCOR,IUHF,1)
C
C     For CCSDT-4 (the World's most important method) replace inter-
C     mediates with integrals on lists 54-59.
C
      IF(CCSDT4)          CALL   CSDT4(ICORE(1),MAXCOR,IUHF)
C
      IF(IUHF.EQ.0.OR.ISKP14) GOTO 1000
C
      DO  110 ISPIN=1,IUHF+1
      NO = NOCCO(ISPIN)
      NV = NVRTO(ISPIN)
C
      IF(NO.LT.3.OR.NV.LT.3) GOTO 110
C
C     TO AVOID THE SILLY FPS500 PROBLEM OF DOUBLE WORD BOUNDARIES,
C     GENUINE INTEGERS ARE ALLOCATED AS THOUGH THEY ARE DOUBLE PRECISION.
C     THIS WILL LEAD TO A WASTE OF AT MOST A FEW THOUSAND WORDS WE
C     THINK. (EXCEPT, OF COURSE, IN ADENOSINE TRIPHOSPHATE PVQZ
C     CALCULATIONS).
C
C     I000   OOOV INTEGRALS
C     I010   EVAL EIGENVALUES
C     I020   IADBLK
C     I030   LENBLK
C     I040   IADT2
C     I050   LENT2
C     I060   IADV
C     I070   LENV
C     I080   OOOVAD
C     I090   LENINT
C     I100   ICORE (THE REST)
C
      I000 =     1
      I010 = I000 + IINTFP * LNOOOV
      I020 = I010 + IINTFP * (NO + NV)
      I030 = I020 + IINTFP * 8
      I040 = I030 + IINTFP * 8
      I050 = I040 + IINTFP * 8
      I060 = I050 + IINTFP * 8
      I070 = I060 + IINTFP * 8
      I080 = I070 + IINTFP * 8
      I090 = I080 + IINTFP * 8
      I100 = I090 + IINTFP * 8
      NEED = I100
C     COMPUTE NUMBER OF FREE INTEGER WORDS
      NFREE = MAXCOR + 1 - NEED
      NALLOC = NEED - 1
      WRITE(6,2001) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
       WRITE(6,2002)
       CALL INSMEM('TRPS',NEED,MAXCOR)
      ENDIF
      CALL TRPS14(ICORE(I000),
     1            ICORE(I010),
     1            ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     1            ICORE(I060),ICORE(I070),ICORE(I080),ICORE(I090),
     1            LNOOOV,ICORE(I100),
     1            NFREE,ISPIN,NO,NV,S1A,S1B)
  110 CONTINUE
 1000 CONTINUE
C
      IF(NOCA .LT.2.OR.NOCB .LT.1.OR.
     1   NVRTA.LT.2.OR.NVRTB.LT.1.OR.ISKP2) GOTO 2000
C
C     Transpose certain OOOV lists.
C
      CALL SETOOOV(ICORE(1),MAXCOR,IUHF,1)
C
C     I000   EVAL
C     I010   SCR1
C     I020   SCR2
C     I030   SCR3
C     I040   T2ABA
C     I050   T2ABBI
C     I060   T2ABBJ
C     I070   T2AAAI
C     I080   T2AAAJ
C     I090   IADT3
C     I100   LENT3
C     I110   IADT2
C     I120   LENT2
C     I130   IADV
C     I140   LENV
C     I150   ICORE (THE REST)
C
      I000 =    1
      I010 = I000 + IINTFP *(NOCA + NVRTA) * 2
      I020 = I010 + IINTFP * SCRDIM
      I030 = I020 + IINTFP * SCRDIM
      I040 = I030 + IINTFP * SCRDIM
      I050 = I040 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I060 = I050 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I070 = I060 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I080 = I070 + IINTFP * MAX1  * MAX(NOCA,NOCB)
      I090 = I080 + IINTFP * MAX1  * MAX(NOCA,NOCB)
      I100 = I090 + IINTFP * 8
      I110 = I100 + IINTFP * 8
      I120 = I110 + IINTFP * 8
      I130 = I120 + IINTFP * 8
      I140 = I130 + IINTFP * 8
      I150 = I140 + IINTFP * 8
      NEED = I150
C
C     Compute number of free integer words.
C
      NFREE = MAXCOR + 1 - NEED
      NALLOC = NEED - 1
      WRITE(6,2011) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
       WRITE(6,2012)
       CALL INSMEM('TRPS',NEED,MAXCOR)
      ENDIF
      CALL TRPS2(ICORE(I000),
     1           ICORE(I010),ICORE(I020),ICORE(I030),
     1           ICORE(I040),ICORE(I050),ICORE(I060),
     1                       ICORE(I070),ICORE(I080),
     1           ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     1           ICORE(I130),ICORE(I140),
     1           ICORE(I150),
     1           NOCA,NOCB,NVRTA,NVRTB,NFREE,IUHF,S1A,S1B)
      CALL SETOOOV(ICORE(1),MAXCOR,IUHF,2)
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
      CALL TSPABCI2(ICORE(1),MAXCOR,1)
      CALL AT2IAB(ICORE(1),0,.FALSE.)
      IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE(1),MAXCOR,IUHF,1)
      ENDIF
C
      I000 =    1
      I010 = I000 + IINTFP *(NOCA + NVRTA) * 2
      I020 = I010 + IINTFP * SCRDIM
      I030 = I020 + IINTFP * SCRDIM
      I040 = I030 + IINTFP * SCRDIM
      I050 = I040 + IINTFP * LNOOOV * 4
      I060 = I050 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I070 = I060 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I080 = I070 + IINTFP * MAX13 * MAX(NOCA,NOCB)
      I090 = I080 + IINTFP * MAX2  * MAX(NOCA,NOCB)
      I100 = I090 + IINTFP * MAX2  * MAX(NOCA,NOCB)
      I110 = I100 + IINTFP * 8
      I120 = I110 + IINTFP * 8
      I130 = I120 + IINTFP * 8
      I140 = I130 + IINTFP * 8
      I150 = I140 + IINTFP * 8
      I160 = I150 + IINTFP * 8
      I170 = I160 + IINTFP * 8             * 4
      I180 = I170 + IINTFP * 8             * 4
      NEED = I180
C
C     Compute number of free integer words.
C
      NFREE = MAXCOR + 1 - NEED
      NALLOC = NEED - 1
      WRITE(6,2021) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
       WRITE(6,2022)
       CALL INSMEM('TRPS',NEED,MAXCOR)
      ENDIF
      CALL TRPS3(ICORE(I000),
     1           ICORE(I010),ICORE(I020),ICORE(I030),
     1           ICORE(I040),
     1           ICORE(I050),ICORE(I060),ICORE(I070),
     1                       ICORE(I080),ICORE(I090),
     1           ICORE(I100),ICORE(I110),ICORE(I120),ICORE(I130),
     1           ICORE(I140),ICORE(I150),ICORE(I160),ICORE(I170),
     1           ICORE(I180),
     1           LNOOOV,NOCA,NOCB,NVRTA,NVRTB,NFREE,IUHF,S1A,S1B)
C
      IF(IUHF.EQ.1)THEN
      CALL TSPABCI2(ICORE(1),MAXCOR,2)
c     CALL TSPABCI(ICORE(1),MAXCOR,2)
      CALL AT2IAB(ICORE(1),1,.FALSE.)
      IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE(1),MAXCOR,IUHF,2)
      ENDIF
C
 3000 CONTINUE
C
C     Put abic/ijka, t2 ab, and d2t2 ab lists back in original form.
C
      IF(IUHF.EQ.1)THEN
c     CALL TSPABCI2(ICORE(1),MAXCOR,2)
      CALL TSPABCI(ICORE(1),MAXCOR,2)
c     CALL AT2IAB(ICORE(1),1,INT2)
c     IF(CCSDT.OR.CCSDT4) CALL SYMABCD(ICORE(1),MAXCOR,IUHF,2)
      ENDIF
C
C     If this is CCSDT or CCSDT-4 reorder MBEJ lists.
C
      IF(CCSDT4.OR.CCSDT) CALL SETVOVO(ICORE(1),MAXCOR,IUHF,2)
C
C     Attempt to compute the RMS error in T3.
C
c     CALL RMST3(ICORE(1),MAXCOR,IUHF)
C
C     Move the triples amplitudes so that the new (ie current) ones will
C     be the old ones in the next iteration.
C
c      IF(T3STOR)THEN
c      CALL MOVLIST(ICORE(1),1,5)
cc      CALL MOVLIST(ICORE(1),2,6)
c      CALL MOVLIST(ICORE(1),3,7)
c      CALL MOVLIST(ICORE(1),4,8)
c      IF(INCORE) CALL PUTALLT3
c      ENDIF
C
C     return if this is an "iterative noniterative" calculation.
C
C     If method is CCSD+T(CCSD), CCSD(T), QCISD(T) then evaluate
C     CCSD/QCISD energy. Note F*T1 contribution must be included
C     if reference is non-Hartree-Fock.
C
      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.21.OR.ICLLVL.EQ.22)THEN
      CALL CCENRG(ICORE(1),IUHF,ECCTOT,NONHF)
      ENDIF
C
C     Summarize results.
C
c      IF(IUHF.EQ.0)THEN
c      E4T  = 2.0D+00 * (E4TAAA + E4TAAB)
c      E5ST = 2.0D+00 * (E5STA  + E5STB )
c      E4ST = 2.0D+00 * (E4STA  + E4STB )
c      ELSE
c      E4T  =            E4TAAA + E4TBBB + E4TAAB + E4TBBA
cC      WRITE(6,*) ' @TRPS-I, the energies ',E4TAAA,E4TBBB,E4TAAB,
cC     1 E4TBBA,E4T
c      E5ST =            E5STA  + E5STB
c      E4ST =            E4STA  + E4STB
c      ENDIF
C
C
c      IF(ICLLVL.EQ.11.OR.ICLLVL.EQ.22)THEN
c      WRITE(6,2033) E4T,E5ST,E4ST,E4T,E5ST,E4ST
c 2033 FORMAT(' @TRPS-I, E4T (T(CCSD)) = ',F15.12,/,
c     1       '          E5ST          = ',F15.12,/,
c     1       '          E4ST          = ',F15.12,/,
c     1       '          In other words  ',/,
c     1       '          <T3[2](Dagger) D3  T3[2]> = ',F15.12,/,
c     1       '          <T1   (Dagger) Q1 WT3[2]> = ',F15.12,/,
c     1       '          <T2   (Dagger) Q2 FT3[2]> = ',F15.12)
c      ENDIF
c      ETCCSD = ECCTOT + E4T
c      EPAREN = ECCTOT + E4T + E5ST
c      WRITE(6,2035) ECCTOT,ETCCSD,EPAREN
c 2035 FORMAT(' @TRPS-I, CCSD           = ',F20.12,/,
c     1       '          CCSD + T(CCSD) = ',F20.12,/,
c     1       '          CCSD(T)        = ',F20.12)
c      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
c      ENDIF
cC
cC     --- CCSD(T) with Non-HF reference ---
cC
c      IF(ICLLVL.EQ.22.AND.NONHF)THEN
c      EPAREN = ECCTOT + E4T + E5ST + E4ST
c      WRITE(6,2036) ECCTOT,EPAREN
c 2036 FORMAT(' @TRPS-I, CCSD           = ',F20.12,/,
c     1       '          CCSD(T)        = ',F20.12)
c      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,EPAREN)
c      ENDIF
C
      RETURN
 1081 FORMAT(' @TRPS-I, Lengths of ijka integrals. #1 ',I10,/,
     1       '                                     #2 ',I10,/,
     1       '                                     #3 ',I10,/,
     1       '                                     #4 ',I10,/,
     1       '                                    Max ',I10)
 2001 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS14 ',I12,/,
     1 '          Number of free words for TRPS14 ',I12)
 2002 FORMAT(' @TRPS-F, Insufficient memory for TRPS14 pre-allocation.')
 2011 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS2 ',I12,/,
     1 '          Number of free words for TRPS2 ',I12)
 2012 FORMAT(' @TRPS-F, Insufficient memory for TRPS2  pre-allocation.')
 2021 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS3 ',I12,/,
     1 '          Number of free words for TRPS3 ',I12)
 2022 FORMAT(' @TRPS-F, Insufficient memory for TRPS3  pre-allocation.')
      END
