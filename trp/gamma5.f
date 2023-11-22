      SUBROUTINE GAMMA5(ICORE,MAXCOR,IUHF)
C
C
C  THIS ROUTINE COMPUTES THE GAMMA INTERMEDIATE 5
C
C  GAMMA 5 IS ONLY REQUIRED FOR FOURTH ORDER MBPT AND    
C  ALL CC METHODS WHICH INVOLVE SINGLES
C
C MBPT(4) :
C
C  G(AB,CI) = SUM M T(M,C) T(MI,AB)
C
C QCISD :
C
C  G(AB,CI) = 1/2 SUM M L(M,C) T(MI,AB) + 1/2 SUM M T(M,C) L(MI,AB)
C
C CCSD :
C
C  G(AB,CI) = 1/2 SUM M L(M,C) TAU(MI,AB) + 1/2 SUM M T(M,C) L(MI,AB)
C
C             - 1/4 SUM M,N SUM F TAU(MN,AB) L(MN,EF) T(I,F)
C
C             - 1/2 P(AB) G(C,A) T(I,B)
C
C             - 1/2 P(AB) SUM N H(NC,IA) T(N,B)
C
C  THE FOLLOWING SPIN CASES HAVE TO BE CONSIDERED
C
C  AAAA    AA    AAAA    (UHF ONLY)
C  BBBB    BB    BBBB    (UHF ONLY)
C  ABBA    BB    BAAB    (UHF ONLY) STORED WITh NEGATIVE SIGN)
C  ABAB    AA    ABAB
C
C  THIS ROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C  CODED AUGUST/90   JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYG,DISSYH,DISSYH1,DISSYH2,POP,VRT
      LOGICAL LAMBDA,MBPT4,TAU
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2(2) 
CJ      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
CJ     &             TRIP2
C
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     1                LWIC15,LWIC16,LWIC17,LWIC18,
     1                LWIC21,LWIC22,LWIC23,LWIC24,
     1                LWIC25,LWIC26,LWIC27,LWIC28,
     1                LWIC31,LWIC32,LWIC33,
     1                LWIC34,LWIC35,LWIC36,
     1                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      DATA ONE,ONEM,HALF,HALFM /1.D0,-1.D0,0.5D0,-0.5D0/
C
      MXCOR=MAXCOR
CJ      LAMBDA=CCSD.OR.QCISD
CJ      TAU=CCSD
CJ      MBPT4=M4SDQ.OR.M4SDTQ
      LAMBDA = .FALSE.
      TAU    = .FALSE.
      TRIP1  = .TRUE.
C
C  ALLOCATE FIRST MEMORY FOR T1-AMPLITUDES
C
C   STRUCTURE OF THE MEMORY ALLOCATION
C
C    I0T1A =       T1 AMPLITUDES AA
C
C    I0T2A =       L1 AMPLITUDES AA  (IF CCSD OR QCISD)
C
C    I0T1B =       T1 AMPLITUDES BB                      UHF ONLY
C
C    I0T2B =       L1 AMPLITUDES BB  (IF CCSD OR QCISD)  UHF ONLY
C
      I0T1A=MXCOR+1-NTAA*IINTFP
      MXCOR=MXCOR-NTAA*IINTFP
      CALL GETLST(ICORE(I0T1A),1,1,1,1,90)  
CJ      IF(LAMBDA) THEN
CJ       I0T2A=I0T1A-NTAA*IINTFP
CJ       MXCOR=MXCOR-NTAA*IINTFP
CJ       CALL GETLST(ICORE(I0T2A),1,1,2,1,190)
CJ      ELSE
       I0T2A=I0T1A
CJ      ENDIF
      IF(IUHF.EQ.0) THEN
       I0T1B=I0T1A
       I0T2B=I0T2A 
      ELSE
       I0T1B=I0T2A-NTBB*IINTFP
       MXCOR=MXCOR-NTBB*IINTFP
       CALL GETLST(ICORE(I0T1B),1,1,1,2,90)
CJ       IF(LAMBDA) THEN
CJ        I0T2B=I0T1B-NTBB*IINTFP
CJ        MXCOR=MXCOR-NTBB*IINTFP
CJ        CALL GETLST(ICORE(I0T2B),1,1,2,2,190)
CJ       ELSE
        I0T2B=I0T1B
CJ       ENDIF
      ENDIF
C
C  NOW PERFORM MULTIPLICATION
C
      DO 1000 ISPIN=1,IUHF+1
C
      IF(ISPIN.EQ.1) THEN
       IOFFT=I0T2A
       I0T1=I0T1A
       I0T2=I0T1B
       NT=NTAA
      ELSE
       IOFFT=I0T2B
       I0T1=I0T1B
       I0T2=I0T1A
       NT=NTBB
      ENDIF
      IF(IUHF.EQ.1) THEN
C
C  FIST PERFORM AAAA MULTIPLICATION
C
C  NOTE THAT THE ACTUAl SIGN MUST HERE THE NEGATIVE OF 
C  THE SIGN GIVEN ABOVE BECAUSE WE ARE SWITCHING INDICES
C
CJ      LISTG=126+ISPIN
      IF(ISPIN.EQ.1)THEN
       LISTG = LWIC25
      ELSE
       LISTG = LWIC26
      ENDIF
CJ      IF(MBPT4) THEN
CJ       LISTT=43+ISPIN
      LISTT = 13 + ISPIN
CJ       FACT=ONE
       FACT=ONE
CJ      ELSE IF(LAMBDA) THEN
CJ       LISTT=43+ISPIN
CJ       LISTL=143+ISPIN
CJ       FACT=HALFM
CJ      ENDIF
CJ      IF(CCSD) THEN
CJ       LISTH=53+ISPIN
CJ       LISTG1=165+ISPIN
CJ       FACTH=HALF
CJ      ENDIF
C
      DO 100 IRREP=1,NIRREP

       NOCCSQ=0
       NVRTSQ=0
       DO 101 IRREPJ=1,NIRREP
        IRREPI=DIRPRD(IRREP,IRREPJ)
        NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*POP(IRREPI,ISPIN)
        NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)*VRT(IRREPI,ISPIN)
101    CONTINUE
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
       NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
C
C  FOR CCSD CALCULATE FIRST THE CONTRIBUTION OF H4 TO GAMMA5
C
CJ       IF(CCSD) THEN
CJ        DISSYH=IRPDPD(IRREP,ISYTYP(1,LISTH))
CJ        NUMSYH=IRPDPD(IRREP,ISYTYP(2,LISTH))
CJ        I001=1
CJ        I002=I001+IINTFP*MAX(NUMSYH*DISSYH,DISSYG*NUMSYG) 
CJ        I003=I002+IINTFP*NVRTSQ*NUMSYG
CJ        I004=I003+3*IINTFP*MAX(NUMSYG,DISSYG,NUMSYH,DISSYH)
CJ        IF(I004.LT.MXCOR) THEN
CJ         CALL H4G5AA(ICORE(I001),ICORE(I002),
CJ     &               ICORE(I001),ICORE(I0T1),FACTH,DISSYH,
CJ     &               DISSYG,NUMSYH,NUMSYG,NVRTSQ,POP(1,ISPIN),
CJ     &               VRT(1,ISPIN),LISTH,LISTG,IRREP,ICORE(I003),
CJ     &               TRIP1,LISTG1)
CJ        ELSE
CJ         CALL INSMEM('H4G5AA',I004,MXCOR)
CJ        ENDIF
CJC
CJC  CALCULATE IN ADDITION THE T1*T2*L2 CONTRIBUTION TO GAMMA5
CJC
CJ        I001=1
CJ        I002=I001+IINTFP*MAX(NUMSYT*DISSYT,NUMSYT*NUMSYG)
CJ        I003=I002+IINTFP*NUMSYT*NVRTSQ
CJ        MAXSIZE=(MXCOR-I003)/IINTFP
CJC
CJ        IF(MIN(NUMSYG,DISSYG,NUMSYT,DISSYT).NE.0) THEN
CJ        IF(MAXSIZE.GT.DISSYG) THEN
CJC
CJ         CALL T1G5AA(ICORE(I001),ICORE(I002),ICORE(I001),ICORE(I002),
CJ     &               ICORE(I003),MAXSIZE,ICORE(I0T1),POP(1,ISPIN),
CJ     &               VRT(1,ISPIN),DISSYT,DISSYT,DISSYG,NUMSYT,
CJ     &               NUMSYT,NUMSYG,LISTL,LISTT,LISTG,IRREP,ISPIN)
CJ        ELSE
CJ         CALL INSMEM('T1G5AA',MAXSIZE,DISSYG)
CJ        ENDIF
CJ        ENDIF
CJC
CJ       ENDIF
C
       I001=1
       I002=I001+IINTFP*MAX(NOCCSQ*DISSYT,NF2(ISPIN))
       I003=I002+IINTFP*NUMSYG*DISSYG
       IF(MIN(NUMSYG,DISSYG).NE.0) THEN
        I004=I003+3*IINTFP*MAX(NUMSYG,NUMSYT,DISSYG,DISSYT)
        IF(MXCOR.GT.I004) THEN
         CALL G5AA(ICORE(I001),ICORE(I002),ICORE(IOFFT),ICORE(I0T1),
     &             ICORE(I001),FACT,LAMBDA,TAU,TRIP1,ISPIN,
     &             POP(1,ISPIN),VRT(1,ISPIN),NT,DISSYT,NUMSYT,
     &             DISSYG,NUMSYG,LISTT,LISTL,LISTG,IRREP,
     &             ICORE(I003))
        ELSE
         CALL INSMEM('G5AA',I004,MXCOR)
        ENDIF
       ENDIF
100   CONTINUE
      ENDIF
C
C  AB SPIN CASE
C
CJ      LISTG=131-ISPIN
      IF(ISPIN.EQ.1)THEN
       LISTG = LWIC28
      ELSE
       LISTG = LWIC27
      ENDIF
CJ      IF(MBPT4) THEN
CJ       LISTT=46
      LISTT = 16
CJ       FACT=ONEM
       FACT=ONEM
CJ      ELSE IF(LAMBDA) THEN
CJ       LISTT=46
CJ       LISTL=146 
CJ       FACT=HALF
CJ      ENDIF
CJ      IF(CCSD) THEN
CJ       LISTH1=55+ISPIN
CJ       LISTH2=57+ISPIN
CJ       LISTG1=170-ISPIN
CJ       FACTH=HALF
CJ      ENDIF
CJC
C  LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
C
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       DISSYG=IRPDPD(IRREP,ISYTYP(1,LISTG))
       NUMSYG=IRPDPD(IRREP,ISYTYP(2,LISTG))
C
C  FOR CCSD CALCULATE FIRST THE CONTRIBUTION OF H4 TO GHAMMA5
C
CJ       IF(CCSD) THEN
CJ        DISSYH1=IRPDPD(IRREP,ISYTYP(1,20+ISPIN))
CJ        DISSYH2=IRPDPD(IRREP,ISYTYP(1,LISTH2))
CJ        NUMSYH1=IRPDPD(IRREP,ISYTYP(2,20+ISPIN))
CJ        NUMSYH2=IRPDPD(IRREP,ISYTYP(2,LISTH2))
CJ        I001=1
CJ        I002=I001+IINTFP*MAX(NUMSYH1*DISSYH1,NUMSYH2*DISSYH2,
CJ     &                       DISSYG*NUMSYG)
CJ        I003=I002+IINTFP*NUMSYG*DISSYG
CJ        I004=I003+3*IINTFP*MAX(NUMSYH1,NUMSYH2,NUMSYG,DISSYH1,
CJ     &                         DISSYH2,DISSYG)
CJ        IF(I004.LT.MXCOR) THEN
CJ         CALL H4G5AB(ICORE(I001),ICORE(I001),ICORE(I002), 
CJ     &               ICORE(I001),ICORE(I0T1),ICORE(I0T2),FACTH,
CJ     &               DISSYH1,DISSYH2,DISSYG,NUMSYH1,NUMSYH2,NUMSYG,
CJ     &               POP(1,ISPIN),POP(1,3-ISPIN),VRT(1,ISPIN),
CJ     &               VRT(1,3-ISPIN),LISTH1,LISTH2,LISTG,ISPIN,
CJ     &               IRREP,IUHF,ICORE(I003),TRIP1,LISTG1)
CJ        ELSE
CJ         CALL INSMEM('H4G5AB',I004,MXCOR)
CJ        ENDIF
CJC
CJC  CALCULATE IN ADDITION THE T1*T2*L2 CONTRIBUTION TO GAMMA5
CJC
CJ        I001=1
CJ        I002=I001+IINTFP*MAX(NUMSYT*DISSYT,NUMSYT*NUMSYG)
CJ        I003=I002+IINTFP*NUMSYT*DISSYT
CJ        I004=I003+3*IINTFP*MAX(NUMSYT,NUMSYG,DISSYT,DISSYG) 
CJ        MAXSIZE=(MXCOR-I004)/IINTFP
CJ        IF(MIN(NUMSYT,DISSYT,NUMSYG,DISSYG).NE.0) THEN
CJ        IF(MAXSIZE.GE.DISSYG) THEN
CJ         CALL T1G5AB(ICORE(I001),ICORE(I002),ICORE(I001),ICORE(I002),
CJ     &               ICORE(I004),MAXSIZE,ICORE(I0T1),ICORE(I0T2),
CJ     &               POP(1,ISPIN),POP(1,3-ISPIN),VRT(1,ISPIN),
CJ     &               VRT(1,3-ISPIN),DISSYT,DISSYT,DISSYG,NUMSYT,
CJ     &               NUMSYT,NUMSYG,LISTL,LISTT,LISTG,IRREP,ISPIN,
CJ     &               ICORE(I003))
CJ        ELSE
CJ         CALL INSMEM('T1G5AB',DISSYG,MAXSIZE)
CJ        ENDIF
CJ        ENDIF
CJ       ENDIF
C
       I001=1
       I002=I001+MAX(IINTFP*NUMSYT*DISSYT,IINTFP*NF2(ISPIN))
       I003=I002+IINTFP*NUMSYG*DISSYG
       IF(MIN(NUMSYG,DISSYG).NE.0) THEN
        I004=I003+3*IINTFP*MAX(NUMSYG,NUMSYT,DISSYG,DISSYT)
        IF(MXCOR.GT.I004) THEN
         CALL G5AB(ICORE(I001),ICORE(I002),ICORE(IOFFT),
     &             ICORE(I0T1A),ICORE(I0T1B),ICORE(I001),
     &             FACT,LAMBDA,TAU,TRIP1,ISPIN,POP(1,ISPIN),
     &             POP(1,3-ISPIN),VRT(1,ISPIN),VRT(1,3-ISPIN),NT,
     &             DISSYT,NUMSYT,DISSYG,NUMSYG,LISTT,LISTL, 
     &             LISTG,IRREP,ICORE(I003))
        ELSE 
         CALL INSMEM('G5AB',I004,MXCOR)
        ENDIF
       ENDIF
200   CONTINUE
1000  CONTINUE
C
C ALL DONE
C
c      CALL CHECKGAM(ICORE,30,130,2.)
c      IF(IUHF.EQ.1) THEN
c      CALL CHECKGAM(ICORE,27,127,2.)
c       CALL CHECKGAM(ICORE,28,128,2.)
c       CALL CHECKGAM(ICORE,29,129,2.)
c      ENDIF
C
      RETURN
      END
