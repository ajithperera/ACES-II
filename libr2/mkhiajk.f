      SUBROUTINE MKHIAJK(ICORE,MAXCOR,IUHF,LZ0,LT20,LT2R,LT2R2,
     &                   LT1U,LT1D,LRAAAA,LRBBBB,LRABAB,LRBABA,
     &                   LRBAAB,LRABBA,LIJKA0,LIJKL0,LFIA,LABCI0,
     &                   IRREPT,IRREPW,TDER,WDER,ANTI,
     &                   TERM1,TERM2,TERM3,TERM4,TERM5,TERM6,TAU,
     &                   ZERLST,CCSD)
C
C CALCULATES HBAR(IAJK) OR ITS DERIVATIVE, DEPENDING ON INPUT PARAMETERS
C
C LZ0      - OFFSET LIST FOR IAJK TARGET [AAAA LIST NUMBER - 1]
C LT20     - OFFSET LIST FOR AB,IJ ORDERED T2
C LT2R     - LIST FOR AI,BJ ORDERED RESORTED T2
C LT2R2    - LIST FOR AJ,BI ORDERED RESORTED T2
C LT1U     - LIST FOR UNDIFFERENTIATED T1 AMPLITUDES
C LT1D     - LIST FOR DIFFERENTIATED T1 AMPLITUDES [SAME AS LT1U
C            IF THIS CALL USES DIFFERENTIATED W]
C LRAAAA   - LIST FOR MBEJ RING INTERMEDIATES [<||>+T2*<||> ONLY.
C LRBBBB   - LIST FOR mbej RING INTERMEDIATES
C LRABAB   - LIST FOR MbEj RING INTERMEDIATES
C LRBABA   - LIST FOR mBeJ RING INTERMEDIATES
C LRBAAB   - LIST FOR mBEj RING INTERMEDIATES
C LRABBA   - LIST FOR MbeJ RING INTERMEDIATES
C LIJKA0   - OFFSET LIST FOR IJKA INTEGRALS
C LIJKL0   - OFFSET LIST FOR IJKL HBAR 
C LFIA     - LIST FOR F(IA) HBAR
C LABCI0   - OFFSET LIST FOR ABCI INTEGRALS
C IRREPT   - OVERALL SYMMETRY OF T QUANTITIES
C IRREPW   - OVERALL SYMMETRY OF W QUANTITIES
C TDER     - LOGICAL WHICH IS .TRUE. IF T'S ARE DIFFERENTIATED
C WDER     - LOGICAL WHICH IS .TRUE. IF W'S ARE DIFFERENTIATED
C ANTI     - LOGICAL WHICH IS .TRUE. FOR IMAGINARY PERTURBATION
C TERM1    - INTEGRALS ARE INCLUDED
C TERM2    - T1*IJKL PART INCLUDED
C TERM3    - F(IA) CONTRIBUTION INCLUDED
C TERM4    - T1*RING CONTRIBUTION INCLUDED
C TERM5    - CONTRACTION FROM HELL INCLUDED
C TERM6    - TAU*T2 TERM INCLUDED
C TAU      - LOGICAL CONTROLLING IF TAU OR T2 IS USED FOR TERM6
C CCSD     - See FORMW5
C Linear CCSD modifications: Making sure Tau=T2 for TERM6, Use 
C bare integral for TERM2
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TAUFAC,FACT,FACTI
      LOGICAL TDER,WDER,ANTI,TAU,ZERLST,CCSD
      LOGICAL TERM1,TERM2,TERM3,TERM4,TERM5,TERM6,HBAR_4LCCSD
      LOGICAL NONHF
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS/IFLAGS2(500)
C
      IRREPX=DIRPRD(IRREPW,IRREPT)

      HBAR_4LCCSD = .FALSE.
      NONHF       = .FALSE.
      HBAR_4LCCSD = (IFLAGS(2)  .EQ. 6 .OR. IFLAGS2(117) .EQ. 7)
      NONHF       = (IFLAGS(38) .NE. 0)
C
C ZERO OUT TARGET LIST IF SPECIFIED.  OTHERWISE, IT IS AUGMENTED!
C
      IF(ZERLST)THEN
       DO 5 ISPIN=4,4-3*IUHF,-1
        CALL DZERSYM(IRREPX,ICORE,LZ0+ISPIN)
5      CONTINUE
      ENDIF
C
      IF(TDER)THEN
       TAUFAC=1.0D0
      ELSE
       TAUFAC=0.5D0
      ENDIF
C
      IF(TERM5)THEN
       FACT=1.0D0
       IF(IUHF.NE.0)THEN
        CALL DHBIAJKA(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                LT2R-4,LIJKA0,LZ0)
        CALL DHBIAJKB(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                LT2R-4,LIJKA0,LZ0)
        CALL DHBIAJKC(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                LT2R-4,LIJKA0,LZ0)
       ELSE 
        CALL DHBIAJK2(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                LT2R,LT2R2,LIJKA0+4,LZ0+4)
       ENDIF
      ELSE
       FACT=0.0D0
      ENDIF
C
C  ADD IN BARE MO INTEGRALS
C
      IF(TERM1)THEN
       FACTI=1.0D0
       IF(ANTI)FACTI=-FACTI
       CALL DHBIAJK0(ICORE,MAXCOR,IUHF,IRREPX,LZ0,LIJKA0,
     &               FACT,FACTI)
      ENDIF

CSSS      TERM4 = .TRUE.
CSSS      TERM2 = .TRUE.
CSSS      TERM3 = .TRUE.
      
      IF(TERM4)THEN
       CALL SORTRING(ICORE,MAXCOR,IUHF,IRREPW,1,LRAAAA,LRBBBB,LRABAB,
     &               LRBABA,LRBAAB,LRABBA,.TRUE.)
       CALL DHBIAJK3(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &               LT1D,LRAAAA,LRBBBB,LRABAB,LRBABA,LRBAAB,LRABBA,
     &               LZ0,ANTI.AND.WDER)
       CALL SORTRING(ICORE,MAXCOR,IUHF,IRREPW,2,LRAAAA,LRBBBB,LRABAB,
     &               LRBABA,LRBAAB,LRABBA,.TRUE.)
      ENDIF
C
C This is the problem child for Linear CCSD. Instead of the W(im,jk)
C intermediate we need the <im||jk>. 

      IF(TERM2)THEN
       CALL DHBIAJK4(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &               LT1D,LIJKL0,LZ0,HBAR_4LCCSD)
      ENDIF
C
C For HF Linear CCSD methods this is zero. Lets generalize this 
C to work for NONHF references. 
C
      IF(TERM3)THEN
       CALL DHBIAJK5(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &               LT20,LFIA,LZ0,HBAR_4LCCSD,NONHF)          
      ENDIF
C
C Make sure that the Tau = t2 for linear CCSD. 05/2015, Ajith Perera
C
      IF(TERM6)THEN
       CALL DHBIAJK6(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &               TAUFAC,LT1U,LT1D,LT20,LABCI0,LZ0,ANTI.AND.WDER,
     &               TAU,CCSD)
      ENDIF
C
      RETURN
      END