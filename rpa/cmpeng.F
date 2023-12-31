C
cjp
cjp in bwcc it computes diagonal elements of the effective hamiltonian,
cjp not the energy, which is obtained from diagonalization of that Heff
cjp after also the offdiagonal elements are computed
cjp
C  DRIVER FOR THE  CALCULATION OF THE CORRELATION ENERGY FOR A GIVEN SET
C  OF AMPLITUDES
C
C  ARGUMENTS :  ICORE ..... ICORE ARRAY
C               MAXCOR .... DIMENSION OF ICORE
C               NLIST2 .... OFFSET OF T2 LIST ON MOINTS (WITH RESPECT TO
C                            TYPE)
C               NLIST1 .... OFFSET OF T1 LISTS ON MOINTS (WITH RESPECT TO
C                            SPIN TYPE)
C               ECORR ..... RETURNS THE CORRELATION ENERGY FOR ALL SPIN CASES
C               ETOT .....  RETURNS TOTAL CORRELATION ENERGY
C               ETOTT2 ...  RETURNS LINEAR CONTRIBUTION TO
C                            THE CORRELATION ENERGY
C               IUHF .....  IUHF FLAG
C               IPRINT ...  PRINT FLAG
C
      SUBROUTINE CMPENG(ICORE,MAXCOR,NLIST2,NLIST1,ECORR,ETOT,ETOTT2,
     &                  IUHF,IPRINT)
      IMPLICIT INTEGER (A-Z)
      LOGICAL TAU,NONHF
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1
      CHARACTER*2 SPCASE(3)
      DOUBLE PRECISION E,ETOT,FACTOR,ECORR(3),ESPIN,ET2,ETOTT2,
     &                 ESING,SDOT,ESING_KEEP, EDBLS_KEEP
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD
      COMMON /NHFREF/NONHF
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DIMENSION I0T(2),I0F(2)
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA SPCASE /'AA','BB','AB'/
cjp
#include "bwcc.com"
C
c Nevin insured that mxcor is even for alignment
c      MXCOR=MAXCOR
      MXCOR=MAXCOR - MOD(MAXCOR,2)

cYAU: The one-particle entities were moved from the end of icore to the front.
      IFREE=1

      TAU=.FALSE.
CSSS      IF((METHOD.GT.9.AND.SING1.AND.METHOD.NE.21.AND.METHOD.NE.23)
CSSS     &    .OR.NONHF.AND.SING1)THEN
C
C   ALLOCATE MEMORY FOR T1 AMPLITUDES
C
cYAU       I0T(1)=MXCOR+1-NT(1)*IINTFP
cYAU       MXCOR=MXCOR-NT(1)*IINTFP
       I0T(1)=IFREE
       IFREE=IFREE+NT(1)*IINTFP
       CALL GETLST(ICORE(I0T(1)),1,1,1,1+NLIST1,90)
       IF(IUHF.EQ.0) THEN
        I0T(2)=I0T(1)
       ELSE
cYAU        I0T(2)=I0T(1)-NT(2)*IINTFP
cYAU        MXCOR=MXCOR-NT(2)*IINTFP
        I0T(2)=IFREE
        IFREE=IFREE+NT(2)*IINTFP
        CALL GETLST(ICORE(I0T(2)),1,1,1,2+NLIST1,90)
       ENDIF
       TAU=.TRUE.
C
C   FOR NON HF REFERENCES ALLOCATE MEMORY FOR f(a,I)
C
       IF(NONHF) THEN
cYAU        I0F(1)=I0T(2)-NT(1)*IINTFP
cYAU        MXCOR=MXCOR-NT(1)*IINTFP
        I0F(1)=IFREE
        IFREE=IFREE+NT(1)*IINTFP
        CALL GETLST(ICORE(I0F(1)),1,1,1,3,93)
        IF(IUHF.EQ.0) THEN
         I0F(2)=I0F(1)
        ELSE
cYAU         I0F(2)=I0F(1)-NT(2)*IINTFP
cYAU         MXCOR=MXCOR-NT(2)*IINTFP
         I0F(2)=IFREE
         IFREE=IFREE+NT(2)*IINTFP
         CALL GETLST(ICORE(I0F(2)),1,1,1,4,93)
        eNDIF
       ENDIF
CSSS      ENDIF
C
      ETOT=0.D0
      ETOTT2=0.D0
      FACTOR=1.D0
      ESING_KEEP = 0.0D0
      EDBLS_KEEP = 0.0D0
      IF(IUHF.EQ.0)FACTOR=2.D0
      DO 10 ISPIN=1,IUHF+1
       LISTT=43+ISPIN
       ESPIN=0.D0
C
C  FOR NON HF REFERENCE FUNCTION ADD HERE THE T(I,A) f(I,A) CONTRIBUTION
C  THERE ARE AA AND BB CONTRIBUTIONS
C
       IF(NONHF.AND.SING1) THEN
        ESING=SDOT(NT(ISPIN),ICORE(I0T(ISPIN)),1,ICORE(I0F(ISPIN)),1)
        ESPIN=ESPIN+ESING
        ETOT=ETOT+FACTOR*ESING
        ETOTT2=ETOTT2+FACTOR*ESING
        ESING_KEEP = ESING_KEEP + FACTOR*ESING
       ENDIF
C
C  THE TAU(IJ,AB) <IJ//AB> CONTRIBUTION TO THE ENERGY (AA AND BB SPIN
C  CASES
C
       DO 100 IRREP=1,NIRREP
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        IF(MIN(NUMSYT,DISSYT).NE.0) THEN
cYAU        I001=1
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &             ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ISPIN,TAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+FACTOR*ET2
        ETOT=ETOT+FACTOR*E
        ESPIN=ESPIN+E
        EDBLS_KEEP = EDBLS_KEEP + FACTOR*ET2
        ESING_KEEP = ESING_KEEP + FACTOR*(E - ET2)
        ENDIF
100    CONTINUE
        ECORR(ISPIN)=ESPIN
       IF(IPRINT.NE.0)WRITE(*,80)SPCASE(ISPIN),ESPIN
80     FORMAT('The ',A2,' contribution to the correlation ',
     &        'energy is: ',F12.7,' a.u.')
10    CONTINUE
C
C  THE TAU(Ij,Ab) <Ij//Ab> CONTRIBUTION TO THE ENERGY (SPIN CASE AB)
C
      ESPIN=0.D0
      DO 200 IRREP=1,NIRREP
       DISSYT=IRPDPD(IRREP,ISYTYP(1,46))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,46))
       IF(MIN(NUMSYT,DISSYT).NE.0) THEN
cYAU        I001=1
        I001=IFREE
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYT*DISSYT
        I004=I003+NUMSYT
        IF(I004.LT.MXCOR) THEN
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &              ICORE(I002),ICORE(I0T(1)),ICORE(I0T(2)),3,TAU,
     &              IRREP,POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &              ICORE(I003))
        ELSE
         CALL INSMEM('CMPENG',I004,MXCOR)
        ENDIF
        ETOTT2=ETOTT2+ET2
        ETOT=ETOT+E
        ESPIN=ESPIN+E
        EDBLS_KEEP = EDBLS_KEEP + ET2
        ESING_KEEP = ESING_KEEP + (E - ET2)
        
       ENDIF
200   CONTINUE
      ECORR(3)=ESPIN
      IF(IPRINT.NE.0)WRITE(*,80)SPCASE(3),ESPIN
      IF(IPRINT.NE.0)WRITE(*,81)ETOT
81    FORMAT('The total correlation energy is ',F15.12,' a.u.')
cjp
cjp add vacuum energy of (iref) and store the result in diagonal element
cjp of the effective hamiltonian
cjp
      if (isbwcc) heff(iref,iref)=etot+hdiagcontr(iref)
      Write(6,*) 
      Write(6, "(a,F12.7,a)") 
     & "The singles contribution to corrrelation energy", ESING_KEEP, 
     & " a.u."
      Write(6, "(a,F12.7,a)") 
     & "The Doubles contribution to corrrelation energy", EDBLS_KEEP,
     & " a.u."
C
      RETURN
      END
