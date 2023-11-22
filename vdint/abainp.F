      SUBROUTINE ABAINP(WORK1,LWORK1,NUCIND1)
C
C READS INPUT FOR ABACUS
C
CEND
C
C  ADAPTED TO CRAPS AND IMPLEMENTED ON CRAY 
C  OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "mxcent.par"
      PARAMETER (LUCMD = 5, LUPRI = 6)
C
      LOGICAL SCF,NONHF,GRAD,INTWRIT
      LOGICAL DCORD,DCORGD,NOORBT,DOPERT
      LOGICAL MOLGRD,MOLHES,DIPDER,POLAR,INPTES,
     *        VIB,RESTAR,DOWALK,GDALL,FOCK,H2MO,FOCK1
C
      CHARACTER NAMEX*6
      CHARACTER WORD*7
C
      DIMENSION WORK1(LWORK1)
C
      COMMON /IPRI/ IPRINT
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK1,
     *                H2MO
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/OPTION/GRAD,INTWRIT,FOCK
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
C
      COMMON/FLAGS/IFLAGS(100)
C
C     Initialize /ABAINF/
C
      IPRDEF=1
      MOLGRD=.TRUE.
C
C DEFAULT IS GRADIENT CALCULATAION ONLY
C
      MOLHES=.FALSE.
      DIPDER=.FALSE.
      POLAR=.FALSE.
      INPTES=.FALSE.
      VIB=.FALSE.
      RESTAR=.FALSE.
      DOWALK=.FALSE.
      GDALL=.FALSE.
      FOCK1=FOCK
      H2MO=.FALSE.
C
C     Initialize /NUCLEI/
C
      nrotco = 0
      NTRACO = 0
      DO 50 I = 1, 3
       IROTCO(I) = 0
       ITRACO(I) = 0
   50 CONTINUE
      DO 60 I = 1, 3*MXCENT
       DOPERT(I) = .TRUE.
   60 CONTINUE
      IPRDEF=1 
C
C IF THE DERIVATIVE LEVEL IS GREATER THAN TWO, EVALUATE ANALYTICAL
C SECOND DERIVAITIVES
C
      IF(IFLAGS(3).GE.2) THEN
         IF(IFLAGS(2).EQ.0) Then 
C
CSSS .AND. IFLAGS(54) .EQ. 1) THEN
C If analytical hessian is turned on the MOLHES can be set to true. 
C The currently limited to SCF (IFLAGS(2) .EQ. 0), but the test for
C VIB=EXACT is redundent and become problematic when VIB=EXACT,GEOM_OPT
C is combined. Ajith Perera, 04/2006.
C
        MOLHES=.TRUE.
        DIPDER=.TRUE.
        POLAR=.TRUE.
       ELSE
        IF(GRAD) THEN
         MOLHES=.FALSE.
         DIPDER=.TRUE.
         POLAR=.TRUE.
        ENDIF
       ENDIF
      ENDIF
C
      WORD='*END OF'
      CALL REAINP(WORD)
      CALL NUCINP(WORD)
      CALL ONEINP(WORD)
      CALL TWOINP(WORD)
      CALL RHSINP(WORD)
      CALL TSTINP(WORD)
      CALL ORTINP(WORD)
      CALL RESINP(WORD)
c      CALL TROINP(WORD)
      CALL DIPINP(WORD)
      CALL RLXINP(WORD)
      CALL CRSINP(WORD)
      CALL CZRINP(WORD)
      CALL CTRINP(WORD)
C
C     Read geometry, orbital spec., etc. (MOLECULE format)
C
C     Logical argument indicates whether information on LUONEL should be
C     written or not.
C
      CALL READIN(.FALSE.,WORK1,LWORK1)
C
      NUCIND1=NUCIND
C
      RETURN
C
      STOP 
      END
