#include "flags.h"
      SUBROUTINE FORMT2(W,DENOM,IUHF)
C
C THIS ROUTINE FORMS THE FIRST ORDER T2 VECTOR IN SYMMETRY 
C  PACKED FORM AND WRITES IT OUT.  ALSO COMPUTE THE CONTRIBUTION
C  TO THE CORRELATION ENERGY.
C
C      ARRAYS : W - USED TO HOLD INTEGRAL ARRAYS.  MUST
C                   BE DIMENSIONED TO LARGEST IRREP OF THIS
C                   ARRAY.
C
C               DENOM - USED FOR DENOMINATOR ARRAYS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W(1),DENOM(1),EAA,EAB,ENERGY,SDOT,ETAA,ESCF
      DOUBLE PRECISION ESING
      LOGICAL NONHF
      LOGICAL POLSCF
      CHARACTER*2 SPCASE(3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /SHIFT/  ISHIFT,NDRGEO
cjp
#include "bwcc.com"
      DATA SPCASE /'AA','BB','AB'/
      ETAA=0.0
cjp avoid printing of senseless numbers in the case of bwcc
      if(isbwcc) write(6,*)
     &  '@FORMT2-I BW-CC: pseudo MBPT(2) energies not printed'
      if(.not.isbwcc) WRITE(6,1010)
1010  FORMAT(/,T3,'@FORMT2-I, Second-order MBPT correlation energies:')
cjp
      NONHF=(IFLAGS(38)+IFLAGS(77).NE.0).or.isbwcc
      IF(NONHF)THEN
       if(.not.isbwcc) WRITE(6,1011)
1011   FORMAT(T3,'@FORMT2-I, Singles contribution will be calculated.')
       IF(IFLAGS(39).EQ.0)THEN
        if(.not.isbwcc) WRITE(6,1012)
1012    FORMAT(T3,'@FORMT2-W, MBPT(2) energies are correct only ',
     &           'for semicanonical orbitals.')
       ENDIF
      ENDIF 
      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
      if(.not.isbwcc) WRITE(6,1005)
      if(.not.isbwcc) WRITE(6,1004)ESCF
cjp here initial t value is computed as integrals.1/denominators
cjp listd .. .denoms; listt ... amplitudes; listw ... integrals
cjp and they are immediatelly used to comp. initial correlation energy
      DO 5 ISPIN=1,1+IUHF
       LISTW=13+ISPIN + ISHIFT
       LISTD=47+ISPIN
       LISTT=43+ISPIN
       EAA=0.0 
       DO 10 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        CALL GETLST(W,1,NUMDIS,2,IRREP,LISTW)
        CALL GETLST(DENOM,1,NUMDIS,1,IRREP,LISTD)
#ifdef _DEBUG_LVLM
       write(6,"(a)") " The AAAA and BBBB denominators and W"
       CALL CHECKSUM("@-INTPRC:W", W, NUMDIS*DISSIZ)
       call output(w,1,DISSIZ,1,NUMDIS,DISSIZ,NUMDIS,1)
#endif 
        DO 20 IELEM=1,DISSIZ*NUMDIS
#ifdef _DEBUG_LVLM
        Write(*,"(2(1x,F15.10))") W(IELEM), DENOM(IELEM)
#endif 
         DENOM(IELEM)=DENOM(IELEM)*W(IELEM)
20      CONTINUE
        CALL PUTLST(DENOM,1,NUMDIS,1,IRREP,LISTT)
#ifdef _DEBUG_LVL0
        CALL CHECKSUM("@-INTPRC:T2", DENOM, NUMDIS*DISSIZ)
#endif 
        EAA=EAA+SDOT(NUMDIS*DISSIZ,W,1,DENOM,1)
10     CONTINUE
       if(.not.isbwcc) WRITE(LUOUT,1000)SPCASE(ISPIN),EAA
       IF(IUHF.EQ.0)EAA=2.0*EAA
       ETAA=ETAA+EAA
5     CONTINUE
      LISTW=16 + ISHIFT
      LISTD=50
      LISTT=46
      EAB=0.0
      DO 110 IRREP=1,NIRREP
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
cjp order of integrals and denominators same, since the integrals were sorted
cjp wrto the true fermi vacuum and the denoms were constructed so that
cjp the order corresponds to the same f vacuum
cjp check it anyway in the dijab and in the integral symmetry packing routine
       CALL GETLST(W,1,NUMDIS,2,IRREP,LISTW)
       CALL GETLST(DENOM,1,NUMDIS,1,IRREP,LISTD)
#ifdef _DEBUG_LVLM
CSSS       write(6,"(a)") " The ABAB denominators"
CSSS       call output(Denom,1,DISSIZ,1,NUMDIS,DISSIZ,NUMDIS,1)
       CALL CHECKSUM("@-INTPRC:D-ABAB ", DENOM, NUMDIS*DISSIZ)
CSSS       write(6,"(a)") " The integrals"
CSSS       call output(w,1,DISSIZ,1,NUMDIS,DISSIZ,NUMDIS,1)
       CALL CHECKSUM("@-INTPRC:W-ABAB ", W, NUMDIS*DISSIZ)
#endif
       DO 120 IELEM=1,DISSIZ*NUMDIS
        DENOM(IELEM)=DENOM(IELEM)*W(IELEM)
120    CONTINUE
       CALL PUTLST(DENOM,1,NUMDIS,1,IRREP,LISTT)
#ifdef _DEBUG_LVL0
       write(6,"(a)") " The T2-ABAB"
CSSS       call output(Denom,1,DISSIZ,1,NUMDIS,DISSIZ,NUMDIS,1)
       CALL CHECKSUM("@-INTPRC:T2-ABAB", DENOM, NUMDIS*DISSIZ)
#endif
       EAB=EAB+SDOT(NUMDIS*DISSIZ,W,1,DENOM,1)
110    CONTINUE
       if(.not.isbwcc) WRITE(LUOUT,1000)SPCASE(3),EAB
       ENERGY=EAB+ETAA
C
C CALCULATE SINGLES CONTRIBUTION TO SECOND-ORDER ENERGY
C  FOR SEMICANONICAL ORBITALS
C
       ESING=0.0
       IF(NONHF)THEN
        DO 130 ISPIN=1,1+IUHF
         I000=1
         I010=I000+NT(ISPIN)
         I020=I010+NT(ISPIN)
         CALL GETLST(W(I000),1,1,1,ISPIN,90)
#ifdef _DEBUG_LVLM1
         Write(6,*) "T1 amps"
         call output(W(I000), 1, NT(ispin), 1, 1, NT(ispin), 1)
#endif
         CALL GETLST(W(I010),1,1,1,ISPIN+2,93)
#ifdef _DEBUG_LVLM1
         Write(6,*) "FAI amps"
         call output(W(I010), 1, NT(ispin), 1, 1, NT(ispin), 1)
#endif
         ESING=ESING+SDOT(NT(ISPIN),W(I000),1,W(I010),1)
         If ((Iflags(h_IFLAGS_calc) .EQ. 45) .OR. 
     &       (Iflags(h_IFLAGS_calc) .EQ. 46)) 
     &       ESING = 2.0D0*ESING 
   
130     CONTINUE
        if(.not.isbwcc) WRITE(LUOUT,1003)ESING
       ENDIF
C
       if(.not.isbwcc) WRITE(LUOUT,1002)ENERGY
       if(.not.isbwcc) WRITE(LUOUT,1001)ENERGY+ESCF+ESING
C
C WRITE TOTAL ENERGY TO JOBARC (THIS IS AT LEAST NEEDED FOR
C OPTIMIZATIONS USING NUMERICALLY EVALUATED GRADIENTS
C
CJDW   8/30/95
C
C      Do not write out the second-order energy in scf analytical hessian
C      polyrate runs.
C Not just for polyrate, but forall analytic SCF Hessian calcualtions
C 12/09, Ajith Perera
C
       POLSCF = .FALSE.
       POLSCF = IFLAGS(2).EQ.0 .AND. IFLAGS(3).EQ.2 
cjp
cjp do this only for reference no. 1; 
       IF(POLSCF.or.isbwcc.and.nref.gt.1)THEN
        if(.not.isbwcc) WRITE(6,1020)
       ELSE
        CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ENERGY+ESCF+ESING)
       ENDIF
cjp       if(.not.isbwcc) print *,'in formt2 energies', ENERGY+ESCF+ESING,ENERGY,ESCF,ESING
C
       WRITE(6,1005)
1000   FORMAT(T15,'E2(',A2,')',T37,'=',F18.12,' a.u.')
1001   FORMAT(T15,'Total MBPT(2) energy',T37,'=',F18.12,' a.u.')
1002   FORMAT(T15,'E2(TOT)',T37,'=',F18.12,' a.u.')
1003   FORMAT(T15,'E2(SINGLE) ',T37,'=',F18.12,' a.u.')
1004   FORMAT(T15,'E(SCF)',T37,'=',F18.12,' a.u.')
1005   FORMAT(T15,46('-'))
cjp
1020   FORMAT(T15,' @FORMT2-I, MR-BWCC or SCF Hessian Polyrate Run.',
     &            ' MBPT(2) energy not written. ')
       RETURN
       END
