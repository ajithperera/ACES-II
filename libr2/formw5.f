      SUBROUTINE FORMW5(ICORE,MAXCOR,IUHF,TERM1,TERM2,TERM3,TERM4,
     &                  TERM5,TERM6,TAU,CCSD)
C
C     This subroutine calculates the ABCI intermediate, referred to in
C     some places as W5 :
C
C     W(abci) = <ab||ci> + t(d,i) <ab||cd> + F(m,c) t(ab,im)
C
C             - P(ab) t(a,m) {\bar {\bar W}}(mb,ci) 
C
C             + P(ab) <am||ce> (t(be,im) - t(b,m) t(e,i))
C
C             - 0.5 W(mn,ci) tau(ab,mn)
C
C     The target lists are 127-130 (130 only for RHF). This routine
C     anticipates the usual MOINTS lists. In particular, lists 7-10
C     must contain plain integrals. W(mn,ci) is generated on lists 7-10
C     in FORMW5 and integrals are recovered at the end. {\bar {\bar W}}}
C     must have been prepared by the caller. W5SCRFIL generates DERGAM
C     and creates list 401. DERGAM is removed by KILLDER at the end of
C     FORMW5.
C
C     This routine now reverts essentially to its original form. In
C     particular, MODT2 is called in CCSDT-3, CCSDT-4, and CCSDT. It
C     is not called in CCSDT-2, of course, since there are no T1 terms
C     in the T3 equation. This leads to potential problems if we want
C     to do CCSDT-2 gradients and want to call FORMW5 to form the CCSD
C     ABCI Hbar elements. We may yet have to have an argument which
C     distinguishes the different reasons why we call FORMW5, i.e. to
C     form intermediates to go in the T3 equation or to form intermed-
C     iates for lambda/vcceh/vee, etc. In addition, we do not any longer
C     use FORMW5 to drive the calculation of the term containing ABCD
C     INTEGRALS. Note how TERM2 is set to .FALSE. below. Basically,
C     FORMW5 now is supposed to compute everything except the <ab||cd>
C     containing term. This term is handled separately by, for example,
C     W5T1ABCD (and later more general routines).
C
C     FORMW5 is called by : TRPINT, VLAMCC.
C
C     CCSD is logical variable which determines whether this routine
C     has been called to compute full abci intermediate regardless of
C     method (i.e. it was called from VLAMCC) or whether it pays att-
C     ention to underlying triples method (i.e. call is from TRPINT).
C     Most significance is for CC3 method, and to some extent CCSDT-2.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL TERM1,TERM2,TERM3,TERM4,TERM5,TERM6,TAU,CCSDT2,CCSD
      LOGICAL bRedundant, CC2, HBAR_4LCCSD, NONHF
      DIMENSION ICORE(MAXCOR)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/IFLAGS2(500)
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      CCSDT2 = METHOD.EQ.15 .AND. .NOT.CCSD
      bRedundant = IFLAGS2(155).EQ.0
      CC2        = IFLAGS(2) .EQ. 47
C
C FOR CCSD SECOND DERIVATIVES, STORE W(AB,CI) ON LISTS 227 TO 230
C
      IOFFLIST=0
      IF(IFLAGS(3).EQ.2) IOFFLIST=100
C
C ZERO OUT THESE LISTS BECAUSE ROUTINES USE FANCYPUT 
C
      CALL ACES_COM_SYMLOC
      CALL ZEROLIST(ICORE,MAXCOR,130+IOFFLIST)
      TERM2=.FALSE.
      IF(IFLAGS(3).EQ.2) TERM2=.TRUE.
      IF(CC2) TERM2= .TRUE.
      HBAR_4LCCSD  = .FALSE. 
      NONHF        = .FALSE.
      HBAR_4LCCSD  = (IFLAGS(2) .EQ. 6 .OR. IFLAGS2(117) .EQ. 7)
      NONHF        = (IFLAGS(38) .NE. 0)
C
C REPLACE IAJK INTEGRALS WITH CORRESPONDING HBAR ELEMENTS
C FOR EVERYTHING EXCEPT ROHF-MBPT AND CCSDT-2
C
      IF(IFLAGS(2).LE.4.AND.IFLAGS(11).EQ.2)THEN
       CONTINUE
      ELSE
       IF(.NOT.CCSDT2 .AND. .NOT. HBAR_4LCCSD) 
     &     CALL MODIAJK(ICORE,MAXCOR,IUHF,ONE)
      ENDIF
C
C EVALUATE THE CONTRIBUTION FROM HELL
C
 
      IF(TERM1.OR.TERM5)THEN

       IF(.NOT.CCSDT2 .AND. .NOT. CC2 .AND. .NOT. HBAR_4LCCSD) THEN
           CALL MODT2(ICORE,MAXCOR,IUHF,ONE,ONE,bRedundant)
       ENDIF

       IF (CC2) THEN
           CALL MODT1T1(ICORE,MAXCOR,IUHF,ONE)
       ENDIF 

       IF(IUHF.EQ.0)THEN
        CALL W5INRHF(ICORE,MAXCOR,IUHF,TERM1,TERM5,IOFFLIST)
       
       ELSE
        CALL ZEROLIST(ICORE,MAXCOR,127+IOFFLIST)
        CALL ZEROLIST(ICORE,MAXCOR,128+IOFFLIST)
        CALL ZEROLIST(ICORE,MAXCOR,129+IOFFLIST)
        IF(TERM1.OR.TERM5)THEN
         CALL W5SCRFIL
         CALL W5INuhf1(ICORE,MAXCOR,IUHF,TERM1,TERM5,IOFFLIST)
         CALL W5INuhf2(ICORE,MAXCOR,IUHF,TERM1,TERM5,IOFFLIST)
         CALL W5INuhf3(ICORE,MAXCOR,IUHF,TERM1,TERM5,IOFFLIST)
        ENDIF

       ENDIF

       IF(.NOT.CCSDT2 .AND. .NOT. CC2 .AND. .NOT. HBAR_4LCCSD) THEN
          CALL MODT2(ICORE,MAXCOR,IUHF,ONEM,ZILCH,bRedundant)
       ENDIF

       IF (CC2) THEN
          CALL RNABIJ(ICORE,MAXCOR,IUHF,'T')
       ENDIF
     
      ENDIF
      
C              
C CALCULATE RING CONTRIBUTION
C
      IF(TERM4)CALL W5RING(ICORE,MAXCOR,IUHF,IOFFLIST)
C
C CALCULATE REMAINING PIECES
C
      IF(TERM2.OR.TERM3.OR.TERM6)THEN
       IF(TERM6)CALL W5TAU(ICORE,MAXCOR,IUHF,TAU,IOFFLIST,CCSD)
       IF(TERM3)CALL W5F(ICORE,MAXCOR,IUHF,IOFFLIST,HBAR_4LCCSD,NONHF)
       IF(IUHF.NE.0.AND.TERM2)THEN
        CALL W5AB2(ICORE,MAXCOR,IUHF,TERM3,TERM2,TERM6,TAU,IOFFLIST)
        CALL W5AB1(ICORE,MAXCOR,IUHF,TERM2,TERM3,TERM6,TAU,IOFFLIST)
        CALL W5AA1(ICORE,MAXCOR,IUHF,TERM2,TERM3,TERM6,TAU,IOFFLIST)
       ELSE
        IF(TERM2)CALL W5AB2RHF(ICORE,MAXCOR,IUHF,TERM3,TERM2,TERM6,
     &                         TAU,IOFFLIST,HBAR_4LCCSD,NONHF)
       ENDIF
      ENDIF
C
C PUT IAJK INTEGRALS BACK ON DISK
C FOR EVERYTHING EXCEPT ROHF-MBPT
C
      IF(IFLAGS(2).LE.4.AND.IFLAGS(11).EQ.2)THEN
       CONTINUE
      ELSE
       IF(.NOT.CCSDT2 .AND. .NOT. HBAR_4LCCSD) 
     &     CALL MODIAJK(ICORE,MAXCOR,IUHF,ONEM)
      ENDIF
C
      CALL ACES_IO_REMOVE(54,'DERGAM')
C
      RETURN
      END
