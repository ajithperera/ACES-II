      SUBROUTINE TPDENS(ICORE,MAXCOR,IUHF,R0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GABCD
      DIMENSION ICORE(MAXCOR)
C
C DRIVER ROUTINE FOR THE CALCULATION OF THE
C TWO-PARTICLE DENSITY MATRIX IN EOM-CCSD
C GRADIENT CALCULATIONS
C
CEND
C
C CODED SEPTEMBER/93 JG
C
      COMMON/STATSYM/IRREPX
      COMMON/ABCD/GABCD
C
      DATA ONE,ONEM,HALF /1.0D0,-1.0D0,0.5D0/

C
C INITIALIZE LISTS FOR TWO-PARTICLE-DENSITY ON GAMLAM
C
      CALL INIGAM(IUHF)
C
C RESTORE PRISTINE ABCI AND IJKA INTEGRAL LISTS
C
      CALL MODAIBC(ICORE,MAXCOR,IUHF,ONEM)
      CALL MODIAJK(ICORE,MAXCOR,IUHF,ONEM)
C
C G(IJ,KL) 
C
      CALL TPDIJKL(ICORE,MAXCOR,IUHF)
C
C G(AB,CD)
C
      IF(GABCD)CALL TPDABCD(ICORE,MAXCOR,IUHF)
C
C G(AI,BC) AND G(AB,CI)
C
      CALL TPDAIBC(ICORE,MAXCOR,IUHF,HALF)
      CALL TPDABCI(ICORE,MAXCOR,IUHF)
C
C G(IA,JK) AND G(IJ,KA)
C
      CALL TPDIAJK(ICORE,MAXCOR,IUHF,HALF)
      CALL TPDIJKA(ICORE,MAXCOR,IUHF)
C
C G(AI,BJ)
C
      CALL TPDAIBJ(ICORE,MAXCOR,IUHF)
C
C G(IJ,AB)
C
      CALL TPDIJAB(ICORE,MAXCOR,IUHF)
C
C CALL ROUTINES TO CALCULATE G(ABCD) CONTRIBUTIONS TO
C I(AB) AND I(AI) IF G(ABCD) IS NOT STORED.
C
      IF(.NOT.GABCD)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@TPDENS-I, Gamma(abcd) contributions to ',
     &           'intermediates calculated directly.')
       CALL IABCONT(ICORE,MAXCOR,IUHF)
       CALL IAICONT(ICORE,MAXCOR,IUHF)
       IONE=1
       CALL PUTREC(20,'JOBARC','STATESYM',IONE,IRREPX)
      ENDIF
C
      IF(IRREPX.EQ.1)CALL NEWAMP(ICORE,MAXCOR,IUHF,R0)
      IF(IUHF.EQ.0)CALL ABTOAA(ICORE,MAXCOR,IUHF)
      CALL RNABIJ(ICORE,MAXCOR,IUHF,'L')
      CALL GFORMG(1,1,144,44,100,ICORE,MAXCOR,0,ONE,ONE,IUHF)
C
      IF(IUHF.EQ.1) THEN
       CALL NEWTYP2(1,56,9,10,.TRUE.)
       CALL NEWTYP2(1,57,10,9,.TRUE.)
      ENDIF
C
      IF(IUHF.NE.0.AND.GABCD)THEN
       CALL S2EXPECT(ICORE,MAXCOR,IUHF,'N',.TRUE.,2)
      ENDIF
      RETURN
      END
