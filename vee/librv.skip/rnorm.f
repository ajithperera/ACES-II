      SUBROUTINE RNORM(IRREPX,NSIZEC,ROOT,SR,TMP,MAXCOR,IUHF,SR0)
C
C THIS ROUTINE PROPERLY NORMALIZES THE RIGHT-HAND EOM-CC EIGENVECTOR
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LAMTHERE, PRINT
      PARAMETER (TOL=1.D-5)
      DIMENSION SR(*),TMP(MAXCOR)
      COMMON/LAMSTATE/LAMTHERE
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
C
      PRINT = IFLAGS(1) .GE. 2
      IF (PRINT) WRITE(6,1000)
1000  FORMAT(T3,'@RNORM-I, Processing right-hand wavefunction.')
C
C THE FOLLOWING CODE ONLY RUNS FOR IRREPX=1
C
      IF(IRREPX.EQ.1)THEN
C
C CALCULATE CONTRIBUTION FROM REFERENCE FUNCTION (S0)
C
C           -1 
C     S0 = W   [ F(ai) C(ai) + W(ab,ij) C(ab,ij) ]
C
C
       CALL LOADVEC1(IRREPX,TMP,MAXCOR,IUHF,93,0,13,
     &               LENTOT,IUHF.EQ.0)
       SR0=SDOT(LENTOT,TMP,1,SR,1)/ROOT
C
C CALCULATE OVERLAP WITH LAMBDA VECTOR [SHOULD BE ZERO]
C
       CALL LOADVEC1(IRREPX,TMP,MAXCOR,IUHF,190,0,143,
     &               LENTOT,IUHF.EQ.0)
       X=SDOT(LENTOT,TMP,1,SR,1)+SR0
       IF(LAMTHERE)THEN
        IF (PRINT) WRITE(6,1001)X
1001    FORMAT(T3,'Overlap with ground bra state : ',F20.10,'.')
        IF(X.GT.TOL)THEN
         IF (PRINT) WRITE(6,1002)
1002     FORMAT(T3,'@RNORM-W, Biorthogonality not satisfied. ',
     &             'Something is wrong!')
        ENDIF
       ENDIF
       Z=ONE+SR0*SR0
       FACT=ONE/DSQRT(Z)
       CALL SSCAL(NSIZEC,FACT,SR,1)
       SR0=SR0*FACT
      ELSE
       SR0=ZILCH
       FACT=ONE
      ENDIF
C
C NOW WRITE THE CONVERGED VECTOR TO THE LISTS AND OBTAIN RESORTED
C ORDERING
C
      CALL PUTLST(SR,1,1,1,1,490)
      IF(IUHF.NE.0)THEN
       IS1AA=1
       IS1BB=IS1AA+IRPDPD(IRREPX,9)
       IS2AB=IS1BB+IRPDPD(IRREPX,10)
       IS2BB=IS2AB+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
       IS2AA=IS2BB+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
       CALL PUTLST(SR(IS1BB),1,1,1,2,490)
       CALL PUTALL(SR(IS2AA),1,IRREPX,444)
       CALL PUTALL(SR(IS2BB),1,IRREPX,445)
      ELSE
       IS2AB=1+IRPDPD(IRREPX,9)
      ENDIF
      CALL PUTALL(SR(IS2AB),1,IRREPX,446)
      CALL RESORT(TMP,MAXCOR,IUHF,IRREPX,444,434)
C
      CALL SSCAL(NSIZEC,ONE/FACT,SR,1)
C
      RETURN
      END
