      SUBROUTINE RCC_RNGPRD(W,T,Q,SCR,MAXCOR,DISSYW,NUMSYW,DISSYT,
     &                      NUMSYT,DISSYQ,NUMSYQ,LISTW,LISTT,IRREP,
     &                      CALMOD,SPINAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT
      LOGICAL SPINAD,bRedundant
      CHARACTER*3 CALMOD
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &        CC2,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      INTEGER DISSYT,DISSYW,DISSYQ
      DIMENSION W(DISSYW,NUMSYW),T(DISSYT,NUMSYT),Q(DISSYQ,NUMSYQ)
      DIMENSION SCR(MAXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(6)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS2/ IFLAGS2(500)
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA ONEM/-1.0/
      DATA HALFM,HALF/-0.5,0.50D0/
      DATA ZILCH/0.0/

      bRedundant = IFLAGS2(155).EQ.0
C
C PICK UP THE RELEVANT T2 AND W PIECES, AND SPIN ADAPT IF NEEDED.
C
      IF(bRedundant) THEN
        IF (CC2) THEN
            CALL ZERO(Q,NUMSYQ*DISSYQ)
            RETURN
        ELSE
            CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT)
        ENDIF
      ELSE
        CALL GETLST_NR(T,SCR,MAXCOR*IINTFP,LISTT, IRREP)
      ENDIF
C
C Spin adaptation is false for drCCD and or rCCD. Nevertheless 
C scale by two to keep consitency with the rest of then code.

      IF(SPINAD)THEN
        IF(bRedundant) THEN
          CALL GETLST(W,1,NUMSYT,1,IRREP,34)
        ELSE
          CALL GETLST_NR(W,SCR,MAXCOR*IINTFP,34, IRREP)
        ENDIF    
       CALL SAXPY(NUMSYT*DISSYT,ONEM,W,1,T,1)
       LISTW=56
      ELSE 
       CALL DSCAL(NUMSYQ*DISSYQ,TWO,T,1)
      ENDIF 
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)
      IF(MIN(DISSYT,NUMSYT,NUMSYW).NE.0)THEN
       CALL XGEMM('N','N',DISSYT,NUMSYW,NUMSYT,ONE,T,DISSYT,
     &            W,DISSYW,ZILCH,Q,DISSYT)
c       CALL MXM(T,DISSYT,W,NUMSYT,Q,NUMSYW)
      ELSE
       CALL ZERO(Q,NUMSYQ*DISSYQ)
      ENDIF

C Multiply by half so that the overall contribution is not twice. 
C Note that in rcc_dwbej and above there is a extra factor of two. 
C It is not clear to me (yet) what this negative sign do here 
C when the symmetry adaptation is done. 
 
      IF(SPINAD) THEN 
        CALL SSCAL(NUMSYQ*DISSYQ,HALFM,Q,1)
      ELSE
       CALL SSCAL(NUMSYQ*DISSYQ,HALF,Q,1)
      ENDIF 

      RETURN
      END
