      SUBROUTINE SAVD(DOOA,DOOB,DVVA,DVVB,DVOA,DVOB,IUHF)
C
C THIS SUBROUTINE SAVES THE RELAXED DENSITY MATRIX ON THE
C GAMLAM FILE. THIS IS ONLY REQUIRED FOR SECOND DERIVATIVE
C CALCULATIONS.
C
C THE LISTS WRITTEN ARE:
C
C  DOOA    1,160
C  DOOB    2,160    UHF AND ROHF ONLY
C  DVVA    3,160
C  DVVB    4,160    UHF AND ROHF ONLY
C  DVOA    5,160  
C  DVOB    6,160    UHF AND ROHF ONLY
C
CEND
C
C CODED APRIL/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL YESNO, CIS,EOM, MRCC, RAMAN
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      INTEGER POP,VRT
      CHARACTER*80 FNAME
      DIMENSION DOOA(100),DOOB(1),DVVA(1),DVVB(1),DVOA(1),DVOB(1)
      COMMON /FLAGS2/ IFLAGS2(500) 
      COMMON/EXCITE/CIS,EOM
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
C For EOM gradients, list 160 was created in vee with three distributions.
C   It does not need to be made here.
C SG 1/96
C
      MRCC  = IFLAGS2(132) .EQ. 3 
      RAMAN = IFLAGS2(151) .EQ. 1
      IF (.NOT. (EOM .OR. MRCC)) THEN
c        CALL GFNAME('GAMLAM  ',FNAME,ILENGTH)
c        INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
c        IF (YESNO) THEN
          IENTER=0
          IOFF=0
c        ELSE
c          IENTER=1
c          IOFF=-1
c        ENDIF
C
C The original updmoi was called with number of distributions set 1
C and it was changed to 3 during the work that lead to pert. 
C density code. As far as the density code is concerned this is
C ineffective. But when both density and vcceh needs to run
C (for Raman intensities) this has to be set to 3, so no 
C inconsitencies arise. Ajith Perera, 12/08


        IF (RAMAN) THEN
           NOF_DIS = 3
        ELSE
           NOF_DIS = 1
        ENDIF
C
        DO 100 ISPIN=1,IUHF+1
C
          CALL UPDMOI(NOF_DIS,NF1(ISPIN),ISPIN,160,IENTER,IOFF)
          IENTER=0
          IOFF=0
          CALL UPDMOI(NOF_DIS,NF2(ISPIN),2+ISPIN,160,IENTER,IOFF)
          CALL UPDMOI(NOF_DIS,NT(ISPIN),4+ISPIN,160,IENTER,IOFF)
C
 100    CONTINUE
      ENDIF
C
      CALL PUTLST(DOOA,1,1,1,1,160)
      CALL PUTLST(DVVA,1,1,1,3,160)
      CALL PUTLST(DVOA,1,1,1,5,160)
C
      IF(IUHF.NE.0) THEN
C
       CALL PUTLST(DOOB,1,1,1,2,160)
       CALL PUTLST(DVVB,1,1,1,4,160)
       CALL PUTLST(DVOB,1,1,1,6,160)
C
      ENDIF
C
      RETURN
      END
