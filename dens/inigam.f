      SUBROUTINE INIGAM(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE INITIALIZES THE LISTS
C FOR THE GAMMA INTERMEDIATES 
C
CEND
C
C  CODED AUGUST/90
C
      IMPLICIT INTEGER(A-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2,
     &        GABCD,RELAXED,TRULY_NONHF
      LOGICAL TDA,EOM
      COMMON/EXCITE/TDA,EOM
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2,GABCD,RELAXED,TRULY_NONHF
C
C  FIRST ALL LISTS REQUIRED FOR MBPT3 AND CCD
C
C  G2 LISTS         
C
      IMODE=0
C
C FOR MBPT2 AND MBPT3 THE GAMLAM FILE HAS TO BE INITIALIZED
C
      IF(TDA)GOTO 3000
      IF (MBPT2.OR.(MBPT3.AND.(.NOT.ROHF))) THEN
         call aces_io_remove(51,'GAMLAM')
      ENDIF
C
      IF(.NOT.GABCD) THEN
       CALL INIPCK(1,13,13,133,IMODE,0,1)
       IF(IUHF.NE.0) THEN
        CALL INIPCK(1,1,1,131,IMODE,0,1)
        CALL INIPCK(1,2,2,132,IMODE,0,1)
       ENDIF
      ENDIF
C
C  G3 LISTS
C
      CALL INIPCK(1,14,14,113,IMODE,0,1)
      IF(IUHF.NE.0) THEN
       CALL INIPCK(1,3,3,111,IMODE,0,1)
       CALL INIPCK(1,4,4,112,IMODE,0,1)
      ENDIF
C
C  G4 LISTS
C
3000  CONTINUE
      CALL INIPCK(1,9,9,123,IMODE,0,1)
      CALL INIPCK(1,9,10,118,IMODE,0,1)
      CALL INIPCK(1,11,11,125,IMODE,0,1)
      IF(IUHF.NE.0) THEN
       CALL INIPCK(1,10,10,124,IMODE,0,1)
       CALL INIPCK(1,10,9,117,IMODE,0,1)
       CALL INIPCK(1,12,12,126,IMODE,0,1)
      ENDIF
      IF(TDA)RETURN
C 
C  G5 AND G6 LISTS
C 
C      IN CASE OF TRIPLE EXCITATIONS, THOSE
C      LISTS ARE ALREADY INITIALIZED,
C      SKIP HERE !
C
      IF((.NOT.(MBPT3.AND.(.NOT.ROHF))).AND.
     &       (.NOT.CCD).AND.(.NOT.TRIP1)) THEN
       CALL INIPCK(1,14,18,110,IMODE,0,1)
       IF(IUHF.EQ.1) THEN
        CALL INIPCK(1,3,16,107,IMODE,0,1)
        CALL INIPCK(1,4,17,108,IMODE,0,1)
        CALl INIPCK(1,14,11,109,IMODE,0,1)
       ENDIF
       CALL INIPCK(1,13,11,130,IMODE,0,1)
       IF(IUHF.EQ.1) THEN
        CALL INIPCK(1,1,9,127,IMODE,0,1)
        CALL INIPCK(1,2,10,128,IMODE,0,1)
        CALL INIPCK(1,13,18,129,IMODE,0,1)
       ENDIF
      ENDIF
C
C  G1 LISTS
C
      IF(.NOT.(MBPT3.AND.(.NOT.ROHF))) THEN
       IF((.NOT.TRIP1).OR.(M4SDTQ)) THEN 
        CALL INIPCK(1,1,3,114,IMODE,0,1)
        CALL INIPCK(1,13,14,116,IMODE,0,1)
        IF(IUHF.EQ.1) THEN
         CALL INIPCK(1,2,4,115,IMODE,0,1)
        ENDIF
       ENDIF
      ENDIF
C
C  ALL DONE
C
      RETURN
      END
