      SUBROUTINE INIGAM(IUHF,icore,maxcor)
C
C THIS ROUTINE INITIALIZES THE LISTS
C FOR THE TWO-PARTICLE DENSITY MATRIX
C
CEND
C
C  CODED SEPTEMBER/93 JG
C
      IMPLICIT INTEGER(A-Z)
      LOGICAL GABCD
      COMMON/ABCD/GABCD
      integer maxcor,icore(maxcor),imode
C G1 Lists
      IMODE=0
      CALL INIPCK(1,1,3,114,IMODE,0,1)
      CALL INIPCK(1,13,14,116,IMODE,0,1)
      IF(IUHF.EQ.1) THEN
       CALL INIPCK(1,2,4,115,IMODE,0,1)
      ENDIF

C
C  G2 LISTS
C
C
      IF(.NOT.GABCD)THEN
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
      CALL INIPCK(1,9,9,123,IMODE,0,1)
      CALL INIPCK(1,9,10,118,IMODE,0,1)
      CALL INIPCK(1,11,11,125,IMODE,0,1)
      IF(IUHF.NE.0) THEN
       CALL INIPCK(1,10,10,124,IMODE,0,1)
       CALL INIPCK(1,10,9,117,IMODE,0,1)
       CALL INIPCK(1,12,12,126,IMODE,0,1)
      ENDIF
C
C  G5 AND G6 LISTS
C
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

C Zero out the lists
      call zerolist(icore,maxcor,133)
      call zerolist(icore,maxcor,113)
      call zerolist(icore,maxcor,123)
      call zerolist(icore,maxcor,118)
      call zerolist(icore,maxcor,125)
      call zerolist(icore,maxcor,110)
      call zerolist(icore,maxcor,130)
      call zerolist(icore,maxcor,114)
      call zerolist(icore,maxcor,116)
      if (iuhf .ne. 0) then
        call zerolist(icore,maxcor,131)
        call zerolist(icore,maxcor,132)
        call zerolist(icore,maxcor,111)
        call zerolist(icore,maxcor,112)
        call zerolist(icore,maxcor,124)
        call zerolist(icore,maxcor,117)
        call zerolist(icore,maxcor,126)
        call zerolist(icore,maxcor,107)
        call zerolist(icore,maxcor,108)
        call zerolist(icore,maxcor,109)
        call zerolist(icore,maxcor,129)
        call zerolist(icore,maxcor,128)
        call zerolist(icore,maxcor,127)
        call zerolist(icore,maxcor,115)
      endif
C
C  ALL DONE
C
      RETURN
      END
