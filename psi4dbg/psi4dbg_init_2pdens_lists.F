      SUBROUTINE PSi4DBG_INIT_2PDENS_LISTS(ICORE,MAXCOR,IUHF)
C
      IMPLICIT INTEGER(A-Z)
      DIMENSION ICORE(MAXCOR)
C
      IMODE=0

C PPPP lists

      CALL INIPCK(1,13,13,133,IMODE,0,1)

C HHHH lists

      CALL INIPCK(1,14,14,113,IMODE,0,1)

C                          -  - 
C HPHP lists (ph|ph,ph|ph,ph|ph)

      CALL INIPCK(1,9,9,123,IMODE,0,1)
      CALL INIPCK(1,9,10,118,IMODE,0,1)
      CALL INIPCK(1,11,11,125,IMODE,0,1)

C H4 lists to avoid overwriting Hbar(mb,ej) lists. Also, note that list 154 is for 
c temporarily keep H4(aj,bi) contributions formed by pccd_quickaa2. 

      CALL INIPCK(1,9,9,154,IMODE,0,1)
      CALL INIPCK(1,9,10,156,IMODE,0,1)
      CALL INIPCK(1,11,11,158,IMODE,0,1)

C HPHH and PPPH lists,list 10 and 30 already initilized 
C for HHHP and PHHPP (noticed that HPHH stored as HHHP)

      CALL INIPCK(1,14,18,110,IMODE,0,1)
      CALL INIPCK(1,13,11,130,IMODE,0,1)

C PPHH lists 
      CALL INIPCK(1,13,14,116,IMODE,0,1)
      CALL INIPCK(1,13,14,216,IMODE,0,1)

      RETURN
      END
