      SUBROUTINE CRAPSO_C4
      LOGICAL YESNO,ISOPN
      CHARACTER*80 FNAME
      CALL GFNAME('MOINTS  ',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO,OPENED=ISOPN)
      IF(YESNO.AND.ISOPN)THEN
       CALL TERMIO_C4
      ENDIF
      CALL DUMPJA_C4('O')
      CALL CHECKOUT_C4
      RETURN
      END
