       SUBROUTINE NONRED_INIT(IUHF)
       IMPLICIT NONE
       INTEGER IUHF

       CALL INIT_NOREDNT_LIST(1,9,9,34)
       CALL INIT_NOREDNT_LIST(1,9,10,37)
       CALL INIT_NOREDNT_LIST(1,11,12,39)
       IF(IUHF.NE.0) THEN
          CALL INIT_NOREDNT_LIST(1,10,9,36)
          CALL INIT_NOREDNT_LIST(1,12,11,38)
          CALL INIT_NOREDNT_LIST(1,10,10,35)
       ENDIF
       RETURN
       END
