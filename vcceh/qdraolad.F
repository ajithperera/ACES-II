      SUBROUTINE QDRAOLAD(ICORE,MAXCOR,IUHF,TAU,TAU2,IRREPX,LISTFLAG,
     &                    LSTMO,LSTMOINC,LSTAO,LSTAOINC)
C 
C Driver for the AO basis algorithm to evaluate contributions from
C Hbar(ab,ef) to quadratic contribution. This is identical to the
C draolad using dt2int2 (in hbarxc) (renamed only for book keeping
C purposes).
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      LOGICAL TAU,TAU2

      COMMON /FLAGS/ IFLAGS(100)
C
      CALL AOLADLST(IUHF,LISTFLAG,IRREPX)
      CALL T2TOAO(ICORE,MAXCOR,IUHF,TAU,TAU2,LSTMO,LSTAO,IRREPX)
      IF(IFLAGS(95).EQ.1)THEN
       CALL AOLAD2(ICORE,MAXCOR,IUHF,TAU,IRREPX,LSTAO,LSTAOINC)
      ELSE
       CALL AOLAD3(ICORE,MAXCOR,IUHF,TAU,IRREPX,LSTAO,LSTAOINC)
      ENDIF
      CALL QZ2TOMO(ICORE,MAXCOR,IUHF,TAU,IRREPX,LSTMO,LSTMOINC,
     &             LSTAOINC,.FALSE.,1)

      RETURN
      END
