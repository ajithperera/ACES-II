      SUBROUTINE NOMEM(CNAME,TRACE,MEMREQ,MEMHAVE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*(*) CNAME,TRACE
      CHARACTER*40 JNAME
      CHARACTER*60 TNAME
      COMMON /FILES/ LUOUT,MOINTS
C
      iln=LEN(CNAME)
      IF(iln.GT.40) iln=40
      JNAME=CNAME(1:iln)
      iln=LEN(TRACE)
      IF(iln.GT.60) iln=60
      TNAME=TRACE(1:iln)
      WRITE(LUOUT,9000)JNAME
      WRITE(LUOUT,9005)TNAME
      WRITE(LUOUT,9010)MEMREQ,MEMHAVE
 9000 FORMAT(T3,'@NOMEM-F, Insufficient memory:  ',A40)
 9005 FORMAT(/,T8,'Traceback:  ',A60)
 9010 FORMAT(/,T8,' Memory required: ',I10,/,
     &       T8,'Memory available: ',I10)
      CALL ERREX
      RETURN
      END
