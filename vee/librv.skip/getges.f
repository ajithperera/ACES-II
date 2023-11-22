      SUBROUTINE GETGES(VEC, NSIZEC, IREAD, IRREPX, IJUNK)
C
C  A NEW GUESS VECTOR IS READ INTO VEC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DOUBLE, NONSTD
      DIMENSION VEC(NSIZEC),IJUNK(100)
C
      COMMON/GUESS/DOUBLE,NONSTD
      COMMON/GUESS3/ZMAP2(10,100,8),IMAP2(10,100,8)
      COMMON/TDALIST/LISTETDA, LISTVTDA
C
      DATA ONE /1.0D0/
C
      CALL ZERO(VEC, NSIZEC)
C
      IF(.NOT.DOUBLE.AND..NOT.NONSTD)THEN
C
C Originally the root was pass in as an argument. It is supposed to
C be an array but it never was. The other problem was that
C it was a double precision in the calling routine(STARTV) and
C integer here. Moreover, the root was not used for anything. So I 
C comment it until someone come up with a use for it. Ajith September, 1997.
C
C        CALL GETLST(ROOT,IREAD,1,1,IRREPX,LISTETDA)
        CALL GETLST(VEC,IREAD,1,1,IRREPX,LISTVTDA)
C         CALL SSCAL (LENSIN,ZFACT,VEC,1)
      ELSEIF(DOUBLE.AND..NOT.NONSTD)THEN
        CALL GETREC(20,'JOBARC','STRTGUES',IREAD,IJUNK)
        VEC(IJUNK(IREAD))=ONE
      ELSEIF(NONSTD)THEN
        DO 100 I=1,10
          VEC(IMAP2(I,IREAD,IRREPX))=ZMAP2(I,IREAD,IRREPX)
  100   CONTINUE
      ENDIF
C
      RETURN
      END
