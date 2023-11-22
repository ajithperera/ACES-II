      SUBROUTINE REOREVL(EVALIN,EVALOT,NBAS,ISPIN)
C
C  This routine reorders the eigenvectors from correlated order
C  to that required for the SCF.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "maxbasfn.par"
C
      CHARACTER*1 ISP(2)
C
      DIMENSION EVALIN(NBAS),EVALOT(NBAS),ireordr(MAXBASFN)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ISP /'A','B'/
C
C ORDER EIGENVECTORS BY IRREP
C
      call getrec(-1,'JOBARC','REORDER'//isp(ispin),nbas,ireordr)
      ITOT=0
      DO 100 I=1,NBAS
        ITOT=ITOT+IREORDR(I)
  100 CONTINUE
      IF(ITOT.EQ.0) THEN
        CALL SCOPY(NBAS,EVALIN,1,EVALOT,1)
        RETURN
      ENDIF
C
      do 10 i=1,nbas
       jtar=ireordr(i)
       evalot(i)=evalin(jtar)
10    continue
      call scopy(nbas,evalot,1,evalin,1)
      return
      end
