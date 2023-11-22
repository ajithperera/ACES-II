      SUBROUTINE GETBLK(WFULL,WBLOCK,NBF,NBAS,IREPS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION WFULL(NBAS,NBAS),WBLOCK(NBF,NBF)
C
      if (ireps-1+nbf.gt.nbas) then
         print *, '@GETBLK: Assertion failed.'
         print *, '         nbas = ',nbas
         print *, '         ireps = ',ireps,'; nbf = ',nbf
         call errex
      end if
C
      DO 10 J=1,NBF
        DO 11 I=1,NBF
          WBLOCK(I,J)=WFULL((IREPS-1)+I,(IREPS-1)+J)
   11   CONTINUE
   10 CONTINUE
C
      RETURN
      END
