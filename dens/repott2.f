      SUBROUTINE REPOTT2(GOLD,G,NBSA0,nda,itype)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GOLD(2),G(2) 
      jipoint(j,i) = (i*i - 3*i)/2 + (j+1) 
c
      IF (NBSA0.LE.1) RETURN   
c---------- drop  vrt-vrt mo -----------
      IOFF = 0
      if (itype.eq.1) go to 111
c---------- drop  occ-occ mo -----------
      IOFF = NDA 
      if (itype.eq.4) go to 111 
      write (6,109) itype
 109  format (//3x,'there must be some error, itype must be ',
     x             ' 1 or 4 for repott2, present itype=',i5)
      call errex
c
c------------------------------------
 111  continue
      NBSA = NBSA0 + NDA
      LENG = (NBSA*(NBSA-1))/2 
       CALL ZERO(G,LENG) 
c
      DO 50 I=2,NBSA0
       JI = I - 1
       DO 40 J=1,JI
        JN = J + IOFF
        IN = I + IOFF 
        G(jipoint(JN,IN)) = GOLD(jipoint(J,I)) 
  40   CONTINUE
  50  CONTINUE 
       CALL ZERO(GOLD,(NBSA0*(NBSA0-1))/2) 
c------------------------------------
      RETURN
      END
