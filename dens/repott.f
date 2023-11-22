      SUBROUTINE REPOTT(GOLD,G,NBSA0,NBSB0,NBSA,NBSB,nda,ndb,itype)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GOLD(NBSA0,NBSB0)
      DIMENSION G(NBSA,NBSB)
c---------- drop  vrt-vrt mo -----------
      IOFF = 0
      JOFF = 0
      if (itype.eq.1) go to 111
c---------- drop  occ-vrt mo -----------
      IOFF = 0
      JOFF = NDA
      if (itype.eq.2) go to 111
c---------- drop  vrt-occ mo -----------
      IOFF = NDB 
      JOFF = 0
      if (itype.eq.3) go to 111 
c---------- drop  occ-occ mo -----------
      IOFF = NDB 
      JOFF = NDA
      if (itype.eq.4) go to 111 
      write (6,109) itype
 109  format (//3x,'there must be some error, itype must be ',
     x             'between 1 and 4,  present itype=',i5)
      call errex
c
c------------------------------------
 111  continue
       call zero(G,NBSA*NBSB) 
c
      DO 50 I=1,NBSB0
       DO 40 J=1,NBSA0
        JN = J + JOFF
        IN = I + IOFF 
        G(JN,IN) = GOLD(J,I) 
  40   CONTINUE
  50  CONTINUE 
       call zero(GOLD,NBSA0*NBSB0) 
c------------------------------------
      RETURN
      END
