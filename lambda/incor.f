
C This routine loads lists into ICORE(I0), increments I0, and decrements
C ICRSIZ. I0 and ICRSIZ do not have to point to the /ISTART/ common
C block, but they must be relative to ICORE(1) in the blank common
C block. The logic is also reentrant, so it will not load lists that are
C already in memory.

      SUBROUTINE INCOR(I0,ICRSIZ,IUHF)
      IMPLICIT NONE
      INTEGER I0,ICRSIZ,IUHF
      INTEGER ICORE(1),IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,IFLAGS(100)
      COMMON / / ICORE
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS

      IF (IFLAGS(35).EQ.0) THEN
c      o NONE
      ELSE IF (IFLAGS(35).EQ.1) THEN
c      o NOABCD
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,  1,130,1)
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,134,200,1)
      ELSE IF (IFLAGS(35).EQ.2) THEN
c      o T
      ELSE IF (IFLAGS(35).EQ.3) THEN
c      o NONE
      ELSE IF (IFLAGS(35).EQ.4) THEN
c      o NOABCI
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,  1, 26,1)
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1, 31,126,1)
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,134,165,1)
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,170,200,1)
      ELSE IF (IFLAGS(35).EQ.5) THEN
c      o NOABIJ
         CALL INCOR_GEN(I0,ICRSIZ,IUHF)
      ELSE IF (IFLAGS(35).EQ.6) THEN
c      o ALL
         CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,1,500,1)
      ELSE
c      o (impossible)
      END IF

      RETURN
      END

