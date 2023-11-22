C
C This routine loads lists into ICORE(I0), increments I0, and decrements
C ICRSIZ. I0 and ICRSIZ do not have to point to the /ISTART/ common
C block, but they must be relative to ICORE(1) in the blank common
C block. The logic is also reentrant, so it will not load lists that are
C already in memory.
C
      SUBROUTINE INCOR(I0,ICRSIZ,IUHF)
      IMPLICIT INTEGER (A-Z)
      COMMON / / ICORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /AUXCACHE/ QUIKGET(10,500)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
cjp
#include "bwcc.com"
C
cjp
cjp because of switching between MO.... files of different references
cjp do not hold anything in core
cjp USE OF INCOR is INCOMPATIBLE with whole filosophy of this MR implementation
cjp the memory buffers would have to be saved to file and restored from
cjp in each iref loop anyway, and the use of cache is thus efficient enough
cjp or we would have to replicate the memory buffers and replicate the common auxcache
cjp
      IF(IFLAGS(2).LE.1) RETURN
      IF(IFLAGS(35).EQ.0)RETURN
       if(isbwcc.and. nref.gt.1) then
        write(6,*)'@INCOR-E: IN CORE REGIME NOT SUPPORTED FOR MR-BW-CC 
     & CALCULATION (with more that one reference per processor)
     &  -> SWITCHED OFF'
        return
        endif
C
      IF(IFLAGS(35).EQ.1)THEN
C
C HOLD EVERYTHING BUT ABCD INTEGRALS IN CORE!
C
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,1,199,1)
C
      ELSEIF(IFLAGS(35).EQ.6)THEN
C
C HOLD ABSOLUTELY EVERYTHING IN CORE
C
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,1,500,1)
C
      ELSEIF(IFLAGS(35).EQ.5)THEN
C
C HOLD IJKL, IJKA AND SINGLES IN CORE
C
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1, 7,13,1)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,51,53,1)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,90,94,1)
C
      ELSEIF(IFLAGS(35).EQ.4)THEN
C
C HOLD EVERYTHING BUT ABCD AND ABCI INTEGRALS IN CORE!
C
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1, 1, 26,1)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,10,1,31,199,1)
C
      ELSEIF(IFLAGS(35).EQ.2)THEN
C
C HOLD T2 AND INCREMENTS IN CORE
C
       IOFF = 2*(1-IUHF)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,NIRREP,1,44+IOFF,46,1)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,NIRREP,1,61+IOFF,63,1)
C
C FOR QCISD AND CCSD KEEP THE SINGLE STUFF IN CORE AS WELL
C
      IF(IFLAGS(2).GE.10) THEN
C
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,1,1+IUHF,1,90,93,1)
       CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,3,3+IUHF,1,90,90,1)
C
       IF(IFLAGS(11).EQ.2.OR.IFLAGS(77).NE.0)THEN
        CALL ACES_AUXCACHE_LDBLK(I0,ICRSIZ,3,3+IUHF,1,91,93,1)
       ENDIF
      ENDIF
      ENDIF
      RETURN
      END
