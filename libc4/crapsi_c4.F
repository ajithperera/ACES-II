      SUBROUTINE CRAPSI_C4(ICORE,I0,MEMALLOC,IUHF,IENTRY)
      IMPLICIT INTEGER (A-Z)
      INTEGER ICORE(1)
      COMMON/RECSIZE/ISCALE
      COMMON/FLAGS/IFLAGS(100),IFLAGS2(100)

#include "machsp.com"

      IONE    = 1
      NFLAGS  = 100 
      NFLAGS2 = 500 
      ISCALE  = 1 

      CALL ACES_INIT_RTE
      CALL DUMPJA_C4('I')
      CALL CHECKIN_C4

      IF(IENTRY.NE.-1) THEN 
         CALL GETREC_C4(20,'JOBARC','IFLAGS  ',NFLAGS,IFLAGS)
         CALL GETREC_C4(20,'JOBARC','IFLAGS2 ',NFLAGS,IFLAGS2)

         IF(IENTRY.EQ.100) THEN
           IENTRY=0
           MEM=IFLAGS(36)
           MEMMIN=MIN(4500000,IFLAGS(36))
           MEMMAX=MAX(MEMMIN,IFLAGS(36)-10000000)
         ELSE
           MEMMAX=IFLAGS(36)
         ENDIF
         ICORE(1) = 0
         CALL GETMEM_A2(ICORE(1),I0,MEMALLOC,MEMMAX,.TRUE.)
         
C         DO WHILE (ICORE(1) .EQ. 0 .AND. (MEMMAX .GT. 1000000))
C            CALL ACES_MALLOC(MEMMAX,ICORE,I0)
C            If (ICORE(1) .EQ. 0) MEMMAX = (MEMMAX - 1000000)
C         ENDDO 
C         IF (ICORE(1) .EQ. 0) THEN
C            WRITE(6,"(2a)") "@-crapsi_c4: Unable to allocate at",
C     &                     " least 1000000 integers of memory"
C            CALL ERREX 
C         ENDIF 
         CALL POPCORN_C4(MEMALLOC)
      ELSE
         I0=1
      ENDIF

CSSS      CALL ACES_IO_INIT_C4(ICORE(I0),I0,MEMALLOC,.TRUE.)
CSSS      CALL ACES_INIT_CHEMSYS_C4

      CALL GETSTF_C4(ICORE(I0),IUHF,IENTRY)

      RETURN
      END
