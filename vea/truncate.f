      SUBROUTINE TRUNCATE(EVAL,EVEC,EVALSEL,EVECSEL,SCR,MAXCOR,NDIMVEC,
     &   NDIMOLD,NDIMNEW,INORM,IRREPX,ISIDE, ISPIN, IUHF, ICALC)
C
C THIS ROUTINE REDUCES THE DIMENSION OF THE CI EXPANSION BASIS
C BY PERFORMING A LINEAR TRANSFORMATION OF THE EXPANSION VECTORS.
C
C  NDIMNEW MAY BE CHANGED IN THIS PROCEDURE. ONLY REAL EIGENVALUES
C  ARE ALLOWED (EVAL < 1000, SAY, IMAGINARY EIGENVALUES WERE PUT TO
C  1.E12 OR SO. IF THERE ARE NOT ENOUGH REAL EIGENVALUES AVAILABLE
C  NDIMNEW IS LOWERED.
C
CEND
      IMPLICIT INTEGER (A-Z)
      PARAMETER (MAXORD=100)
      LOGICAL INORM,PRINT
      DOUBLE PRECISION EVEC(NDIMOLD,NDIMOLD),EVECSEL(NDIMOLD,NDIMNEW),
     $   SCR(MAXCOR), EVALSEL(NDIMNEW), EVAL(NDIMOLD)
      DOUBLE PRECISION RESID, XCLOSE, ROOT, EIGVAL
C
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/ROOTS/EIGVAL(100,8,3), OSCSTR(100,8,3)
      COMMON/EXTINF2/ROOT
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      DATA ONE,ZILCH /1.0D0,0.0D0/
C
      PRINT=.FALSE.
      IF(PRINT)THEN
       WRITE(6,2000) NDIMNEW
2000   FORMAT(T3,' Expansion vector space truncated to ',
     &           'dimension ',I3,'.')
      ENDIF
C
C  SELECT EIGENVECTORS THAT WILL BE RETAINED IN TRUNCATED SPACE
C  FIRST FIND CURRENT EIGENVECTOR AND ALREADY FOUND VECTORS
C
      CALL SCOPY(NDIMOLD,EVAL,1,EVALSEL,1)
      ICOUNT = 0
      CALL FNDCLOSE(NDIMOLD,EVAL,ROOT, XCLOSE,ICLOSE)
      ICOUNT = ICOUNT + 1
      CALL SCOPY(NDIMOLD,EVEC(1,ICLOSE),1,EVECSEL(1,ICOUNT),1)
      EVALSEL(ICOUNT) = XCLOSE
      EVAL(ICLOSE) = 1.0D30
         DO 10 I = MIN(NDIMNEW-1,IROOT),1, -1
            CALL FNDCLOSE(NDIMOLD,EVAL,
     $         EIGVAL(I,IRREPX,ISPIN+1-ITOTALS),XCLOSE,ICLOSE)
            ICOUNT = ICOUNT + 1
            CALL SCOPY(NDIMOLD,EVEC(1,ICLOSE),1,EVECSEL(1,ICOUNT),1)
            EVALSEL(ICOUNT) = XCLOSE
            EVAL(ICLOSE) = 1.0D30
 10      CONTINUE
C
C  NEXT RETAIN EIGENVALUES CLOSE TO ROOT
C
      DO 20 I = 1, NDIMNEW - IROOT - 1
         CALL FNDCLOSE(NDIMOLD,EVAL,ROOT, XCLOSE,ICLOSE)
         ICOUNT = ICOUNT + 1
         CALL SCOPY(NDIMOLD,EVEC(1,ICLOSE),1,EVECSEL(1,ICOUNT),1)
         EVALSEL(ICOUNT) = XCLOSE
         EVAL(ICLOSE) = 1.0D30
 20      CONTINUE
C
C  CHECK IF ALL RETAINED EIGENVALUES ARE MEANINGFUL (REAL ORIGINALLY)
C
 25      IF (EVALSEL(NDIMNEW) .GT. 1000.0) THEN
            NDIMNEW = NDIMNEW - 1
            IF (NDIMNEW .GE. 1) THEN
               GO TO 25
            ELSE
               WRITE(6,*) 'SOMETHING WRONG IN TRUNCATE '
               WRITE(6,*) ' SELECTED EIGENVALUES '
               DO 26 I = 1, NDIMOLD
                  WRITE(6,*) EVALSEL(I)
 26            CONTINUE
            ENDIF
         ENDIF
C
         CALL SCOPY(NDIMNEW,EVALSEL,1,EVAL,1)
C
C SCHMIDT ORTHOGONALIZE EXPANSION VECTORS
C
      DO 30 I=2,NDIMNEW
       CALL GSCHMIDT(EVECSEL(1,I),EVECSEL,NDIMOLD,I-1,SCR,RESID)
 30   CONTINUE
C
C
C CARRY OUT TRANSFORMATION TO NEW SPACE
C
      CALL TRNSPACE(NDIMVEC,NDIMNEW,NDIMOLD,SCR,EVECSEL,MAXCOR,IRREPX,
     &   ISIDE,LISTC,.TRUE.,MAXORD,IOLDEST,PRINT,IUHF,ICALC)
      CALL TRNSPACE(NDIMVEC,NDIMNEW,NDIMOLD,SCR,EVECSEL,MAXCOR,IRREPX,
     &   ISIDE,LISTHC,.FALSE.,MAXORD,IOLDEST,PRINT,IUHF,ICALC)
C
C TRANSFORM LASTVECT TO NEW BASIS
C
      IF(INORM)THEN
       I000=1
       I010=I000+100
       I020=I010+100
       CALL ZERO(SCR(I010),100)
       CALL GETREC(20,'JOBARC','LASTVCTR',100*IINTFP,SCR)
       CALL XGEMM('N','N',1,NDIMNEW,NDIMOLD,ONE,SCR(2),1,EVECSEL,
     $    NDIMOLD, ZILCH,SCR(I010+1),1)
       CALL PUTREC(20,'JOBARC','LASTVCTR',100*IINTFP,SCR(I010))
      ENDIF
C 
      RETURN
1000  FORMAT(T3,'@TRUNCATE-F, Insufficient core to truncate expansion ',
     &          'space.')
      END