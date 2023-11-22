      SUBROUTINE MATEXP2(IRREPX,WPACK,WFULL,NBAS)
C
C THIS ROUTINE TAKES THE COMPRESSED REPRESENTATION OF A ONE-ELECTRON
C MATRIX IN THE SYMMETRY ADAPTED ATOMIC ORBITAL BASIS AND EXPANDS IT 
C TO THE FULL MATRIX.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      DIMENSION WPACK(*),WFULL(NBAS,NBAS)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
C
      CALL ZERO(WFULL,NBAS*NBAS)
C
      IF(IRREPX.EQ.1)THEN
       IOFF=1
       DO 10 IRREP=1,NIRREP
        NDIM=IAOPOP(IRREP)
        IPOS=IOFFAO(IRREP)
        CALL BLKCPY(WPACK(IOFF),NDIM,NDIM,WFULL,NBAS,NBAS,IPOS,IPOS)
        IOFF=IOFF+NDIM*NDIM
10     CONTINUE 
C
      ELSE
C
       IOFF=1
       DO 20 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NDIML=IAOPOP(IRREPL)
        NDIMR=IAOPOP(IRREPR)
        IPOSL=IOFFAO(IRREPL)
        IPOSR=IOFFAO(IRREPR)
        CALL BLKCPY(WPACK(IOFF),NDIML,NDIMR,WFULL,NBAS,NBAS,IPOSL,IPOSR)
        IOFF=IOFF+NDIML*NDIMR
20     CONTINUE 
C
      ENDIF
C
      RETURN
      END