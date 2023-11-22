      SUBROUTINE PROCEVEC(ICORE,MAXCOR,IUHF,NBASIS,NCOMP)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /AOOFST/ INDOCC(8,2)
C
      NBAS2=NBASIS*NCOMP
      I000=1
      I010=I000+IINTFP*NBAS2
      I020=I010+IINTFP*NBAS2
      CALL GETREC(20,'JOBARC','SCFEVECA',NBAS2*IINTFP,ICORE(I000))
C
C CALL ROUTINE WHICH FORMS THE SYMMETRY-ORDERED EIGENVECTOR MATRIX
C AND THE SYMMETRY OFFSETS (STORED IN COMMON BLOCK AOOFST).
C
      IOFF0=0
      CALL SYPKEV(ICORE(I000),ICORE(I010),NBASIS,TOTLENA,IOFF0,1)
      I010=I000+IINTFP*TOTLENA
C
C  NOW REORDER BETA EIGENVECTORS (UHF ONLY)
C
      IF(IUHF.NE.0)THEN
       IOFF0=IOFF0+TOTLENA
       CALL GETREC(20,'JOBARC','SCFEVECB',NBAS2*IINTFP,ICORE(I010))
       I030=I020+IINTFP*NBAS2
       CALL SYPKEV(ICORE(I010),ICORE(I020),NBASIS,TOTLENB,IOFF0,2)
       I020=I010+IINTFP*TOTLENB
      ELSE
       CALL ICOPY(NIRREP,INDOCC(1,1),1,INDOCC(1,2),1)
       I020=I010
      ENDIF
      RETURN
      END
