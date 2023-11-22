      SUBROUTINE QUIKAB(ICORE,MAXCOR,IUHF)
C
C THIS SUBROUTINE FORMS THE ALL ALPHA W(mbej) INTERMEDIATE
C  FROM THE ABAB AND ABBA PIECES FOR RHF REFERENCE FUNCTIONS.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONEM
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONEM /-1.0/
      IF(IUHF.NE.0)THEN
       WRITE(6,100)
100    FORMAT(T3,'@QUIKAB, Fatal logic error.  QUIKAB called for UHF.')
      ENDIF
      TOTSIZ=ISYMSZ(ISYTYP(1,56),ISYTYP(1,56))
      I000=1
      I010=I000+TOTSIZ*IINTFP
      I020=I010+TOTSIZ*IINTFP 
      IF(I020.GT.MAXCOR)CALL INSMEM('QUIKAA',I020,MAXCOR)
      CALL GETALL(ICORE(I000),TOTSIZ,1,54)
      CALL GETALL(ICORE(I010),TOTSIZ,1,58)
      CALL VADD(ICORE(I000),ICORE(I010),ICORE(I000),TOTSIZ,ONEM)
      CALL PUTALL(ICORE,TOTSIZ,1,56)
      RETURN
      END
