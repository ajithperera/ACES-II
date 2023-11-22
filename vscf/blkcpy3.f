      SUBROUTINE BLKCPY3(MATFRM,NROWFRM,NCOLFRM,MATTAR,NROWTAR,NCOLTAR,
     &                  IROWFRM,ICOLFRM)
C
C THIS ROUTINE COPIES AN NROWTAR BY NCOLTAR MATRIX INTO AN NROWTAR BY 
C NCOLTAR TARGET BLOCK SUBMATRIX (IN MATFRM) SUCH THAT THE (1,1)
C ELEMENT IN MATTAR BECOMES THE IROWFRM,ICOLFRM ELEMENT IN MATFRM.
C
C THE PHYSICAL DIMENSION OF MATFRM IS (NROWFRM,NCOLFRM)
C THE PHYSICAL DIMENSION OF MATTAR IS THE SAME AS THE SUBMATRIX SIZE
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION MATTAR(NROWTAR,NCOLTAR),MATFRM(NROWFRM,NCOLFRM)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DO 10 ICOL=1,NCOLTAR
       ICOL0=ICOL+ICOLFRM-1
       IROW0=IROWFRM
       CALL SCOPY(NROWTAR,MATTAR(1,ICOL),1,MATFRM(IROW0,ICOL0),1)
10    CONTINUE
      RETURN
      END