         SUBROUTINE  MAT__C_EQ_UNIT_FLOAT
     +
     +                    ( DDROW,DDCOL,
     +                      ROW,COL,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_UNIT_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a two-dimensional unit matrix C
C                of floating point type.
C
C                The unit matrix need not to be a square matrix.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDROW,DDCOL
         INTEGER    I,J
         INTEGER    NDIAG
         INTEGER    ROW,COL

         DOUBLE PRECISION   C (1:DDROW,1:DDCOL)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  ( ROW .GT. DDROW   .OR.
     +         COL .GT. DDCOL         )  THEN

               WRITE (1,*) ' Dimensions of matrix C too small: '
               WRITE (1,*) ' mat__c_eq_unit_float '
               WRITE (1,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               WRITE (*,*) ' Dimensions of matrix C too small: '
               WRITE (*,*) ' mat__c_eq_unit_float '
               WRITE (*,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               STOP

         END IF
C
C
C             ...initialize matrix to zero.
C
C
         DO  10  J = 1,COL
         DO  10  I = 1,ROW  
             C (I,J) = 0.D0
   10    CONTINUE
C
C
C             ...set all diagonal matrix elements to 1.
C
C
         NDIAG  =  MIN0 ( ROW,COL )

         DO  20  I = 1,NDIAG  
             C (I,I) = 1.D0
   20    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
