         SUBROUTINE  MAT__C_EQ_ZERO_INTEGER
     +
     +                    ( DDROW,DDCOL,
     +                      ROW,COL,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_ZERO_INTEGER
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a two-dimensional zero matrix C
C                of integer type.
C
C                The zero matrix need not to be a square matrix.
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
         INTEGER    ROW,COL

         INTEGER    C (1:DDROW,1:DDCOL)
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
               WRITE (1,*) ' mat__c_eq_zero_integer '
               WRITE (1,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               WRITE (*,*) ' Dimensions of matrix C too small: '
               WRITE (*,*) ' mat__c_eq_zero_integer '
               WRITE (*,*) ' DDROW,DDCOL,ROW,COL = ',
     +                       DDROW,DDCOL,ROW,COL

               STOP

         END IF
C
C
C             ...build zero matrix.
C
C
         DO  10  J = 1,COL
         DO  10  I = 1,ROW  
             C (I,J) = 0
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
