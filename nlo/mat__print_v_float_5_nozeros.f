         SUBROUTINE  MAT__PRINT_V_FLOAT_5_NOZEROS
     +
     +                    ( UNITID,
     +                      TITLE,
     +                      DDVEC,
     +                      VEC,
     +                      V )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__PRINT_V_FLOAT_5_NOZEROS
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation prints a vector V to the unit specified
C                by UNITID in floating point format F11.5 . Values
C                below 5.0D-6 are not printed.
C
C                The vector is printed in one column.                    
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDVEC
         INTEGER    I
         INTEGER    UNITID
         INTEGER    VEC

         CHARACTER*11       S
         CHARACTER*(*)      TITLE

         DOUBLE PRECISION   LIMIT

         DOUBLE PRECISION   V (1:DDVEC)

         DATA LIMIT  /5.0D-6/
C
C
C------------------------------------------------------------------------
C
C
C             ...printout title.
C
C
         WRITE (UNITID,6000)
         WRITE (UNITID, *  ) TITLE
         WRITE (UNITID,7000)
C
C
C             ...immediate return if length is zero.
C
C
         IF  ( VEC.EQ.0 )  RETURN
C
C
C             ...check dimension.
C
C
         IF  ( VEC .GT. DDVEC )  THEN

               WRITE (1,*) ' Dimension of vector V too small: '
               WRITE (1,*) ' mat__print_v_float_5_nozeros '
               WRITE (1,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               WRITE (*,*) ' Dimensions of vector V too small: '
               WRITE (*,*) ' mat__print_v_float_5_nozeros '
               WRITE (*,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               STOP

         END IF
C
C
C             ...print vector.
C
C
         DO  10  I = 1,VEC

             IF  ( ABS ( V(I) ) .LT. LIMIT )  THEN
                   S = ' '
             ELSE 
                   WRITE (S,8000) V(I)
             END IF

             WRITE (UNITID,9000) I,S

   10    CONTINUE
C
C
C             ...formats of printing.
C
C
 6000    FORMAT  (/)
 7000    FORMAT  (/)
 8000    FORMAT  (F11.5)
 9000    FORMAT  (1X,I8,A11)
C
C
C             ...ready!
C
C
         RETURN
         END
