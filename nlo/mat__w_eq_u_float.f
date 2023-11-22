         SUBROUTINE  MAT__W_EQ_U_FLOAT
     +
     +                    ( DDVECU,
     +                      DDVECW,
     +                      VEC,
     +                      U,
     +
     +                              W )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__W_EQ_U_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation copies a vector U of floating point
C                type to vector W.
C
C                The dimensions of both vectors need not be the same.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDVECU
         INTEGER    DDVECW
         INTEGER    I
         INTEGER    VEC

         DOUBLE PRECISION   U (1:DDVECU)
         DOUBLE PRECISION   W (1:DDVECW)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  ( VEC .GT. DDVECU )  THEN

               WRITE (1,*) ' Dimension of vector U too small: '
               WRITE (1,*) ' mat__w_eq_u_float '
               WRITE (1,*) ' DDVECU,VEC = ',
     +                       DDVECU,VEC

               WRITE (*,*) ' Dimension of vector U too small: '
               WRITE (*,*) ' mat__w_eq_u_float '
               WRITE (*,*) ' DDVECU,VEC = ',
     +                       DDVECU,VEC

               STOP

         END IF

         IF  ( VEC .GT. DDVECW )  THEN

               WRITE (1,*) ' Dimension of vector W too small: '
               WRITE (1,*) ' mat__w_eq_u_float '
               WRITE (1,*) ' DDVECW,VEC = ',
     +                       DDVECW,VEC

               WRITE (*,*) ' Dimension of vector W too small: '
               WRITE (*,*) ' mat__w_eq_u_float '
               WRITE (*,*) ' DDVECW,VEC = ',
     +                       DDVECW,VEC

               STOP

         END IF
C
C
C             ...loop over all vector elements.
C
C
         DO  10  I = 1,VEC  
             W (I) = U (I)
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
