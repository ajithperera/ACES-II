         SUBROUTINE  MAT__W_EQ_ZERO_FLOAT
     +
     +                    ( DDVEC,
     +                      VEC,
     +
     +                              W )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__W_EQ_ZERO_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a zero vector W of floating
C                point type.
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
         INTEGER    VEC

         DOUBLE PRECISION   W (1:DDVEC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimension.
C
C
         IF  ( VEC .GT. DDVEC )  THEN

               WRITE (1,*) ' Dimension of vector W too small: '
               WRITE (1,*) ' mat__w_eq_zero_float '
               WRITE (1,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               WRITE (*,*) ' Dimension of vector W too small: '
               WRITE (*,*) ' mat__w_eq_zero_float '
               WRITE (*,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               STOP

         END IF
C
C
C             ...build zero vector.
C
C
         DO  10  I = 1,VEC  
             W (I) = 0.D0
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
