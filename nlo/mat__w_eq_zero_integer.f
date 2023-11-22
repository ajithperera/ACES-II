         SUBROUTINE  MAT__W_EQ_ZERO_INTEGER
     +
     +                    ( DDVEC,
     +                      VEC,
     +
     +                              W )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__W_EQ_ZERO_INTEGER
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation builds a zero vector W of integer type.
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

         INTEGER    W (1:DDVEC)
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
               WRITE (1,*) ' mat__w_eq_zero_integer '
               WRITE (1,*) ' DDVEC,VEC = ',
     +                       DDVEC,VEC

               WRITE (*,*) ' Dimension of vector W too small: '
               WRITE (*,*) ' mat__w_eq_zero_integer '
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
             W (I) = 0
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
