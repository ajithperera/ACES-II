         SUBROUTINE  MAT__C_EQ_C_TRANSPOSED_FLOAT
     +
     +                    ( DDROWC,DDCOLC,
     +                      ROW,COL,
     +                      C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_C_TRANSPOSED_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation transposes a 2-dimensional matrix C of
C                floating point type.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDROWC,DDCOLC
         INTEGER    I,J,K
         INTEGER    ROW,COL

         DOUBLE PRECISION   TEMP

         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  ( ROW .GT. DDCOLC   .OR.
     +         COL .GT. DDROWC         )  THEN

               WRITE (1,*) ' Dimensions of matrix C too small: '
               WRITE (1,*) ' mat__c_eq_c_transposed_float '
               WRITE (1,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                       DDROWC,DDCOLC,ROW,COL

               WRITE (*,*) ' Dimensions of matrix C too small: '
               WRITE (*,*) ' mat__c_eq_c_transposed_float '
               WRITE (*,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                       DDROWC,DDCOLC,ROW,COL

               STOP

         END IF
C
C
C             ...transpose common upper left square.
C
C
         K = MIN0 (ROW,COL)
         DO 10 I = 1,K
         DO 10 J = I,K
            TEMP = C (I,J)
            C (I,J) = C (J,I)
            C (J,I) = TEMP
   10    CONTINUE
C
C
C             ...transpose rest of C.
C
C
         IF (ROW.GT.COL) THEN
             DO 20 J = 1,COL
             DO 20 I = COL+1,ROW
                C (J,I) = C (I,J)
   20        CONTINUE
         ELSE IF (ROW.LT.COL) THEN
             DO 30 I = 1,ROW
             DO 30 J = ROW+1,COL
                C (J,I) = C (I,J)
   30        CONTINUE
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
