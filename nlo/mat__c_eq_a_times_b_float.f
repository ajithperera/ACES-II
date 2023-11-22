         SUBROUTINE  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      ROW,SUM,COL,
     +                      A,B,
     +
     +                              C )
     +
C-----------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_A_TIMES_B_FLOAT
C  MODULE      : MATRIX
C  MODULE-ID   : MAT
C  DESCRIPTION : This operation multiplies the two-dimensional matrix A
C                by the two-dimensional matrix B using the middle 
C                product algorithm, which is the fastest algorithm for
C                matrix multiplication on scalar as well as vector
C                machines:
C
C                                C = A * B
C
C                Matrices A and B need not be distinct, but matrix C
C                must be distinct from A and B and all matrices must
C                be of floating point type.
C
C  AUTHOR      : Norbert Flocke
C-----------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         INTEGER    DDROWA,DDCOLA
         INTEGER    DDROWB,DDCOLB
         INTEGER    DDROWC,DDCOLC
         INTEGER    I,J,K
         INTEGER    ROW,SUM,COL

         DOUBLE PRECISION   S

         DOUBLE PRECISION   A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION   B (1:DDROWB,1:DDCOLB)
         DOUBLE PRECISION   C (1:DDROWC,1:DDCOLC)

C         INTEGER    DDROWA,DDCOLA
C         INTEGER    DDROWB,DDCOLB
C         INTEGER    DDROWC,DDCOLC
C         INTEGER    I,J,K
C         INTEGER    INJ
C         INTEGER    NJ,NJB,NJC,NK
C         INTEGER    ROW,SUM,COL
C
C         DOUBLE PRECISION   S
C
C         DOUBLE PRECISION   A (1:ROW*SUM)
C         DOUBLE PRECISION   B (1:SUM*COL)
C         DOUBLE PRECISION   C (1:ROW*COL)
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions.
C
C
         IF  ( ROW .GT. DDROWA   .OR.
     +         SUM .GT. DDCOLA         )  THEN
               WRITE (1,*) ' Dimensions of matrix A too small: '
               WRITE (1,*) ' mat__c_eq_a_times_b_float '
               WRITE (1,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                       DDROWA,DDCOLA,ROW,SUM
               WRITE (*,*) ' Dimensions of matrix A too small: '
               WRITE (*,*) ' mat__c_eq_a_times_b_float '
               WRITE (*,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                       DDROWA,DDCOLA,ROW,SUM
               STOP
         END IF

         IF  ( SUM .GT. DDROWB   .OR.
     +         COL .GT. DDCOLB         )  THEN
               WRITE (1,*) ' Dimensions of matrix B too small: '
               WRITE (1,*) ' mat__c_eq_a_times_b_float '
               WRITE (1,*) ' DDROWB,DDCOLB,SUM,COL = ',
     +                       DDROWB,DDCOLB,SUM,COL
               WRITE (*,*) ' Dimensions of matrix B too small: '
               WRITE (*,*) ' mat__c_eq_a_times_b_float '
               WRITE (*,*) ' DDROWB,DDCOLB,SUM,COL = ',
     +                       DDROWB,DDCOLB,SUM,COL
               STOP
         END IF

         IF  ( ROW .GT. DDROWC   .OR.
     +         COL .GT. DDCOLC         )  THEN
               WRITE (1,*) ' Dimensions of matrix C too small: '
               WRITE (1,*) ' mat__c_eq_a_times_b_float '
               WRITE (1,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                       DDROWC,DDCOLC,ROW,COL
               WRITE (*,*) ' Dimensions of matrix C too small: '
               WRITE (*,*) ' mat__c_eq_a_times_b_float '
               WRITE (*,*) ' DDROWC,DDCOLC,ROW,COL = ',
     +                       DDROWC,DDCOLC,ROW,COL
               STOP
         END IF
C
C
C             ...middle product algorithm for small dimensions:
C                   1) initialize elements of matrix C.
C                   2) build the product such that first index
C                      varies most.
C
C
C         DO 10 J = 1,COL
C            S = B (1,J)
C            DO 20 I = 1,ROW  
C               C (I,J) = A (I,1) * S
C   20       CONTINUE
C   10    CONTINUE
C
C         DO 30 K = 2,SUM
C         DO 30 J = 1,COL
C            S = B (K,J)
C            DO 40 I = 1,ROW
C               C (I,J) = C (I,J) + A (I,K) * S
C   40       CONTINUE
C   30    CONTINUE
C
C
C             ...highly optimized routine for large dimensions.
C
C
         CALL  MAT__DGEMM_C_EQ_AXB_FLOAT
     +
     +              ( DDROWA,DDCOLA,
     +                DDROWB,DDCOLB,
     +                DDROWC,DDCOLC,
     +                ROW,SUM,COL,
     +                A,B,
     +
     +                       C )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
