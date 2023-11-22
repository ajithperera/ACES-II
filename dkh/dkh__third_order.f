         SUBROUTINE  DKH__THIRD_ORDER
     +
     +                   ( KNOW,NBAS,
     +                     NAI,PVP,
     +                     SCR1,SCR2,
     +                     E1,PE1P,
     +                     AA,RR,TT,EP,
     +                     AVHA,APVHPA,
     +                     W1W1,W1E1W1,
     +
     +                               SCR )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : DKH__THIRD_ORDER
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__MATMUL_NN
C                DKH__MATMULD
C                DKH__MATMULM
C                DKH__MATRIXADD
C                DKH__SQ_PRINT
C
C  DESCRIPTION : Perform the third order Douglas-Kroll-Hess
C                transformation.
C
C
C                  Input:
C
C                    KNOW         =  Level of printing, rather, how much
C                                    the user wants to "KNOW"
C                    NBAS         =  number of atomic orbital basis
C                                    functions
C                    NAI          =  transformed nuclear attraction integrals
C                    PVP          =  transformed pVp integrals
C                    SCRx         =  scratch matrix for transformation
C                    PE1P         =  transformed pVp potential
C                    E1           =  first order DKH potential
C                    AA           =  A matrix in DKH routine
C                    RR           =  R matrix in DKH routine
C                    TT           =  kinetic energy eigenvalues
C                    EP           =  relativistic kinetic energies
C                    AVHA         =  AVA divided by energy denominator
C                    APVHPA       =  APVPA divided by energy denominator
C                    W1W1         =  DKH specific matrix
C                    W1E1W1       =  DKH specific matrix
C
C                  Output:
C
C                    SCR          =  third order DKH potential
C
C
C  AUTHOR      : Thomas Watson Jr.
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     KNOW,NBAS
         INTEGER     I,J

         DOUBLE PRECISION  ZERO,ONE,ONEM
         DOUBLE PRECISION  HALF,HALFM

         DOUBLE PRECISION  AA (1:NBAS)
         DOUBLE PRECISION  RR (1:NBAS)
         DOUBLE PRECISION  TT (1:NBAS)
         DOUBLE PRECISION  EP (1:NBAS)

         DOUBLE PRECISION  E1  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  NAI (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR (1:NBAS,1:NBAS)

         DOUBLE PRECISION  SCR1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVHA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1W1 (1:NBAS,1:NBAS)

         DOUBLE PRECISION  PE1P   (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1E1W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  APVHPA (1:NBAS,1:NBAS)

         CHARACTER    STRING*50

         PARAMETER  ( ZERO  =  0.0D0 )
         PARAMETER  ( ONE   =  1.0D0 )
         PARAMETER  ( ONEM  = -1.0D0 )
         PARAMETER  ( HALF  =  0.5D0 )
         PARAMETER  ( HALFM = -0.5D0 )
C
C
C------------------------------------------------------------------------
C
C
C                             ~
C             ...create the 'AVA' array!
C
C
         DO I = 1, NBAS
         DO J = 1, NBAS
            AVHA (I,J) = AA (I) * NAI (I,J) * AA (J)
         END DO
         END DO

         DO I = 1, NBAS
         DO J = 1, NBAS
            AVHA (I,J) = AVHA (I,J) / (EP (I) + EP (J))
         END DO
         END DO
C
C
C                          ~
C             ...create 'APVPA' array!
C
C
         DO I = 1, NBAS
         DO J = 1, NBAS
            APVHPA (I,J) = AA (I) * RR (I) * PVP (I,J) * RR (J) * AA (J)
         END DO
         END DO
      
         DO I = 1, NBAS
         DO J = 1, NBAS
            APVHPA (I,J) = APVHPA (I,J) / (EP (I) + EP (J))
         END DO
         END DO
C
C
C             ...create 'W1W1' array!
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,AVHA,W1W1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,APVHPA,W1W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,AVHA,AVHA,W1W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,APVHPA,W1W1,ONE,ONE )
     +
     +
C
C
C             ...create 'W1E1W1' array!
C
C
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,PE1P,SCR1,SCR,ONE,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,PE1P,SCR2,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR1,AVHA,W1E1W1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR1,APVHPA,W1E1W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR2,AVHA,W1E1W1,ONEM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR2,APVHPA,W1E1W1,SCR,ONE,ONE,TT,RR )
     +
     +
C
C
C             ...DKH3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,E1,SCR1,HALF,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,W1W1,SCR1,HALF,ONE )
     +
     +
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SCR1,ONEM,W1E1W1,SCR )
     +
     +
C
C
C             ...print the potential if the user wants!
C
C
         IF (KNOW .GE. 10) THEN
             STRING = "*** Third order potential (p-space) ***"
             CALL  DKH__SQ_PRINT (NBAS,SCR,STRING)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END

