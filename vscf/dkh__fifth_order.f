         SUBROUTINE  DKH__FIFTH_ORDER
     +
     +                   ( KNOW,NBAS,
     +                     NAI,PVP,
     +                     SCRH1,SCRH2,SCRH3,
     +                     SCRH4,SCRH5,SCRH6,
     +                     SCRH7,SCRH8,SCRHP5,
     +                     SCRHP6,SCRHP7,SCRHP8,
     +                     E1,PE1P,AA,RR,TT,EP,
     +                     AVA,AVHA,APVPA,
     +                     APVHPA,W1W1,W1O1,
     +                     O1W1,W1W2,W2W1,W2W2,
     +                     W2O1,O1W2,W1W12,W2W13,
     +                     W13W2,W2E1W2,SUM1,
     +                     SUM2,SUM3,SUM4,
     +
     +                                         SCR )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : DKH__FIFTH_ORDER
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__MATMUL_NN
C                DKH__MATMULD
C                DKH__MATMULM
C                DKH__MATDIVE
C                DKH__MATRIXADD
C                DKH__SQ_PRINT
C                DKH__MATRHTE
C                DKH__MATRHTE
C
C  DESCRIPTION : Perform the fifth order Douglas-Kroll-Hess
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
C                    E1           =  first order DKH potential
C                    PE1P         =  transformed pVp potential
C                    AA           =  A matrix in DKH routine
C                    RR           =  R matrix in DKH routine
C                    TT           =  kinetic energy eigenvalues
C                    EP           =  relativistic kinetic energies
C                    AVA          =  transformed NAI (p-space) array
C                    APVPA        =  transformed pVp (p-space) array
C                    AVHA         =  AVA divided by energy denominator
C                    APVHPA       =  APVPA divided by energy denominator
C                    WxWx         =  DKH specific matricies
C                    WxOx         =  DKH specific matricies
C                    SUMx         =  scratch matrices used in matrix additions
C
C                  Output:
C
C                    SCR          =  fifth order DKH potential
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

         DOUBLE PRECISION  A3,ZERO,ONE,ONEM
         DOUBLE PRECISION  HALF,HALFM,EIGHTH,QUART

         DOUBLE PRECISION  AA (1:NBAS)
         DOUBLE PRECISION  RR (1:NBAS)
         DOUBLE PRECISION  TT (1:NBAS)
         DOUBLE PRECISION  EP (1:NBAS)

         DOUBLE PRECISION  E1  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  NAI (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR (1:NBAS,1:NBAS)

         DOUBLE PRECISION  SUM1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SUM2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SUM3 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SUM4 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVHA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1O1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  O1W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1W2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W2W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W2W2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W2O1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  O1W2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PE1P (1:NBAS,1:NBAS)

         DOUBLE PRECISION  SCRH1  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH2  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH3  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH4  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH5  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH6  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH7  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH8  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRHP5 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRHP6 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRHP7 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRHP8 (1:NBAS,1:NBAS)

         DOUBLE PRECISION  APVPA  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W13W2  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W2W13  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1W12  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  APVHPA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W2E1W2 (1:NBAS,1:NBAS)

         CHARACTER    STRING*50

         PARAMETER  ( A3 = 0.14644661D0 )
         PARAMETER  ( ZERO   =  0.0D0   )
         PARAMETER  ( ONE    =  1.0D0   )
         PARAMETER  ( ONEM   = -1.0D0   )
         PARAMETER  ( HALF   =  0.50D0  )
         PARAMETER  ( HALFM  = -0.50D0  )
         PARAMETER  ( EIGHTH =  0.125D0 )
         PARAMETER  ( QUART  =  0.250D0 )
C
C
C------------------------------------------------------------------------
C
C                                       ~
C             ...create the 'AVA' and 'AVA' arrays!
C
C
         DO I = 1, NBAS
         DO J = 1, NBAS
            AVA (I,J) = AA (I) * NAI (I,J) * AA (J)
         END DO
         END DO

         DO I = 1, NBAS
         DO J = 1, NBAS
            AVHA (I,J) = AVA (I,J) / (EP (I) + EP (J))
         END DO
         END DO
C
C
C                                      ~
C             ...create 'APVPA' and 'APVPA arrays!
C
C        
         DO I = 1, NBAS
         DO J = 1, NBAS 
            APVPA (I,J) = AA (I) * RR (I) * PVP (I,J) * RR (J) * AA (J)
         END DO
         END DO

         DO I = 1, NBAS
         DO J = 1, NBAS
            APVHPA (I,J) = APVPA (I,J) / (EP (I) + EP (J))
         END DO
         END DO
C
C
C             ...calculate scratch-matrices SCRH_i (i=1,..,8):  
C
C
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,PE1P,SCRH1,SCR,ONE,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,PE1P,SCRH2,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,APVHPA,SCRH3,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,AVHA,SCRH4,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,E1,SCRH5,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,E1,SCRH6,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,PE1P,AVHA,SCRH7,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,APVHPA,SCRH8,SCR,ONE,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATDIVE ( NBAS,SCRH1,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH2,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH3,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH4,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH5,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH6,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH7,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH8,EP )
C
C
C     2.  Calculation of W1^2, W1*W2, W2*W1, W2^2, (W1^3)*W2, W2*(W1^3)         
C
C
C             a)  Calculation of W1^2:
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
C             b)  Calculation of W1^4:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,W1W1,W1W12,ONE,ZERO )
     +
     +
C
C
C             ...create 'W1O1' array.
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,AVA,W1O1,ONEM,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,APVPA,W1O1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,AVHA,AVA,W1O1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,APVPA,W1O1,ONEM,ONE )
     +
     +
C
C
C             ...create 'O1W1' array.
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVPA,AVHA,O1W1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVPA,APVHPA,O1W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,AVA,AVHA,O1W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVA,APVHPA,O1W1,ONE,ONE )
     +
     +
C
C
C             e)  Calculation of W1*W2:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,SCRH5,W1W2,ONEM,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,SCRH6,W1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,AVHA,SCRH5,W1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,SCRH6,W1W2,ONEM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,SCRH7,W1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVHPA,SCRH8,W1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,SCRH7,W1W2,ONEM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,SCRH8,W1W2,ONE,ONE )
     +
     +
C
C
C             f)  Calculation of W2*W1:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH1,AVHA,W2W1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,APVHPA,W2W1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH2,AVHA,W2W1,ONEM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,APVHPA,W2W1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH3,AVHA,W2W1,ONEM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,APVHPA,W2W1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,SCRH4,AVHA,W2W1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,APVHPA,W2W1,ONEM,ONE )
     +
     +
C
C
C             g)  Calculation of W2^2:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH1,SCRH5,W2W2,ONEM,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRH6,W2W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH2,SCRH5,W2W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRH6,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRH7,W2W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRH8,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRH7,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRH8,W2W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH3,SCRH5,W2W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRH6,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,SCRH4,SCRH5,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRH6,W2W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRH7,W2W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRH8,W2W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRH7,W2W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRH8,W2W2,ONEM,ONE )
     +
     +
C
C
C             h)  Calculation of (W1^3)*W2:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,W1W2,W13W2,ONE,ZERO )
     +
     +
C      
C
C             i)  Calculation of W2*(W1^3):
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W2W1,W1W1,W2W13,ONE,ZERO )
     +
     +
C
C
C             j)  Calculation of W2*O1:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH1,AVA,W2O1,ONEM,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,APVPA,W2O1,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH2,AVA,W2O1,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,APVPA,W2O1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH3,AVA,W2O1,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,APVPA,W2O1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,SCRH4,AVA,W2O1,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,APVPA,W2O1,ONE,ONE )
     +
     +
C
C
C             k)  Calculation of O1*W2:
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVPA,SCRH5,O1W2,ONEM,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVPA,SCRH6,O1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,AVA,SCRH5,O1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVA,SCRH6,O1W2,ONEM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVPA,SCRH7,O1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,APVPA,SCRH8,O1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVA,SCRH7,O1W2,ONEM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVA,SCRH8,O1W2,ONE,ONE )
     +
     +
C
C
C             l)  Calculation of W2*E1*W2:
C
C             ...calculate new scratch-matrices SCRHP_i (i=5,..,8):
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,PE1P,SCRH5,SCRHP5,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,SCRH6,SCRHP6,SCR,ONE,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,SCRH7,SCRHP7,SCR,ONE,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,SCRH8,SCRHP8,SCR,ONE,ZERO,TT,RR )
     +
     +
C 
C
C             ...calculate W2E1W2:
C
C
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRHP5,W2E1W2,SCR,ONEM,ZERO,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRHP6,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRHP5,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRHP6,W2E1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRHP7,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCRHP8,W2E1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRHP7,W2E1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCRHP8,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRHP5,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRHP6,W2E1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRHP5,W2E1W2,ONEM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRHP6,W2E1W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRHP7,W2E1W2,SCR,ONEM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCRHP8,W2E1W2,SCR,ONE,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRHP7,W2E1W2,ONE,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCRHP8,W2E1W2,ONEM,ONE )
     +
     +
C
C
C     3. sum_1 = 1/2 [W2,[W2,E1]] = 1/2 (W2^2)E1 - W2E1W2 + 1/2 E1 (W2^2)
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W2W2,E1,SUM1,HALF,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,W2W2,SUM1,HALF,ONE )
     +
     +
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SUM1,ONEM,W2E1W2,SUM1 )
     +
     +
C
C
C     4. sum_2  =  1/2 [W2,[W1,[W1,O1]]]  +  1/2 [W2,W1O1W1]   =
C
C               =  1/2 W2(W1^2)O1  -  1/2 W2W1O1W1  +  1/2 W2O1(W1^2) 
C
C                  -  1/2 (W1^2)O1W2  +  1/2 W1O1W1W2  -  1/2 O1(W1^2)W2
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W2W1,W1O1,SUM2,HALF,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W2W1,O1W1,SUM2,HALFM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W2O1,W1W1,SUM2,HALF,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,O1W2,SUM2,HALFM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1O1,W1W2,SUM2,HALF,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,O1W1,W1W2,SUM2,HALFM,ONE )
     +
     +
C
C
C     5. sum_3  =  - 1/8 [W1^2,[W1^2,E1]]  = 
C
C               =  - 1/8 (W1^4)E1 + 1/4 (W1^2)E1(W1^2) - 1/8 E1(W1^4)
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W12,E1,SUM3,-EIGHTH,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,W1W12,SUM3,-EIGHTH,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,E1,SCRH1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH1,W1W1,SUM3,QUART,ONE )
     +
     +
C
C
C     6. sum_4  =  a3 [[W2,W1^3],E0]  =  
C
C               =  a3 * ( W2(W1^3)E - EW2(W1^3) - (W1^3)W2E + E(W1^3)W2 )
C
C
         CALL  DKH__MATRHTE
     +
     +              ( NBAS,A3,ZERO,W2W13,SUM4,SCR,EP )
     +
     +
         CALL  DKH__MATLFTE
     +
     +              ( NBAS,-A3,ONE,W2W13,SUM4,SCR,EP )
     +
     +
         CALL  DKH__MATRHTE
     +
     +              ( NBAS,-A3,ONE,W13W2,SUM4,SCR,EP )
     +
     +
         CALL  DKH__MATLFTE
     +
     +              ( NBAS,A3,ONE,W13W2,SUM4,SCR,EP )
     +
     +
C
C
C     7. SCR =  ev5 = sum_1 + sum_2 + sum_3 + sum_4
C
C
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SUM1,ONE,SUM2,SCRH1 )
     +
     +
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SCRH1,ONE,SUM3,SCRH2 )
     +
     +
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SCRH2,ONE,SUM4,SCR )
     +
     +
C
C
C             ...print the potential if the user wants!
C
C
         IF (KNOW .GE. 10) THEN
             STRING = "*** Fifth order potential (p-space) ***"
             CALL  DKH__SQ_PRINT (NBAS,SCR,STRING)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END
