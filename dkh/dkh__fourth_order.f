         SUBROUTINE  DKH__FOURTH_ORDER
     +
     +                   ( KNOW,NBAS,
     +                     NAI,PVP,
     +                     SCR1,SCR2,
     +                     SCR3,SCR4,
     +                     SCRH1,SCRH2,
     +                     SCRH3,SCRH4,
     +                     E1,PE1P,
     +                     AA,RR,
     +                     TT,EP,
     +                     AVA,AVHA,
     +                     APVPA,APVHPA,
     +                     W1W1,W1O1,
     +                     O1W1,SUM1,
     +
     +                                 SCR )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : DKH__FOURTH_ORDER
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__MATMUL_NN
C                DKH__MATMULD
C                DKH__MATMULM
C                DKH__MATDIVE
C                DKH__MATRIXADD
C                DKH__SQ_PRINT
C
C  DESCRIPTION : Perform the fourth order Douglas-Kroll-Hess
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
C                    W1W1         =  DKH specific matrix
C                    W1O1         =  DKH specific matrix
C                    O1W1         =  DKH specific matrix
C                    SUM1         =  scratch matrix used in matrix addition
C
C                  Output:
C
C                    SCR          =  fourth order DKH potential
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
         DOUBLE PRECISION  HALF,HALFM,EIGHTH

         DOUBLE PRECISION  AA (1:NBAS)
         DOUBLE PRECISION  RR (1:NBAS)
         DOUBLE PRECISION  TT (1:NBAS)
         DOUBLE PRECISION  EP (1:NBAS)

         DOUBLE PRECISION  E1  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  NAI (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR (1:NBAS,1:NBAS)

         DOUBLE PRECISION  SCR1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR3 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR4 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SUM1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVHA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1O1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  O1W1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PE1P (1:NBAS,1:NBAS)

         DOUBLE PRECISION  SCRH1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH3 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCRH4 (1:NBAS,1:NBAS)

         DOUBLE PRECISION  APVPA  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  APVHPA (1:NBAS,1:NBAS)

         CHARACTER    STRING*50

         PARAMETER  ( ZERO   =  0.0D0  )
         PARAMETER  ( ONE    =  1.0D0  )
         PARAMETER  ( ONEM   = -1.0D0  )
         PARAMETER  ( HALF   =  0.50D0 )
         PARAMETER  ( HALFM  = -0.50D0 )
         PARAMETER  ( EIGHTH = 0.125D0 )
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
C             ...create 'W1O1' array!
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
C             ...Create 'O1W1' array!
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
C             ...SUM1 = 1/2 [W2,[W1,E1]] = 1/2(W2W1E1-W2E1W1-W1E1W2+E1W1W2)
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,E1,SCR1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,E1,SCR2,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,PE1P,AVHA,SCR3,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,APVHPA,SCR4,SCR,ONE,ZERO,TT,RR )
     +
     +

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
         CALL  DKH__MATDIVE ( NBAS,SCRH1,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH2,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH3,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH4,EP )
C
C
C             ...SUM1 = 1/2 W2W1E1 ([1] - [8])
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH1,SCR1,SUM1,HALF,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCR2,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH2,SCR1,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCR2,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH3,SCR1,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCR2,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,SCRH4,SCR1,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCR2,SUM1,HALFM,ONE )
     +
     +
C
C
C             ...2b)  sum_1 = - 1/2 W2E1W1 (+ sum_1)   ( [9]-[16] )
C
C
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCR3,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH1,SCR4,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCR3,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH2,SCR4,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCR3,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCRH3,SCR4,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCR3,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCRH4,SCR4,SUM1,HALF,ONE  )
     +
     +
C
C
C               ...scr_i & scrh_i  for steps 2c (W1E1W2)
C                  and 2d (E1W1W2):
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
     +              ( NBAS,E1,APVHPA,SCR3,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,E1,AVHA,SCR4,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,AVHA,E1,SCRH1,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,APVHPA,E1,SCRH2,ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,PE1P,AVHA,SCRH3,ONE,ZERO )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,PE1P,APVHPA,SCRH4,SCR,ONE,ZERO,TT,RR )
     +
     +

         CALL  DKH__MATDIVE ( NBAS,SCRH1,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH2,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH3,EP )
         CALL  DKH__MATDIVE ( NBAS,SCRH4,EP )
C
C
C             ...2c)  sum_1 = - 1/2 W1E1W2 (+ sum_1)   ( [17]-[24] )
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR1,SCRH1,SUM1,HALF,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR1,SCRH2,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR2,SCRH1,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR2,SCRH2,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR1,SCRH3,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR1,SCRH4,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR2,SCRH3,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR2,SCRH4,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
C
C
C             ...2d)  sum_1 = 1/2 E1W1W2 (+ sum_1)     ( [25]-[32] )
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR3,SCRH1,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR3,SCRH2,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULM
     +
     +              ( NBAS,SCR4,SCRH1,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR4,SCRH2,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR3,SCRH3,SUM1,SCR,HALF,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMULD
     +
     +              ( NBAS,SCR3,SCRH4,SUM1,SCR,HALFM,ONE,TT,RR )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR4,SCRH3,SUM1,HALFM,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,SCR4,SCRH4,SUM1,HALF,ONE )
     +
     +
C
C
C             3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =  
C
C                  = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,W1O1,SCR1,EIGHTH,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1W1,O1W1,SCR1,-3*EIGHTH,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,W1O1,W1W1,SCR1,3*EIGHTH,ONE )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,O1W1,W1W1,SCR1,-EIGHTH,ONE )
     +
     +
C
C
C             4.  SCR = sum1 + sum2  
C
C
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,ONE,SUM1,ONE,SCR1,SCR )
     +
     +
C
C
C             ...print the potential if the user wants!
C
C
         IF (KNOW .GE. 10) THEN
             STRING = "*** Fourth order potential (p-space) ***"
             CALL  DKH__SQ_PRINT (NBAS,SCR,STRING)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END
