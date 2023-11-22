         SUBROUTINE  DKH__SECOND_ORDER
     +
     +                   ( KNOW,NBAS,
     +                     NAI,PVP,AA,
     +                     RR,TT,EP,
     +                     AVA,APVPA,
     +                     AVHA,APVHPA,
     +                     W1O1,O1W1,
     +
     +                               SCR )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : DKH__SECOND_ORDER
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__MATMUL_NN
C                DKH__MATMULD
C                DKH__MATMULM
C                DKH__MATRIXADD
C                DKH__SQ_PRINT
C
C  DESCRIPTION : Perform the second order Douglas-Kroll-Hess
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
C                    AA           =  A matrix in DKH routine
C                    RR           =  R matrix in DKH routine
C                    TT           =  kinetic energy eigenvalues
C                    EP           =  relativistic kinetic energies
C                    AVA          =  transformed NAI (p-space) array
C                    APVPA        =  transformed pVp (p-space) array
C                    AVHA         =  AVA divided by energy denominator
C                    APVHPA       =  APVPA divided by energy denominator
C                    W1O1         =  DKH specific matrix
C                    O1W1         =  DKH specific matrix
C
C                  Output:
C
C                    SCR          =  second order DKH potential
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

         DOUBLE PRECISION  NAI (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR (1:NBAS,1:NBAS)
         DOUBLE PRECISION  AVA (1:NBAS,1:NBAS)

         DOUBLE PRECISION  AVHA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  W1O1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  O1W1 (1:NBAS,1:NBAS)

         DOUBLE PRECISION  APVPA  (1:NBAS,1:NBAS)
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
C             ...commpute 1/2 [W1, O1] and squeeze
C                final answer into 'SCR'!
C
C
         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,HALF,W1O1,HALFM,O1W1,SCR )
     +
     +
C
C
C             ...print the potential if the user wants!
C
C
         IF (KNOW .GE. 10) THEN
             STRING = "*** Second order potential (p-space) ***"
             CALL  DKH__SQ_PRINT (NBAS,SCR,STRING)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END

