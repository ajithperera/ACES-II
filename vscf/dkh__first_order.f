         SUBROUTINE  DKH__FIRST_ORDER
     +
     +                   ( KNOW,NBAS,
     +                     KIN,NAI,PVP,OVL,
     +                     SCR1,SCR2,
     +                     SHF,OMGA,EW,
     +
     +                             EP,TT,AA,
     +                             RR,TRN,
     +                             PE1P,E1,POT )
     +
     +
C------------------------------------------------------------------------
C  OPERATION   : DKH__FIRST_ORDER
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DCOPY                    !!! ACES3 Specific !!!
C                DIAG                     !!! ACES3 Specific !!!
C                DKH__MATMUL_NN
C                DKH__MATMUL_NT
C                DKH__MATMUL_TN
C                DKH__SQ_PRINT
C
C  DESCRIPTION : Perform the first order Douglas-Kroll-Hess
C                transformation.
C
C
C                  Input:
C
C                    KNOW         =  Level of printing, rather, how much
C                                    the user wants to "KNOW"
C                    NBAS         =  number of atomic orbital basis
C                                    functions
C                    KIN          =  kinetic energy integrals in AO basis
C                    NAI          =  nuclear attraction integrals in AO
C                                    basis
C                    PVP          =  second derivative nuclear attraction
C                                    integrals in the AO basis
C                    OVL          =  overlap integrals in AO basis
C                    SCRx         =  scratch matrix for transformation
C                    SHF          =  overlap matrix raised to the minus 
C                                    one half
C                    OMGA         =  kinetic energy eigenvectors
C                    EW           =  shifted relativistic energies
C
C                  Output:
C
C                    KIN          =  relativistic kinetic energy integrals
C                    NAI          =  transformed nuclear attraction integrals
C                    PVP          =  transformed pVp integrals
C                    OVL          =  eigenvalues of overlap matrix
C                    EP           =  relativistic kinetic energies
C                    TT           =  kinetic energy eigenvalues
C                    AA           =  A matrix in DKH routine
C                    RR           =  R matrix in DKH routine
C                    TRN          =  transformation matrix out of p-space
C                    PE1P         =  transformed pVp potential
C                    E1           =  first order DKH potential
C                    POT          =  first order DKH potential
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

         INTEGER     KNOW
         INTEGER     NBAS,JUNK
         INTEGER     I,J,K

         DOUBLE PRECISION  ZERO,ONE,TWO
         DOUBLE PRECISION  THREE,FOUR
         DOUBLE PRECISION  RATIO,THRESH,VELIT
         DOUBLE PRECISION  CSQRD,CONST,INVC2
         DOUBLE PRECISION  TV1,TV2,TV3,TV4

         PARAMETER  ( ZERO   = 0.000000000 D0 )
         PARAMETER  ( ONE    = 1.000000000 D0 )
         PARAMETER  ( TWO    = 2.000000000 D0 )
         PARAMETER  ( THREE  = 3.000000000 D0 )
         PARAMETER  ( FOUR   = 4.000000000 D0 )
         PARAMETER  ( VELIT  = 137.0359895 D0 )
         PARAMETER  ( THRESH = 0.020000000 D0 )

         CHARACTER    STRING*50

         DOUBLE PRECISION  EP (1:NBAS)
         DOUBLE PRECISION  EW (1:NBAS)
         DOUBLE PRECISION  AA (1:NBAS)
         DOUBLE PRECISION  RR (1:NBAS)
         DOUBLE PRECISION  TT (1:NBAS)

         DOUBLE PRECISION  KIN  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  NAI  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  OVL  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SHF  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  TRN  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  OMGA (1:NBAS,1:NBAS)
         DOUBLE PRECISION  E1   (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PE1P (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR1 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  SCR2 (1:NBAS,1:NBAS)
         DOUBLE PRECISION  POT  (1:NBAS,1:NBAS)
C
C
C------------------------------------------------------------------------
C
C
C                                                     -1/2
C             ...diagonalize the overlap matrix, form S
C
C
         CALL  DCOPY
     +
     +              ( NBAS*NBAS,OVL,1,
     +
     +                           SCR1,1 )
     +
     +
         CALL  EIG
     +
     +              ( SCR1,OMGA,
     +                JUNK,NBAS, 1 )
     +
     +
         DO I = 1, NBAS
            SCR1 (I,I) = ONE / DSQRT (SCR1 (I,I))
         ENDDO

         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                OMGA,SCR1,SCR2,
     +                ONE,ZERO       )
     +
     +
         CALL  DKH__MATMUL_NT
     +
     +              ( NBAS,
     +                SCR2,OMGA,SHF,
     +                ONE,ZERO       )
     +
     +
C
C
C             ...rotate normalize the kinetic energy matrix
C             ...diagonalize the kinetic energy matrix and
C                retain the eigenvectors for transformations
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                KIN,SHF,SCR1,
     +                ONE,ZERO       )
     +
     +
         CALL  DKH__MATMUL_TN
     +
     +              ( NBAS,
     +                SHF,SCR1,SCR2,
     +                ONE,ZERO       )
     +
     +
         CALL  EIG
     +
     +              ( SCR2,OMGA,
     +                JUNK,NBAS,-1 )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                SHF,OMGA,TRN,
     +                ONE,ZERO      )
     +
     +
C
C
C             ...form the relativistic energies!
C
C
         INVC2 = ONE / (VELIT * VELIT)
         CONST = INVC2 + INVC2
         CSQRD = ONE / INVC2
         DO I = 1, NBAS

            RATIO = SCR2 (I,I) / VELIT
            IF (RATIO .LT. THRESH) THEN
                TV1 = SCR2 (I,I)
                TV2 = -TV1 * SCR2 (I,I) * INVC2 / TWO
                TV3 = -TV2 * SCR2 (I,I) * INVC2
                TV4 = -TV3 * SCR2 (I,I) * INVC2 * 1.25D0
                EW (I) = TV1 + TV2 + TV3 + TV4
            ELSE
                EW (I) = CSQRD * (DSQRT (ONE + CONST*SCR2 (I,I)) - ONE)
            ENDIF

            EP (I) = EW (I) + CSQRD
            AA (I) = DSQRT ((CSQRD + EP (I)) / (TWO * EP (I)))
            RR (I) = VELIT / (CSQRD + EP (I))
            TT (I) = SCR2 (I,I)

         ENDDO
C
C
C             ...transform the nuclear attraction and 
C                second derivative nuclear attraction 
C                AO integrals!
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                NAI,TRN,SCR1,
     +                ONE,ZERO       )
     +
     +
         CALL  DKH__MATMUL_TN
     +
     +              ( NBAS,
     +                TRN,SCR1,NAI,
     +                ONE,ZERO       )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                PVP,TRN,SCR1,
     +                ONE,ZERO       )
     +
     +
         CALL  DKH__MATMUL_TN
     +
     +              ( NBAS,
     +                TRN,SCR1,PVP,
     +                ONE,ZERO       )
     +
     +
C
C
C             ...form the DKH potentials!
C
C
         DO I = 1, NBAS
         DO J = 1, NBAS
            E1   (I,J) = AA (I) * NAI (I,J) * AA (J)
            PE1P (I,J) =       FOUR * TT (I) * RR (I) * RR (I)
     +                   * E1 (I,J) * TT (J) * RR (J) * RR (J)
         ENDDO
         ENDDO

         DO I = 1, NBAS
         DO J = 1, NBAS
            E1   (I,J) =   E1  (I,J) + AA (I) * RR (I)
     +                   * PVP (I,J) * RR (J) * AA (J)
            PE1P (I,J) =  PE1P (I,J) + AA (I) * RR (I)
     +                   * PVP (I,J) * RR (J) * AA (J)
         ENDDO
         ENDDO
C
C
C             ...compute the reverse transformation.
C                currently, TRN contains the inverse square
C                root of the overlap matrix times OMGA, so that
C                times the overlap matrix is the square root overlap
C                times OMGA!
C
C
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,
     +                OVL,TRN,SCR1,
     +                ONE,ZERO       )
     +
     +
         CALL  DCOPY
     +   
     +              ( NBAS*NBAS,SCR1,1,
     +   
     +                           TRN,1 )
     +   
     +   
C
C
C             ...backtransform the kinetic energy!
C
C
         DO I = 1, NBAS
         DO J = 1, NBAS
            KIN (I,J) = 0.0D0
            DO K = 1, NBAS
               KIN (I,J) =   KIN (I,J)
     +                     + TRN (I,K) * TRN (J,K) * EW (K)
            END DO
         END DO
         END DO
C
C
C             ...copy the E1 potential into POT array!
C
C
         CALL  DCOPY
     +
     +              ( NBAS*NBAS,E1,1,
     +
     +                           POT,1 )
     +
     +
C
C
C             ...print information if user wants it!
C
C
         IF (KNOW .GE. 10) THEN
             CALL  DKH__SQ_PRINT (NBAS,NAI,"*** NAI (p-space) ***")
             CALL  DKH__SQ_PRINT (NBAS,PVP,"*** pVp (p-space) ***")
             CALL  DKH__SQ_PRINT
     +             (NBAS,KIN,"*** Relativistic Kinetic Energy ***")

             STRING = "*** First order potential (p-space) ***"
             CALL  DKH__SQ_PRINT (NBAS,E1,STRING)
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END
