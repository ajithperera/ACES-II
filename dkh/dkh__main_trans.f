         SUBROUTINE  DKH__MAIN_TRANS
     +
     +                    ( ORDER,NBAS,
     +                      KNOW,ZMAX,
     +                      ZCORE,
     +                      OVL,NAI,PVP,
     +
     +                                ERR,
     +                                KIN)
C------------------------------------------------------------------------
C  OPERATION   : DKH__START_TRANS
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__SQ_PRINT
C                DKH__FIRST_ORDER
C                DKH__SECOND_ORDER
C                DKH__THIRD_ORDER
C                DKH__FOURTH_ORDER
C                DKH__FIFTH_ORDER
C                DKH__MATRIXADD
C                DKH__MATMUL_NT            !!! ACES3 specific !!!
C                DKH__MATMUL_NN            !!! ACES3 specific !!!
C
C  DESCRIPTION : Main operation that drives the Douglas-Kroll-Hess
C                transformation of the one electron hamiltonian to a
C                given order (ORDER).
C
C                This routine will transform the kinetic energy
C                integrals and the nuclear attraction integrals,
C                which should be added together to form the DKH
C                Hamiltonian.  However, it is formed in the present
C                routine and stored, on exit, in ZCORE (1).
C
C
C                  Input:
C
C                    ORDER        =  order of the transformation (1-5)
C                    NBAS         =  number of atomic orbital basis
C                                    functions
C                    KNOW         =  level of printing
C                    ZMAX         =  maximum floating point memory
C                    ZCORE        =  double precision scratch space
C                    OVL          =  overlap integrals in AO basis,
C                                    modified on exit
C                    PVP          =  second derivative nuclear attraction
C                                    integrals in the AO basis,
C                                    modified on exit
C                    KIN          =  kinetic energy integrals in AO basis,
C                                    modified on exit
C                    NAI          =  nuclear attraction integrals in AO
C                                    basis, modified on exit
C
C                  Output:
C
C                    ERR          =  error status; zero if everything
C                                    executed normally, one if not
C                    KIN          =  relativistic kinetic energy integrals
C                    NAI          =  transformed nuclear attraction integrals
C                    ZCORE(1)     =  beginning at ZCORE(1), the full
C                                    Douglas-Kroll-Hess one electron
C                                    Hamiltonian is stored;
C                                    Length = NBAS * NBAS
C
C
C
C                            !!! IMPORTANT !!!
C
C                The input arrays containing the overlap, kinetic,
C                nuclear attraction, second derivative nuclear
C                attraction atomic orbital integrals are all modified
C                on exit.
C
C                          OVL will be unchanged on exit
C                          KIN will contain relativistic KIN integrals
C                          NAI will contain the DKH potential
C                          PVP will be in p-space
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

         INTEGER     ORDER,NBAS
         INTEGER     ZMAX,ERR,KNOW

         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  KIN  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  NAI  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  PVP  (1:NBAS,1:NBAS)
         DOUBLE PRECISION  OVL  (1:NBAS,1:NBAS)

         DOUBLE PRECISION  ZCORE (1:ZMAX)

         INTEGER     ZEP,ZEW,ZAA,ZRR,
     +               ZTT,ZSHF,ZPOT,ZE1,
     +               ZPE1P,ZGG,ZPP,ZTRN,
     +               ZOMGA,ZAVA,ZAPVPA,
     +               ZAVHA,ZAPVHPA,ZSCR,
     +               ZW1O1,ZO1W1,ZW1W1,
     +               ZSCR1,ZSCR2,ZSCR3,
     +               ZSCR4,ZSCRH1,ZSCRH2,
     +               ZSCRH3,ZSCRH4,ZSCRH5,
     +               ZSCRH6,ZSCRH7,ZSCRH8,
     +               ZSCRHP5,ZSCRHP6,ZSCRHP7,
     +               ZSCRHP8,ZW1W2,ZW2W1,ZW2W2,
     +               ZW2O1,ZO1W2,ZW1W12,ZW2W13,
     +               ZW13W2,ZW1E1W1,ZW2E1W2,
     +               ZSUM1,ZSUM2,ZSUM3,ZSUM4,ZEND

         PARAMETER  ( ZERO = 0.0D0 )
         PARAMETER  ( ONE  = 1.0D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...set up the ERR variable and print information
C                depending on what the user wants!
C
C
         ERR = 0

         IF (KNOW .GE. 1) THEN
            WRITE (*,*) ""
            WRITE (*,'(X,A29,A24,I2)') "Performing Douglas-Kroll-Hess",
     +                                 "transformation to order",ORDER
         ENDIF

         IF (KNOW .GE. 20) THEN

             CALL DKH__SQ_PRINT
     +
     +               ( NBAS,KIN,"*** Kinetic Energy Integrals ***" )
     +
     +
             CALL DKH__SQ_PRINT
     +
     +               ( NBAS,NAI,"*** Nuclear Attraction Integrals ***" )
     +
     +
             CALL DKH__SQ_PRINT
     +
     +               ( NBAS,OVL,"*** Overlap Integrals ***" )
     +
     +
             CALL DKH__SQ_PRINT
     +
     +               ( NBAS,PVP,"*** pVp Integrals ***" )
     +
     +
         ENDIF
C
C
C             ...partition the memory for the
C                first order transformation
C
C
         ZEP   =  1
         ZTT   =  ZEP   +  NBAS
         ZAA   =  ZTT   +  NBAS
         ZRR   =  ZAA   +  NBAS
         ZE1   =  ZRR   +  NBAS
         ZPE1P =  ZE1   +  NBAS * NBAS
         ZTRN  =  ZPE1P +  NBAS * NBAS
         ZPOT  =  ZTRN  +  NBAS * NBAS

         ZSCR1 =  ZPOT  +  NBAS * NBAS
         ZSCR2 =  ZSCR1 +  NBAS * NBAS
         ZSHF  =  ZSCR2 +  NBAS * NBAS
         ZOMGA =  ZSHF  +  NBAS * NBAS
         ZEW   =  ZOMGA +  NBAS * NBAS
         ZEND  =  ZEW   +  NBAS

         IF ( (ZEND - ZEP) .GT. ZMAX ) THEN
            WRITE (*,*) "Error in allocating memory in ",
     +                  "DKH__MAIN_TRANS in DKH-1:"
            WRITE (*,*) "  Available: ", ZMAX
            WRITE (*,*) "  Needed:    ", ZEND - ZEP
            ERR = 1
            RETURN
         ENDIF
C
C
C             ...perform the first order transformation.
C
C
         IF (KNOW .GE. 1) THEN
            WRITE (*,'(/,A38)') "Performing first order transformation."
         ENDIF

         CALL  DKH__FIRST_ORDER
     +
     +              ( KNOW,NBAS,
     +                KIN,NAI,PVP,OVL,
     +                ZCORE (ZSCR1),ZCORE (ZSCR2),
     +                ZCORE (ZSHF), ZCORE (ZOMGA),
     +                ZCORE (ZEW),
     +
     +                        ZCORE (ZEP),  ZCORE (ZTT),
     +                        ZCORE (ZAA),  ZCORE (ZRR),
     +                        ZCORE (ZTRN), ZCORE (ZPE1P),
     +                        ZCORE (ZE1),  ZCORE (ZPOT)    )
     +
     +
C
C
C             ...Allocate memory for second order transformation.
C                Then perform the transformation!
C
C
         IF (ORDER .GE. 2) THEN

             ZAVA    =  ZPOT    +  NBAS * NBAS
             ZAPVPA  =  ZAVA    +  NBAS * NBAS
             ZAVHA   =  ZAPVPA  +  NBAS * NBAS
             ZAPVHPA =  ZAVHA   +  NBAS * NBAS
             ZW1O1   =  ZAPVHPA +  NBAS * NBAS
             ZO1W1   =  ZW1O1   +  NBAS * NBAS
             ZSCR    =  ZO1W1   +  NBAS * NBAS
             ZEND    =  ZSCR    +  NBAS * NBAS

             IF ( (ZEND - ZEP) .GT. ZMAX ) THEN
                 WRITE (*,*) "Error in allocating memory in ",
     +                       "DKH__MAIN_TRANS for DKH-2:"
                 WRITE (*,*) "  Available: ", ZMAX
                 WRITE (*,*) "  Needed:    ", ZEND - ZEP
                 ERR = 1
                 RETURN
             ENDIF

             IF (KNOW .GE. 1) THEN
                 WRITE (*,*) "Performing second order transformation."
             ENDIF

             CALL  DKH__SECOND_ORDER
     +
     +                  ( KNOW,NBAS,
     +                    NAI,PVP,
     +                    ZCORE (ZAA),  ZCORE (ZRR),
     +                    ZCORE (ZTT),  ZCORE (ZEP),
     +                    ZCORE (ZAVA), ZCORE (ZAPVPA),
     +                    ZCORE (ZAVHA),ZCORE (ZAPVHPA),
     +                    ZCORE (ZW1O1),ZCORE (ZO1W1),
     +
     +                                         ZCORE (ZSCR) )
     +
     +
             CALL  DKH__MATRIXADD
     +
     +                  ( NBAS*NBAS,ONE,ZCORE (ZPOT),
     +                    ONE,ZCORE (ZSCR),
     +
     +                                       ZCORE (ZPOT) )
     +
     +
         ENDIF
C        
C           
C             ...Allocate memory for third order transformation.
C                Then perform the transformation!
C    
C    
         IF (ORDER .GE. 3) THEN

             ZSCR1   =  ZPOT    +  NBAS * NBAS
             ZSCR2   =  ZSCR1   +  NBAS * NBAS
             ZAVHA   =  ZSCR2   +  NBAS * NBAS
             ZAPVHPA =  ZAVHA   +  NBAS * NBAS
             ZW1W1   =  ZAPVHPA +  NBAS * NBAS
             ZW1E1W1 =  ZW1W1   +  NBAS * NBAS
             ZSCR    =  ZW1E1W1 +  NBAS * NBAS
             ZEND    =  ZSCR    +  NBAS * NBAS

             IF ( (ZEND - ZEP) .GT. ZMAX ) THEN
                 WRITE (*,*) "Error in allocating memory in ",
     +                       "DKH__MAIN_TRANS for DKH-3:"
                 WRITE (*,*) "  Available: ", ZMAX
                 WRITE (*,*) "  Needed:    ", ZEND - ZEP
                 ERR = 1
                 RETURN 
             ENDIF

             IF (KNOW .GE. 1) THEN
                 WRITE (*,*) "Performing third order transformation."
             ENDIF

             CALL  DKH__THIRD_ORDER
     +       
     +                  ( KNOW,NBAS,
     +                    NAI,PVP,
     +                    ZCORE (ZSCR1),ZCORE (ZSCR2),
     +                    ZCORE (ZE1),  ZCORE (ZPE1P),
     +                    ZCORE (ZAA),  ZCORE (ZRR),
     +                    ZCORE (ZTT),  ZCORE (ZEP),
     +                    ZCORE (ZAVHA),ZCORE (ZAPVHPA),
     +                    ZCORE (ZW1W1),ZCORE (ZW1E1W1),
     +           
     +                                       ZCORE (ZSCR) )
     +           
     +           
             CALL  DKH__MATRIXADD
     +                    
     +                  ( NBAS*NBAS,ONE,ZCORE (ZPOT),
     +                    ONE,ZCORE (ZSCR),
     +                    
     +                                      ZCORE (ZPOT) )
     +                    
     +                    
         ENDIF
C
C           
C             ...Allocate memory for fourth order transformation.
C                Then perform the transformation!
C    
C    
         IF (ORDER .GE. 4) THEN

             ZSCR1   =  ZPOT    +  NBAS * NBAS
             ZSCR2   =  ZSCR1   +  NBAS * NBAS
             ZSCR3   =  ZSCR2   +  NBAS * NBAS
             ZSCR4   =  ZSCR3   +  NBAS * NBAS
             ZSCRH1  =  ZSCR4   +  NBAS * NBAS
             ZSCRH2  =  ZSCRH1  +  NBAS * NBAS
             ZSCRH3  =  ZSCRH2  +  NBAS * NBAS      
             ZSCRH4  =  ZSCRH3  +  NBAS * NBAS
             ZAVA    =  ZSCRH4  +  NBAS * NBAS
             ZAVHA   =  ZAVA    +  NBAS * NBAS
             ZAPVPA  =  ZAVHA   +  NBAS * NBAS
             ZAPVHPA =  ZAPVPA  +  NBAS * NBAS
             ZW1W1   =  ZAPVHPA +  NBAS * NBAS
             ZW1O1   =  ZW1W1   +  NBAS * NBAS
             ZO1W1   =  ZW1O1   +  NBAS * NBAS
             ZSUM1   =  ZO1W1   +  NBAS * NBAS
             ZSCR    =  ZSUM1   +  NBAS * NBAS
             ZEND    =  ZSCR    +  NBAS * NBAS

             IF ( (ZEND - ZEP) .GT. ZMAX ) THEN
                 WRITE (*,*) "Error in allocating memory in ",
     +                       "DKH__MAIN_TRANS for DKH-4:"
                 WRITE (*,*) "  Available: ", ZMAX
                 WRITE (*,*) "  Needed:    ", ZEND - ZEP
                 ERR = 1
                 RETURN
             ENDIF

             IF (KNOW .GE. 1) THEN
                 WRITE (*,*) "Performing fourth order transformation."
             ENDIF

             CALL  DKH__FOURTH_ORDER
     +      
     +                  ( KNOW,NBAS,
     +                    NAI,PVP,
     +                    ZCORE (ZSCR1),  ZCORE (ZSCR2),
     +                    ZCORE (ZSCR3),  ZCORE (ZSCR4),
     +                    ZCORE (ZSCRH1), ZCORE (ZSCRH2),
     +                    ZCORE (ZSCRH3), ZCORE (ZSCRH4),
     +                    ZCORE (ZE1),    ZCORE (ZPE1P),
     +                    ZCORE (ZAA),    ZCORE (ZRR),
     +                    ZCORE (ZTT),    ZCORE (ZEP),
     +                    ZCORE (ZAVA),   ZCORE (ZAVHA),
     +                    ZCORE (ZAPVPA), ZCORE (ZAPVHPA),
     +                    ZCORE (ZW1W1),  ZCORE (ZW1O1),
     +                    ZCORE (ZO1W1),  ZCORE (ZSUM1),
     +
     +                                       ZCORE (ZSCR) )
     +
     +
             CALL  DKH__MATRIXADD
     +
     +                  ( NBAS*NBAS,ONE,ZCORE (ZPOT),
     +                    ONE,ZCORE (ZSCR),
     +
     +                                       ZCORE (ZPOT) )
     +
     +
         ENDIF
C
C           
C             ...Allocate memory for fifth order transformation.
C                Then perform the transformation!
C    
C        
         IF (ORDER .GE. 5) THEN

             ZSCRH1  =  ZPOT    +  NBAS * NBAS
             ZSCRH2  =  ZSCRH1  +  NBAS * NBAS
             ZSCRH3  =  ZSCRH2  +  NBAS * NBAS
             ZSCRH4  =  ZSCRH3  +  NBAS * NBAS
             ZSCRH5  =  ZSCRH4  +  NBAS * NBAS
             ZSCRH6  =  ZSCRH5  +  NBAS * NBAS
             ZSCRH7  =  ZSCRH6  +  NBAS * NBAS
             ZSCRH8  =  ZSCRH7  +  NBAS * NBAS
             ZSCRHP5 =  ZSCRH8  +  NBAS * NBAS
             ZSCRHP6 =  ZSCRHP5 +  NBAS * NBAS
             ZSCRHP7 =  ZSCRHP6 +  NBAS * NBAS
             ZSCRHP8 =  ZSCRHP7 +  NBAS * NBAS
             ZAVA    =  ZSCRHP8 +  NBAS * NBAS
             ZAVHA   =  ZAVA    +  NBAS * NBAS
             ZAPVPA  =  ZAVHA   +  NBAS * NBAS
             ZAPVHPA =  ZAPVPA  +  NBAS * NBAS
             ZW1W1   =  ZAPVHPA +  NBAS * NBAS
             ZW1O1   =  ZW1W1   +  NBAS * NBAS
             ZO1W1   =  ZW1O1   +  NBAS * NBAS
             ZW1W2   =  ZO1W1   +  NBAS * NBAS
             ZW2W1   =  ZW1W2   +  NBAS * NBAS
             ZW2W2   =  ZW2W1   +  NBAS * NBAS
             ZW2O1   =  ZW2W2   +  NBAS * NBAS
             ZO1W2   =  ZW2O1   +  NBAS * NBAS
             ZW1W12  =  ZO1W2   +  NBAS * NBAS
             ZW2W13  =  ZW1W12  +  NBAS * NBAS
             ZW13W2  =  ZW2W13  +  NBAS * NBAS
             ZW2E1W2 =  ZW13W2  +  NBAS * NBAS
             ZSUM1   =  ZW2E1W2 +  NBAS * NBAS
             ZSUM2   =  ZSUM1   +  NBAS * NBAS
             ZSUM3   =  ZSUM2   +  NBAS * NBAS
             ZSUM4   =  ZSUM3   +  NBAS * NBAS
             ZSCR    =  ZSUM4   +  NBAS * NBAS
             ZEND    =  ZSCR    +  NBAS * NBAS

             IF ( (ZEND - ZEP) .GT. ZMAX ) THEN
                 WRITE (*,*) "Error in allocating memory in ",
     +                       "DKH__MAIN_TRANS for DKH-5:"
                 WRITE (*,*) "  Available: ", ZMAX
                 WRITE (*,*) "  Needed:    ", ZEND - ZEP
                 ERR = 1
                 RETURN 
             ENDIF

             IF (KNOW .GE. 1) THEN
                 WRITE (*,*) "Performing fifth order transformation."
             ENDIF

             CALL  DKH__FIFTH_ORDER
     +       
     +                  ( KNOW,NBAS,
     +                    NAI,PVP, 
     +                    ZCORE (ZSCRH1),
     +                    ZCORE (ZSCRH2), ZCORE (ZSCRH3),
     +                    ZCORE (ZSCRH4), ZCORE (ZSCRH5),
     +                    ZCORE (ZSCRH6), ZCORE (ZSCRH7),
     +                    ZCORE (ZSCRH8), ZCORE (ZSCRHP5),
     +                    ZCORE (ZSCRHP6),ZCORE (ZSCRHP7),
     +                    ZCORE (ZSCRHP8),ZCORE (ZE1),
     +                    ZCORE (ZPE1P),  ZCORE (ZAA),
     +                    ZCORE (ZRR),    ZCORE (ZTT),
     +                    ZCORE (ZEP),    ZCORE (ZAVA),
     +                    ZCORE (ZAVHA),  ZCORE (ZAPVPA),
     +                    ZCORE (ZAPVHPA),ZCORE (ZW1W1),
     +                    ZCORE (ZW1O1),  ZCORE (ZO1W1),
     +                    ZCORE (ZW1W2),  ZCORE (ZW2W1),
     +                    ZCORE (ZW2W2),  ZCORE (ZW2O1),
     +                    ZCORE (ZO1W2),  ZCORE (ZW1W12),
     +                    ZCORE (ZW2W13), ZCORE (ZW13W2),
     +                    ZCORE (ZW2E1W2),ZCORE (ZSUM1),
     +                    ZCORE (ZSUM2),  ZCORE (ZSUM3),
     +                    ZCORE (ZSUM4),
     +                    
     +                                            ZCORE (ZSCR) )
     +                                               
     +
             CALL  DKH__MATRIXADD
     +
     +                  ( NBAS*NBAS,ONE,ZCORE (ZPOT),
     +                    ONE,ZCORE (ZSCR),
     +
     +                                            ZCORE (ZPOT) )
     +
     +
         ENDIF
C
C
C             ...transform the potential back to the original
C                basis and then add it and the kinetic energy
C                to finally form the DKH transformed one
C                electron Hamiltonian!
C
C
         ZSCR = ZPOT + NBAS * NBAS

         CALL  DKH__MATMUL_NT
     +
     +              ( NBAS,ZCORE (ZPOT),ZCORE (ZTRN),
     +
     +                              ZCORE (ZSCR),ONE,ZERO )
     +
     +
         CALL  DKH__MATMUL_NN
     +
     +              ( NBAS,ZCORE (ZTRN),ZCORE (ZSCR),
     +
     +                                    NAI,ONE,ZERO )
     +
     +
C
C
C             ...store the DKH Hamiltonian in case that is what
C                the user would use instead of the KIN and NAI
C                transformed integrals!
C
C

         CALL  DKH__MATRIXADD
     +
     +              ( NBAS*NBAS,
     +                ONE,KIN,
     +                ONE,NAI,
     +
     +                       ZCORE (1) )
     +
     +
C
C
C             ...print the DKH Hamiltonian, if the user wants!
C
C
         IF (KNOW .GE. 5) THEN

             CALL DKH__SQ_PRINT
     +
     +               ( NBAS,ZCORE (1),"*** DKH Hamiltonian ***" )
     +
     +
             WRITE (*,*) ""
         ENDIF
C
C
C             ...ready!
C
C
         RETURN
         END

