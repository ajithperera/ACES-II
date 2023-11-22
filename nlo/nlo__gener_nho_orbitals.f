         SUBROUTINE  NLO__GENER_NHO_ORBITALS
     +
     +                    ( NBAS,NATOM,N2CEN,
     +                      MXNHBA,
     +                      MAXOCC,BONDSIZE,
     +                      NHB,NHA,
     +                      ATORD,AT2CEN,ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                      RYD2HYB,
     +                      COLMAP,
     +                      NHYB,
     +                      NWSTEP,
     +                      P,
     +                      SAH,
     +                      PH,PHSUB,PHDEP,
     +                      BOMAT,
     +                      WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      IVEC,XVEC,
     +                      XMAT,
     +
     +                             NBOND,
     +                             NLPAIR,NEPAIR,
     +                             MXNCEN,
     +                             WBOND,WSTAR,
     +                             NHOLP,
     +                             NHOEP,
     +                             NHORYD,
     +                             BDNCEN,
     +                             BDCEN,
     +                             BDBAS,
     +                             BDOCC,
     +                             H,
     +                             W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NHO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the sets of natural hybrid
C                orbitals (NHO's) as a set of coefficients in terms
C                of the NHB Hybrid NAOs. The procedure by which all
C                NHB Hybrid NAOs are transformed into NHOs is as
C                follows:
C
C
C                   1) Calculate the NHB x NHB occupation matrix PH
C                      in the Hybrid NAO basis. This matrix will
C                      have the following structure:
C
C
C                             A     B     C     D
C                           -----------------------
C                          |     |     |     |     |
C                        A |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        B |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        C |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C                          |     |     |     |     |
C                        D |     |     |     |     |
C                          |     |     |     |     |
C                           -----------------------
C
C
C                      where A,B,C,D denote atoms.
C
C                   2) Perform the following loop, until the sets of
C                      NHOs found exhaust the NHB space:
C
C                          # of centers NHCEN = 1
C
C                          loop over all NHCEN atomic combinations
C
C                             for each NHCEN atomic combination:
C
C                               i) set up submatrix PHSUB of PH
C                                  by taking appropriate subblocks
C                                  of PH corresponding to the NHCEN
C                                  atomic combination
C
C                              ii) Diagonalize PHSUB and extract
C                                  those eigenfunctions which have
C                                  eigenvalues > Threshold value
C
C                             iii) Deplete PHSUB from those
C                                  eigenfunctions
C
C                              iv) Replace PHSUB section of PH
C                                  by PHSUB matrix
C
C                             if # of NHOs complete, exit
C
C                          continue
C
C                          remove all atomic sites which have by
C                          now a complete set of NHOs.
C
C                          if # of NHOs incomplete, NHCEN = NHCEN + 1
C                          and return to loop beginning
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    N2CEN        =  # of atomic center pairs that
C                                    will be considered for 2-center
C                                    bond construction.
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NHB          =  total # of Hybrid NAOs
C                    NHA          =  total # of hybrid atoms
C                    ATORD (I)    =  I-th atomic index of optimum
C                                    atomic index ordering array of
C                                    all hybrid atoms to be used for
C                                    bond/antibond formation search.
C                    AT2CEN (1,N) =  1st atomic center label of N-th
C                                    2-center pair.
C                    AT2CEN (2,N) =  2nd atomic center label of N-th
C                                    2-center pair.
C                    ATHCEN (I)   =  will contain hybrid atomic labels
C                                    (indices) for a subset of all NHA
C                                    atomic hybrid centers
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHVAL (A)   =  # of Valence NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid NAOs
C                                    on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th NAO, based
C                                    on atomic order, in the NHB/NCB/NRB
C                                    order.
C                    NHYB (A)     =  will contain total # of pre-NHOs
C                                    on each hybrid atom A at each stage
C                                    of the calculation.
C                    NWSTEP (x)   =  will contain the # of steps for
C                                    decreasing/increasing the weight
C                                    acceptance criterion for
C                                    bond/antibond formation for bonds
C                                    over x centers.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    SAH          =  will contain the atomic hybrid
C                                    overlap matrices between pre-NHOs.
C                    PH           =  will contain the (possibly depleted)
C                                    NHB x NHB hybrid occupation matrix
C                    PHSUB        =  will contain the submatrices of
C                                    the occupation matrix PH for
C                                    finding the NHOs.
C                    PHDEP        =  will be used for accumulation of
C                                    the  depleted occupation matrix PH.
C                    BOMAT (A,B)  =  simplified bond order matrix
C                                    containing info if atoms A and B
C                                    are considered to be bonded or not
C                    WBDMIN       =  initial minimum accepted weight
C                                    for bond construction.
C                    WSTMAX       =  initial maximum accepted weight
C                                    for antibond construction.
C                    WBDCRT       =  critical weight limit for bond
C                                    construction.
C                    WSTEP        =  weight stepping size to decrease/
C                                    increase the bond/antibond weight
C                                    limits.
C                    MXCHOOSE     =  maximum # of bonds selected to
C                                    be chosen. The maximum is build
C                                    from all # of chosen bonds for all
C                                    bondsizes.
C                    NCHOOSE (B)  =  # of bonds to be chosen for bonds
C                                    of size B. Four cases:
C                                    1) = 9999 => skip search for bonds
C                                       of size B.
C                                    2) = 0 => complete search for all
C                                       possible bonds of size B will
C                                       be performed.
C                                    3) = -n => only n bonds of size B
C                                       will be searched between those
C                                       atomic indices as provided by
C                                       the CHOOSE array.
C                                    4) = +n => same as case 3) but
C                                       followed by a complete search
C                                       for all possible remaining
C                                       bonds of size B.
C                                    Priority level: 1) > 3) > 4) > 2).
C                    CHOOSE       =  element CHOOSE (I,N,B) contains
C                                    the I-th atomic index of the N-th
C                                    chosen bond of size B. The order
C                                    of the atomic indices is arbitrary.
C                    IVEC,XVEC    =  int/flp scratch array of vector
C                                    type.
C                    XMAT         =  flp scratch array of matrix type.
C                    W            =  NAO weight vector in atomic order.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in atomic order.
C
C
C                  Output:
C
C                    NBOND        =  total # of bonds found.
C                    NLPAIR       =  total # of Lone-pair type bonds
C                                    found.
C                    NEPAIR       =  total # of Empty-pair type bonds
C                                    found.
C                    MXNCEN       =  maximum # of atomic hybrid centers
C                                    per bond.
C                    WBOND (x)    =  will contain the weight acceptance
C                                    criterion for x-center bond
C                                    formation during the NHO search
C                                    process. The final values when
C                                    exiting this routine will be the
C                                    lowest for each x-center bond
C                                    that lead to the final complete
C                                    set of NHO.
C                    WSTAR (x)    =  will contain the weight acceptance
C                                    criterion for x-center antibond
C                                    formation during the NHO search
C                                    process. The final values when
C                                    exiting this routine will be the
C                                    highest for each x-center antibond
C                                    that lead to the final complete
C                                    set of NHO.
C                    NHOLP (A)    =  total # of Lone-pair NHOs found
C                                    on each atom A.
C                    NHOEP (A)    =  total # of Empty-pair NHOs found
C                                    on each atom A.
C                    NHORYD (A)   =  total # of Rydberg NHOs found on
C                                    each atom A.
C                    BDNCEN (I)   =  # of atomic centers for I-th bond.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond.
C                    BDOCC (I)    =  # of occupied levels for I-th bond.
C                    H (I,J)      =  MXNHBA x NHB matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB NHOs,
C                                    with all NHOs belonging to a
C                                    specific atomic center being
C                                    grouped together.
C                    W            =  NAO weight vector in NHB/NCB/NRB
C                                    order.
C                    C            =  NAO coefficient matrix in AO basis
C                                    with columns in NHB/NCB/NRB order.
C
C
C                  Note that the concept of a bond is here synonymous
C                  with atomic center connectivity. Although a specific
C                  atomic center connection might lead to many occupied
C                  'bonds' in the usual chemists point of view, all
C                  those occupied 'bonds' are treated here as one bond
C                  and its different occupied levels.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         CHARACTER*16  STAGE

         LOGICAL     FAILED
         LOGICAL     LTRG,UTRG
         LOGICAL     MORE
         LOGICAL     RYD2HYB
         LOGICAL     SAVEP

         INTEGER     ATOM
         INTEGER     BONDSIZE
         INTEGER     I
         INTEGER     INCSIZE
         INTEGER     MXCHOOSE
         INTEGER     MXNHBA,MXNCEN
         INTEGER     NATOM,NHATOM
         INTEGER     NBAS
         INTEGER     NBOND
         INTEGER     NLPAIR,NEPAIR
         INTEGER     NHB,NHA
         INTEGER     N2CEN,NHCEN
         INTEGER     NSTEP,STEP

         INTEGER     ATNHB   (1:NHA   )
         INTEGER     ATHCEN  (1:NHA   )
         INTEGER     ATHOFF  (1:NHA   )
         INTEGER     ATHVAL  (1:NHA   )
         INTEGER     ATORD   (1:NHA   )
         INTEGER     BDNCEN  (1:NHB   )
         INTEGER     BDOCC   (1:NHB   )
         INTEGER     COLMAP  (1:NBAS  )
         INTEGER     IVEC    (1:NBAS+2*NHB+2*NHA)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NHOLP   (1:NHA   )
         INTEGER     NHOEP   (1:NHA   )
         INTEGER     NHORYD  (1:NHA   )
         INTEGER     NHYB    (1:NHA   )
         INTEGER     NWSTEP  (1:BONDSIZE)

         INTEGER     AT2CEN  (1:2  ,1:N2CEN)
         INTEGER     BDBAS   (1:NHA,1:NHB  )
         INTEGER     BDCEN   (1:NHA,1:NHB  )

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  W     (1:NBAS)
         DOUBLE PRECISION  WBOND (1:BONDSIZE)
         DOUBLE PRECISION  WSTAR (1:BONDSIZE)
         DOUBLE PRECISION  XVEC  (1:3*NBAS)

         DOUBLE PRECISION  BOMAT (1:NATOM ,1:NATOM )
         DOUBLE PRECISION  C     (1:NBAS  ,1:NBAS  )
         DOUBLE PRECISION  P     (1:NBAS  ,1:NBAS  )
         DOUBLE PRECISION  H     (1:MXNHBA,1:NHB   )
         DOUBLE PRECISION  PH    (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHDEP (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  PHSUB (1:NHB   ,1:NHB   )
         DOUBLE PRECISION  SAH   (1:MXNHBA,1:MXNHBA)
         DOUBLE PRECISION  XMAT  (1:NBAS  ,1:NHB   )

         DATA  ZERO   /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...reorder the NAO weights and coefficient matrix from
C                atomic to NHB/NCB/NRB order.
C
C
C         CALL    MAT__PRINT_V_INTEGER_NOZEROS
C     +
C     +                ( 1,
C     +                  ' COLMAP ',
C     +                  NBAS,
C     +                  NBAS,
C     +                  COLMAP )
C     +
C     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' Hybrid atomic ordering vector ',
     +                  1,NHA,
     +                  1,NHA,
     +                  ATORD )
     +
     +
         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,
     +                NBAS,NBAS,
     +                COLMAP,
     +                IVEC,
     +                XVEC,
     +
     +                        C )
     +
     +
C
C
C             ...initialize the # of bond/antibond weight limit steps.
C
C
         DO 10 NHCEN = 1,BONDSIZE
            NWSTEP (NHCEN) = 1
   10    CONTINUE

         INCSIZE = 0
C
C
C             ...start iterates on NHO search. Increase # of
C                bond/antibond weight limit steps for appropriate
C                bond size by + 1.
C
C
 9000    IF (INCSIZE.GT.BONDSIZE) THEN
             INCSIZE = 1
             NWSTEP (1) = NWSTEP (1) + 1
         ELSE IF (INCSIZE.GT.0) THEN
             NWSTEP (INCSIZE) = NWSTEP (INCSIZE) + 1
         END IF

         STAGE = 'start'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
         LTRG  = .TRUE.
         UTRG  = .FALSE.
         SAVEP = .TRUE.

         CALL  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +              ( NBAS,NHB,
     +                NBAS,NBAS,
     +                NBAS,NHB,
     +                NHB,
     +                NHB,NBAS,
     +                0,
     +                SAVEP,LTRG,UTRG,
     +                C,
     +                P,
     +                XVEC,
     +
     +                        XMAT )
     +
     +
         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NHB,
     +                NHB,NHB,
     +                NHB,NHB,
     +                XMAT,
     +
     +                        PH )
     +
     +
C         CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
C     +
C     +              ( NBAS,NHB,NHB,
C     +                NHB,
C     +                .TRUE.,
C     +
C     +                        XVEC,
C     +                        XMAT )
C     +
C     +
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' PH matrix ',
C     +                  NHB,NHB,
C     +                  NHB,0,
C     +                  PH )
C     +
C     +
C         CALL    MAT__PRINT_V_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Eigenvalues of PH matrix ',
C     +                  NHB,
C     +                  NHB,
C     +                  XVEC )
C     +
C     +
C
C
C             ...initialize search for pre-NHOs.
C
C
         CALL  MAT__W_EQ_ZERO_INTEGER  (NHA,NHA,NHOLP)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NHA,NHA,NHOEP)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NHA,NHA,NHORYD)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NHB,NHB,BDNCEN)
         CALL  MAT__W_EQ_ZERO_INTEGER  (NHB,NHB,BDOCC)
         CALL  MAT__C_EQ_ZERO_INTEGER  (NHA,NHB,NHA,NHB,BDBAS)
         CALL  MAT__C_EQ_ZERO_INTEGER  (NHA,NHB,NHA,NHB,BDCEN)
         CALL  MAT__C_EQ_ZERO_FLOAT    (MXNHBA,NHB,MXNHBA,NHB,H)
C
C
C             ...perform the loop over all allowed atomic hybrid
C                center combinations. Start with all hybrid atoms
C                in their optimum order using the optimum ordering
C                array. Copy initial PH matrix to the matrix that
C                will accumulate all depletions during a specific
C                bondsize run (if any). Keep going in reducing
C                the weight limits for each x-centered bond type
C                and check for pre-NHOs as long as that is still
C                possible.
C
C
         NBOND = 0
         MXNCEN = 0

         NHATOM = 0
         DO 100 I = 1,NHA
            ATOM = ATORD (I)
            NHATOM = NHATOM + 1
            ATHCEN (NHATOM) = ATOM
            NHYB (ATOM) = 0
  100    CONTINUE

         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NHB,NHB,
     +                NHB,NHB,
     +                NHB,NHB,
     +                PH,
     +
     +                        PHDEP )
     +
     +
         DO 1000 NHCEN = 1,BONDSIZE

            NSTEP = NWSTEP (NHCEN)
            WBOND (NHCEN) = WBDMIN
            WSTAR (NHCEN) = WSTMAX

            DO 1100 STEP = 1,NSTEP

               IF (NHATOM.GT.0 .AND. NHCEN.LE.NHATOM) THEN

                   CALL  NLO__CHECK_NHO_ALL_X_CENTERS
     +
     +                        ( N2CEN,NATOM,
     +                          BONDSIZE,
     +                          NHCEN,NHATOM,
     +                          MXNHBA,
     +                          AT2CEN,ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                          RYD2HYB,
     +                          IVEC (1),IVEC (NHA+1),
     +                          NHB,NHA,
     +                          WBOND (NHCEN),WSTAR (NHCEN),
     +                          SAH,
     +                          PH,PHSUB,
     +                          BOMAT,
     +                          MXCHOOSE,NCHOOSE (NHCEN),
     +                          CHOOSE (1,1,NHCEN),
     +                          IVEC (2*NHA+1),
     +                          XVEC,XMAT,
     +
     +                                   FAILED,
     +                                   NBOND,
     +                                   MXNCEN,
     +                                   NHYB,
     +                                   NHOLP,
     +                                   BDNCEN,
     +                                   BDCEN,
     +                                   BDBAS,
     +                                   BDOCC,
     +                                   W,H,
     +                                   PHDEP )
     +
     +
                   IF (FAILED) THEN
                       STAGE = 'failure'
                       CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                       ( STAGE,
     +                         BONDSIZE,
     +                         MAXOCC,
     +                         NHCEN,
     +                         NWSTEP,
     +                         WBDMIN,WSTMAX,WSTEP,
     +                         WBOND,WSTAR )
     +
     +
                       INCSIZE = INCSIZE + 1
                       GOTO 9000
                   END IF
C
C
C             ...check, if any hybrid atoms have a nonexhaustive bond
C                forming pre-NHO set. Determine their number and check
C                for more bond forming pre-NHO's. If there are more
C                allowed, update the PH matrix with all depletions
C                accumulated during the present bondsize run and
C                continue.
C
C
                   NHATOM = 0
                   DO 110 I = 1,NHA
                      ATOM = ATORD (I)
                      IF (NHYB (ATOM) .LT. ATHVAL (ATOM)) THEN
                          NHATOM = NHATOM + 1
                          ATHCEN (NHATOM) = ATOM
                      END IF
  110              CONTINUE

                   CALL  MAT__C_EQ_A_FLOAT
     +
     +                        ( NHB,NHB,
     +                          NHB,NHB,
     +                          NHB,NHB,
     +                          PHDEP,
     +
     +                               PH )
     +
     +
               END IF

               IF (STEP.NE.NSTEP) THEN
                   WBOND (NHCEN) = WBOND (NHCEN) - WSTEP
                   WSTAR (NHCEN) = WSTAR (NHCEN) + WSTEP
               END IF

 1100       CONTINUE

            IF (NHCEN.EQ.1) THEN
                NLPAIR = NBOND
            END IF

 1000    CONTINUE
C
C
C             ...try to complete the present sets of bond forming
C                pre-NHOs on each atom with eventual Empty-pair
C                pre-NHOs and possibly Rydberg pre-NHOs. If the
C                routine is unable to do so, we have to retry the
C                the bond forming pre-NHO search with lower bond
C                weight limits. If the routine is successful, the
C                complete set of orthogonalized NHOs on each atom
C                have been obtained.
C
C
 2000    STAGE = 'complete'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
         CALL  NLO__COMPLETE_NHO_SPACE
     +
     +              ( NBAS,NHATOM,
     +                MXNHBA,
     +                ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                IVEC (1),
     +                NHB,NHA,
     +                NHYB,
     +                WBDCRT,WSTAR (1),
     +                SAH,
     +                P,PH,PHSUB,
     +                W,C,
     +                IVEC (NHATOM+1),
     +                XVEC,XMAT,
     +
     +                        FAILED,
     +                        MORE,
     +                        NBOND,
     +                        NEPAIR,
     +                        NHOEP,
     +                        NHORYD,
     +                        BDNCEN,
     +                        BDCEN,
     +                        BDBAS,
     +                        BDOCC,
     +                        H )
     +
     +
         IF (FAILED) THEN
             STAGE = 'failure complete'
             CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                  ( STAGE,
     +                    BONDSIZE,
     +                    MAXOCC,
     +                    NHCEN,
     +                    NWSTEP,
     +                    WBDMIN,WSTMAX,WSTEP,
     +                    WBOND,WSTAR )
     +
     +
             INCSIZE = INCSIZE + 1
             GOTO 9000
         END IF

         IF (MORE) THEN

             STAGE = 'more'

             CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                  ( STAGE,
     +                    BONDSIZE,
     +                    MAXOCC,
     +                    NHCEN,
     +                    NWSTEP,
     +                    WBDMIN,WSTMAX,WSTEP,
     +                    WBOND,WSTAR )
     +
     +
             DO 2200 NHCEN = 1,BONDSIZE

                IF (NHCEN.EQ.1) THEN
                    NLPAIR = NLPAIR - NBOND
                END IF

                WBOND (NHCEN) = WBOND (NHCEN) - WSTEP
                WSTAR (NHCEN) = WSTAR (NHCEN) + WSTEP

                IF (NHATOM.GT.0 .AND. NHCEN.LE.NHATOM) THEN

                    CALL  NLO__CHECK_NHO_ALL_X_CENTERS
     +
     +                         ( N2CEN,NATOM,
     +                           BONDSIZE,
     +                           NHCEN,NHATOM,
     +                           MXNHBA,
     +                           AT2CEN,ATHCEN,ATNHB,ATHVAL,ATHOFF,
     +                           RYD2HYB,
     +                           IVEC (1),IVEC (NHA+1),
     +                           NHB,NHA,
     +                           WBOND (NHCEN),WSTAR (NHCEN),
     +                           SAH,
     +                           PH,PHSUB,
     +                           BOMAT,
     +                           MXCHOOSE,NCHOOSE (NHCEN),
     +                           CHOOSE (1,1,NHCEN),
     +                           IVEC (2*NHA+1),
     +                           XVEC,XMAT,
     +
     +                                    FAILED,
     +                                    NBOND,
     +                                    MXNCEN,
     +                                    NHYB,
     +                                    NHOLP,
     +                                    BDNCEN,
     +                                    BDCEN,
     +                                    BDBAS,
     +                                    BDOCC,
     +                                    W,H,
     +                                    PHDEP )
     +
     +
                    IF (FAILED) THEN
                        STAGE = 'failure'
                        CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +                        ( STAGE,
     +                          BONDSIZE,
     +                          MAXOCC,
     +                          NHCEN,
     +                          NWSTEP,
     +                          WBDMIN,WSTMAX,WSTEP,
     +                          WBOND,WSTAR )
     +
     +
                        INCSIZE = INCSIZE + 1
                        GOTO 9000
                    END IF

                    NHATOM = 0
                    DO 200 I = 1,NHA
                       ATOM = ATORD (I)
                       IF (NHYB (ATOM) .LT. ATHVAL (ATOM)) THEN
                           NHATOM = NHATOM + 1
                           ATHCEN (NHATOM) = ATOM
                       END IF
  200               CONTINUE

                    CALL  MAT__C_EQ_A_FLOAT
     +
     +                         ( NHB,NHB,
     +                           NHB,NHB,
     +                           NHB,NHB,
     +                           PHDEP,
     +
     +                                PH )
     +
     +
                END IF

                IF (NHCEN.EQ.1) THEN
                    NLPAIR = NLPAIR + NBOND
                END IF

 2200       CONTINUE
            GOTO 2000
         END IF

         STAGE = 'success'

         CALL  NLO__PRINT_NHO_SEARCH_INFO
     +
     +              ( STAGE,
     +                BONDSIZE,
     +                MAXOCC,
     +                NHCEN,
     +                NWSTEP,
     +                WBDMIN,WSTMAX,WSTEP,
     +                WBOND,WSTAR )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
