         SUBROUTINE  NLO__GENER_NBO_ORBITALS
     +
     +                    ( NBAS,NBOND,
     +                      MXNHBA,MXNCEN,
     +                      BONDSIZE,
     +                      NHB,NCB,NHA,
     +                      ATNHB,ATHOFF,
     +                      BDNCEN,BDCEN,BDBAS,BDOCC,
     +                      COLMAP,
     +                      WBOND,WSTAR,
     +                      WRYD,WOCC,
     +                      P,H,
     +                      PH,
     +                      IVEC,
     +                      XVEC,XMAT,
     +
     +                              NLB,NBB,NEB,NAB,NYB,
     +                              NBOBD,
     +                              B,
     +                              W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NBO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the sets of natural bond
C                orbitals (NBOs) by combining the final set of NHOs
C                in the appropriate way.
C
C
C                  Input:
C
C                    NBAS         =  total # of AOs in AO basis
C                    NBOND        =  total # of NHO bonds found
C                    MXNHBA       =  maximum # of Hybrid NAOs per atom.
C                    MXNCEN       =  maximum # of atomic hybrid centers
C                                    per bond.
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NHB          =  total # of Hybrid NAOs
C                    NCB          =  total # of Core NAOs
C                    NHA          =  total # of hybrid atoms
C                    ATNHB (A)    =  # of Hybrid NAOs on hybrid atom A.
C                    ATHOFF (A)   =  index offset for Hybrid NAOs for
C                                    hybrid atom A. This index is equal
C                                    to the total number of Hybrid
C                                    NAOs on all hybrid atoms preceeding
C                                    hybrid atom A.
C                    BDNCEN (I)   =  # of atomic centers for I-th bond.
C                    BDCEN (I,J)  =  I-th atomic center index for
C                                    J-th bond.
C                    BDBAS (I,J)  =  I-th global basis (NHO) index for
C                                    J-th bond.
C                    BDOCC (I)    =  # of occupied levels for I-th bond.
C                    COLMAP (I)   =  will contain the column map in the
C                                    active sense for the pure NBO
C                                    ordering NHB -> NLB/NBB/NEB/NAB/NYB
C                                    and the subsequent mixed NBO/NAO
C                                    orderng NHB/NCB -> NCB/NHB to bring
C                                    the Core NAOs up front. COLMAP (I)
C                                    contains the position index of the
C                                    I-th NBO in the final respective
C                                    order.
C                    WBOND (x)    =  lowest weight criterion for each
C                                    x-center bond.
C                    WSTAR (x)    =  highest weight criterion for each
C                                    x-center antibond.
C                    WRYD         =  highest weight below which a one-
C                                    center bond is classified as an
C                                    atomic Rydberg NBO.
C                    WOCC         =  the weight limit to decide when
C                                    a NBO is considered occupied or
C                                    virtual.
C                    P            =  full NBAS x NBAS occupation matrix.
C                    H (I,J)      =  MXNHBA x NHB matrix containing the
C                                    atomic NHOs. I is the local atomic
C                                    index labeling the atomic Hybrid
C                                    NAOs from which the NHOs are
C                                    constructed. J is the global NHO
C                                    index running over all NHB NHOs,
C                                    with all NHOs belonging to a
C                                    specific atomic center being
C                                    grouped together.
C                    PH           =  will contain the NHB x NHB valence
C                                    occupation matrix in NAO basis
C                    IVEC         =  int scratch array of vector type.
C                    XVEC         =  flp scratch array of vector type.
C                    XMAT         =  flp scratch array of matrix type.
C                    W            =  weight vector in NHB/NCB order. The
C                                    NCB part contains the original
C                                    Core NAO weights. The NHB part was
C                                    used during NHO generation and will
C                                    be replaced by the NBO Lone-pair,
C                                    Bond, Empty-pair, Antibond and
C                                    Rydberg weights.
C                    C            =  NBAS x (NHB+NCB) part of the NAO
C                                    coefficient matrix in AO basis
C                                    with columns in NHB/NCB order.
C
C
C                  Output:
C
C                    NxB          =  total # of Lone-pair, Bond,
C                                    Empty-pair, Antibond and Rydberg
C                                    NBOs found (x=L,B,E,A,Y).
C                    NBOBD (I)    =  contains the NHO bond index number
C                                    number for the I-th NBO in the NHB
C                                    part. The NBOBD elements are in
C                                    NLB/NBB/NEB/NAB/NYB order.
C                    B (I,J)      =  MXNCEN x NHB matrix containing
C                                    the J-th NBO expansion coefficients
C                                    in terms of the I-th atomic NHOs
C                                    forming the J-th NBO. The NBO
C                                    column index is in NLB/NBB/NEB/
C                                    NAB/NYB order.
C                    W            =  NBO weight vector with elements in
C                                    NCB/NLB/NBB/NEB/NAB/NYB order.
C                    C            =  NBO coefficient matrix in AO basis
C                                    with columns in NCB/NLB/NBB/NEB/
C                                    NAB/NYB order.
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

         LOGICAL     LTRG,UTRG
         LOGICAL     SAVEP

         INTEGER     BOND
         INTEGER     BONDSIZE
         INTEGER     I
         INTEGER     INDEX
         INTEGER     MXNHBA,MXNCEN
         INTEGER     NBAS
         INTEGER     NBO,NBOOLD
         INTEGER     NBOND
         INTEGER     NHB,NCB,NHA
         INTEGER     NHCEN
         INTEGER     NLB,NBB,NEB,NAB,NYB
         INTEGER     NMOVE
         INTEGER     NOCC
         INTEGER     OFFLP,OFFBD,OFFEP,OFFAB,OFFRY

         INTEGER     ATNHB   (1:NHA    )
         INTEGER     ATHOFF  (1:NHA    )
         INTEGER     BDNCEN  (1:NHB    )
         INTEGER     BDOCC   (1:NHB    )
         INTEGER     COLMAP  (1:NCB+NHB)
         INTEGER     IVEC    (1:NCB+NHB)
         INTEGER     NBOBD   (1:NHB    )

         INTEGER     BDBAS   (1:NHA,1:NHB)
         INTEGER     BDCEN   (1:NHA,1:NHB)

         DOUBLE PRECISION  WRYD,WOCC

         DOUBLE PRECISION  W     (1:NCB+NHB)
         DOUBLE PRECISION  WBOND (1:BONDSIZE)
         DOUBLE PRECISION  WSTAR (1:BONDSIZE)
         DOUBLE PRECISION  XVEC  (1:NBAS   )

         DOUBLE PRECISION  B     (1:MXNCEN,1:NHB    )
         DOUBLE PRECISION  C     (1:NBAS  ,1:NCB+NHB)
         DOUBLE PRECISION  P     (1:NBAS  ,1:NBAS   )
         DOUBLE PRECISION  H     (1:MXNHBA,1:NHB    )
         DOUBLE PRECISION  PH    (1:NHB   ,1:NHB    )
         DOUBLE PRECISION  XMAT  (1:NBAS  ,1:NHB    )
C
C
C------------------------------------------------------------------------
C
C
C             ...form the full NHB x NHB occupation matrix PH.
C
C
         LTRG  = .TRUE.
         UTRG  = .TRUE.
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
C
C
C             ...loop over all NHO bonds found and accumulate the NBO
C                weights, the coefficient matrix, the NHO bond index
C                numbers and the local column NBO ordering map
C                NHB -> NLB/NBB/NEB/NAB/NYB with added offsets for
C                temporary protection. The difference between the
C                offsets is set to NHB, since this ensures no overlap
C                of indices between different types of NBOs.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( MXNCEN,NHB,
     +                MXNCEN,NHB,
     +
     +                        B )
     +
     +
         NBO = 0
         NLB = 0
         NBB = 0
         NEB = 0
         NAB = 0
         NYB = 0
         NBOOLD = 0

         OFFLP = 0
         OFFBD = OFFLP + NHB
         OFFEP = OFFBD + NHB
         OFFAB = OFFEP + NHB
         OFFRY = OFFAB + NHB

         DO 100 BOND = 1,NBOND

            NOCC = BDOCC (BOND)
            NHCEN = BDNCEN (BOND)

            CALL  NLO__FORM_NBO_X_CENTERS
     +
     +              ( NBAS,NHCEN,NOCC,
     +                MXNHBA,MXNCEN,
     +                BONDSIZE,
     +                ATNHB,ATHOFF,
     +                NHB,NHA,
     +                OFFLP,OFFBD,OFFEP,OFFAB,OFFRY,
     +                BDCEN (1,BOND),
     +                BDBAS (1,BOND),
     +                WBOND (NHCEN),
     +                WSTAR (NHCEN),
     +                WRYD,WOCC,
     +                C,
     +                PH,
     +                H,
     +                IVEC,
     +
     +                        NBO,
     +                        NLB,NBB,NEB,NAB,NYB,
     +                        COLMAP,
     +                        B (1,NBOOLD+1),
     +                        W (NBOOLD+1),
     +                        XMAT (1,NBO+1) )
     +
     +
            DO 110 I = NBOOLD+1,NBO
               NBOBD (I) = BOND
  110       CONTINUE

            NBOOLD = NBO
  100    CONTINUE
C
C
C             ...check, if the # of Lone-pairs, Bonds, Empty-pairs,
C                Antibonds and Rydbergs add up ok.
C
C
         IF ((NLB+NBB+NEB+NAB+NYB).NE.NHB) THEN
             WRITE (*,*) ' Problems in finding correct # of NBOs! '
             WRITE (*,*) ' NLB,NBB,NEB,NAB,NYB,NHB = ',
     +                     NLB,NBB,NEB,NAB,NYB,NHB
             WRITE (*,*) ' nlo__gener_nbo_orbitals '
             WRITE (1,*) ' Problems in finding correct # of NBOs! '
             WRITE (1,*) ' NLB,NBB,NEB,NAB,NYB,NHB = ',
     +                     NLB,NBB,NEB,NAB,NYB,NHB
             WRITE (1,*) ' nlo__gener_nbo_orbitals '
             STOP
         END IF
C
C
C             ...remove the protective offsets on the local column
C                NBO ordering NHB -> NLB/NBB/NEB/NAB/NYB and reorder
C                the NBO -> NHO expansion coefficient matrix as well
C                as the vector containing the NHO bond index numbers.
C                Use the obsolete PH matrix for scratch.
C
C
         DO 200 I = 1,NHB
            INDEX = COLMAP (I)
            IF (INDEX.GT.OFFRY) THEN
                INDEX = INDEX - OFFRY + NLB + NBB + NEB + NAB
            ELSE IF (INDEX.GT.OFFAB) THEN
                INDEX = INDEX - OFFAB + NLB + NBB + NEB
            ELSE IF (INDEX.GT.OFFEP) THEN
                INDEX = INDEX - OFFEP + NLB + NBB
            ELSE IF (INDEX.GT.OFFBD) THEN
                INDEX = INDEX - OFFBD + NLB
            ELSE IF (INDEX.GT.OFFLP) THEN
                INDEX = INDEX - OFFLP
            END IF
            COLMAP (I) = INDEX
  200    CONTINUE

         CALL  MAT__REORDER_VECTOR_INTEGER
     +
     +              ( NHB,
     +                NHB,
     +                NHB,
     +                NHB,
     +                COLMAP,
     +                IVEC,
     +
     +                        NBOBD )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( MXNCEN,NHB,
     +                NHB,
     +                NHB,
     +                MXNCEN,
     +                MXNCEN,NHB,
     +                COLMAP,
     +                IVEC,
     +                PH,
     +
     +                        B )
     +
     +
C
C
C             ...determine global mixed column NBO/NAO ordering NHB/NCB
C                -> NCB/NHB. Copy the NBO coefficient accumulation
C                matrix to final NBO coefficient matrix and reorder
C                the latter from NHB/NCB order to NCB/NHB order, such
C                that the Core NAOs come first. Reorder also the
C                NHB/NCB part of the weight vector.
C
C
         DO 300 I = 1,NHB
            COLMAP (I) = COLMAP (I) + NCB
  300    CONTINUE

         DO 310 I = 1,NCB
            COLMAP (NHB+I) = I
  310    CONTINUE

         CALL  MAT__C_EQ_A_FLOAT
     +
     +              ( NBAS,NHB,
     +                NBAS,NBAS,
     +                NBAS,NHB,
     +                XMAT,
     +
     +                        C )
     +
     +
         NMOVE = NCB + NHB

         CALL  MAT__REORDER_VECTOR_FLOAT
     +
     +              ( NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                COLMAP,
     +                XVEC,
     +
     +                        W )
     +
     +
         CALL  MAT__REORDER_MATRIX_COLUMNS
     +
     +              ( NBAS,NMOVE,
     +                NMOVE,
     +                NMOVE,
     +                NBAS,
     +                NBAS,NMOVE,
     +                COLMAP,
     +                IVEC,
     +                XMAT,
     +
     +                        C )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
