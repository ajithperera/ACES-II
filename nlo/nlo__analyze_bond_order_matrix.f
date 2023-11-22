         SUBROUTINE  NLO__ANALYZE_BOND_ORDER_MATRIX
     +
     +                    ( NATOM,
     +                      MX2CEN,
     +                      NO2CEN,
     +                      INDEX,
     +                      X,
     +
     +                            N2CEN,
     +                            AT2CEN,
     +                            ATORD,
     +                            BOMAT )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ANALYZE_BOND_ORDER_MATRIX
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine analyzes the atomic bond order matrix
C                and extracts those pairs of atomic indices which are
C                most likely to give 2 center bonds. The routine
C                searches through the lower triangle of the bond order
C                matrix and orders the atomic pairs such that those
C                having the largest bond orders come first. This will
C                be the order in which the search for 2-center bonds
C                will be performed. Also at this stage the bond order
C                matrix will be transformed to a more simple form:
C                If a 2-center bond between atoms A and B has been
C                found, reset BOMAT (A,B) = 1.0. If no 2-center bond
C                is present then BOMAT (A,B) = zero. This makes
C                searches of the bond order matrix at a later stage
C                much easier and avoids passing the 2-center bond
C                formation limit NO2CEN to subsequent routines.
C
C
C                  Input:
C
C                    NATOM        =  total # of atomic centers
C                    MX2CEN       =  total # of atomic center pairs
C                    NO2CEN       =  2-center bond formation criterion
C                                    for analysis of the atomic bond
C                                    order matrix.
C                    INDEX        =  will hold reordering indices
C                    X            =  will hold lower triangle elements
C                                    of atomic bond order matrix
C                    BOMAT        =  initial atomic bond order matrix
C                                    containing the true bond orders
C
C
C                  Output:
C
C                    N2CEN        =  # of atomic center pairs that
C                                    will be considered for 2-center
C                                    bond construction.
C                    AT2CEN (1,N) =  1st atomic center label of N-th
C                                    pair.
C                    AT2CEN (2,N) =  2nd atomic center label of N-th
C                                    pair. Always the 2nd label is >
C                                    than the 1st label.
C                    ATORD        =  ordering of atomic indices to
C                                    be used for NBO construction.
C                    BOMAT        =  simplified atomic bond order
C                                    matrix with ones and zeros.
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

         LOGICAL     ABSOLUT
         LOGICAL     INCRESE

         INTEGER     I,J,M,N
         INTEGER     MX2CEN
         INTEGER     N2CEN
         INTEGER     NATOM
         INTEGER     NSET

         INTEGER     ATORD (1:NATOM)
         INTEGER     INDEX (1:MX2CEN)

         INTEGER     AT2CEN  (1:2,1:MX2CEN)

         DOUBLE PRECISION  BOVAL
         DOUBLE PRECISION  DIFF
         DOUBLE PRECISION  NO2CEN
         DOUBLE PRECISION  SYMTRY
         DOUBLE PRECISION  ZERO,ONE,TWO,EIGHT

         DOUBLE PRECISION  X (1:MX2CEN)

         DOUBLE PRECISION  BOMAT (1:NATOM,1:NATOM)

         DATA  ZERO   /0.D0/
         DATA  ONE    /1.D0/
         DATA  TWO    /2.D0/
         DATA  EIGHT  /8.D0/
         DATA  SYMTRY /1.D-5/
C
C
C------------------------------------------------------------------------
C
C
C             ...copy lower triangle of bond order matrix rowwise
C                into array X, starting at first row. The diagonal
C                elements of the bond order matrix are not considered
C                and a zero will be placed in the corresponding place
C                in X.
C
C
         N = 0
         DO 100 J = 1,NATOM
            DO 110 I = 1,J-1
               N = N + 1
               X (N) = BOMAT (I,J)
  110       CONTINUE
            N = N + 1
            X (N) = ZERO
  100    CONTINUE
C
C
C             ...order array X values and get corresponding index
C                vector.
C
C
         ABSOLUT = .TRUE.
         INCRESE = .FALSE.

         CALL  NLO__SORT_FLP_VECTOR_ELEMENTS
     +
     +              ( N,N,
     +                1,N,
     +                ABSOLUT,INCRESE,
     +                0,
     +                X,
     +
     +                         INDEX )
     +
     +
C
C
C             ...determine the # of relevant 2-center pairs.
C
C
         N2CEN = 0
         DO 200 I = 1,N
            BOVAL = DABS (X (INDEX (I)))
            IF (BOVAL.GT.NO2CEN) THEN
                N2CEN = N2CEN + 1
            END IF
  200    CONTINUE
C
C
C             ...reconstruct the atomic indices for all the 2-center
C                pairs determined.
C
C
         DO 300 N = 1,N2CEN
            M = INDEX (N)
            J = INT ( (ONE + DSQRT (EIGHT * DFLOAT (M))) / TWO )
            I = M - J*(J-1)/2
            AT2CEN (1,N) = I
            AT2CEN (2,N) = J
  300    CONTINUE
C
C
C             ...group the atomic indices into sets of equal bond
C                orders (within a certain tolerance) and reorder
C                the atomic indices within each set to improve
C                their connectivity. This is necessary to avoid
C                random order of symmetry equivalent atomic indices
C                leading in some cases to disconnected sets which
C                during the 2 center bond search result in polyradical
C                structures. As an example consider D6h benzene:
C
C                                      1
C                                     / \
C                                    /   \
C                                   6     2
C                                   |     |
C                                   |     |
C                                   5     3
C                                    \   /
C                                     \ /
C                                      4
C
C                Due to its symmetry, all bond orders:
C
C                             1  4  2  5  3  1
C                             2  5  3  6  4  6
C
C                will be large and of equal magnitude. But a search
C                for the 2 center bonds employing that particular
C                order will place pi-bonds between 1 and 2 and between
C                4 and 5, leaving atoms 3 and 6 disconnected from the
C                rest with subsequent 'empty' pairs on 3 and 6 with
C                occupation numbers near 1/2 the pi-bond occupation
C                number (i.e. a diradical).
C
C
         N = 1

 4000    IF (N.LE.N2CEN) THEN
             NSET = 1
             I = AT2CEN (1,N)
             J = AT2CEN (2,N)
             BOVAL = DABS (BOMAT (I,J))
             DO 400 M = N+1,N2CEN
                I = AT2CEN (1,M)
                J = AT2CEN (2,M)
                DIFF = DABS (BOMAT (I,J)) - BOVAL
                IF (DABS (DIFF).LT.SYMTRY) THEN
                    NSET = NSET + 1
                END IF
  400        CONTINUE

             CALL  NLO__IMPROVE_2CEN_CONNECTIONS
     +
     +                  ( NSET,
     +                    INDEX,
     +
     +                           AT2CEN (1,N) )
     +
     +
             N = N + NSET
             GOTO 4000
         END IF
C
C
C             ...determine new atomic order from the ordered list
C                of atomic index pairs.
C
C
         CALL  NLO__OPTIMUM_ATOM_INDEX_ORDER
     +
     +              ( NATOM,N2CEN,
     +                INDEX,
     +                AT2CEN,
     +
     +                        ATORD )
     +
     +
         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' Atomic ordering vector ',
     +                  1,NATOM,
     +                  1,NATOM,
     +                  ATORD )
     +
     +
C
C
C             ...replace bond order matrix with simplified version.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NATOM,NATOM,
     +                NATOM,NATOM,
     +
     +                        BOMAT )
     +
     +
         DO 500 N = 1,N2CEN
            I = AT2CEN (1,N)
            J = AT2CEN (2,N)
            BOMAT (I,J) = ONE
            BOMAT (J,I) = ONE
  500    CONTINUE

         CALL    MAT__PRINT_A_INTEGER_NOZEROS
     +
     +                ( 1,
     +                  ' Atomic 2-center labels for NBO search ',
     +                  2,MX2CEN,
     +                  2,N2CEN,
     +                  AT2CEN )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
