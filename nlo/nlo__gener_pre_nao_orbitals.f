         SUBROUTINE  NLO__GENER_PRE_NAO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      MXSHELL,MXNAL,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      LSIZE,
     +                      P,S,
     +                      SAL,CAL,WAL,
     +                      WRYDAT,
     +                      MJUMP,
     +
     +                              NMB,NRB,
     +                              COLMAP,
     +                              NRYDAL,
     +                              W,C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_PRE_NAO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine forms all the pre-NAO's by diagonalizing
C                the systems of equations P(AL)C(AL) = S(AL)C(AL)W(AL)
C                for all (AL) spaces. Analyze the symmetry average
C                weights W(AL) and decompose the pre-NAO space into
C                the minimal and Rydberg pre-NAO spaces.
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A.
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    LSIZE (I)    =  I-th l-shell size
C                    P            =  full NBAS x NBAS occupation matrix.
C                    S            =  full NBAS x NBAS overlap matrix.
C                    SAL,CAL,WAL  =  submatrices S(AL) and C(AL) for
C                                    overlap and density/pre-NAO coeffs
C                                    and subvector W(AL) for weights.
C                    WRYDAT (A)   =  the pre-NAO weight threshold below
C                                    which a pre-NAO will be considered
C                                    of Rydberg type for atom A.
C                    MJUMP        =  is .true., if the m values in the
C                                    m-space are ordered such that the
C                                    same m values are separated. This
C                                    keyword is necessary because some
C                                    AO basis functions are m-ordered
C                                    differently within each l-shell.
C                                    It invokes different types of
C                                    m-averaging algorithms.
C
C
C                  Output:
C
C                    NMB,NRB      =  # of minimal and Rydberg pre-NAO's
C                                    found.
C                    COLMAP (I)   =  column map in the active sense,
C                                    such that COLMAP (I) contains the
C                                    position index of I-th pre-NAO,
C                                    based on atomic clustering, in the
C                                    final NMB and NRB bundled pre-NAO
C                                    coefficient array C.
C                    NRYDAL (I,A) =  # of Rydberg pre-NAO's found for
C                                    the I-th atomic l-shell for atom A.
C                                    Needed later on for updating the
C                                    Rydberg pre-NAO's when generating
C                                    the final NAO's.
C                    W            =  pre-NAO weight vector.
C                    C            =  pre-NAO coefficient matrix in AO
C                                    basis.
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

         LOGICAL     LTRG
         LOGICAL     MJUMP
         LOGICAL     REVERS
         LOGICAL     ZEROC

         INTEGER     ATOM
         INTEGER     BASNR
         INTEGER     INDEX
         INTEGER     J,L,M,N
         INTEGER     LDIM,LTOT,LTYPE
         INTEGER     MXSHELL,MXNAL
         INTEGER     NBAS,NATOM
         INTEGER     NAL
         INTEGER     NLTYPE
         INTEGER     NMB,NRB,NRBOLD

         INTEGER     COLMAP  (1:NBAS)
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NSHELLS (1:NATOM)

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     NRYDAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)

         DOUBLE PRECISION  WEIGHT,WTHRSH
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  W      (1:NBAS)
         DOUBLE PRECISION  WAL    (1:MXNAL)
         DOUBLE PRECISION  WRYDAT (1:NATOM)

         DOUBLE PRECISION  C   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  S   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  P   (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  SAL (1:MXNAL,1:MXNAL)
         DOUBLE PRECISION  CAL (1:MXNAL,1:MXNAL)

         PARAMETER    (ZERO    =  0.D0 )
         PARAMETER    (ONE     =  1.D0 )
         PARAMETER    (REVERS  = .TRUE.)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loop over all atoms, inner loop over all l-shell
C                types for each atom.
C
C
         LTRG  = .TRUE.
         ZEROC = .FALSE.

         NMB = 0
         NRB = 0
         BASNR = 0

         DO 1000 ATOM = 1,NATOM
            WTHRSH = WRYDAT (ATOM)
            NLTYPE = NSHELLS (ATOM)
            DO 2000 N = 1,NLTYPE
C
C
C             ...weighted average in m-symmetry for present shell.
C                Extracts lower triangular parts of relevant diagonal
C                subblocks from the occupation matrix P and the overlap
C                matrix S and places them in tiny matrices CAL and SAL,
C                respectively. Subblock occupation matrix CAL will be
C                used for storing eigenvectors later on and will hence
C                be destroyed.
C
C
               NAL = NBASAL (N,ATOM)
               LTYPE = SHELLS (N,ATOM)
               LDIM = LSIZE (LTYPE)
               LTOT = NAL * LDIM

               CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                    ( NBAS,LTOT,
     +                      MXNAL,MXNAL,
     +                      LTOT,LDIM,
     +                      BASNR,
     +                      MJUMP,LTRG,
     +                      P (1,BASNR+1),
     +
     +                              CAL )
     +
     +
               CALL  NLO__FORM_M_AVERAGED_MATRIX
     +
     +                    ( NBAS,LTOT,
     +                      MXNAL,MXNAL,
     +                      LTOT,LDIM,
     +                      BASNR,
     +                      MJUMP,LTRG,
     +                      S (1,BASNR+1),
     +
     +                              SAL )
     +
     +
C
C
C             ...solve the generalized eigenvalue problem for
C                the m-symmetry averaged P(AL) and S(AL) matrices.
C
C
               CALL  MAT__GEN_EIGSYS_REAL_SYMMETRIC
     +
     +                    ( MXNAL,MXNAL,
     +                      MXNAL,MXNAL,
     +                      MXNAL,
     +                      NAL,
     +                      REVERS,
     +
     +                              WAL,
     +                              SAL,
     +                              CAL )
     +
     +
C
C
C             ...zero appropriate section of the pre-NAO coefficient
C                matrix and form the pre-NAO coefficients for that
C                section.
C
C
               CALL  MAT__W_EQ_ZERO_FLOAT
     +
     +                    ( LTOT*NBAS,
     +                      LTOT*NBAS,
     +
     +                              C (1,BASNR+1) )
     +
     +
               CALL  NLO__FORM_M_EXPANDED_MATRIX
     +
     +                    ( MXNAL,MXNAL,
     +                      NBAS,LTOT,
     +                      NAL,LTOT,
     +                      BASNR,
     +                      MJUMP,.FALSE.,ZEROC,
     +                      CAL,
     +
     +                              C (1,BASNR+1) )
     +
C
C
C             ...determine decomposition into minimal and Rydberg
C                spaces and generate the pre-NAO coefficients in
C                the AO basis together with the corresponding average
C                weights.
C
C                The decomposition into minimal and Rydberg spaces
C                for each atom A is governed by the critical value
C                WRYDAT (A), which has the following meaning:
C
C                          weight  >  WRYDAT (A)    minimal
C                          weight =<  WRYDAT (A)    Rydberg
C
C                As we move along, determining each atomic section of
C                the pre-NAO coefficients, we also need to keep track
C                of which column indices will belong to the NMB and
C                NRB spaces. A column map COLMAP is established, which
C                will contain the mapping of the columns from the
C                original atomic structure to the new NMB and NRB
C                clustered structure in the active sense. Pictorially,
C                this means a columnwise restructuring of the pre-NAO
C                coefficient matrix:
C
C
C                 at 1    at 2   ...             NMB          NRB
C               --------------------         ------------------------
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              | N | N | N | N | ...        | a | a |... | a | a |...
C              | M | R | M | R | ...    ->  | t | t |... | t | t |...
C              | B | B | B | B | ...        | 1 | 2 |... | 1 | 2 |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C              |   |   |   |   | ...        |   |   |... |   |   |...
C
C
C
C                The active column map has the following info:
C
C                  COLMAP (atomic based index) = NMB/NRB based index
C
C                Note, that while it is easy to determine the NMB
C                index, the NRB index has to be determined after the
C                we know the size of the NMB space. Hence the initial
C                COLMAP values for the NRB space must be protected
C                from identification at a later stage, which is easy
C                to do by adding the constant NBAS, since the maximum
C                value that NMB can take is NBAS.
C
C                The whole idea behind creating the column map is the
C                fact that the NMB and NRB spaces will be manipulated
C                separately, hence for efficiency reasons in handling
C                matrix operations on these spaces we will temporarily
C                bundle each kind together in later routines.
C
C
               J = BASNR
               NRBOLD = NRB
C hughes
C               write(*,*) 'atom = ',atom
C               write(*,*) 'nal  = ',nal
               DO 200 L = 1,NAL

                  WEIGHT = WAL (L)
C                  write(*,*) 'weight = ',weight

                  IF (WEIGHT.GT.WTHRSH) THEN
                      DO 210 M = 1,LDIM
                         J = J + 1
                         NMB = NMB + 1
                         COLMAP (J) = NMB
                         W (J) = WEIGHT
  210                 CONTINUE
                  ELSE
                      DO 220 M = 1,LDIM
                         J = J + 1
                         NRB = NRB + 1
                         COLMAP (J) = NRB + NBAS
                         W (J) = WEIGHT
  220                 CONTINUE
                  END IF

  200          CONTINUE

               NRYDAL (N,ATOM) = NRB - NRBOLD
C
C
C             ...next shell type / atomic site.
C
C
               BASNR = BASNR + LTOT

 2000       CONTINUE
 1000    CONTINUE
C
C
C             ...the complete pre-NAO coefficient matrix is ready in
C                atomic order. Update the NRB section of the column
C                map.
C
C
         DO 300 J = 1,NBAS
            INDEX = COLMAP (J)
            IF (INDEX.GT.NBAS) THEN
                INDEX = INDEX - NBAS + NMB
                COLMAP (J) = INDEX
            END IF
  300    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
