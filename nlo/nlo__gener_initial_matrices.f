         SUBROUTINE  NLO__GENER_INITIAL_MATRICES
     +
     +                    ( NBAS,NATOM,
     +                      BASBEG,BASEND,
     +                      S,
     +                      XMAT,
     +
     +                            BOMAT,
     +                            D )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_INITIAL_MATRICES
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine generates the atomic bond order matrix
C                and the occupation matrix from the density matrix
C                (i.e. the expansion coefficients of the density rho(r)
C                in terms of basis function products) and the overlap
C                matrix. The occupation matrix is defined as a matrix
C                representation of the density rho(r) in terms of the
C                basis functions.
C
C                Procedure:
C
C                The definition of the density matrix with elements
C                D (i,j) in terms of the basis functions is:
C
C                       rho (r)  =  sum  D (i,j) * chi (i) * chi (j)
C                                    ij
C
C                The occupation matrix element P (k,l) is thus:
C
C                       P (k,l)  =  < chi (k) | rho (r) | chi (l) >
C
C                                =  sum  S (k,i) * D (i,j) * S (j,l)
C                                    ij
C
C                where S is the overlap matrix of the basis functions.
C                The atomic bond order matrix is defined as summations
C                of elements of D*S between atomic subblocks a and b:
C
C                    BOMAT (a,b) =  sum     sum    DS (i,j) * DS (j,i)
C                                 i in a  j in b
C
C                Note that the product D*S is needed for evaluation
C                of both matrices, hence the first step is the
C                generation of D*S.
C
C                Since the density matrix D is not needed later on,
C                the array D will be overwritten with the occupation
C                matrix.
C
C
C                  Input:
C
C                    NBAS         =  total # of basis functions
C                    NATOM        =  total # of atomic centers
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    S            =  overlap matrix
C                    XMAT         =  flp scratch matrix
C                    D            =  density matrix
C
C
C                  Output:
C
C                    BOMAT        =  atomic bond order matrix
C                    D            =  occupation matrix in AO basis
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

         INTEGER     A,B,I,J
         INTEGER     IBEG,IEND,JBEG,JEND
         INTEGER     NATOM
         INTEGER     NBAS

         INTEGER     BASBEG  (1:NATOM)
         INTEGER     BASEND  (1:NATOM)

         DOUBLE PRECISION  XVEC  (1:NBAS)

         DOUBLE PRECISION  BOMAT (1:NATOM,1:NATOM)
         DOUBLE PRECISION  D     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  S     (1:NBAS ,1:NBAS )
         DOUBLE PRECISION  XMAT  (1:NBAS ,1:NBAS )
C
C
C------------------------------------------------------------------------
C
C
C             ...form D*S first.
C
C
         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,NBAS,
     +                D,S,
     +
     +                         XMAT )
     +
     +
C
C
C             ...evaluate the atomic bond order matrix.
C
C
         CALL  MAT__C_EQ_ZERO_FLOAT
     +
     +              ( NATOM,NATOM,
     +                NATOM,NATOM,
     +
     +                         BOMAT )
     +
     +
         DO 100 B = 1,NATOM
            JBEG = BASBEG (B)
            JEND = BASEND (B)
            DO 110 A = 1,NATOM
               IBEG = BASBEG (A)
               IEND = BASEND (A)

               DO 120 I = IBEG,IEND
               DO 120 J = JBEG,JEND
                  BOMAT (A,B) = BOMAT (A,B) + XMAT (I,J) * XMAT (J,I)
  120          CONTINUE

  110       CONTINUE
  100    CONTINUE

         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +                ( 1,
     +                  ' Atomic bond order matrix ',
     +                  NATOM,NATOM,
     +                  NATOM,NATOM,
     +                  BOMAT )
     +
     +
C
C
C             ...form occupation matrix P = S * DS overwriting D.
C
C
         CALL  MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +              ( NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,
     +                NBAS,NBAS,NBAS,
     +                S,XMAT,
     +
     +                         D )
     +
     +
         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +                ( 1,
     +                  'printing SPS with P the HF',
     +                  NBAS,NBAS,
     +                  NBAS,NBAS,
     +                  D )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
