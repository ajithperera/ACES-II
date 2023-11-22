         SUBROUTINE  NLO__ANALYZE_NBO_ORBITALS
     +
     +                    ( NBAS,NATOM,
     +                      ATIDX,
     +                      BASBEG,BASEND,
     +                      C,
     +                      SHALF,
     +
     +                              LOCAL )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__ANALYZE_NBO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine analyzes the NBO orbitals obtained.
C                An important quantity that is determined here is
C                the atomic localization, which measures the degree
C                of atomic orbital involvment in each NBO. We can
C                extract the content of the correponding AOs in the
C                NBO expansion via the following analysis:
C
C                In order to get rid of the AO nonorthogonality
C                problem, we express each NBO in terms of the
C                orthonormal symmetrized AO basis (SAO), which
C                resembles most closely the original AO basis. The
C                new NBO expansion coefficients c' in terms of SAOs
C                are then given by:
C
C                                 c' = S**(1/2) * c
C
C                where c is the coefficient vector for the NBO in AO
C                basis. Suppose that we want to calculate the atomic
C                localization due to atom A form the NBO. Then we form
C                the quantity:
C
C                                 x = c'(T)*c
C
C                where we restrict the summation over elements of c'
C                only to those corresponding to atom A. This in turn
C                means that we only take those rows of S**(1/2)
C                corresponding to atom A in evaluating c'. The obtained
C                number x we identify as the atomic localization
C                corresponding to atom A in the NBO.
C
C                The routine evaluates the complete atomic localization
C                map LOCAL (ATOM,NBO) for each atom and NBO.
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atoms
C                    ATIDX (A)    =  will contain atomic index for
C                                    atom A.
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    C            =  NBAS x NBAS NBO coefficient matrix
C                                    in AO basis.
C                    SHALF        =  full NBAS x NBAS square root of
C                                    the overlap matrix in AO basis.
C
C
C                  Output:
C
C                    LOCAL (A,I)  =  atomic localization content of
C                                    atom A for I-th NBO.
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

         LOGICAL     TOTAL

         INTEGER     ATOM
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NBO

         INTEGER     ATIDX  (1:NATOM)
         INTEGER     BASBEG (1:NATOM)
         INTEGER     BASEND (1:NATOM)

         DOUBLE PRECISION  C     (1:NBAS ,1:NBAS)
         DOUBLE PRECISION  LOCAL (1:NATOM,1:NBAS)
         DOUBLE PRECISION  SHALF (1:NBAS ,1:NBAS)
C
C
C------------------------------------------------------------------------
C
C
C             ...form general atomic index vector.
C
C
         DO 10 ATOM = 1,NATOM
            ATIDX (ATOM) = ATOM
   10    CONTINUE
C
C
C             ...calculate the individual atomic localization contents
C                for each NBO.
C
C
         TOTAL = .FALSE.

         DO 100 NBO = 1,NBAS
            CALL  NLO__EXTRACT_ATOMIC_CONTENT
     +
     +                 ( NBAS,NATOM,
     +                   NATOM,
     +                   ATIDX,
     +                   BASBEG,BASEND,
     +                   TOTAL,
     +                   C (1,NBO),
     +                   SHALF,
     +
     +                           LOCAL (1,NBO) )
     +
     +
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
