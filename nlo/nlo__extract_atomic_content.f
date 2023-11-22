         SUBROUTINE  NLO__EXTRACT_ATOMIC_CONTENT
     +
     +                    ( NBAS,NATOM,
     +                      NCEN,
     +                      ATIDX,
     +                      BASBEG,BASEND,
     +                      TOTAL,
     +                      C,
     +                      SHALF,
     +
     +                              X )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__EXTRACT_ATOMIC_CONTENT
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine extracts from a given MO coefficient
C                vector C in AO basis its atomic center content
C                in terms of the symmetrized AO basis corresponding
C                to NCEN atoms. The content X is simply the number
C                which results in the multiplication:
C
C                                X = C'(T) * C'
C
C                with C' being the MO coefficent expansion in terms
C                of the symmetrized AO basis:
C
C                             C' = S**(1/2) * C
C
C                and restricting the summation to only those rows
C                of S**(1/2) belonging to the NCEN atoms.
C
C                The routine allows for two options:
C
C                    1) the total atomic content is evaluated with
C                       summation over all NCEN atoms and adding up
C                       the individual atomic contents. The result
C                       is placed in X (1).
C
C                    2) the individual atomic contents are evaluated
C                       and the results placed into the NCEN first
C                       places of array X.
C
C                Which option is performed is controlled by the
C                keyword TOTAL, which is true, if the first option
C                is wanted and false for option 2).
C
C
C                  Input:
C
C                    NBAS         =  total # of AO's in AO basis.
C                    NATOM        =  total # of atoms.
C                    NCEN         =  # of atomic centers for content
C                                    evaluation.
C                    ATIDX (A)    =  atomic index for atom A for
C                                    content evaluation.
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
C                    TOTAL        =  is true, if the total atomic
C                                    content is wanted, and false, if
C                                    the individual ones are wanted.
C                    C            =  MO coefficient vector in AO basis.
C                    SHALF        =  full NBAS x NBAS square root of
C                                    the overlap matrix in AO basis.
C
C
C                  Output:
C
C                    X            =  atomic content of MO in C.
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
         INTEGER     I,J,N
         INTEGER     JBEG,JEND
         INTEGER     NATOM
         INTEGER     NBAS
         INTEGER     NCEN

         INTEGER     ATIDX  (1:NCEN )
         INTEGER     BASBEG (1:NATOM)
         INTEGER     BASEND (1:NATOM)

         DOUBLE PRECISION  CS
         DOUBLE PRECISION  XATOM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  C (1:NBAS)
         DOUBLE PRECISION  X (1:NCEN)

         DOUBLE PRECISION  SHALF (1:NBAS,1:NBAS)

         DATA  ZERO  /0.D0/
C
C
C------------------------------------------------------------------------
C
C
C             ...form X = C'(T)*C using only atoms in ATIDX.
C
C
         X (1) = ZERO

         DO 100 N = 1,NCEN
            ATOM = ATIDX (N)
            JBEG = BASBEG (ATOM)
            JEND = BASEND (ATOM)

            XATOM = ZERO
            DO 110 J = JBEG,JEND
               CS = ZERO
               DO 120 I = 1,NBAS
                  CS = CS + C (I) * SHALF (I,J)
  120          CONTINUE
               XATOM = XATOM + CS * CS
  110       CONTINUE

            IF (TOTAL) THEN
                X (1) = X (1) + XATOM
            ELSE
                X (N) = XATOM
            END IF

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
