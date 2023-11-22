         SUBROUTINE  NLO__IMPROVE_2CEN_CONNECTIONS
     +
     +                    ( NSET,
     +                      INDEX,
     +
     +                             AT2CEN )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__IMPROVE_2CEN_CONNECTIONS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Given a set of atomic pair (2 center) indices, this
C                routine reorders them in such a way that it improves
C                connectivity between them, i.e. that each index
C                pair is connected to the preceeding one. By connected
C                we mean that two consecutive atomic index pairs have
C                one atom index in common.
C
C                Example:
C
C                   Consider the atomic index pairs corresponding
C                   to naphtalene:
C
C
C                                   10     2
C                                   /\     /\
C                                  /  \   /  \
C                                9/    \1/    \3
C                                |      |      |
C                                |      |      |
C                                |      |      |
C                                8\    /6\    /4
C                                  \  /   \  /
C                                   \/     \/
C                                   7      5
C
C
C                   given in some random order:
C
C                          9  3  7  1  4  5  1  8  6  1  2
C                         10  4  8  2  5  6 10  9  7  6  3
C
C                This order has only connections at two places: one
C                involving atom 5 and the other one involving atom 6.
C                After exiting this routine the index pairs will be
C                rearranged the following way:
C
C                          9  8  7  6  5  1  1  1  3  2  4
C                         10  9  8  7  6  6  2 10  4  3  5
C
C                by the current implemented algorithm. Note that
C                connectivity is improved but not maximized. There
C                are two disconnected sets of index pairs containing
C                the atoms {1,2,5,6,7,8,9,10} and {2,3,4,5}, which is
C                a vast improvement over the original set but which
C                obviously could be merged to give an overall connected
C                set. The outcome of this routine depends very much
C                on the original order of the index pairs.
C
C                Algorithm:
C
C                   1) Start with the first atom index (9) of the
C                      first pair {9,10} and collect from the list
C                      all those index pairs involving atom 9:
C
C                                      9  8
C                                     10  9
C
C                   2) The lastest new atom index not involving 9
C                      is atom index 8.
C
C                   3) Add all those index pairs from the list
C                      involving atom 8:
C
C                                      9  8  7
C                                     10  9  8
C
C                   4) Identify 7 as the latest new atom index added.
C                      Go back to step 3) but adding all index pairs
C                      with atom 7.
C
C                   5) Keep on cycling between steps 3) and 4) until
C                      no more left index pairs contain the latest
C                      new atom index.
C
C                   6) If more index pairs are still left from the
C                      original list, start a new search for them
C                      going back to 1).
C
C
C                  Input:
C
C                    NSET        =  # of atomic index pairs that
C                                    will be considered for improving
C                                    connectivity.
C                    INDEX        =  will hold reordering indices
C                    AT2CEN (1,N) =  original 1st atomic index of N-th
C                                    pair.
C                    AT2CEN (2,N) =  original 2nd atomic index of N-th
C                                    pair.
C
C
C                  Output:
C
C                    AT2CEN (1,N) =  new improved 1st atomic index
C                                    of N-th pair.
C                    AT2CEN (2,N) =  new improved 2nd atomic index
C                                    of N-th pair.
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

         INTEGER     BASE
         INTEGER     CYCLE
         INTEGER     I,N
         INTEGER     NEW
         INTEGER     NOLD
         INTEGER     NSET
         INTEGER     TEST

         INTEGER     INDEX (1:NSET)

         INTEGER     AT2CEN  (1:2,1:NSET)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine base value to be used for compressing
C                both indices and set index array to zero.
C
C
         BASE = 0
         DO 10 I = 1,NSET
            INDEX (I) = 0
            BASE = MAX0 (BASE,AT2CEN (2,I))
   10    CONTINUE
         BASE = BASE + 1
C
C
C             ...start search.
C
C
         N = 1
         NOLD = 1
         INDEX (1) = 1
         TEST = AT2CEN  (1,1)

         DO 100 CYCLE = 1,NSET

            DO 110 I = 1,NSET
               IF (INDEX (I).EQ.0) THEN
                   IF (AT2CEN (1,I).EQ.TEST) THEN
                       N = N + 1
                       INDEX (I) = N
                       NEW = AT2CEN (2,I)
                   ELSE IF (AT2CEN (2,I).EQ.TEST) THEN
                       N = N + 1
                       INDEX (I) = N
                       NEW = AT2CEN (1,I)
                   END IF
               END IF
  110       CONTINUE

            IF (N.EQ.NSET) THEN
C
C
C             ...search complete. Reorder the atomic index pairs by
C                compressing and decompressing the pairs and using
C                one part of the index pair arrays for reordering.
C
C
                DO 120 I = 1,NSET
                   AT2CEN (1,I) = AT2CEN (1,I) * BASE + AT2CEN (2,I)
  120           CONTINUE

                DO 130 I = 1,NSET
                   AT2CEN (2,INDEX (I)) = AT2CEN (1,I)
  130           CONTINUE

                DO 140 I = 1,NSET
                   AT2CEN (1,I) = AT2CEN (2,I) / BASE
                   AT2CEN (2,I) = AT2CEN (2,I) - BASE * AT2CEN (1,I)
  140           CONTINUE

                RETURN

            ELSE IF (N.EQ.NOLD) THEN
C
C
C             ...one connectivity cycle complete. Check for remaining
C                index pairs and enter new cycle if any index pairs
C                left.
C
C
                DO 150 I = 1,NSET
                   IF (INDEX (I).EQ.0) THEN
                       N = N + 1
                       NOLD = N
                       TEST = AT2CEN (1,I)
                       INDEX (I) = N
                       GOTO 100
                   END IF
  150           CONTINUE
            ELSE
C
C
C             ...cycle not complete. Search for more index pairs
C                to connect.
C
C
                NOLD = N
                TEST = NEW
            END IF

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
