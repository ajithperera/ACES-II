         SUBROUTINE  NLO__PARTIAL_SCHMIDT
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWD,DDCOLD,
     +                      DDROWS,DDCOLS,
     +                      ROW,COLC,COLD,
     +                      S,
     +                      SAVES,SAVEC,SAVED,
     +                      C,
     +                      XVEC,
     +                      XMAT,
     +
     +                              D )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__PARTIAL_SCHMIDT
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine orthogonalizes a set of MO vectors to
C                an already orthogonalized set of MO vectors using
C                the Schmidt orthogonalization method. The overlaps
C                needed between the two sets are evaluated from the
C                MO expansion coefficients in terms of AO's and the
C                AO overlap matrix.
C
C                The Schmidt orthogonalization leads to the following
C                expression of the new orthogonal D vectors:
C
C
C                                       M1
C                        D(i) = D(i) - sum S(k,i) C(k)
C                                       k
C
C                where i refers to the initially non-orthogonal D
C                functions and S (k,i) denotes the MO overlap matrix
C                element between the k-th orthogonal C function and
C                the i-th D function.
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

         LOGICAL     SAVES,SAVEC,SAVED

         INTEGER     DDROWC,DDCOLC,DDROWD,DDCOLD,DDROWS,DDCOLS
         INTEGER     I,J,K
         INTEGER     ROW,COLC,COLD

         DOUBLE PRECISION  X

         DOUBLE PRECISION  XVEC  (1:MAX0(COLC,COLD))

         DOUBLE PRECISION  C     (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  D     (1:DDROWD,1:DDCOLD)
         DOUBLE PRECISION  S     (1:DDROWS,1:DDCOLS)
         DOUBLE PRECISION  XMAT  (1:ROW,1:MIN0(COLC,COLD))
C
C
C------------------------------------------------------------------------
C
C
C             ...check dimensions of C,D and S matrices supplied.
C
C
         IF (ROW.GT.DDROWC .OR. COLC.GT.DDCOLC) THEN
             WRITE (*,*) ' Dimensions of matrix C too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWC,DDCOLC,ROW,COLC = ',
     +                     DDROWC,DDCOLC,ROW,COLC
             WRITE (1,*) ' Dimensions of matrix C too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWC,DDCOLC,ROW,COLC = ',
     +                     DDROWC,DDCOLC,ROW,COLC
             STOP
         END IF

         IF (ROW.GT.DDROWD .OR. COLD.GT.DDCOLD) THEN
             WRITE (*,*) ' Dimensions of matrix D too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWD,DDCOLD,ROW,COLD = ',
     +                     DDROWD,DDCOLD,ROW,COLD
             WRITE (1,*) ' Dimensions of matrix D too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWD,DDCOLD,ROW,COLD = ',
     +                     DDROWD,DDCOLD,ROW,COLD
             STOP
         END IF

         IF (ROW.GT.DDROWS .OR. ROW.GT.DDCOLS) THEN
             WRITE (*,*) ' Dimensions of matrix S too small! '
             WRITE (*,*) ' nlo__partial_schmidt '
             WRITE (*,*) ' DDROWS,DDCOLS,ROW = ',DDROWS,DDCOLS,ROW
             WRITE (1,*) ' Dimensions of matrix S too small! '
             WRITE (1,*) ' nlo__partial_schmidt '
             WRITE (1,*) ' DDROWS,DDCOLS,ROW = ',DDROWS,DDCOLS,ROW
             STOP
         END IF
C
C
C             ...proceed according to sizes of C and D. If the size of
C                C is smaller than D, then we generate the COLD x COLC
C                molecular orbital overlap matrix, otherwise we generate
C                its transpose of COLC x COLD size.
C
C
         IF (COLC.LT.COLD) THEN

             CALL    MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                    ( DDROWD,DDCOLD,
     +                      DDROWS,DDCOLS,
     +                      DDROWC,DDCOLC,
     +                      ROW,COLC,
     +                      COLD,
     +                      COLD,ROW,COLC,
     +                      0,0,
     +                      SAVED,SAVES,SAVEC,
     +                      D,S,C,
     +                      XVEC,
     +
     +                             XMAT )
     +
     +
             DO 100 I = 1,COLC
             DO 100 J = 1,COLD
                X = XMAT (J,I)
                DO 110 K = 1,ROW
                   D (K,J) = D (K,J) - X * C (K,I)
  110           CONTINUE
  100        CONTINUE

         ELSE

             CALL    MAT__C_EQ_ORTHOTRAN_OFFDIAG
     +
     +                    ( DDROWC,DDCOLC,
     +                      DDROWS,DDCOLS,
     +                      DDROWD,DDCOLD,
     +                      ROW,COLD,
     +                      COLC,
     +                      COLC,ROW,COLD,
     +                      0,0,
     +                      SAVEC,SAVES,SAVED,
     +                      C,S,D,
     +                      XVEC,
     +
     +                             XMAT )
     +
     +
             DO 200 J = 1,COLD
             DO 200 I = 1,COLC
                X = XMAT (I,J)
                DO 210 K = 1,ROW
                   D (K,J) = D (K,J) - X * C (K,I)
  210           CONTINUE
  200        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
