         SUBROUTINE  NLO__SYMMETRIC_JACOBI
     +
     +                    ( DDROWX,DDCOLX,
     +                      DDROWC,DDCOLC,
     +                      DDVECD,
     +                      N,
     +                      ZROW,ZCOL,
     +                      REVERS,
     +                      LUSED,
     +                      X,
     +
     +                              NITER,
     +                              D,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__SYMMETRIC_JACOBI
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine does a symmetric Jacobi transformation
C                on the input matrix X. This procedure differs from
C                the standard Jacobi procedure in that here a bunch of
C                identical 2x2 submatrices are rotated all at once
C                by a symmetrized Jacobi rotation matrix instead of each
C                2x2 submatrix at a time. The reason for doing this
C                is to ensure any symmetry preservation that might
C                be present in the original matrix X. A standard serial
C                Jacobi diagonalization, or one based on the Householder
C                tridiagonalization procedure, would result in mixing
C                of eigenvectors corresponding to different symmetries
C                if their eigenvalues are degenerate to the precision
C                used.
C
C                A particular symmetrized Jacobi rotation performed
C                on an already partly transformed matrix X has the
C                following essential steps:
C
C
C                 1) Find largest element in magnitude |Xij| of X.
C                    If |Xij| below the transformation threshold
C                    => stop.
C
C                 2) Find all elements Xkl of X for which |Xkl| = |Xij|,
C                    Xkk = Xii and Xll = Xjj to within the symmetry
C                    recognition accuracy limits. All these 2x2
C                    subblocks of X are related by symmetry and
C                    should be treated all at once in order not to
C                    break the symmetry.
C
C                 3) Set up the symmetrized Jacobi rotation matrix.
C
C                 4) Perform a symmetrized Jacobi rotation on all the
C                    relevant rows and columns of X.
C
C                 5) Return to 1).
C
C
C                The routine allows for partial transformation of X
C                as well, using the transformation delimiters ZROW
C                and ZCOL, their meaning being best explained by the
C                following picture showing the final X on exit from
C                this routine:
C
C                                              
C                              *
C                              * *
C                              * * *
C                              * * * *
C                    ZROW  ->  0 0 0 0 *
C                              0 0 0 0 0 *
C                              0 0 0 0 0 0 *
C                              0 0 0 0 0 0 0 *
C                              0 0 0 0 0 0 0 * *
C                              0 0 0 0 0 0 0 * * *
C                              0 0 0 0 0 0 0 * * * *
C
C                                          |
C
C                                         ZCOL
C
C                Several special transformation cases arise from
C                different values of these delimiters:
C
C                  ZROW=1,ZCOL=1 : Only the first column of X
C                                  below the first element X(1,1)
C                                  will be zero.
C
C                  ZROW=N,ZCOL=N : Only the last row of X to the
C                                  left of the last element X(N,N)
C                                  will be zero.
C
C                  ZROW=1,ZCOL=N : This corresponds to a complete
C                                  diagonalization of X.
C
C
C                  Input:
C
C                    DDROWX,DDCOLX  =  declared dimensions of matrix X
C                                      to be transformation.
C                    DDROWC,DDCOLC  =  declared dimensions of
C                                      transformation matrix C.
C                    DDVECD         =  declared dimensions of vector D.
C                    N              =  actual dimensions of N x N
C                                      matrices X,C and vector D.
C                    ZROW,ZCOL      =  defines the section to be
C                                      transformed.
C                    REVERS         =  if true => possible eigenvalues
C                                      and associated eigenvectors are
C                                      in decreasing sequence, if
C                                      false in ascending sequence.
C                    LUSED          =  will hold info about used rows/
C                                      columns of X.
C                    X              =  N x N matrix X to be transformed.
C
C
C                  Output:
C
C                    NITER          =  # of iterations performed.
C                    D              =  diagonal elements of transformed
C                                      X. The section of the D vector
C                                      corresponding to eigenvalues
C                                      will be in the specified order.
C                    C              =  transformation matrix of X.
C                                      The section of columns of C
C                                      corresponding to eigenvectors
C                                      of X will be in the specified
C                                      order.
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

         CHARACTER*40  SYMINFO

         LOGICAL     CASE1,CASE2
         LOGICAL     KKTRUE,KLTRUE,LLTRUE
         LOGICAL     KUSED
         LOGICAL     MAXSYM
         LOGICAL     REVERS
         LOGICAL     SKIPK,SKIPL

         INTEGER     DDROWX,DDCOLX,DDROWC,DDCOLC,DDVECD
         INTEGER     I,J,K,L,M,N
         INTEGER     IBEG,IEND,KBEG,KEND,LBEG
         INTEGER     IMAX,JMAX,KMAX,LMAX
         INTEGER     IX,JX
         INTEGER     KADD,LADD
         INTEGER     KBASE,KSUB
         INTEGER     KR,KW,KX
         INTEGER     L1CACHE
         INTEGER     LTEST
         INTEGER     MXNROT,MXBLCK,MXITER
         INTEGER     MXNSUB
         INTEGER     NITER
         INTEGER     NKQ,NLQ
         INTEGER     NSUB
         INTEGER     NROT
         INTEGER     ZROW,ZCOL

         INTEGER     LUSED (1:N)

         DOUBLE PRECISION  ABSXIJ,ABSXKL
         DOUBLE PRECISION  AVGEKK,AVGEKL,AVGELL
         DOUBLE PRECISION  COSINE,SINE
         DOUBLE PRECISION  E,Q,S,T
         DOUBLE PRECISION  SQ2INV
         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  SYMDIA,SYMOFF,CONVGE
         DOUBLE PRECISION  USED
         DOUBLE PRECISION  XII,XIJ,XJJ,XKK,XLL,XKL
         DOUBLE PRECISION  XIJMAX,XIJMIN,XKLMAX
         DOUBLE PRECISION  Z11,Z21,Z22
C         DOUBLE PRECISION  XIIMIN, XJJMIN
         DOUBLE PRECISION  XIIMIN, XJJMIN
         DOUBLE PRECISION  ZERO,QUART,HALF,ONE,TWO,VSMALL

         DOUBLE PRECISION  TRACE

         DOUBLE PRECISION  D (1:DDVECD)

         DOUBLE PRECISION  C (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  X (1:DDROWX,1:DDCOLX)

         DATA  SYMDIA /5.D-6/
         DATA  SYMOFF /1.D-6/
C         DATA  SYMDIA /5.D-3/
C         DATA  SYMOFF /1.D-3/
C         DATA  CONVGE /1.D-12/
C         DATA  VSMALL /1.D-12/
         DATA  CONVGE /1.D-16/
         DATA  VSMALL /1.D-16/
         DATA  TWO    /2.D0/
         DATA  ONE    /1.D0/
         DATA  USED   /-1.D0/
         DATA  SQ2INV /0.70710678118655D0/
         DATA  HALF   /0.5D0/
         DATA  QUART  /0.25D0/
         DATA  ZERO   /0.D0/
C hughes         DATA  MXITER /1000000/
         DATA  MXITER /10000000/

         DATA  SYMINFO /'High symmetry in matrix X! Proceeding...'/
C
C
C             ...hardwired arrays.
C
C
         PARAMETER  (L1CACHE = 16384)
         PARAMETER  (MXNROT  = 20)
         PARAMETER  (MXBLCK  = ((3*L1CACHE/4)-3*MXNROT*MXNROT)/MXNROT)
         PARAMETER  (MXNSUB  = (MXNROT*MXNROT)/4)

         INTEGER     KLQ (1:2*MXNROT)
         INTEGER     KP  (1:MXNSUB)
         INTEGER     LP  (1:MXNSUB)

         DOUBLE PRECISION  P (1:MXNSUB)

         DOUBLE PRECISION  ROT (1:MXNROT,1:MXNROT)
         DOUBLE PRECISION  Y   (1:MXNROT,1:MXNROT)
         DOUBLE PRECISION  Z   (1:MXNROT,1:MXNROT)
         DOUBLE PRECISION  W   (1:MXBLCK,1:MXNROT)
C
C
C------------------------------------------------------------------------
C
C
C             ...check passed dimensions of matrices X,C and vector D.
C
C
         IF (N.GT.DDROWX .OR. N.GT.DDCOLX) THEN
             WRITE (1,*) ' Dimension of matrix X too small: '
             WRITE (1,*) ' nlo__symmetric_jacobi '
             WRITE (1,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             WRITE (*,*) ' Dimension of matrix X too small: '
             WRITE (*,*) ' nlo__symmetric_jacobi '
             WRITE (*,*) ' DDROWX,DDCOLX,N = ',DDROWX,DDCOLX,N
             STOP
         END IF

         IF (N.GT.DDROWC .OR. N.GT.DDCOLC) THEN
             WRITE (1,*) ' Dimension of matrix C too small: '
             WRITE (1,*) ' nlo__symmetric_jacobi '
             WRITE (1,*) ' DDROWC,DDCOLC,N = ',DDROWC,DDCOLC,N
             WRITE (*,*) ' Dimension of matrix C too small: '
             WRITE (*,*) ' nlo__symmetric_jacobi '
             WRITE (*,*) ' DDROWC,DDCOLC,N = ',DDROWC,DDCOLC,N
             STOP
         END IF

         IF (N.GT.DDVECD) THEN
             WRITE (1,*) ' Dimension of vector D too small: '
             WRITE (1,*) ' nlo__symmetric_jacobi '
             WRITE (1,*) ' DDVECD,N = ',DDVECD,N
             WRITE (*,*) ' Dimension of vector D too small: '
             WRITE (*,*) ' nlo__symmetric_jacobi '
             WRITE (*,*) ' DDVECD,N = ',DDVECD,N
             STOP
         END IF
C
C
C             ...check, if meaningful diagonalization section
C                is specified.
C
C
         IF (ZROW.LT.1.OR.ZROW.GT.N.OR.ZCOL.LT.1.OR.ZCOL.GT.N) THEN
             WRITE (1,*) ' Diagonalization section out of range! '
             WRITE (1,*) ' nlo__symmetric_jacobi '
             WRITE (1,*) ' ZROW,ZCOL,N = ',ZROW,ZCOL,N
             WRITE (*,*) ' Diagonalization section out of range! '
             WRITE (*,*) ' nlo__symmetric_jacobi '
             WRITE (*,*) ' ZROW,ZCOL,N = ',ZROW,ZCOL,N
             STOP
         END IF
C
C
C             ...handle special case, if order of matrix is 1.
C
C
         IF (N.EQ.1) THEN
             D (1) = X (1,1)
             C (1,1) = ONE
             RETURN
         END IF
C
C
C             ...preliminary steps before starting the iterations.
C                Set initial eigenvector matrix equal to a unit
C                matrix. Find initial 2x2 subblock Xii,Xjj,Xij with
C                index relation i < j containing the largest absolute
C                offdiagonal element in X within the specified
C                diagonalization range and symmetrize X.
C                
C
C
         DO 10 J = 1,N
            DO 20 I = 1,N
               C (I,J) = ZERO
   20       CONTINUE
            C (J,J) = ONE
   10    CONTINUE

         XIJMAX = ZERO

         M = MIN0 (ZROW-1,ZCOL)

         DO 30 I = 1,M
            D (I) = X (I,I)
            DO 32 J = I+1,ZROW-1
               X (I,J) = X (J,I)
   32       CONTINUE
            DO 34 J = ZROW,N
               ABSXIJ = DABS (X (J,I))
               IF (ABSXIJ.GT.XIJMAX) THEN
                   IMAX = I
                   JMAX = J
                   XIJMAX = ABSXIJ
               END IF
               X (I,J) = X (J,I)
   34       CONTINUE
   30    CONTINUE

         DO 40 I = ZROW,ZCOL
            D (I) = X (I,I)
            DO 42 J = I+1,N
               ABSXIJ = DABS (X (J,I))
               IF (ABSXIJ.GT.XIJMAX) THEN
                   IMAX = I
                   JMAX = J
                   XIJMAX = ABSXIJ
               END IF
               X (I,J) = X (J,I)
   42       CONTINUE
   40    CONTINUE

         DO 50 I = ZCOL+1,N
            D (I) = X (I,I)
            DO 52 J = I+1,N
               X (I,J) = X (J,I)
   52       CONTINUE
   50    CONTINUE

         XII = X (IMAX,IMAX)
         XIJ = X (JMAX,IMAX)
         XJJ = X (JMAX,JMAX)

C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Initial X matrix ',
C     +                  DDROWX,DDCOLX,
C     +                  N,0,
C     +                  X )
C     +
C     +
C         WRITE (1,*) ' IMAX,JMAX = ',IMAX,JMAX
C         WRITE (1,*) ' XIJMAX = ',XIJMAX
C         WRITE (1,*) ' MXNROT = ',MXNROT
C
C
C             ...start the symmetrized Jacobi iterations.
C
C
         DO 1000 NITER = 1,MXITER
C
C
C             ...if convergence was achieved, order eigenvalues
C                and eigenvectors in desired order (if any) and
C                return.
C
C
            IF (XIJMAX.LT.CONVGE) THEN

CSSS                WRITE (*,*) ' # of iterations = ',NITER

                DO 100 I = ZROW,ZCOL
                   K = I
                   E = D (I)
                   IF (REVERS) THEN
                       DO 110 J = I+1,ZCOL
                          IF (D (J).GT.E) THEN
                              K = J
                              E = D (J)
                          END IF
  110                  CONTINUE
                   ELSE
                       DO 120 J = I+1,ZCOL
                          IF (D (J).LT.E) THEN
                              K = J
                              E = D (J)
                          END IF
  120                  CONTINUE
                   END IF

                   IF (K.NE.I) THEN
                       D (K) = D (I)
                       D (I) = E
                       DO 130 J = 1,N
                          E = C (J,I)
                          C (J,I) = C (J,K)
                          C (J,K) = E
  130                  CONTINUE
                   END IF
  100           CONTINUE

                RETURN
            END IF
C
C
C             ...determine the set {Xkk,Xll,Xkl;k<l} of all identical
C                2x2 subblocks, including the first 2x2 subblock
C                {Xii,Xjj,Xij;i<j}, within symmetry defining limits
C                and within the specified diagonalization range.
C                The set {Xkk,Xll,Xkl;k<l} should be such that both
C                sets {k} and {l} of obtained k and l indices are
C                nonoverlapping.
C
C                Determine also the following information:
C
C                 a) number of 2x2 subblocks in the set: NSUB
C                 b) unique lists of {k} and {l} without reps: KQ,LQ
C                    stored into common array KLQ. First KQ is stored
C                    starting at KLQ (1) and LQ is stored starting at
C                    KLQ (MXNROT+1). Then KLQ is compressed for later
C                    use.
C                 c) dimensions of the {k} and {l} sets: NKQ,NLQ
C                 d) positions of NSUB k/l indices in KQ,LQ: KP,LP
C                 e) largest absolute offdiagonal element XKLMAX in X
C                    excluding the rows and columns corresponding
C                    to the found {k} and {l} sets to be changed
C                    by rotation.
C
C                Step e) is done here to avoid a complete search in X
C                for the largest absolute offdiagonal element of X
C                after! the rotation has been done. After having the
C                info of step e), the search for the largest absolute
C                offdiagonal element of X boils down to a search
C                between the rotated rows and columns and the maximum
C                found in step e).
C
C                The eigenvector D will be used to accumulate the
C                info necessary to evaluate step e). Hence the diagonal
C                elements of X must be restored into D after step e)
C                is complete.
C
C                If there are more symmetry related elements than
C                the rotation matrix is able to accomodate, the
C                routine proceeds with those symmetry elements already
C                found, leaving the extra ones aside. Since this might
C                (rarely) result in symmetry loss, a message is issued
C                informing the user.
C
C
            NKQ = 0
            NLQ = 0
            NSUB = 0
            MAXSYM = .FALSE.

            XIJMIN = (ONE - SYMOFF) * XIJMAX
C            XIIMIN = (ONE - SYMDIA) * XII
C            XJJMIN = (ONE - SYMDIA) * XJJ
            XIIMIN = (ONE - SYMDIA) * XII
            XJJMIN = (ONE - SYMDIA) * XJJ

            DO 200 L = 1,N
               D (L) = ZERO
               LUSED (L) = 0
  200       CONTINUE

            DO 210 K = 1,ZCOL
               SKIPK = .FALSE.

               IF (LUSED (K).EQ.0) THEN
                   LBEG = MAX0 (K+1,ZROW)

                   IF (.NOT.MAXSYM) THEN
                       KUSED = .FALSE.
                       XKK = X (K,K)
                       DO 230 L = LBEG,N
                          SKIPL = .FALSE.
                          XKL = X (L,K)
                          XLL = X (L,L)
                          KLTRUE = DABS (XKL) .GE. XIJMIN
                          KKTRUE = DABS (XKK) .GE. XIIMIN
C                          KKTRUE = DABS (XKK - XII) .LE. SYMDIA
                          LLTRUE = DABS (XLL) .GE. XJJMIN
C                          LLTRUE = DABS (XLL - XJJ) .LE. SYMDIA

                          IF (KKTRUE .AND. KLTRUE .AND. LLTRUE) THEN
                              KADD = 0
                              LADD = 0

                              IF (.NOT.KUSED) THEN
                                  KADD = 1
                                  NKQ = NKQ + 1
                                  KLQ (NKQ) = K
                                  KUSED = .TRUE.
                              END IF

                              LTEST = LUSED (L)
                              IF (LTEST.EQ.0) THEN
                                  LADD = 1
                                  NLQ = NLQ + 1
                                  KLQ (MXNROT+NLQ) = L
                                  LUSED (L) = NLQ
                              END IF

                              IF ((NKQ+NLQ).GT.MXNROT) THEN
                                  WRITE (*,*) SYMINFO
                                  WRITE (*,*) ' nlo__symmetric_jacobi '
                                  WRITE (1,*) SYMINFO
                                  WRITE (1,*) ' nlo__symmetric_jacobi '
                                  NKQ = NKQ - KADD
                                  NLQ = NLQ - LADD
                                  LUSED (L) = LTEST
                                  MAXSYM = .TRUE.
                                  GOTO 232
                              ELSE
                                  NSUB = NSUB + 1
                                  P (NSUB) = DSIGN (ONE,XKL*(XLL-XKK))
                                  KP (NSUB) = NKQ
                                  LP (NSUB) = LUSED (L)
                                  SKIPK = .TRUE.
                                  SKIPL = .TRUE.
                              END IF

                          END IF

                          IF (SKIPL) THEN
                              D (K) = USED
                              D (L) = USED
                          END IF

  230                  CONTINUE
                   END IF
               ELSE
                   SKIPK = .TRUE.
               END IF

  232          IF (.NOT.SKIPK) THEN
                   DO 220 L = LBEG,N
                      S = D (L)
                      IF (S.GE.ZERO) THEN
                          D (L) = DMAX1 (DABS(X(L,K)),S)
                      END IF
  220              CONTINUE
               END IF

  210       CONTINUE
C
C
C             ...check, if the set {Xkk,Xll,Xkl;k<l} is empty.
C                If it is, it signals a bug someplace, since the
C                first set {Xii,Xjj,Xij;i<j} must be recovered.
C
C
            IF (NSUB.EQ.0) THEN
                WRITE (*,*) ' Problems finding 2x2 subblocks of X! '
                WRITE (*,*) ' NSUB = ',NSUB
                WRITE (*,*) ' nlo__symmetric_jacobi '
                WRITE (1,*) ' Problems finding 2x2 subblocks of X! '
                WRITE (1,*) ' NSUB = ',NSUB
                WRITE (1,*) ' nlo__symmetric_jacobi '
                STOP
            END IF
C
C
C             ...extract largest absolute offdiagonal element XKLMAX
C                in X excluding the rotation rows and columns using
C                the info present in array D. Reset all diagonal
C                elements of X into D.
C
C
            XKLMAX = ZERO
            DO 240 K = ZROW,N
               IF (D (K).GT.XKLMAX) THEN
                   KMAX = K
                   XKLMAX = D (K)
               END IF
  240       CONTINUE

            DO 250 L = 1,ZCOL
               SKIPL = D (L) .EQ. USED
               IF (DABS (X(L,KMAX)).EQ.XKLMAX .AND. .NOT.SKIPL) THEN
                   LMAX = L
               END IF
  250       CONTINUE

            DO 260 I = 1,N
               D (I) = X (I,I)
  260       CONTINUE
C
C
C             ...calculate cosine and absolute value of sine of
C                the Jacobi rotation(s).
C
C
            S = DABS (XJJ - XII)

            IF (S.LT.VSMALL*XIJMAX) THEN
                COSINE = SQ2INV
                SINE = SQ2INV
            ELSE
                T = XIJMAX / DABS (XJJ - XII)
                Q = QUART / DSQRT (QUART + T*T)
                COSINE = DSQRT (HALF + Q)
                SINE = TWO*T*Q / COSINE
            END IF

C            WRITE (1,*) ' T,Q = ',T,Q
C            WRITE (1,*) ' COS,SIN = ',COSINE,SINE
C
C
C             ...handle special case of 2x2 rotation matrix (simple
C                Jacobi case) separately. For comments on the
C                sequence of steps consult the general case below.
C
C
            IF (NSUB.EQ.1) THEN

                SINE = DSIGN (SINE,P (1))

                IX = KLQ (1)
                JX = KLQ (MXNROT+1)

                CASE1 = IX.GE.ZROW
                CASE2 = JX.LE.ZCOL

                DO 9000 KBASE = 0,N-1,MXBLCK
                   KSUB = MIN0 (MXBLCK,N-KBASE)
                   DO 9010 K = 1,KSUB
                      W (K,1) = COSINE * C (KBASE+K,IX)
                      W (K,2) =   SINE * C (KBASE+K,IX)
 9010              CONTINUE
                   DO 9020 K = 1,KSUB
                      W (K,1)        = W (K,1) -   SINE * C (KBASE+K,JX)
                      C (KBASE+K,JX) = W (K,2) + COSINE * C (KBASE+K,JX)
 9020              CONTINUE
                   DO 9030 K = 1,KSUB
                      C (KBASE+K,IX) = W (K,1)
 9030              CONTINUE
 9000           CONTINUE

                Z11 = ZERO
                Z21 = ZERO
                Z22 = ZERO

                DO 9100 KBASE = 0,N-1,MXBLCK
                   KSUB = MIN0 (MXBLCK,N-KBASE)

                   DO 9110 K = 1,KSUB
                      W (K,1) = COSINE * X (KBASE+K,IX)
                      W (K,2) =   SINE * X (KBASE+K,IX)
 9110              CONTINUE
                   DO 9120 K = 1,KSUB
                      W (K,1)        = W (K,1) -   SINE * X (KBASE+K,JX)
                      W (K,2)        = W (K,2) + COSINE * X (KBASE+K,JX)
                      X (KBASE+K,JX) = W (K,2)
 9120              CONTINUE
                   DO 9130 K = 1,KSUB
                      X (KBASE+K,IX) = W (K,1)
 9130              CONTINUE

                   K = IX - KBASE
                   L = JX - KBASE
                   KKTRUE = K.GT.0 .AND. K.LE.KSUB
                   LLTRUE = L.GT.0 .AND. L.LE.KSUB

                   IF (KKTRUE .AND. LLTRUE) THEN
                       Z11 = Z11 + COSINE * W (K,1) -   SINE * W (L,1)
                       Z21 = Z21 +   SINE * W (K,1) + COSINE * W (L,1)
                       W (K,1) = ZERO
                       W (L,1) = ZERO
                       Z22 = Z22 +   SINE * W (K,2) + COSINE * W (L,2)
                       W (K,2) = ZERO
                       W (L,2) = ZERO
                   ELSE IF (KKTRUE) THEN
                       Z11 = Z11 + COSINE * W (K,1)
                       Z21 = Z21 +   SINE * W (K,1)
                       W (K,1) = ZERO
                       Z22 = Z22 +   SINE * W (K,2)
                       W (K,2) = ZERO
                   ELSE IF (LLTRUE) THEN
                       Z11 = Z11 -   SINE * W (L,1)
                       Z21 = Z21 + COSINE * W (L,1)
                       W (L,1) = ZERO
                       Z22 = Z22 + COSINE * W (L,2)
                       W (L,2) = ZERO
                   END IF

                   IF (CASE1) THEN
                       KBEG = 1
                   ELSE
                       KBEG = MAX0 (1,ZROW-KBASE)
                   END IF

                   DO 9140 K = KBEG,KSUB
                      ABSXKL = DABS (W (K,1))
                      IF (ABSXKL.GT.XKLMAX) THEN
                          KMAX = KBASE + K
                          LMAX = IX
                          XKLMAX = ABSXKL
                      END IF
 9140              CONTINUE

                   IF (CASE2) THEN
                       KEND = KSUB
                   ELSE
                       KEND = MIN0 (KSUB,ZCOL-KBASE)
                   END IF

                   DO 9150 K = 1,KEND
                      ABSXKL = DABS (W (K,2))
                      IF (ABSXKL.GT.XKLMAX) THEN
                          KMAX = KBASE + K
                          LMAX = JX
                          XKLMAX = ABSXKL
                      END IF
 9150              CONTINUE

 9100           CONTINUE

                D (IX) = Z11
                D (JX) = Z22
                X (IX,IX) = Z11
                X (JX,IX) = Z21
                X (IX,JX) = Z21
                X (JX,JX) = Z22

                ABSXKL = DABS (Z21)
                IF (ABSXKL.GT.XKLMAX) THEN
                    KMAX = IX
                    LMAX = JX
                    XKLMAX = ABSXKL
                END IF

                DO 9200 K = 1,N
                   X (IX,K) = X (K,IX)
                   X (JX,K) = X (K,JX)
 9200           CONTINUE

            ELSE
C
C
C
C
C             ...the general case of n x n rotation matrix with n > 2.
C                Update LP array and compress {k} and {l} index sets
C                in array KLQ.
C
C
                DO 270 I = 1,NSUB
                   LP (I) = LP (I) + NKQ
  270           CONTINUE

                DO 280 L = 1,NLQ
                   KLQ (NKQ+L) = KLQ (MXNROT+L)
  280           CONTINUE
C
C
C             ...find the symmetrized Jacobi rotation matrix.
C                Note, that its dimension will be NROT x NROT,
C                where NROT is the sum of all distinct rows and
C                columns involved.
C
C                Calculate initial Jacobi rotation matrix.
C
C
                NROT = NKQ + NLQ

                DO 300 J = 1,NROT
                   DO 310 I = 1,NROT
                      ROT (I,J) = ZERO
  310              CONTINUE
                   ROT (J,J) = ONE
  300           CONTINUE

                K = KP (1)
                L = LP (1)
                S = DSIGN (SINE,P (1))

                ROT (K,K) = COSINE
                ROT (L,K) = - S
                ROT (K,L) = S
                ROT (L,L) = COSINE

                DO 320 J = 2,NSUB
                   K = KP (J)
                   L = LP (J)
                   S = DSIGN (SINE,P (J))
                   DO 330 I = 1,NROT
                      T = ROT (I,K)
                      ROT (I,K) = COSINE * T - S * ROT (I,L)
                      ROT (I,L) = COSINE * ROT (I,L) + S * T
  330              CONTINUE
  320           CONTINUE
C
C
C             ...average diagonal KK- and LL-block elements and
C                offdiagonal KL-block elements, and reset rotation
C                matrix elements with averages. The resulting
C                averaged rotation matrix is placed in working array
C                Y and is skew symmetric, that is Y (I,J) = - Y (J,I),
C                and is not yet orthonormal. For efficiency, the
C                transpose of the nonorthonormal averaged Y is first
C                produced first, since this makes the Lowdin
C                orthonormalization procedure more efficient.
C
C
                AVGEKK = ROT (1,1)
                DO 400 K = 2,NKQ
                   AVGEKK = AVGEKK + ROT (K,K)
  400           CONTINUE

                AVGELL = ROT (NKQ+1,NKQ+1)
                DO 410 L = NKQ+2,NROT
                   AVGELL = AVGELL + ROT (L,L)
  410           CONTINUE

                AVGEKL = ZERO
                DO 420 I = 1,NSUB
                   K = KP (I)
                   L = LP (I)
                   AVGEKL = AVGEKL + DABS (ROT(K,L)) + DABS (ROT(L,K))
  420           CONTINUE
                AVGEKL = AVGEKL / DFLOAT (2*NSUB)

                DO 430 K = 1,NKQ
                   DO 440 I = 1,NROT
                      Y (I,K) = ZERO
                      Z (I,K) = ZERO
                      ROT (I,K) = ZERO
  440              CONTINUE
                   Y (K,K) = AVGEKK
  430           CONTINUE

                DO 450 L = NKQ+1,NROT
                   DO 460 I = 1,NROT
                      Y (I,L) = ZERO
                      Z (I,L) = ZERO
                      ROT (I,L) = ZERO
  460              CONTINUE
                   Y (L,L) = AVGELL
  450           CONTINUE

                DO 470 I = 1,NSUB
                   K = KP (I)
                   L = LP (I)
                   S = DSIGN (AVGEKL,P (I))
                   Y (K,L) = - S
                   Y (L,K) =   S
  470           CONTINUE
C
C
C             ...Lowdin orthonormalize columns of nonorthonormal
C                average rotation matrix to ensure orthogonal
C                rotation matrix. Note that Matrix ROT is free to use
C                and, as well as matrix Z, has been set to zero
C                to accumulate the overlap and the inverse square
C                root of the overlap matrix.
C
C
                DO 500 K = 1,NROT
                DO 500 J = 1,NROT
                   S = Y (J,K)
                   DO 510 I = J,NROT
                      ROT (I,J) = ROT (I,J) + S * Y (I,K)
  510              CONTINUE
  500           CONTINUE

                CALL  MAT__DIAGONALIZE_REAL_SYMMETRIC
     +
     +                     ( MXNROT,MXNROT,MXNROT,
     +                       NROT,
     +                       REVERS,
     +
     +                              P,
     +                              ROT )
     +
     +
                DO 520 K = 1,NROT
                   S = ONE / DSQRT (P (K))
                   DO 530 J = 1,NROT
                      T = S * ROT (J,K)
                      DO 540 I = 1,NROT
                         Z (I,J) = Z (I,J) + T * ROT (I,K)
  540                 CONTINUE
  530              CONTINUE
  520           CONTINUE

                DO 550 J = 1,NROT
                   DO 560 I = 1,NROT
                      S = ZERO
                      DO 570 K = 1,NROT
                         S = S + Y (K,I) * Z (K,J)
  570                 CONTINUE
                      ROT (I,J) = S
  560              CONTINUE
                   DO 580 K = 1,NROT
                      Z (K,J) = ZERO
  580              CONTINUE
  550           CONTINUE
C
C
C             ...orthonormal average rotation matrix ROT is ready
C                and scratch array Z has been reset to zero for further
C                use. Update the eigenvector matrix C first by
C                performing a row-blocked matrix multiplication:
C                C = C*ROT with the rotation matrix.
C
C
                DO 600 KBASE = 0,N-1,MXBLCK
                   KSUB = MIN0 (MXBLCK,N-KBASE)

                   IX = KLQ (1)
                   DO 610 J = 1,NROT
                      S = ROT (1,J)
                      DO 620 K = 1,KSUB
                         W (K,J) = S * C (KBASE+K,IX)
  620                 CONTINUE
  610              CONTINUE

                   DO 630 I = 2,NROT
                      IX = KLQ (I)
                      DO 640 J = 1,NROT
                         S = ROT (I,J)
                         DO 650 K = 1,KSUB
                            W (K,J) = W (K,J) + S * C (KBASE+K,IX)
  650                    CONTINUE
  640                 CONTINUE
  630              CONTINUE

                   DO 660 J = 1,NROT
                      JX = KLQ (J)
                      DO 670 K = 1,KSUB
                         C (KBASE+K,JX) = W (K,J)
  670                 CONTINUE
  660              CONTINUE

  600           CONTINUE
C
C
C             ...update the X matrix performing the matrix
C                multiplications:
C
C                          X = ROT (transposed) * X * ROT
C
C                Procedure:
C
C                Loop over block of rows K.
C                   For each K block:
C                     i) Perform W = X * ROT for all rows in block
C                        by taking all columns of X to be transformed
C                        and store result into sequential columns
C                        of W. Determine also which row indices of X
C                        will need the left multiplication with the
C                        transposed of ROT to get the final 'diagonal'
C                        rotated elements.
C                    ii) Update the 'diagonal' rotation elements in
C                        scratch array Z and determine overall maximum
C                        offdiagonal element of X by searching all
C                        relevant rotated columns of X using the
C                        current W elements.
C                   iii) Replace all W into current row section of 
C                        X at all transformed columns of X.
C                continue
C
C                iv) Scatter all final 'diagonal' elements sitting
C                    in Z into appropriate places in X and search for
C                    overall maximum offdiagonal element of X using Z.
C
C                 v) Transpose all transformed columns of X to the
C                    corresponding rows in X (this is the most cache
C                    inefficient step!).
C
C
                DO 700 KBASE = 0,N-1,MXBLCK
                   KSUB = MIN0 (MXBLCK,N-KBASE)
C
C
C             ...step i): W = X * ROT. Arrays KP and LP will contain
C                the (local!) row indices of X on which ROT * W must
C                be performed.
C
C
                   M = 0
                   IX = KLQ (1)

                   DO 710 J = 1,NROT
                      JX = KLQ (J)
                      K = JX - KBASE
                      IF (K.GT.0 .AND. K.LE.KSUB) THEN
                          M = M + 1
                          KP (M) = K
                          LP (M) = J
                      END IF
                      S = ROT (1,J)
                      DO 712 K = 1,KSUB
                         W (K,J) = S * X (KBASE+K,IX)
  712                 CONTINUE
  710              CONTINUE

                   DO 714 I = 2,NROT
                      IX = KLQ (I)
                      DO 716 J = 1,NROT
                         S = ROT (I,J)
                         DO 718 K = 1,KSUB
                            W (K,J) = W (K,J) + S * X (KBASE+K,IX)
  718                    CONTINUE
  716                 CONTINUE
  714              CONTINUE
C
C
C             ...step ii): Z = ROT * W. Search for overall largest
C                offdiagonal X element (if needed) using elements of
C                W after the 'diagonal' positions in W have been set
C                to zero.
C
C
                   DO 720 J = 1,NROT
                      DO 722 I = J,NROT
                         DO 724 K = 1,M
                            KW = KP (K)
                            KR = LP (K)
                            Z (I,J) = Z (I,J) + ROT (KR,I) * W (KW,J)
  724                    CONTINUE
                         Z (J,I) = Z (I,J)
  722                 CONTINUE

                      JX = KLQ (J)

                      CASE1 = JX.GE.ZROW
                      CASE2 = JX.LE.ZCOL

                      IF (CASE1 .AND. CASE2) THEN
                          KBEG = 1
                          KEND = KSUB
                      ELSE IF (CASE1) THEN
                          KBEG = 1
                          KEND = MIN0 (KSUB,ZCOL-KBASE)
                      ELSE IF (CASE2) THEN
                          KBEG = MAX0 (1,ZROW-KBASE)
                          KEND = KSUB
                      END IF

                      IF (KBEG.LE.KEND) THEN
                          DO 726 K = 1,M
                             W (KP(K),J) = ZERO
  726                     CONTINUE
                          DO 728 K = KBEG,KEND
                             ABSXKL = DABS (W (K,J))
                             IF (ABSXKL.GT.XKLMAX) THEN
                                 KMAX = KBASE + K
                                 LMAX = JX
                                 XKLMAX = ABSXKL
                             END IF
  728                     CONTINUE
                      END IF
  720              CONTINUE
C
C
C             ...step iii): scatter W into appropriate columns of X.
C
C
                   DO 730 J = 1,NROT
                      JX = KLQ (J)
                      DO 732 K = 1,KSUB
                         X (KBASE+K,JX) = W (K,J)
  732                 CONTINUE
  730              CONTINUE

  700           CONTINUE
C
C
C             ...steps iv and v): Scatter 'diagonals' in Z into X and
C                search for largest offdiagonal X element (if needed)
C                using Z. Transpose all transformed columns of X to the
C                corresponding rows in X.
C
C
                DO 740 J = 1,NROT
                   JX = KLQ (J)
                   D (JX) = Z (J,J)
                   CASE1 = JX.LE.ZCOL

                   IF (CASE1) THEN
                       DO 742 I = 1,NROT
                          IX  = KLQ (I)
                          X (IX,JX) = Z (I,J)
                          IF (I.GT.J .AND. IX.GE.ZROW) THEN
                              ABSXKL = DABS (Z (I,J))
                              IF (ABSXKL.GT.XKLMAX) THEN
                                  KMAX = IX
                                  LMAX = JX
                                  XKLMAX = ABSXKL
                              END IF
                          END IF
  742                  CONTINUE
                   ELSE
                       DO 744 I = 1,NROT
                          IX  = KLQ (I)
                          X (IX,JX) = Z (I,J)
  744                  CONTINUE
                   END IF

                   DO 750 K = 1,N
                      X (JX,K) = X (K,JX)
  750              CONTINUE

  740           CONTINUE

            END IF
C
C
C             ...set next largest offdiagonal element of transformed
C                X found and proceed with the next iteration.
C
C
            K = KMAX
            KMAX = MIN0 (K,LMAX)
            LMAX = MAX0 (K,LMAX)

            XII = X (KMAX,KMAX)
            XIJ = X (LMAX,KMAX)
            XJJ = X (LMAX,LMAX)
            XIJMAX = XKLMAX

C            IF (N.EQ.46) THEN
C            CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +                   ( 1,
C     +                     ' C matrix ',
C     +                     DDROWC,DDCOLC,
C     +                     N,N,
C     +                     C )
C     +
C     +
C            CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
C     +
C     +                   ( 1,
C     +                     ' X matrix ',
C     +                     DDROWX,DDCOLX,
C     +                     N,0,
C     +                     X )
C     +
C     +
C            WRITE (1,*) ' KMAX,LMAX = ',KMAX,LMAX
C            WRITE (1,*) ' XKLMAX = ',XKLMAX
C            END IF

 1000    CONTINUE

         WRITE (*,*) ' Symmetric Jacobi procedure did not converge! '
         WRITE (*,*) ' Maximum # of iterations performed = ',MXITER
         WRITE (*,*) ' nlo__symmetric_jacobi '
         WRITE (1,*) ' Symmetric Jacobi procedure did not converge! '
         WRITE (1,*) ' Maximum # of iterations performed = ',MXITER
         WRITE (1,*) ' nlo__symmetric_jacobi '

         STOP
C
C
C             ...ready!
C
C
         RETURN
         END
