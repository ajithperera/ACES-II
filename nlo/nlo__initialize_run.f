         SUBROUTINE  NLO__INITIALIZE_RUN
     +
     +                    ( NATOM,
     +                      MXSHELL,
     +                      MAXOCC,BONDSIZE,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      SPHERIC,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      INDEX,
     +
     +                            NO2CEN,
     +                            WRYDAT,
     +                            WNAOVAL,WNAORYD,WNBOOCC,
     +                            WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                            LSIZE,
     +                            BASBEG,BASEND )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__INITIALIZE_RUN
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine initializes the present run. It generates
C                all data, which is needed several times later on.
C                Also transmitted data is checked for consistency and
C                the program is aborted, if any meaningless data is
C                found.
C
C                  Input:
C
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A.
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    SPHERIC      =  is true, if the l-shells are
C                                    spherical, false if they are
C                                    cartesian.
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
C                    INDEX        =  will hold info about used atoms
C
C                  Output:
C
C                    NO2CEN       =  2-center bond formation criterion
C                                    for analysis of the atomic bond
C                                    order matrix.
C                    WRYDAT       =  initial Rydberg weight criterion
C                                    for pre-NAO construction.
C                    WNAOVAL      =  the Valence NAO weight threshold.
C                    WNAORYD      =  the Rydberg NAO weight threshold.
C                    WNBOOCC      =  the weight limit to decide when
C                                    a NBO is considered occupied or
C                                    virtual.
C                    WBDMIN       =  initial minimum accepted weight
C                                    for bond construction.
C                    WSTMAX       =  initial maximum accepted weight
C                                    for antibond construction.
C                    WBDCRT       =  critical weight limit for bond
C                                    construction.
C                    WSTEP        =  weight stepping size to decrease/
C                                    increase the bond/antibond weight
C                                    limits.
C                    LSIZE (I)    =  I-th l-shell size
C                    BASBEG (A)   =  first basis index number for
C                                    atom A.
C                    BASEND (A)   =  last basis index number for
C                                    atom A.
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

         LOGICAL     SPHERIC

         INTEGER     A,B,L,N
         INTEGER     ATOM
         INTEGER     BASNR
         INTEGER     BONDSIZE
         INTEGER     MXCHOOSE
         INTEGER     MXSHELL
         INTEGER     NATOM
         INTEGER     NSHELL

         INTEGER     BASBEG  (1:NATOM)
         INTEGER     BASEND  (1:NATOM)
         INTEGER     INDEX   (1:NATOM)
         INTEGER     LSIZE   (0:MXSHELL)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NSHELLS (1:NATOM)

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION  MAXOCC
         DOUBLE PRECISION  NO2CEN
         DOUBLE PRECISION  WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION  WNAORYD,WNAOVAL,WNBOOCC

         DOUBLE PRECISION  WRYDAT (1:NATOM)
C
C
C------------------------------------------------------------------------
C
C
C             ...set the limits.
C
C
         NO2CEN = 0.025D0 * MAXOCC

         DO 10 N = 1,NATOM
            WRYDAT (N) = 0.25D0 * MAXOCC
   10    CONTINUE

         WNAORYD = 0.025D0 * MAXOCC
         WNAOVAL = 0.99D0  * MAXOCC
         WNBOOCC = 0.500D0 * MAXOCC
CC HUGHES        WBDMIN  = 0.95D0  * MAXOCC
         WBDMIN  = 0.975D0  * MAXOCC
         WSTMAX  = 0.05D0  * MAXOCC
         WBDCRT  = 0.75D0  * MAXOCC
         WSTEP   = 0.01D0  * MAXOCC
C
C
C             ...predetermine l-shell sizes needed.
C
C
         LSIZE (0) = 1
         IF (SPHERIC) THEN
             DO 20 L = 1,MXSHELL
                LSIZE (L) = LSIZE (L-1) + 2
   20        CONTINUE
         ELSE
             DO 30 L = 1,MXSHELL
                LSIZE (L) = LSIZE (L-1) + L + 1
   30        CONTINUE
         END IF
C
C
C             ...determine atomic basis indices for further use.
C
C
         BASNR = 0
         DO 40 A = 1,NATOM
            BASBEG (A) = BASNR + 1
            NSHELL = NSHELLS (A)
            DO 50 N = 1,NSHELL
               BASNR = BASNR + NBASAL (N,A) * LSIZE (SHELLS (N,A))
   50       CONTINUE
            BASEND (A) = BASNR
   40    CONTINUE
C
C
C             ...check the selected atomic indices for bond formation
C                (if any). Stop, if any inconsistencies are found.
C
C
         DO 100 B = 1,BONDSIZE

            N = IABS (NCHOOSE (B))

            IF (N.EQ.9999) GOTO 100

            IF (N.GT.MXCHOOSE) THEN
                WRITE (*,*) ' Dimensions for CHOOSE array too small! '
                WRITE (*,*) ' B,N,MXCHOOSE = ',B,N,MXCHOOSE
                WRITE (*,*) ' nlo__initialize_run '
                WRITE (1,*) ' Dimensions for CHOOSE array too small! '
                WRITE (1,*) ' B,N,MXCHOOSE = ',B,N,MXCHOOSE
                WRITE (1,*) ' nlo__initialize_run '
                STOP
            ELSE IF (N.GT.0) THEN

                DO 110 L = 1,N
                   DO 120 ATOM = 1,NATOM
                      INDEX (ATOM) = 0
  120              CONTINUE
                   DO 130 A = 1,B
                      ATOM = CHOOSE (A,L,B)
                      INDEX (ATOM) = INDEX (ATOM) + 1
  130              CONTINUE
                   DO 140 ATOM = 1,NATOM
                      IF (INDEX (ATOM).GT.1) THEN
                          WRITE (*,*) ' Atom idx used more than once! '
                          WRITE (*,*) ' ATOM,L,B = ',ATOM,L,B
                          WRITE (*,*) ' nlo__initialize_run '
                          WRITE (1,*) ' Atom idx used more than once! '
                          WRITE (1,*) ' ATOM,L,B = ',ATOM,L,B
                          WRITE (1,*) ' nlo__initialize_run '
                          STOP
                      END IF
  140              CONTINUE
  110           CONTINUE

            END IF
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
