         SUBROUTINE  NLO__GENER_NLO_ORBITALS
     +
     +                    ( first,IMAX,ZMAX,
     +                      NBAS,NATOM,
     +                      MXSHELL,MXNAL,
     +                      MAXOCC,BONDSIZE,RYD2HYB,
     +                      NSHELLS,SHELLS,NBASAL,
     +                      ZATOM,
     +                      DENSITY,OVERLAP,
     +                      MJUMP,
     +                      SPHERIC,
     +                      ONLYNAO,ONLYNHO,ONLYNBO,ONLYNLMO,
     +                      MXCHOOSE,NCHOOSE,CHOOSE,
     +                      ICORE,ZCORE,
     +
     +                              COEFFS )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__GENER_NLO_ORBITALS
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : Master routine for finding natural localized orbitals
C                within a molecular environment using the first order
C                density matrix and the overlap matrix.
C
C                The density matrix can be of pure spin alpha or beta
C                in case of UHF functions or of mixed spin in case of
C                RHF functions. To determine the separate alpha and
C                beta natural localized orbitals simply run the routine
C                twice with the corresponding alpha and beta first order
C                density matrix.
C
C
C                  Input:
C
C                    IMAX,ZMAX    =  maximum integer + flp memory
C                    NBAS         =  total # of AO's in AO basis
C                    NATOM        =  total # of atomic centers
C                    MXSHELL      =  largest l-shell value
C                    MXNAL        =  maximum size of atomic l-shell
C                                    space. The atomic l-shell space
C                                    is the total # of contractions for
C                                    an atomic l-shell.
C                    MAXOCC       =  maximum orbital occupancy number
C                                    (can be only 1 or 2).
C                    BONDSIZE     =  maximum # of atomic centers that
C                                    are allowed to form a bond.
C                    RYD2HYB      =  is true, if the Rydberg NAO space
C                                    should be included next to the
C                                    Valence NAO space for natural 
C                                    hybrid orbital (NHO) construction.
C                                    If false, only the Valence NAO
C                                    space will be used.
C                    NSHELLS (A)  =  # of l-shells for atom A.
C                    SHELLS (I,A) =  I-th l-shell type (s=0,p=1,etc...)
C                                    for atom A.
C                    NBASAL (I,A) =  size of I-th atomic l-shell space
C                                    for atom A.
C                    ZATOM (A)    =  atomic number for atom A.
C                    DENSITY      =  full NBAS x NBAS matrix containing
C                                    the expansion coeffs of the density
C                                    in terms of products of the AO
C                                    basis functions.
C                    OVERLAP      =  full NBAS x NBAS overlap matrix.
C                    MJUMP        =  is .true., if the m values in the
C                                    m-space are ordered such that the
C                                    same m values are separated. This
C                                    keyword is necessary because some
C                                    AO basis functions are m-ordered
C                                    differently within each l-shell.
C                    SPHERIC      =  is true, if the l-shells are
C                                    spherical, false if they are
C                                    cartesian.
C                    ONLYx        =  is true, if only x=NAO,NHO,NBO or
C                                    NLMO orbitals are wanted. Note,
C                                    that if any of the ONLYx is true,
C                                    then all those ONLYx preceeding
C                                    that particular ONLYx are forced
C                                    to be false.
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
C                    ICORE,ZCORE  =  integer/flp scratch space
C
C
C                  Output:
C
C                    COEFFS       =  Natural localized orbitals
C                                    coefficient matrix in AO basis.
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

         LOGICAL     CYCLE
         LOGICAL     MJUMP
         LOGICAL     ONLYNAO,ONLYNHO,ONLYNBO,ONLYNLMO
         LOGICAL     RYD2HYB
         LOGICAL     SPHERIC

         INTEGER     BONDSIZE
         INTEGER     IAT2CEN,ICOLMAP,INRYDAL,ILSIZE,INAOTYP,ICSHELL,
     +               IVSHELL,IRSHELL,IHSHELL,IATNCB,IATCIDX,IATCOFF,
     +               IATNVB,IATVIDX,IATVOFF,IATNRB,IATRIDX,IATROFF,
     +               IATNHB,IATHVAL,IATHIDX,IATHOFF,IATHCEN,IATORD,
     +               ILOCMAP,IVALIDX,IRYDIDX,INHYB,IBDNCEN,IBDCEN,
     +               IBDBAS,IBDOCC,INBOBD,INHOLP,INHOEP,INHORYD,INWSTEP,
     +               IBASBEG,IBASEND,IWORK1,IWORK2
         INTEGER     IMAX,ZMAX
         INTEGER     MXCHOOSE
         INTEGER     MXNCBA,MXNVBA,MXNRBA,MXNHBA
         INTEGER     MXSHELL,MXNAL,MXNCEN,MX2CEN
         INTEGER     NBAS,NATOM
         INTEGER     NBOND,N2CEN
         INTEGER     NCA,NVA,NRA,NHA
         INTEGER     NCYCLE
         INTEGER     NLPAIR,NEPAIR
         INTEGER     NOCC
         INTEGER     NMB,NCB,NVB,NRB,NHB,NLB,NBB,NEB,NAB,NYB
         INTEGER     UNITID,UNITDB
         INTEGER     ZCAL,ZSAL,ZWAL,ZWRYDAT,ZWPRE,ZWNAO,ZWORK1,ZWORK2,
     +               ZPOPCOR,ZPOPVAL,ZPOPRYD,ZPOPSUM,ZH,ZSAH,ZPH,ZPHSUB,
     +               ZPHDEP,ZW,ZB,ZBOMAT,ZWBOND,ZWSTAR,ZLOCAL,ZDENORM

         INTEGER     ICORE   (1:IMAX)
         INTEGER     NCHOOSE (1:BONDSIZE)
         INTEGER     NSHELLS (1:NATOM)
         INTEGER     ZATOM   (1:NATOM)

         INTEGER     NBASAL  (1:MXSHELL+1,1:NATOM)
         INTEGER     SHELLS  (1:MXSHELL+1,1:NATOM)

         INTEGER     CHOOSE  (1:BONDSIZE,1:MXCHOOSE,1:BONDSIZE)

         DOUBLE PRECISION     MAXOCC
         DOUBLE PRECISION     NO2CEN
         DOUBLE PRECISION     WBDMIN,WSTMAX,WBDCRT,WSTEP
         DOUBLE PRECISION     WNAOVAL,WNAORYD,WNBOOCC

         DOUBLE PRECISION     ZCORE   (1:ZMAX)

         DOUBLE PRECISION     COEFFS  (1:NBAS,1:NBAS)
         DOUBLE PRECISION     DENSITY (1:NBAS,1:NBAS)
         DOUBLE PRECISION     OVERLAP (1:NBAS,1:NBAS)

         integer first
C
C
C------------------------------------------------------------------------
C
C
C             ...open general printout file and diagnostic/debug file.
C
C

         if (first.eq.1) then
            CALL  NLO__OPEN_FILES
     +
     +              ( 'NLO-results-printout-alpha',
     +                'NLO-diagnostic-debug-alpha',
     +                'UNKNOWN',
     +
     +                           UNITID,
     +                           UNITDB )
     +
     +
         end if
         if (first.eq.0) then
            CALL  NLO__OPEN_FILES
     +
     +              ( 'NLO-results-printout-beta',
     +                'NLO-diagnostic-debug-beta',
     +                'UNKNOWN',
     +
     +                           UNITID,
     +                           UNITDB )
     +
     +
         end if
C
C
C             ...preliminary steps.
C
C
         IF (ONLYNLMO) THEN
             ONLYNAO = .FALSE.
             ONLYNHO = .FALSE.
             ONLYNBO = .FALSE.
         END IF

         IF (ONLYNBO) THEN
             ONLYNAO = .FALSE.
             ONLYNHO = .FALSE.
         END IF

         IF (ONLYNHO) THEN
             ONLYNAO = .FALSE.
         END IF

         MX2CEN  = NATOM * (NATOM + 1) / 2

         IBASBEG = 1
         IBASEND = IBASBEG + NATOM
         IATORD  = IBASEND + NATOM
         IAT2CEN = IATORD  + NATOM
         ILSIZE  = IAT2CEN + 2*MX2CEN
         IWORK1  = ILSIZE  + (MXSHELL + 1)
         IWORK2  = IWORK1  + NBAS

         ZDENORM = 1
         ZBOMAT  = ZDENORM + NBAS
         ZWRYDAT = ZBOMAT  + NATOM * NATOM
         ZWORK1  = ZWRYDAT + NATOM
         ZWORK2  = ZWORK1  + NBAS

         CALL  NLO__INITIALIZE_RUN
     +
     +              ( NATOM,
     +                MXSHELL,
     +                MAXOCC,BONDSIZE,
     +                NSHELLS,SHELLS,NBASAL,
     +                SPHERIC,
     +                MXCHOOSE,NCHOOSE,CHOOSE,
     +                ICORE (IWORK1),
     +
     +                          NO2CEN,
     +                          ZCORE (ZWRYDAT),
     +                          WNAOVAL,WNAORYD,WNBOOCC,
     +                          WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                          ICORE (ILSIZE),
     +                          ICORE (IBASBEG),
     +                          ICORE (IBASEND) )
     +
     +
         CALL  NLO__NORMALIZE_OVERLAP_DENSITY
     +
     +              ( NBAS,
     +                ZCORE (ZWORK1),
     +                ICORE (IWORK1),
     +
     +                          ZCORE (ZDENORM),
     +                          DENSITY,OVERLAP )
     +
     +
         CALL  NLO__GENER_INITIAL_MATRICES
     +
     +              ( NBAS,NATOM,
     +                ICORE (IBASBEG),
     +                ICORE (IBASEND),
     +                OVERLAP,
     +                ZCORE (ZWORK2),
     +
     +                          ZCORE (ZBOMAT),
     +                          DENSITY )
     +
     +
         CALL  NLO__ANALYZE_BOND_ORDER_MATRIX
     +
     +              ( NATOM,
     +                MX2CEN,
     +                NO2CEN,
     +                ICORE (IWORK1),
     +                ZCORE (ZWORK1),
     +
     +                          N2CEN,
     +                          ICORE (IAT2CEN),
     +                          ICORE (IATORD),
     +                          ZCORE (ZBOMAT) )
     +
     +
C
C
C             ...the NAO generation section follows.
C
C
         ICOLMAP = ILSIZE  + (MXSHELL + 1)
         INRYDAL = ICOLMAP + NBAS
         INAOTYP = INRYDAL + NATOM * (MXSHELL + 1)
         ICSHELL = INAOTYP + NBAS
         IVSHELL = ICSHELL + NBAS
         IRSHELL = IVSHELL + NBAS
         IATNCB  = IRSHELL + NBAS
         IATCIDX = IATNCB  + NATOM
         IATCOFF = IATCIDX + NATOM
         IATNVB  = IATCOFF + NATOM
         IATVIDX = IATNVB  + NATOM
         IATVOFF = IATVIDX + NATOM
         IATNRB  = IATVOFF + NATOM
         IATRIDX = IATNRB  + NATOM
         IATROFF = IATRIDX + NATOM
         IWORK1  = IATROFF + NATOM

         ZWNAO   = ZWRYDAT + NATOM
         ZWPRE   = ZWNAO   + NBAS
         ZSAL    = ZWPRE   + NBAS
         ZCAL    = ZSAL    + MXNAL * MXNAL
         ZWAL    = ZCAL    + MXNAL * MXNAL
         ZPOPCOR = ZWAL    + MXNAL
         ZPOPVAL = ZPOPCOR + NATOM
         ZPOPRYD = ZPOPVAL + NATOM
         ZPOPSUM = ZPOPRYD + NATOM
         ZWORK1  = ZPOPSUM + NATOM
         ZWORK2  = ZWORK1  + NBAS + NBAS

         NCYCLE = 0

 1000    NCYCLE = NCYCLE + 1
         WRITE (*,*) ' pre-NAO/NAO cycle # ',NCYCLE

         CALL  NLO__GENER_PRE_NAO_ORBITALS
     +
     +              ( NBAS,NATOM,
     +                MXSHELL,MXNAL,
     +                NSHELLS,SHELLS,NBASAL,
     +                ICORE (ILSIZE),
     +                DENSITY,OVERLAP,
     +                ZCORE (ZSAL),ZCORE (ZCAL),ZCORE (ZWAL),
     +                ZCORE (ZWRYDAT),
     +                MJUMP,
     +
     +                         NMB,NRB,
     +                         ICORE (ICOLMAP),
     +                         ICORE (INRYDAL),
     +                         ZCORE (ZWPRE),
     +                         COEFFS )
     +
     +
         CALL  NLO__GENER_NAO_ORBITALS
     +
     +              ( NBAS,NATOM,
     +                MXSHELL,MXNAL,
     +                NSHELLS,SHELLS,NBASAL,
     +                DENSITY,OVERLAP,
     +                ZCORE (ZSAL),ZCORE (ZCAL),ZCORE (ZWAL),
     +                MJUMP,
     +                NMB,NRB,
     +                ICORE (ILSIZE),
     +                ICORE (ICOLMAP),
     +                ICORE (INRYDAL),
     +                ICORE (IWORK1),
     +                ZCORE (ZWORK1), ZCORE (ZWORK2),
     +                ZCORE (ZWPRE),
     +
     +                         ZCORE (ZWNAO),
     +                         COEFFS )
     +
     +
         CALL  NLO__ANALYZE_NAO_ORBITALS
     +
     +              ( NBAS,NATOM,
     +                MXSHELL,MXNAL,
     +                MAXOCC,
     +                NSHELLS,SHELLS,NBASAL,
     +                NMB,
     +                ICORE (ILSIZE),
     +                WNAOVAL,WNAORYD,
     +                ZCORE (ZWPRE), ZCORE (ZWNAO),
     +
     +                         NCB,NVB,NRB,
     +                         ZCORE (ZWRYDAT),
     +                         ZCORE (ZPOPCOR),
     +                         ZCORE (ZPOPVAL),
     +                         ZCORE (ZPOPRYD),
     +                         ZCORE (ZPOPSUM),
     +                         ICORE (INAOTYP),
     +                         ICORE (ICOLMAP),
     +                         ICORE (ICSHELL),
     +                         ICORE (IVSHELL),
     +                         ICORE (IRSHELL),
     +                         NCA,
     +                         MXNCBA,
     +                         ICORE (IATNCB),
     +                         ICORE (IATCIDX),
     +                         ICORE (IATCOFF),
     +                         NVA,
     +                         MXNVBA,
     +                         ICORE (IATNVB),
     +                         ICORE (IATVIDX),
     +                         ICORE (IATVOFF),
     +                         NRA,
     +                         MXNRBA,
     +                         ICORE (IATNRB),
     +                         ICORE (IATRIDX),
     +                         ICORE (IATROFF),
     +                         CYCLE )
     +
     +
         IF (CYCLE) GOTO 1000

         CALL  NLO__PRINT_NAO_RESULTS
     +
     +              ( UNITID,
     +                NBAS,NATOM,
     +                MXSHELL,
     +                NSHELLS,SHELLS,NBASAL,
     +                ZATOM,
     +                ICORE (ILSIZE),
     +                ZCORE (ZPOPCOR),
     +                ZCORE (ZPOPVAL),
     +                ZCORE (ZPOPRYD),
     +                ZCORE (ZPOPSUM),
     +                ICORE (INAOTYP),
     +                ZCORE (ZWNAO) )
     +
     +
         IF (ONLYNAO) THEN
             CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                  ( NBAS,
     +                    ZCORE (ZDENORM),
     +
     +                             COEFFS )
     +
     +
             RETURN
         END IF
C
C
C             ...enter the NHO generation.
C
C
         IHSHELL = IATROFF + NATOM
         IATNHB  = IHSHELL + NBAS
         IATHVAL = IATNHB  + NATOM
         IATHIDX = IATHVAL + NATOM
         IATHOFF = IATHIDX + NATOM
         ILOCMAP = IATHOFF + NATOM
         IVALIDX = ILOCMAP + NBAS
         IRYDIDX = IVALIDX + NATOM
         IWORK1  = IRYDIDX + NATOM

         CALL  NLO__CHARACTERIZE_NHO_SPACE
     +
     +              ( NBAS,NATOM,
     +                NVB,NVA,
     +                MXNVBA,
     +                ICORE (IVSHELL),
     +                ICORE (IATNVB),
     +                ICORE (IATVIDX),
     +                ICORE (IATVOFF),
     +                RYD2HYB,
     +                ICORE (ILOCMAP),
     +                ICORE (IVALIDX),
     +                ICORE (IRYDIDX),
     +                ICORE (IWORK1),
     +
     +                         ICORE (IATORD),
     +                         NHB,NCB,NRB,
     +                         ICORE (ICOLMAP),
     +                         ICORE (IHSHELL),
     +                         ICORE (ICSHELL),
     +                         ICORE (IRSHELL),
     +                         NHA,
     +                         MXNHBA,
     +                         ICORE (IATNHB),
     +                         ICORE (IATHVAL),
     +                         ICORE (IATHIDX),
     +                         ICORE (IATHOFF),
     +                         NCA,
     +                         MXNCBA,
     +                         ICORE (IATNCB),
     +                         ICORE (IATCIDX),
     +                         ICORE (IATCOFF),
     +                         NRA,
     +                         MXNRBA,
     +                         ICORE (IATNRB),
     +                         ICORE (IATRIDX),
     +                         ICORE (IATROFF) )
     +
     +
         IATHCEN = IATHOFF + NATOM
         INHYB   = IATHCEN + NATOM
         INHOLP  = INHYB   + NATOM
         INHOEP  = INHOLP  + NATOM
         INHORYD = INHOEP  + NATOM
         IBDNCEN = INHORYD + NATOM
         IBDOCC  = IBDNCEN + NHB
         IBDCEN  = IBDOCC  + NHB
         IBDBAS  = IBDCEN  + NHA * NHB
         INWSTEP = IBDBAS  + NHA * NHB
         IWORK1  = INWSTEP + BONDSIZE

         ZW      = ZWNAO
         ZH      = ZW + NBAS
         ZWBOND  = ZH + MXNHBA * NHB
         ZWSTAR  = ZWBOND + BONDSIZE
         ZSAH    = ZWSTAR + BONDSIZE
         ZPH     = ZSAH + MXNHBA * MXNHBA
         ZPHSUB  = ZPH + NHB * NHB
         ZPHDEP  = ZPHSUB + NHB * NHB
         ZWORK1  = ZPHDEP + NHB * NHB
         ZWORK2  = ZWORK1 + 3 * NBAS

         CALL  NLO__GENER_NHO_ORBITALS
     +
     +              ( NBAS,NATOM,N2CEN,
     +                MXNHBA,
     +                MAXOCC,BONDSIZE,
     +                NHB,NHA,
     +                ICORE (IATORD),
     +                ICORE (IAT2CEN),
     +                ICORE (IATHCEN),
     +                ICORE (IATNHB),
     +                ICORE (IATHVAL),
     +                ICORE (IATHOFF),
     +                RYD2HYB,
     +                ICORE (ICOLMAP),
     +                ICORE (INHYB),
     +                ICORE (INWSTEP),
     +                DENSITY,
     +                ZCORE (ZSAH),
     +                ZCORE (ZPH),
     +                ZCORE (ZPHSUB),
     +                ZCORE (ZPHDEP),
     +                ZCORE (ZBOMAT),
     +                WBDMIN,WSTMAX,WBDCRT,WSTEP,
     +                MXCHOOSE,NCHOOSE,CHOOSE,
     +                ICORE (IWORK1),
     +                ZCORE (ZWORK1),
     +                ZCORE (ZWORK2),
     +
     +                         NBOND,
     +                         NLPAIR,NEPAIR,
     +                         MXNCEN,
     +                         ZCORE (ZWBOND),
     +                         ZCORE (ZWSTAR),
     +                         ICORE (INHOLP),
     +                         ICORE (INHOEP),
     +                         ICORE (INHORYD),
     +                         ICORE (IBDNCEN),
     +                         ICORE (IBDCEN),
     +                         ICORE (IBDBAS),
     +                         ICORE (IBDOCC),
     +                         ZCORE (ZH),
     +                         ZCORE (ZW),
     +                         COEFFS )
     +
     +
         CALL  NLO__PRINT_NHO_RESULTS
     +
     +              ( UNITID,
     +                NATOM,NBOND,
     +                MXNHBA,
     +                NHB,NHA,
     +                ZATOM,
     +                ICORE (IATNHB),
     +                ICORE (IATHIDX),
     +                ICORE (IATHOFF),
     +                ICORE (IHSHELL),
     +                ICORE (IBDNCEN),
     +                ICORE (IBDCEN),
     +                ICORE (IBDBAS),
     +                ICORE (IBDOCC),
     +                ZCORE (ZH) )
     +
     +
         IF (ONLYNHO) THEN
             RETURN
         END IF
C
C
C             ...generate the square root of the overlap matrix
C                needed for NBO and NLMO locality analysis. The
C                original overlap matrix is overwritten.
C
C
         ZWORK1  = ZWSTAR + BONDSIZE
         ZWORK2  = ZWORK1 + NBAS

         CALL  NLO__GENER_SQROOT_OVERLAP
     +
     +              ( NBAS,
     +                ZCORE (ZWORK1),
     +                ZCORE (ZWORK2),
     +
     +                         OVERLAP )
     +
     +
C
C
C             ...enter the NBO generation.
C
C
         INBOBD  = IBDBAS  + NHA * NHB
         IWORK1  = INBOBD  + NHB

         ZLOCAL  = ZWSTAR + BONDSIZE
         ZB      = ZLOCAL + NATOM * NBAS
         ZPH     = ZB + MXNCEN * NHB
         ZWORK1  = ZPH + NHB * NHB
         ZWORK2  = ZWORK1 + NBAS

         CALL  NLO__GENER_NBO_ORBITALS
     +
     +              ( NBAS,NBOND,
     +                MXNHBA,MXNCEN,
     +                BONDSIZE,
     +                NHB,NCB,NHA,
     +                ICORE (IATNHB),
     +                ICORE (IATHOFF),
     +                ICORE (IBDNCEN),
     +                ICORE (IBDCEN),
     +                ICORE (IBDBAS),
     +                ICORE (IBDOCC),
     +                ICORE (ICOLMAP),
     +                ZCORE (ZWBOND),
     +                ZCORE (ZWSTAR),
     +                WNAORYD,WNBOOCC,
     +                DENSITY,
     +                ZCORE (ZH),
     +                ZCORE (ZPH),
     +                ICORE (IWORK1),
     +                ZCORE (ZWORK1),
     +                ZCORE (ZWORK2),
     +
     +                         NLB,NBB,NEB,NAB,NYB,
     +                         ICORE (INBOBD),
     +                         ZCORE (ZB),
     +                         ZCORE (ZW),
     +                         COEFFS )
     +
     +
         CALL  NLO__ANALYZE_NBO_ORBITALS
     +
     +              ( NBAS,NATOM,
     +                ICORE (IWORK1),
     +                ICORE (IBASBEG),
     +                ICORE (IBASEND),
     +                COEFFS,
     +                OVERLAP,
     +
     +                         ZCORE (ZLOCAL) )
     +
     +
         CALL  NLO__PRINT_NBO_RESULTS
     +
     +              ( UNITID,
     +                NBAS,NATOM,
     +                MXNCEN,
     +                NHB,NHA,
     +                NCA,NRA,
     +                NCB,NLB,NBB,NEB,NAB,NYB,NRB,
     +                ZATOM,
     +                ICORE (IATNCB),
     +                ICORE (IATCIDX),
     +                ICORE (IATCOFF),
     +                ICORE (IATNRB),
     +                ICORE (IATRIDX),
     +                ICORE (IATROFF),
     +                ICORE (ICSHELL),
     +                ICORE (IRSHELL),
     +                ICORE (IBDNCEN),
     +                ICORE (IBDCEN),
     +                ICORE (IBDBAS),
     +                ICORE (INBOBD),
     +                ZCORE (ZLOCAL),
     +                ZCORE (ZW),
     +                ZCORE (ZB) )
     +
     +
         IF (ONLYNBO) THEN
             CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +                  ( NBAS,
     +                    ZCORE (ZDENORM),
     +
     +                             COEFFS )
     +
     +
             RETURN
         END IF
C
C
C             ...enter the NLMO generation.
C
C
         IWORK1  = IBASEND + NATOM

         ZW      = ZDENORM + NBAS
         ZLOCAL  = ZW + NBAS
         ZWORK1  = ZLOCAL + NATOM * NBAS

         NOCC = NCB + NLB + NBB

         Write(6,*) "Generating NLMOs"
         CALL  NLO__GENER_NLMO_ORBITALS
     +
     +              ( NBAS,NOCC,
     +                ICORE (IWORK1),
     +                DENSITY,
     +                ZCORE (ZWORK1),
     +
     +                         ZCORE (ZW),
     +                         COEFFS )
     +
     +
         CALL  NLO__ANALYZE_NLMO_ORBITALS
     +
     +              ( NBAS,NATOM,
     +                ICORE (IWORK1),
     +                ICORE (IBASBEG),
     +                ICORE (IBASEND),
     +                COEFFS,
     +                OVERLAP,
     +
     +                         ZCORE (ZLOCAL) )
     +
     +
         CALL  NLO__PRINT_NLMO_RESULTS
     +
     +              ( UNITID,
     +                NBAS,NATOM,
     +                ZCORE (ZLOCAL) )
     +
     +
         CALL  NLO__DENORMALIZE_COEFF_MATRIX
     +
     +              ( NBAS,
     +                ZCORE (ZDENORM),
     +
     +                         COEFFS )
     +
     +
         CALL    MAT__PRINT_A_FLOAT_12_NOZEROS
     +
     +                ( 1,
     +                  'C(NLMO) AO bas de-norm',
     +                  NBAS,NBAS,
     +                  NBAS,NBAS,
     +                  COEFFS )
     +
     +
C
C
C             ...close general printout file and diagnostic/debug file.
C
C
         CALL  NLO__CLOSE_FILES
     +
     +              ( UNITID,UNITDB )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
