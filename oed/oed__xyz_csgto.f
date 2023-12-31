C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  OED__XYZ_CSGTO
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      L1CACHE,TILE,NCTROW,
     +                      SPHERIC,SCREEN,
     +                      ICORE,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      XC,YC,ZC,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_CSGTO
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__XYZ_SET_AB
C                OED__OVL_SET_IJ_PAIRS
C                OED__XYZ_E0_DEF_BLOCKS
C                OED__OVL_PREPARE_CTR
C                OED__OVL_E0_PCGTO_BLOCK
C                OED__CTR_2INDEX_BLOCK
C                OED__CTR_RS_EXPAND
C                OED__CTR_2INDEX_REORDER
C                OED__TRANSPOSE_BATCH
C                OED__XYZ_TO_RY_AB
C                OED__CARTESIAN_NORMS
C                OED__HRR_MATRIX
C                OED__HRR_TRANSFORM
C                OED__SPHERICAL_TRANSFORM
C                OED__NORMALIZE_CARTESIAN
C                OED__MOVE_RY
C  DESCRIPTION : This operation calculates a batch of contracted
C                overlap integrals on up to two different centers
C                between spherical or cartesian gaussian type shells.
C
C
C                  Input (x = 1 and 2):
C
C                    IMAX,ZMAX    =  maximum integer,flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2 and C
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for csh
C                                    1 and 2 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1 and 2 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    TILE         =  Number of rows and columns in
C                                    units of 8 Byte of level 1 cache
C                                    square tile array used for
C                                    performing optimum matrix
C                                    transpositions
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C                    SCREEN       =  is true, if screening will be
C                                    done at primitive integral level
C                    ICORE        =  integer scratch space
C                    MOMENTx      =  tells which power of X, Y, or Z
C                                    integrals to compute
C                    ZCORE (part) =  flp scratch space
C
C                  Output:
C
C                    NBATCH       =  # of integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    integral
C                    ZCORE        =  full batch of contracted (1|2)
C                                    overlap integrals over cartesian
C                                    or spherical gaussians starting
C                                    at ZCORE (NFIRST)
C
C
C
C              --- NOTES ABOUT THE OVERALL OVLP INTEGRAL PREFACTOR ---
C
C                The overal overlap integral prefactor is defined here
C                as follows. Consider the normalization factors for a
C                primitive cartesian GTO and for a spherical GTO
C                belonging to the angular momentum L = l+m+n:
C
C
C                    lmn                        l m n           2
C                 GTO  (x,y,z) = N (l,m,n,a) * x y z  * exp (-ar )
C                    a
C
C
C                    LM                     L    LM                2
C                 GTO  (r,t,p) = N (L,a) * r  * Y  (t,p) * exp (-ar )
C                    a
C
C
C                where a = alpha exponent, t = theta and p = phi and
C                N (l,m,n,a) and N (L,a) denote the respective
C                cartesian and spherical normalization factors such
C                that:
C
C
C                              lmn            lmn
C                 integral {GTO  (x,y,z) * GTO   (x,y,z) dx dy dz} = 1
C                              a              a
C
C
C                              LM             LM
C                 integral {GTO  (r,t,p) * GTO  (r,t,p) dr dt dp} = 1
C                              a              a
C
C
C                The normalization constants have then the following
C                values, assuming the spherical harmonics are
C                normalized:
C
C                              _____________________________________
C                             /      2^(2L+1+1/2) * a^((2L+3)/2)
C            N (l,m,n,a) =   / ----------------------------------------
C                          \/ (2l-1)!!(2m-1)!!(2n-1)!! * pi * sqrt (pi)
C
C
C                                   ____________________________
C                                  / 2^(2L+3+1/2) * a^((2L+3)/2)
C                     N (L,a) =   / -----------------------------
C                               \/     (2L+1)!! * sqrt (pi)
C
C
C                Note, that the extra pi under the square root in
C                N (l,m,n,a) belongs to the normalization of the
C                spherical harmonic functions and therefore does not
C                appear in the expression for N (L,a). The common
C                L-,l-,m-,n- and a-independent part of the cartesian
C                norm is a scalar quantity needed for all integrals
C                no matter what L-,l-,m-,n- and a-values they have:
C
C                                       _____________
C                                      /  2^(1+1/2)
C                     N (0,0,0,0) =   / --------------
C                                   \/  pi * sqrt (pi)
C
C
C                Also every overlap integral has a factor of pi**(3/2)
C                associated with it, hence the overall common factor
C                for all overlap integrals will be N(0,0,0,0)**2 times
C                pi**(3/2), which is equal to:
C
C                                          ___
C                            PREFACT  =  \/ 8 
C
C
C                and is set as a parameter inside the present routine.
C                The alpha exponent dependent part of the norms:
C
C                                   ______________
C                                 \/ a^((2L+3)/2)
C
C                will be calculated separately (see below) and their
C                inclusion in evaluating the primitive cartesian
C                overlap [E|0] integrals will be essential for numerical
C                stability during contraction.
C
C
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED package
C
C  MODIFIED    : Thomas Watson Jr.           p  q  r
C                  - Modified OED to handle X, Y, Z integrals
C              : Ajith Perere                p    q    r
C                  - Modified to handle (X-C)(Y-C)(Z-C)
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC
         LOGICAL     BLOCKED
         LOGICAL     EMPTY
         LOGICAL     EQUALAB
         LOGICAL     MEMORY
         LOGICAL     REORDER
         LOGICAL     SCREEN
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAPRS

         INTEGER     IHNROW,IHROW,IHSCR
         INTEGER     IMAX,ZMAX
         INTEGER     IN,OUT
         INTEGER     INDEXA,INDEXB
         INTEGER     INDEXR,INDEXS
         INTEGER     IPRIMA,IPRIMB
         INTEGER     IPUSED,IPSAVE,IPPAIR
         INTEGER     ISNROWA,ISNROWB
         INTEGER     ISROWA,ISROWB
         INTEGER     IUSED,ZUSED
         INTEGER     L1CACHE,TILE,NCTROW
         INTEGER     LCC1,LCC2
         INTEGER     LCCA,LCCB
         INTEGER     LCCSEGA,LCCSEGB
         INTEGER     LEXP1,LEXP2
         INTEGER     LEXPA,LEXPB
         INTEGER     MIJ
         INTEGER     MOVE,NOTMOVE
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NABCOOR
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB,NCGTOAB
         INTEGER     NCGTOR,NCGTOS
         INTEGER     NCOLHRR,NROWHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NIJ,NIJBLK,NIJBEG,NIJEND
         INTEGER     NINT1D
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB,NPGTOAB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NROTA,NROTB
         INTEGER     NROWA,NROWB
         INTEGER     NRYA,NRYB
         INTEGER     NXYZA,NXYZB
         INTEGER     NXYZET,NXYZP
         INTEGER     POS1,POS2
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB
         INTEGER     SHELLP
         INTEGER     TEMP
         INTEGER     MOMENTX,MOMENTY,MOMENTZ

         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,
     +               ZBASE,ZCNORM,
     +               ZRHOAB,ZXP,ZYP,ZZP,
     +               ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZINT1DX,ZINT1DY,ZINT1DZ,
     +               ZSCRMX,ZSCRMY,ZSCRMZ,
     +               ZXMINUS,ZYMINUS,ZZMINUS,
     +               ZSROTA,ZSROTB,
     +               ZHROT

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)
         INTEGER     IXOFF (1:2)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RNABSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,XC,YC,ZC
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

         PARAMETER  (PREFACT = 2.828427124746190D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...fix the A,B labels from the 1,2 ones. Calculate
C                the relevant data for the A,B batch of overlap
C                integrals.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LCC1  = 1
         LCC2  = LCC1 + NPGTO1 * NCGTO1

         CALL  OED__XYZ_SET_AB
     +
     +              ( NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                CC (LCC1),CC (LCC2),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,
     +                            NPGTOA,NPGTOB,
     +                            SHELLA,SHELLB,SHELLP,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,
     +                            ATOMIC,EQUALAB,
     +                            ABX,ABY,ABZ,NABCOOR,RNABSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZET,NXYZP,
     +                            NRYA,NRYB,
     +                            INDEXA,INDEXB,
     +                            SWAP12,SWAPRS,
     +                            LEXPA,LEXPB,
     +                            LCCA,LCCB,
     +                            LCCSEGA,LCCSEGB,
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
     +
C         WRITE (*,*) ' Finished set ab '
C         WRITE (*,*) ' Index A,B = ',INDEXA,INDEXB

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...enter the cartesian contracted (e|0) overlap batch
C                generation. Set the i and j primitive exponent
C                sets and the corresponding exponential prefactors.
C
C
         IF (EQUALAB) THEN
             NPGTOAB = (NPGTOA*(NPGTOA+1))/2
             NCGTOAB = (NCGTOA*(NCGTOA+1))/2
         ELSE
             NPGTOAB = NPGTOA * NPGTOB
             NCGTOAB = NCGTOA * NCGTOB
         END IF

         IPRIMA = 1
         IPRIMB = IPRIMA + NPGTOAB

         CALL  OED__OVL_SET_IJ_PAIRS
     +
     +              ( NPGTOA,NPGTOB,NPGTOAB,
     +                ATOMIC,EQUALAB,
     +                SWAPRS,
     +                RNABSQ,
     +                PREFACT,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                SCREEN,
     +
     +                         EMPTY,
     +                         NIJ,
     +                         ICORE (IPRIMA),ICORE (IPRIMB),
     +                         ZCORE (1) )
     +
     +
         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...decide on the primitive [e|0] block size and
C                return array sizes and pointers for the primitive
C                [e|0] generation. Perform also some preparation
C                steps for contraction.
C
C
         MEMORY = .FALSE.

         CALL  OED__XYZ_E0_DEF_BLOCKS
     +
     +              ( ZMAX,
     +                MOMENTX,MOMENTY,MOMENTZ,
     +                NPGTOA,NPGTOB,
     +                SHELLP,
     +                NIJ,NCGTOAB,
     +                NXYZET,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                        NIJBLK,
     +                        NPSIZE,NCSIZE,NWSIZE,
     +                        NINT1D,
     +                        MXPRIM,MNPRIM,
     +                        ZCBATCH,ZPBATCH,ZWORK,
     +                        ZNORMA,ZNORMB,
     +                        ZRHOAB,ZXP,ZYP,ZZP,
     +                        ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                        ZINT1DX,ZINT1DY,ZINT1DZ,
     +                        ZSCRMX,ZSCRMY,ZSCRMZ,
     +                        ZXMINUS,ZYMINUS,ZZMINUS )
     +
     +
         BLOCKED = NIJBLK .LT. NIJ

         CALL  OED__OVL_PREPARE_CTR
     +
     +              ( NCSIZE,
     +                NIJ,
     +                NPGTOA,NPGTOB,
     +                SHELLA,SHELLB,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                PREFACT,SPNORM,
     +                EQUALAB,
     +                BLOCKED,
     +                ZCORE (1),
     +
     +                          ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                          ZCORE (ZRHOAB),
     +                          ZCORE (ZCBATCH) )
     +
     +
         IPUSED = IPRIMB + NPGTOAB
         IPSAVE = IPUSED + MNPRIM
         IPPAIR = IPSAVE + MXPRIM
C
C
C             ...evaluate unnormalized rescaled [e|0] moment integrals
C                in blocks over ij pairs and add to final contracted
C                (e|0) overlap integrals. The keyword REORDER indicates,
C                if the primitive [e|0] moment integrals need to be
C                transposed before being contracted.
C
C
         REORDER = .TRUE.

         DO 1000 NIJBEG = 1,NIJ,NIJBLK
            NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
            MIJ = NIJEND - NIJBEG + 1

            CALL  OED__XYZ_E0_PCGTO_BLOCK
     +
     +                 ( NPSIZE,NINT1D,
     +                   ATOMIC,
     +                   MIJ,NIJ,NIJBEG,NIJEND,
     +                   NPGTOA,NPGTOB,
     +                   NXYZET,NXYZP,
     +                   SHELLA,SHELLP,
     +                   XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,
     +                   ABX,ABY,ABZ,
     +                   MOMENTX,MOMENTY,MOMENTZ,
     +                   ALPHA (LEXPA),ALPHA (LEXPB),
     +                   ICORE (IPRIMA+NIJBEG-1),
     +                   ICORE (IPRIMB+NIJBEG-1),
     +                   ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                   ZCORE (ZRHOAB),
     +                   ZCORE (ZXP),ZCORE (ZYP),ZCORE (ZZP),
     +                   ZCORE (ZPAX),ZCORE (ZPAY),ZCORE (ZPAZ),
     +                   ZCORE (ZPINVHF),ZCORE (ZSCALE),
     +                   ZCORE (ZINT1DX),ZCORE (ZSCRMX),ZCORE (ZXMINUS),
     +                   ZCORE (ZINT1DY),ZCORE (ZSCRMY),ZCORE (ZYMINUS),
     +                   ZCORE (ZINT1DZ),ZCORE (ZSCRMZ),ZCORE (ZZMINUS),
     +
     +                             ZCORE (ZPBATCH) )
     +
     +
C            WRITE (*,*) ' Finished e0 moment pcgto block '

            CALL  OED__CTR_2INDEX_BLOCK
     +
     +                 ( NPSIZE,NCSIZE,NWSIZE,
     +                   NXYZET,
     +                   MIJ,NCGTOAB,
     +                   NPGTOA,NPGTOB,
     +                   NCGTOA,NCGTOB,
     +                   MXPRIM,MNPRIM,
     +                   CC (LCCA),CC (LCCB),
     +                   CCBEG (LCCSEGA),CCBEG (LCCSEGB),
     +                   CCEND (LCCSEGA),CCEND (LCCSEGB),
     +                   ICORE (IPRIMA+NIJBEG-1),
     +                   ICORE (IPRIMB+NIJBEG-1),
     +                   L1CACHE,TILE,NCTROW,
     +                   EQUALAB,
     +                   SWAPRS,
     +                   REORDER,
     +                   BLOCKED,
     +                   ICORE (IPUSED),
     +                   ICORE (IPSAVE),
     +                   ICORE (IPPAIR),
     +                   ZCORE (ZPBATCH),
     +                   ZCORE (ZWORK),
     +
     +                             ZCORE (ZCBATCH) )
     +
     +
C            WRITE (*,*) ' Finished 2 index ctr block '

 1000    CONTINUE
C
C
C             ...the unnormalized cartesian (e|0) contracted moment
C                batch is ready. Expand the contraction indices
C                (if necessary):
C
C                   batch (nxyzet,r>=s) --> batch (nxyzet,r,s)
C
C                and reorder the contraction index part (if necessary):
C
C                   batch (nxyzet,r,s) --> batch (nxyzet,1,2)
C
C                The array IXOFF (x) indicates the total # of indices
C                to the left of x without including the nxyzet-part.
C                For the left most IXOFF value it is convenient to
C                set it equal to 1 instead of 0. Note, that the IXOFF
C                array indicates the true # of indices to the left
C                after! the batch has been transposed (see below) and
C                can be used as initial values when moving the
C                ry-components later on during the HRR and cartesian ->
C                spherical transformation procedure.
C                
C                For efficient application of the HRR contraction
C                scheme we need the batch elements ordered as:
C
C                            batch (1,2,nxyzet)
C
C                hence we transpose the batch after the reordering.
C
C                The space partitioning of the flp array for all of
C                these steps will be as follows:
C
C
C                          |  Zone 1  |  Zone 2  |
C
C                in which Zone 1 and 2 are 2 batches of HRR maximum
C                size. This can be done because we have always the
C                following dimension inequality:
C
C                              NXYZET =< NXYZHRR
C
C
         IXOFF (1) = 1
         IXOFF (2) = NCGTO1

         NCTR = NCGTO1 * NCGTO2
         MXSIZE = NCTR * NXYZHRR

         IN = ZCBATCH
         OUT = IN + MXSIZE

         IF (EQUALAB .AND. NCGTOAB.GT.1) THEN
             CALL  OED__CTR_RS_EXPAND
     +
     +                  ( NXYZET,NCGTOAB,
     +                    NCGTOA,NCGTOB,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr rs expansion '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         REORDER = SWAP12 .NEQV. SWAPRS

         IF (REORDER .AND. NCTR.GT.1) THEN

             IF (SWAPRS) THEN
                 INDEXR = INDEXB
                 INDEXS = INDEXA
                 NCGTOR = NCGTOB
                 NCGTOS = NCGTOA
             ELSE
                 INDEXR = INDEXA
                 INDEXS = INDEXB
                 NCGTOR = NCGTOA
                 NCGTOS = NCGTOB
             END IF

             CALL  OED__CTR_2INDEX_REORDER
     +
     +                  ( NXYZET,NCTR,
     +                    NCGTOR,NCGTOS,
     +                    IXOFF (INDEXR),IXOFF (INDEXS),
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr reorder '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         IF (NXYZET.GT.1 .AND. NCTR.GT.1) THEN

             CALL  OED__TRANSPOSE_BATCH
     +
     +                  ( NXYZET,NCTR,
     +                    TILE,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished transpose batch '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF
C
C
C             ...enter the HRR contraction and cartesian -> spherical
C                transformation / cartesian normalization section.
C                The sequence of events is to apply the HRR followed
C                by cartesian -> spherical transformations or cartesian
C                normalizations and immediate final positioning of
C                the finished parts to correspond with the contraction
C                indices i,j. The sequence of events is (where ' means
C                spherical or cartesian normalization and [] means the
C                indices are in correspondence):
C
C                       batch (ij,e0) --> batch (ij,ab)
C                       batch (ij,a,b) --> batch (ij,a,b')
C                       batch (ij,a,b') --> batch (ij[b'],a)
C                       batch (ij[b'],a) --> batch (ij[b'],a')
C                       batch (ij[b'],a') --> batch (ij[a'b'])
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
C
C
C                 Zone 1 and 2:  2 batches of MXSIZE maximum size
C                                (set previously)
C
C                       Zone 3:  cart -> spher transformation data
C                                             or
C                                cartesian normalization factors
C
C                       Zone 4:  HRR contraction data
C
C
C                Determine memory allocation offsets for the entire HRR
C                procedure + cartesian -> spherical transformations or
C                cartesian normalizations and generate the transformation
C                matrices + associated data for those shells > p-shell.
C                The offsets are as follows (x=A,B):
C
C                    IN = offset for input HRR batch
C                   OUT = offset for output HRR batch
C
C                ZSROTx = offset for x-part transformation matrix
C               ISNROWx = offset for # of non-zero XYZ contribution row
C                         labels for x-part transformation matrix
C                ISROWx = offset for non-zero XYZ contribution row
C                         labels for x-part transformation matrix
C
C                 ZHROT = offset for HRR transformation matrix
C                IHNROW = offset for # of nonzero row labels for
C                         each HRR matrix column
C                 IHROW = offset for nonzero row labels for the HRR
C                         matrix
C                 IHSCR = integer scratch space for HRR matrix
C                         assembly
C
C
C                In case of s- or p-shells no transformation matrix is
C                generated, hence if we have s- and/or p-shells, then
c                no call to the cartesian -> spherical transformation
C                or cartesian normalization routines needs to be done.
C                All integrals have already been multiplied by a factor
C                SPNORM, which has the following value for each s- and
C                p-shell:
C
C                       For s-type shell  =  1
C
C                       For p-type shell  =  2 * norm for s-type
C
C                This factor was introduced together with the overall
C                prefactor during evaluation of the primitive integrals
C                in order to save multiplications.
C
C
         ZBASE = MAX (IN,OUT) + MXSIZE

         IF (SPHERIC) THEN
             IF (MXSHELL.GT.1) THEN
                 CALL  OED__XYZ_TO_RY_AB
     +
     +                      ( NXYZA,NXYZB,
     +                        NRYA,NRYB,
     +                        SHELLA,SHELLB,
     +                        1,ZBASE,
     +
     +                                 NROWA,NROWB,
     +                                 NROTA,NROTB,
     +                                 ZSROTA,ZSROTB,
     +                                 ISNROWA,ISNROWB,
     +                                 ISROWA,ISROWB,
     +                                 IUSED,ZUSED,
     +                                 ICORE,ZCORE )
     +
     +
C                 WRITE (*,*) ' Finished xyz to ry ab '
             ELSE
                 IUSED = 0
                 ZUSED = 0
             END IF
         ELSE
C
C
             IF (MXSHELL.GT.1) THEN
                 ZCNORM = ZBASE
                 CALL  OED__CARTESIAN_NORMS
     +
     +                      ( MXSHELL,
     +
     +                                 ZCORE (ZCNORM))
     +
     +
C                 WRITE (*,*) ' Finished cartesian partial norms '
                 IUSED = 0
                 ZUSED = MXSHELL + 1
             ELSE
                 IUSED = 0
                 ZUSED = 0
             END IF
         END IF

         IHNROW = 1 + IUSED
         IHROW  = IHNROW + NCOLHRR + NCOLHRR
         IHSCR  = IHROW + NROTHRR + NROTHRR
         ZHROT  = ZBASE + ZUSED
C
C
C             ...do HRR step (if b > 0):
C
C                          (ij,e0) --> (ij,ab)
C
C
         IF (SHELLB.NE.0) THEN

             CALL  OED__HRR_MATRIX
     +
     +                  ( NROTHRR,NCOLHRR,
     +                    NXYZET,NXYZA,NXYZP,
     +                    SHELLA,SHELLB,SHELLP,
     +                    NABCOOR,
     +                    ABX,ABY,ABZ,
     +                    ICORE (IHSCR),
     +
     +                             POS1,POS2,
     +                             NROWHRR,
     +                             ICORE (IHNROW),
     +                             ICORE (IHROW),
     +                             ZCORE (ZHROT) )
     +
     +
C             WRITE (*,*) ' Finished HRR e0 matrix '

             CALL  OED__HRR_TRANSFORM
     +
     +                  ( NCTR,
     +                    NROWHRR,NXYZET,NXYZA*NXYZB,
     +                    NXYZA,NXYZB,
     +                    ICORE (IHNROW+POS1-1),
     +                    ICORE (IHROW+POS2-1),
     +                    ZCORE (ZHROT+POS2-1),
     +                    ZCORE (IN),
     +
     +                             ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished HRR e0 '
             TEMP = IN
             IN = OUT
             OUT = TEMP
C
C
C             ...do cart -> spher transformation / cart normalization
C                (if b > 1):
C
C                        (ij,a,b) --> (ij,a,b')
C
C
             IF (SHELLB.GT.1) THEN
                 IF (SPHERIC) THEN
                     CALL  OED__SPHERICAL_TRANSFORM
     +
     +                          ( NCTR*NXYZA,
     +                            NROWB,NXYZB,NRYB,
     +                            ICORE (ISNROWB),
     +                            ICORE (ISROWB),
     +                            ZCORE (ZSROTB),
     +                            ZCORE (IN),
     +
     +                                    ZCORE (OUT) )
     +
     +
C                     WRITE (*,*) ' Finished sph quart b '
                     TEMP = IN
                     IN = OUT
                     OUT = TEMP
                 ELSE
                     CALL  OED__NORMALIZE_CARTESIAN
     +
     +                          ( NCTR*NXYZA,
     +                            NXYZB,
     +                            SHELLB,
     +                            ZCORE (ZCNORM),
     +
     +                                    ZCORE (IN) )
     +
     +
C                     WRITE (*,*) ' Finished normalized cart b '
                 END IF
             END IF
         END IF
C
C
C             ...move transformed b shell (if size > 1):
C
C                        (ij,a,b') --> (ij[b'],a)
C
C
         IF (NRYB.GT.1) THEN
             NBATCH = NCTR * NXYZA * NRYB
             NOTMOVE = IXOFF (INDEXB)
             MOVE = NBATCH / (NOTMOVE * NRYB)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,2,
     +                        NOTMOVE,MOVE,NRYB,
     +                        INDEXB,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry b '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF
C
C
C             ...do cart -> spher / cart normalization (if a > 1):
C
C                        (ij[b'],a) --> (ij[b'],a')
C
C
         IF (SHELLA.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  OED__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NRYB,
     +                        NROWA,NXYZA,NRYA,
     +                        ICORE (ISNROWA),
     +                        ICORE (ISROWA),
     +                        ZCORE (ZSROTA),
     +                        ZCORE (IN),
     +
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished sph quart a '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             ELSE
                 CALL  OED__NORMALIZE_CARTESIAN
     +
     +                      ( NCTR*NRYB,
     +                        NXYZA,
     +                        SHELLA,
     +                        ZCORE (ZCNORM),
     +
     +                                ZCORE (IN) )
     +
     +
C                 WRITE (*,*) ' Finished normalized cart a '
             END IF
         END IF
C
C
C             ...move transformed a shell (if size > 1):
C
C                        (ij[b'],a') --> (ij[a'b'])
C
C
         NBATCH = NCTR * NRYB * NRYA

         IF (NRYA.GT.1) THEN
             NOTMOVE = IXOFF (INDEXA)
             MOVE = NBATCH / (NOTMOVE * NRYA)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,2,
     +                        NOTMOVE,MOVE,NRYA,
     +                        INDEXA,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry a '

                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF
C
C
C             ...set final pointer to integrals in ZCORE array.
C
C
         NFIRST = IN
C
C
C             ...ready!
C
C
         RETURN
         END
