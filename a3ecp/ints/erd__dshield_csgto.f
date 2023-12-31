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
         SUBROUTINE  ERD__DSHIELD_CSGTO
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      L1CACHE,TILE,NCTROW,
     +                      SPHERIC,SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DSHIELD_CSGTO
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__SET_ABCD
C                ERD__SET_IJ_KL_PAIRS
C                ERD__E0F0_DEF_DSHIELD_BLOCKS
C                ERD__PREPARE_CTR
C                ERD__E0F0_PCGTO_DSHIELD_BLOCK
C                ERD__CTR_4INDEX_BLOCK
C                ERD__CTR_RS_EXPAND
C                ERD__CTR_TU_EXPAND
C                ERD__CTR_4INDEX_REORDER
C                ERD__TRANSPOSE_BATCH
C                ERD__XYZ_TO_RY_ABCD
C                ERD__CARTESIAN_NORMS
C                ERD__HRR_MATRIX
C                ERD__HRR_TRANSFORM
C                ERD__SPHERICAL_TRANSFORM
C                ERD__NORMALIZE_CARTESIAN
C                ERD__MOVE_RY
C  DESCRIPTION : This operation calculates a batch of contracted
C                electron repulsion integrals on up to four different
C                centers between spherical or cartesian gaussian type
C                shells.
C
C
C                  Input (x = 1,2,3 and 4):
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
C                                    y = 1,2,3 and 4
C                    XX,XY,XZ     =  identify which tensor component is
C                    YX,YY,YZ     =  is computed. The one that is being
C                    ZX,ZY,ZZ     =  computed carries 1 and the others
C                                    are zero
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for csh
C                                    1,2,3,4 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1,2,3,4 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
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
C                    ZCORE (part) =  flp scratch space
C
C                  Output:
C
C                    NBATCH       =  # of integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    integral
C                    ZCORE        =  full batch of contracted (12|34)
C                                    integrals over cartesian or
C                                    spherical gaussians starting at
C                                    ZCORE (NFIRST)
C
C
C
C              --- NOTES ABOUT THE OVERALL INTEGRAL PREFACTOR ---
C
C                The overal integral prefactor is defined here as
C                follows. Consider the normalization factors for a
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
C                Also every ERI integral has a factor of 2*pi**(5/2)
C                associated with it, hence the overall common factor
C                for all integrals will be N(0,0,0,0)**4 times
C                2*pi**(5/2), which is equal to:
C
C
C                            PREFACT = 16 / sqrt(pi)
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
C                [E0|F0] integrals will be essential for numerical
C                stability during contraction.
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

         LOGICAL     ATOMIC,ATOMAB,ATOMCD
         LOGICAL     BLOCKED
         LOGICAL     EMPTY
         LOGICAL     EQUALAB,EQUALCD
         LOGICAL     MEMORY
         LOGICAL     REORDER
         LOGICAL     SCREEN
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAP34,SWAPRS,SWAPTU
         LOGICAL     TR1234

         INTEGER     IHNROW,IHROW,IHSCR
         INTEGER     IMAX,ZMAX
         INTEGER     IN,OUT
         INTEGER     INDEXA,INDEXB,INDEXC,INDEXD
         INTEGER     INDEXR,INDEXS,INDEXT,INDEXU
         INTEGER     IPRIMA,IPRIMB,IPRIMC,IPRIMD
         INTEGER     IPUSED,IPSAVE,IPPAIR
         INTEGER     ISNROWA,ISNROWB,ISNROWC,ISNROWD
         INTEGER     ISROWA,ISROWB,ISROWC,ISROWD
         INTEGER     IUSED,ZUSED
         INTEGER     L1CACHE,TILE,NCTROW
         INTEGER     LCC1,LCC2,LCC3,LCC4
         INTEGER     LCCA,LCCB,LCCC,LCCD
         INTEGER     LCCSEGA,LCCSEGB,LCCSEGC,LCCSEGD
         INTEGER     LEXP1,LEXP2,LEXP3,LEXP4
         INTEGER     LEXPA,LEXPB,LEXPC,LEXPD
         INTEGER     MIJ,MKL,MIJKL,MGQIJKL
         INTEGER     MGRID,NGRID
         INTEGER     MOVE,NOTMOVE
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NABCOOR,NCDCOOR
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NCGTOA,NCGTOB,NCGTOC,NCGTOD,NCGTOAB,NCGTOCD
         INTEGER     NCGTOR,NCGTOS,NCGTOT,NCGTOU
         INTEGER     NCOLHRR,NROWHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NGQP,NMOM,NGQSCR
         INTEGER     NIJ,NKL
         INTEGER     NIJBLK,NKLBLK,NIJBEG,NKLBEG,NIJEND,NKLEND
         INTEGER     NINT2D
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD,NPGTOAB,NPGTOCD
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NROTA,NROTB,NROTC,NROTD
         INTEGER     NROWA,NROWB,NROWC,NROWD
         INTEGER     NRYA,NRYB,NRYC,NRYD
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER     NXYZP,NXYZQ,NXYZT,NXYZET,NXYZFT
         INTEGER     POS1,POS2
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     SHELLP,SHELLQ,SHELLT
         INTEGER     TEMP

         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +               ZBASE,ZCNORM,
     +               ZRHOAB,ZRHOCD,
     +               ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +               ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +               ZB00,ZB01,ZB10,
     +               ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +               ZINT2DX,ZINT2DY,ZINT2DZ,
     +               ZINT2DX1,ZINT2DY1,ZINT2DZ1,
     +               ZSROTA,ZSROTB,ZSROTC,ZSROTD,
     +               ZHROT

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)
         INTEGER     IXOFF (1:4)
         INTEGER     XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ

         DOUBLE PRECISION  ABX,ABY,ABZ,CDX,CDY,CDZ
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RNABSQ,RNCDSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

         DOUBLE PRECISION  FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (PREFACT = 9.027033336764101D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...fix the A,B,C,D labels from the 1,2,3,4 ones.
C                Calculate the relevant data for the A,B,C,D batch of
C                integrals.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2
         LEXP4 = LEXP3 + NPGTO3

         LCC1 = 1
         LCC2 = LCC1 + NPGTO1 * NCGTO1
         LCC3 = LCC2 + NPGTO2 * NCGTO2
         LCC4 = LCC3 + NPGTO3 * NCGTO3

         CALL  ERD__SET_ABCD
     +
     +              ( NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                SHELL1,SHELL2,SHELL3,SHELL4,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                ALPHA (LEXP3),ALPHA (LEXP4),
     +                CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,NCGTOD,
     +                            NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                            SHELLA,SHELLB,SHELLC,SHELLD,
     +                            SHELLP,SHELLQ,SHELLT,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                            ATOMIC,ATOMAB,ATOMCD,
     +                            EQUALAB,EQUALCD,
     +                            ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                            NABCOOR,NCDCOOR,
     +                            RNABSQ,RNCDSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZC,NXYZD,
     +                            NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                            NRYA,NRYB,NRYC,NRYD,
     +                            INDEXA,INDEXB,INDEXC,INDEXD,
     +                            SWAP12,SWAP34,SWAPRS,SWAPTU,TR1234,
     +                            LEXPA,LEXPB,LEXPC,LEXPD,
     +                            LCCA,LCCB,LCCC,LCCD,
     +                            LCCSEGA,LCCSEGB,LCCSEGC,LCCSEGD,
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
     +
C         WRITE (*,*) ' Finished set abcd '
C         WRITE (*,*) ' Index A,B,C,D = ',INDEXA,INDEXB,INDEXC,INDEXD
C         WRITE (*,*) ' EQUALAB = ',EQUALAB
C         WRITE (*,*) ' EQUALCD = ',EQUALCD

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...enter the cartesian contracted (e0|f0) batch
C                generation. Set the ij and kl primitive exponent
C                pairs and the corresponding exponential prefactors.
C
C
         IF (EQUALAB) THEN
             NPGTOAB = (NPGTOA*(NPGTOA+1))/2
             NCGTOAB = (NCGTOA*(NCGTOA+1))/2
         ELSE
             NPGTOAB = NPGTOA * NPGTOB
             NCGTOAB = NCGTOA * NCGTOB
         END IF

         IF (EQUALCD) THEN
             NPGTOCD = (NPGTOC*(NPGTOC+1))/2
             NCGTOCD = (NCGTOC*(NCGTOC+1))/2
         ELSE
             NPGTOCD = NPGTOC * NPGTOD
             NCGTOCD = NCGTOC * NCGTOD
         END IF

         NCTR = NCGTOAB * NCGTOCD
         NXYZT = NXYZET * NXYZFT

         IPRIMA = 1
         IPRIMB = IPRIMA + NPGTOAB
         IPRIMC = IPRIMB + NPGTOAB
         IPRIMD = IPRIMC + NPGTOCD

         Write(*,*) "ENTER ERD__SET_IJ_KL_PAIRS"
         CALL  ERD__SET_IJ_KL_PAIRS
     +
     +              ( NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                NPGTOAB,NPGTOCD,
     +                ATOMAB,ATOMCD,
     +                EQUALAB,EQUALCD,
     +                SWAPRS,SWAPTU,
     +                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                RNABSQ,RNCDSQ,
     +                PREFACT,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                ALPHA (LEXPC),ALPHA (LEXPD),
     +                FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                SCREEN,
     +
     +                         EMPTY,
     +                         NIJ,NKL,
     +                         ICORE (IPRIMA),ICORE (IPRIMB),
     +                         ICORE (IPRIMC),ICORE (IPRIMD),
     +                         ZCORE (1) )
     +
     +
         Write(*,*) "END ERD__SET_IJ_KL_PAIRS"

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...decide on the primitive [e0|f0] block size and
C                return array sizes and pointers for the primitive
C                [e0|f0] generation. Perform also some preparation
C                steps for contraction.
C
C
         NGQP = 2 + SHELLT / 2
         NMOM = 2 * NGQP - 1
         NGQSCR = 5 * NMOM + 2 * NGQP - 2

         MEMORY = .FALSE.

         Write(*,*) "The number of grid pts", NGQP
         Write(*,*) "ENTER ERD__E0F0_DEF_DSHIELD_BLOCKS"
         CALL  ERD__E0F0_DEF_DSHIELD_BLOCKS
     +
     +              ( ZMAX,
     +                NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                SHELLP,SHELLQ,
     +                NIJ,NKL,
     +                NCGTOAB,NCGTOCD,NCTR,
     +                NGQP,NGQSCR,
     +                NXYZT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                     NIJBLK,NKLBLK,
     +                     NPSIZE,NCSIZE,NWSIZE,
     +                     NINT2D,
     +                     MXPRIM,MNPRIM,
     +                     ZCBATCH,ZPBATCH,ZWORK,
     +                     ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +                     ZRHOAB,ZRHOCD,
     +                     ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +                     ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +                     ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +                     ZB00,ZB01,ZB10,
     +                     ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +                     ZINT2DX,ZINT2DY,ZINT2DZ,
     +                     ZINT2DX1,ZINT2DY1,ZINT2DZ1 )
     +
     +
         BLOCKED = (NIJBLK.LT.NIJ) .OR. (NKLBLK.LT.NKL)

         Write(*,*) "END ERD__E0F0_DEF_DSHIELD_BLOCKS"
         CALL  ERD__PREPARE_CTR
     +
     +              ( NCSIZE,
     +                NIJ,NKL,
     +                NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                SHELLA,SHELLB,SHELLC,SHELLD,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                ALPHA (LEXPC),ALPHA (LEXPD),
     +                PREFACT,SPNORM,
     +                EQUALAB,EQUALCD,
     +                BLOCKED,
     +                ZCORE (1),
     +
     +                          ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                          ZCORE (ZNORMC),ZCORE (ZNORMD),
     +                          ZCORE (ZRHOAB),ZCORE (ZRHOCD),
     +                          ZCORE (ZCBATCH) )
     +
     +
         IPUSED = IPRIMD + NPGTOCD
         IPSAVE = IPUSED + MNPRIM
         IPPAIR = IPSAVE + MXPRIM
C
C
C             ...evaluate unnormalized rescaled [e0|f0] in blocks
C                over ij and kl pairs and add to final contracted
C                (e0|f0). The keyword REORDER indicates, if the
C                primitive [e0|f0] blocks need to be transposed
C                before being contracted.
C
C
         REORDER = .TRUE.

         DO 1000 NIJBEG = 1,NIJ,NIJBLK
            NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
            MIJ = NIJEND - NIJBEG + 1
            DO 1100 NKLBEG = 1,NKL,NKLBLK
               NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
               MKL = NKLEND - NKLBEG + 1

               MIJKL = MIJ * MKL
               MGQIJKL = NGQP * MIJKL

               Write(*,*) "ENTER ERD__E0F0_PCGTO_DSHIELD_BLOCKS"
               CALL  ERD__E0F0_PCGTO_DSHIELD_BLOCK
     +
     +                    ( NPSIZE,NINT2D,
     +                      ATOMIC,ATOMAB,ATOMCD,
     +                      MIJ,MKL,MIJKL,
     +                      NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                      NGQP,NMOM,NGQSCR,MGQIJKL,
     +                      NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                      SHELLA,SHELLP,SHELLC,SHELLQ,
     +                      XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                      ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                      XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ,
     +                      ALPHA (LEXPA),ALPHA (LEXPB),
     +                      ALPHA (LEXPC),ALPHA (LEXPD),
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      ICORE (IPRIMA+NIJBEG-1),
     +                      ICORE (IPRIMB+NIJBEG-1),
     +                      ICORE (IPRIMC+NKLBEG-1),
     +                      ICORE (IPRIMD+NKLBEG-1),
     +                      ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                      ZCORE (ZNORMC),ZCORE (ZNORMD),
     +                      ZCORE (ZRHOAB),ZCORE (ZRHOCD),
     +                      ZCORE (ZP),
     +                      ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                      ZCORE (ZPAX),ZCORE (ZPAY),ZCORE (ZPAZ),
     +                      ZCORE (ZPINVHF),ZCORE (ZSCPK2),
     +                      ZCORE (ZQ),
     +                      ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                      ZCORE (ZQCX),ZCORE (ZQCY),ZCORE (ZQCZ),
     +                      ZCORE (ZQINVHF),ZCORE (ZSCQK2),
     +                      ZCORE (ZRTS),ZCORE (ZWTS),ZCORE (ZGQSCR),
     +                      ZCORE (ZTVAL),ZCORE (ZPQPINV),
     +                      ZCORE (ZSCPQK4),
     +                      ZCORE (ZB00),ZCORE (ZB01),ZCORE (ZB10),
     +                      ZCORE (ZC00X),ZCORE (ZC00Y),ZCORE (ZC00Z),
     +                      ZCORE (ZD00X),ZCORE (ZD00Y),ZCORE (ZD00Z),
     +                      ZCORE (ZINT2DX),
     +                      ZCORE (ZINT2DY),
     +                      ZCORE (ZINT2DZ),
     +                      ZCORE (ZINT2DX1),
     +                      ZCORE (ZINT2DY1),
     +                      ZCORE (ZINT2DZ1),
     +
     +                                ZCORE (ZPBATCH) )
     +
     +
               WRITE (*,*) ' Finished e0f0 pcgto block '

               CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                    ( NPSIZE,NCSIZE,NWSIZE,
     +                      NXYZT,MIJKL,
     +                      MIJ,MKL,NCGTOAB,NCGTOCD,
     +                      NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      NCGTOA,NCGTOB,NCGTOC,NCGTOD,
     +                      MXPRIM,MNPRIM,
     +                      CC (LCCA),CC (LCCB),CC (LCCC),CC (LCCD),
     +                      CCBEG (LCCSEGA),CCBEG (LCCSEGB),
     +                      CCBEG (LCCSEGC),CCBEG (LCCSEGD),
     +                      CCEND (LCCSEGA),CCEND (LCCSEGB),
     +                      CCEND (LCCSEGC),CCEND (LCCSEGD),
     +                      ICORE (IPRIMA+NIJBEG-1),
     +                      ICORE (IPRIMB+NIJBEG-1),
     +                      ICORE (IPRIMC+NKLBEG-1),
     +                      ICORE (IPRIMD+NKLBEG-1),
     +                      L1CACHE,TILE,NCTROW,
     +                      EQUALAB,EQUALCD,
     +                      SWAPRS,SWAPTU,
     +                      REORDER,
     +                      BLOCKED,
     +                      ICORE (IPUSED),
     +                      ICORE (IPSAVE),
     +                      ICORE (IPPAIR),
     +                      ZCORE (ZPBATCH),
     +                      ZCORE (ZWORK),
     +
     +                                ZCORE (ZCBATCH) )
     +
     +
C               WRITE (*,*) ' Finished 4 index ctr block '

 1100       CONTINUE
 1000    CONTINUE
C
C
C             ...the unnormalized cartesian (e0|f0) contracted batch is
C                ready. Expand the contraction indices (if necessary):
C
C                   batch (nxyzt,r>=s,t>=u) --> batch (nxyzt,r,s,t,u)
C
C                and reorder the contraction index part (if necessary):
C
C                   batch (nxyzt,r,s,t,u) --> batch (nxyzt,1,2,3,4)
C
C                The array IXOFF (x) indicates the total # of indices
C                to the left of x without including the nxyzt-part.
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
C                          batch (1,2,3,4,nxyzt)
C
C                hence we transpose the batch after the reordering.
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
C                              NXYZT =< NXYZHRR
C
C
         IXOFF (1) = 1
         IXOFF (2) = NCGTO1
         IXOFF (3) = NCGTO1 * NCGTO2
         IXOFF (4) = IXOFF (3) * NCGTO3

         NCTR = IXOFF (4) * NCGTO4
         MXSIZE = NCTR * NXYZHRR

         IN = ZCBATCH
         OUT = IN + MXSIZE

         IF (EQUALAB .AND. NCGTOAB.GT.1) THEN
             CALL  ERD__CTR_RS_EXPAND
     +
     +                  ( NXYZT,NCGTOAB,NCGTOCD,
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

         IF (EQUALCD .AND. NCGTOCD.GT.1) THEN
             CALL  ERD__CTR_TU_EXPAND
     +
     +                  ( NXYZT*NCGTOA*NCGTOB,NCGTOCD,
     +                    NCGTOC,NCGTOD,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr tu expansion '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         REORDER = TR1234 .OR. (SWAP12.NEQV.SWAPRS)
     +                    .OR. (SWAP34.NEQV.SWAPTU)

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

             IF (SWAPTU) THEN
                 INDEXT = INDEXD
                 INDEXU = INDEXC
                 NCGTOT = NCGTOD
                 NCGTOU = NCGTOC
             ELSE
                 INDEXT = INDEXC
                 INDEXU = INDEXD
                 NCGTOT = NCGTOC
                 NCGTOU = NCGTOD
             END IF

             CALL  ERD__CTR_4INDEX_REORDER
     +
     +                  ( NXYZT,NCTR,
     +                    NCGTOR,NCGTOS,NCGTOT,NCGTOU,
     +                    IXOFF (INDEXR),IXOFF (INDEXS),
     +                    IXOFF (INDEXT),IXOFF (INDEXU),
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

         IF (NXYZT.GT.1 .AND. NCTR.GT.1) THEN

             CALL  ERD__TRANSPOSE_BATCH
     +
     +                  ( NXYZT,NCTR,
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

C         CALL  ERD__PRINT_BATCH
C     +
C     +              ( NCTR,NXYZT,1,1,  ZCORE (IN) )
C     +
C         STOP
C
C
C             ...enter the HRR contraction and cartesian -> spherical
C                transformation / cartesian normalization section.
C                The sequence of events is to apply the HRR followed
C                by cartesian -> spherical transformations or cartesian
C                normalizations and immediate final positioning of
C                the finished parts to correspond with the contraction
C                indices i,j,k,l. First we do the f0-part followed
C                by the e0-part, which hence gives rise to the sequence
C                (where ' means spherical or cartesian normalization
C                and [] means the indices are in correspondence):
C
C                   batch (ijkl,e0,f0) --> batch (ijkl,e0,cd)
C                   batch (ijkl,e0,c,d) --> batch (ijkl,e0,c,d')
C                   batch (ijkl,e0,c,d') --> batch (ijkl[d'],e0,c)
C                   batch (ijkl[d'],e0,c) --> batch (ijkl[d'],e0,c')
C                   batch (ijkl[d'],e0,c') --> batch (ijkl[c'd'],e0)
C
C                   batch (ijkl[c'd'],e0) --> batch (ijkl[c'd'],ab)
C                   batch (ijkl[c'd'],a,b) --> batch (ijkl[c'd'],a,b')
C                   batch (ijkl[c'd'],a,b') --> batch (ijkl[b'c'd'],a)
C                   batch (ijkl[b'c'd'],a) --> batch (ijkl[b'c'd'],a')
C                   batch (ijkl[b'c'd'],a') --> batch (ijkl[a'b'c'd'])
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                  |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
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
C                The offsets are as follows (x=A,B,C,D):
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
         ZBASE = MAX0 (IN,OUT) + MXSIZE

         IF (SPHERIC) THEN
             IF (MXSHELL.GT.1) THEN
                 CALL  ERD__XYZ_TO_RY_ABCD
     +
     +                      ( NXYZA,NXYZB,NXYZC,NXYZD,
     +                        NRYA,NRYB,NRYC,NRYD,
     +                        SHELLA,SHELLB,SHELLC,SHELLD,
     +                        1,ZBASE,
     +
     +                                 NROWA,NROWB,NROWC,NROWD,
     +                                 NROTA,NROTB,NROTC,NROTD,
     +                                 ZSROTA,ZSROTB,ZSROTC,ZSROTD,
     +                                 ISNROWA,ISNROWB,ISNROWC,ISNROWD,
     +                                 ISROWA,ISROWB,ISROWC,ISROWD,
     +                                 IUSED,ZUSED,
     +                                 ICORE,ZCORE )
     +
     +
C                 WRITE (*,*) ' Finished xyz to ry abcd '
             ELSE
                 IUSED = 0
                 ZUSED = 0
             END IF
         ELSE
             IF (MXSHELL.GT.1) THEN
                 ZCNORM = ZBASE
                 CALL  ERD__CARTESIAN_NORMS
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
C             ...do the first stage of processing the integrals:
C
C                   batch (ijkl,e0,f0) --> batch (ijkl,e0,cd)
C                   batch (ijkl,e0,c,d) --> batch (ijkl,e0,c,d')
C                   batch (ijkl,e0,c,d') --> batch (ijkl[d'],e0,c)
C                   batch (ijkl[d'],e0,c) --> batch (ijkl[d'],e0,c')
C                   batch (ijkl[d'],e0,c') --> batch (ijkl[c'd'],e0)
C
C
         IF (SHELLD.NE.0) THEN

             CALL  ERD__HRR_MATRIX
     +
     +                  ( NROTHRR,NCOLHRR,
     +                    NXYZFT,NXYZC,NXYZQ,
     +                    SHELLC,SHELLD,SHELLQ,
     +                    NCDCOOR,
     +                    CDX,CDY,CDZ,
     +                    ICORE (IHSCR),
     +
     +                             POS1,POS2,
     +                             NROWHRR,
     +                             ICORE (IHNROW),
     +                             ICORE (IHROW),
     +                             ZCORE (ZHROT) )
     +
     +
C             WRITE (*,*) ' Finished HRR f0 matrix '

             CALL  ERD__HRR_TRANSFORM
     +
     +                  ( NCTR*NXYZET,
     +                    NROWHRR,NXYZFT,NXYZC*NXYZD,
     +                    NXYZC,NXYZD,
     +                    ICORE (IHNROW+POS1-1),
     +                    ICORE (IHROW+POS2-1),
     +                    ZCORE (ZHROT+POS2-1),
     +                    ZCORE (IN),
     +
     +                             ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished HRR f0 '

             TEMP = IN
             IN = OUT
             OUT = TEMP

             IF (SHELLD.GT.1) THEN
                 IF (SPHERIC) THEN
                     CALL  ERD__SPHERICAL_TRANSFORM
     +
     +                          ( NCTR*NXYZET*NXYZC,
     +                            NROWD,NXYZD,NRYD,
     +                            ICORE (ISNROWD),
     +                            ICORE (ISROWD),
     +                            ZCORE (ZSROTD),
     +                            ZCORE (IN),
     +
     +                                    ZCORE (OUT) )
     +
     +
C                     WRITE (*,*) ' Finished sph quart d '
                     TEMP = IN
                     IN = OUT
                     OUT = TEMP
                 ELSE
                     CALL  ERD__NORMALIZE_CARTESIAN
     +
     +                          ( NCTR*NXYZET*NXYZC,
     +                            NXYZD,
     +                            SHELLD,
     +                            ZCORE (ZCNORM),
     +
     +                                    ZCORE (IN) )
     +
     +
C                     WRITE (*,*) ' Finished normalize cart d '
                 END IF
             END IF
         END IF

         IF (NRYD.GT.1) THEN
             NBATCH = NCTR * NXYZET * NXYZC * NRYD
             NOTMOVE = IXOFF (INDEXD)
             MOVE = NBATCH / (NOTMOVE * NRYD)
             IF (MOVE.GT.1) THEN

                 CALL  ERD__MOVE_RY
     +
     +                      ( NBATCH,4,
     +                        NOTMOVE,MOVE,NRYD,
     +                        INDEXD,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry d '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF

         IF (SHELLC.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  ERD__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NXYZET*NRYD,
     +                        NROWC,NXYZC,NRYC,
     +                        ICORE (ISNROWC),
     +                        ICORE (ISROWC),
     +                        ZCORE (ZSROTC),
     +                        ZCORE (IN),
     +
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished sph quart c '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             ELSE
                 CALL  ERD__NORMALIZE_CARTESIAN
     +
     +                      ( NCTR*NXYZET*NRYD,
     +                        NXYZC,
     +                        SHELLC,
     +                        ZCORE (ZCNORM),
     +
     +                                ZCORE (IN) )
     +
     +
C                 WRITE (*,*) ' Finished normalize cart c '
             END IF
         END IF

         IF (NRYC.GT.1) THEN
             NBATCH = NCTR * NRYD * NXYZET * NRYC
             NOTMOVE = IXOFF (INDEXC)
             MOVE = NBATCH / (NOTMOVE * NRYC)
             IF (MOVE.GT.1) THEN

                 CALL  ERD__MOVE_RY
     +
     +                      ( NBATCH,4,
     +                        NOTMOVE,MOVE,NRYC,
     +                        INDEXC,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry c '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF
C
C
C             ...do the second stage of processing the integrals:
C
C                   batch (ijkl[c'd'],e0) --> batch (ijkl[c'd'],ab)
C                   batch (ijkl[c'd'],a,b) --> batch (ijkl[c'd'],a,b')
C                   batch (ijkl[c'd'],a,b') --> batch (ijkl[b'c'd'],a)
C                   batch (ijkl[b'c'd'],a) --> batch (ijkl[b'c'd'],a')
C                   batch (ijkl[b'c'd'],a') --> batch (ijkl[a'b'c'd'])
C
C
         IF (SHELLB.NE.0) THEN

             CALL  ERD__HRR_MATRIX
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

             CALL  ERD__HRR_TRANSFORM
     +
     +                  ( NCTR*NRYC*NRYD,
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

             IF (SHELLB.GT.1) THEN
                 IF (SPHERIC) THEN
                     CALL  ERD__SPHERICAL_TRANSFORM
     +
     +                          ( NCTR*NRYC*NRYD*NXYZA,
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
                     CALL  ERD__NORMALIZE_CARTESIAN
     +
     +                          ( NCTR*NRYC*NRYD*NXYZA,
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

         IF (NRYB.GT.1) THEN
             NBATCH = NCTR * NRYC * NRYD * NXYZA * NRYB
             NOTMOVE = IXOFF (INDEXB)
             MOVE = NBATCH / (NOTMOVE * NRYB)
             IF (MOVE.GT.1) THEN

                 CALL  ERD__MOVE_RY
     +
     +                      ( NBATCH,4,
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

         IF (SHELLA.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  ERD__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NRYB*NRYC*NRYD,
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
                 CALL  ERD__NORMALIZE_CARTESIAN
     +
     +                      ( NCTR*NRYB*NRYC*NRYD,
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

         NBATCH = NCTR * NRYB * NRYC * NRYD * NRYA

         IF (NRYA.GT.1) THEN
             NOTMOVE = IXOFF (INDEXA)
             MOVE = NBATCH / (NOTMOVE * NRYA)
             IF (MOVE.GT.1) THEN

                 CALL  ERD__MOVE_RY
     +
     +                      ( NBATCH,4,
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
