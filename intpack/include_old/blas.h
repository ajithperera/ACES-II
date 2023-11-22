
#ifndef _BLAS_H_
#define _BLAS_H_

c The following definitions define the blas routines to call.  It is done
c to distinguish between single and double precision routines.

c Instead of calling any blas routine by name, call it using the appropriate
c B_ macro below.  For the following lines:
c     call scopy(w,x,incx,y,incy)
c     call dcopy(w,x,incx,y,incy)
c should both be replaced with:
c     call B_COPY(w,x,incx,y,incy)

c All macros are of the form B_NAME where NAME is the name of the routine
c (minus the leading S or D).  The only exception is that the functions
c isamax and idamax have the macro B_AMAX.

#ifdef USE_SP_BLAS

c Blas 1
#  define B_ROTG   srotg
#  define B_ROTMG  srotmg
#  define B_ROT    srot
#  define B_ROTM   srotm
#  define B_SWAP   sswap
#  define B_SCAL   sscal
#  define B_COPY   scopy
#  define B_AXPY   saxpy
#  define B_DOT    sdot
#  define B_WRM2   swrm2
#  define B_ASUM   sasum
#  define B_AMAX   isamax

c Blas 2
#  define B_GEMV   sgemv
#  define B_GBMV   sgbmv
#  define B_SYMV   ssymv
#  define B_SBMV   ssbmv
#  define B_SPMV   sspmv
#  define B_TRMV   strmv
#  define B_TBMV   stbmv
#  define B_TPMV   stpmv
#  define B_TRSV   strsv
#  define B_TBSV   stbsv
#  define B_TPSV   stpsv
#  define B_GER    sger
#  define B_SYR    ssyr
#  define B_SPR    sspr
#  define B_SYR2   ssyr2
#  define B_SPR2   sspr2

c Blas 3
#  define B_GEMM   sgemm
#  define B_SYMM   ssymm
#  define B_SYRK   ssyrk
#  define B_SYR2K  ssyr2k
#  define B_TRMM   strmm
#  define B_TRSM   strsm

c Lapack linear equations
#  define B_GESV   sgesv
#  define B_GBSV   sgbsv
#  define B_GTSV   sgtsv
#  define B_POSV   sposv
#  define B_PPSV   sppsv
#  define B_PBSV   spbsv
#  define B_PTSV   sptsv
#  define B_SYSV   ssysv
#  define B_SPSV   sspsv

c Lapack least squares
#  define B_GELS   sgels
#  define B_GELSS  sgelss

c Lapack eigenvalue
#  define B_SYEV   ssyev
#  define B_SPEV   sspev
#  define B_SBEV   ssbev
#  define B_STEV   sstev
#  define B_GEES   sgees
#  define B_GEEV   sgeev
#  define B_GESVD  sgesvd

c Lapack generalized eigenvalue
#  define B_SYGV   ssygv
#  define B_SPGV   sspgv

c Lapack linear equations (expert)
#  define B_GESVX   sgesvx
#  define B_GBSVX   sgbsvx
#  define B_GTSVX   sgtsvx
#  define B_POSVX   sposvx
#  define B_PPSVX   sppsvx
#  define B_PBSVX   spbsvx
#  define B_PTSVX   sptsvx
#  define B_SYSVX   ssysvx
#  define B_SPSVX   sspsvx

c Lapack least squares (expert)
#  define B_GELSX   sgelsx

c Lapack eigenvalue (expert)
#  define B_SYEVX   ssyevx
#  define B_SPEVX   sspevx
#  define B_SBEVX   ssbevx
#  define B_STEVX   sstevx
#  define B_GEESX   sgeesx
#  define B_GEEVX   sgeevx

c Linpack (minv is stupid)
#  define B_GECO    sgeco
#  define B_GEDI    sgedi
#  define B_GESL    sgesl

#else

c Blas 1
#  define B_ROTG   drotg
#  define B_ROTMG  drotmg
#  define B_ROT    drot
#  define B_ROTM   drotm
#  define B_SWAP   dswap
#  define B_SCAL   dscal
#  define B_COPY   dcopy
#  define B_AXPY   daxpy
#  define B_DOT    ddot
#  define B_WRM2   dwrm2
#  define B_ASUM   dasum
#  define B_AMAX   idamax

c Blas 2
#  define B_GEMV   dgemv
#  define B_GBMV   dgbmv
#  define B_SYMV   dsymv
#  define B_SBMV   dsbmv
#  define B_SPMV   dspmv
#  define B_TRMV   dtrmv
#  define B_TBMV   dtbmv
#  define B_TPMV   dtpmv
#  define B_TRSV   dtrsv
#  define B_TBSV   dtbsv
#  define B_TPSV   dtpsv
#  define B_GER    dger
#  define B_SYR    dsyr
#  define B_SPR    dspr
#  define B_SYR2   dsyr2
#  define B_SPR2   dspr2

c Blas 3
#  define B_GEMM   dgemm
#  define B_SYMM   dsymm
#  define B_SYRK   dsyrk
#  define B_SYR2K  dsyr2k
#  define B_TRMM   dtrmm
#  define B_TRSM   dtrsm

c Lapack linear equations
#  define B_GESV   dgesv
#  define B_GBSV   dgbsv
#  define B_GTSV   dgtsv
#  define B_POSV   dposv
#  define B_PPSV   dppsv
#  define B_PBSV   dpbsv
#  define B_PTSV   dptsv
#  define B_SYSV   dsysv
#  define B_SPSV   dspsv

c Lapack least squares
#  define B_GELS   dgels
#  define B_GELSS  dgelss

c Lapack eigenvalue
#  define B_SYEV   dsyev
#  define B_SPEV   dspev
#  define B_SBEV   dsbev
#  define B_STEV   dstev
#  define B_GEES   dgees
#  define B_GEEV   dgeev
#  define B_GESVD  dgesvd

c Lapack generalized eigenvalue
#  define B_SYGV   dsygv
#  define B_SPGV   dspgv

c Lapack linear equations (expert)
#  define B_GESVX   dgesvx
#  define B_GBSVX   dgbsvx
#  define B_GTSVX   dgtsvx
#  define B_POSVX   dposvx
#  define B_PPSVX   dppsvx
#  define B_PBSVX   dpbsvx
#  define B_PTSVX   dptsvx
#  define B_SYSVX   dsysvx
#  define B_SPSVX   dspsvx

c Lapack least squares (expert)
#  define B_GELSX   dgelsx

c Lapack eigenvalue (expert)
#  define B_SYEVX   dsyevx
#  define B_SPEVX   dspevx
#  define B_SBEVX   dsbevx
#  define B_STEVX   dstevx
#  define B_GEESX   dgeesx
#  define B_GEEVX   dgeevx

c Linpack (minv is stupid)
#  define B_GECO    dgeco
#  define B_GEDI    dgedi
#  define B_GESL    dgesl

#endif

#endif /* _BLAS_H_ */

