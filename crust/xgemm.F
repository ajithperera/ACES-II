
cYAU - WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING ! WARNING
c
c    This is the original xgemm, which filtered lda, ldb, and ldc through
c max(1,ld?). This subtly changes the behavior of dgemm depending on the
c other arguments. After wrapping 20 or so xgemm calls (checking for
c proper dimensions), I decided just to let their code go to hell if it wants.
c Have fun.

      subroutine xgemm(transa,transb,
     &                 m,n,k,alpha,a,lda,
     &                             b,ldb,
     &                       beta, c,ldc)

#ifdef USE_SP_BLAS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* USE_SP_BLAS */

      character*1 transa, transb
      integer m, n, k, lda, ldb, ldc, max
      double precision alpha, beta
c      double precision a(lda,*), b(ldb,*), c(ldc,*)
      double precision a, b, c

c ----------------------------------------------------------------------

      if (m.eq.0) return
      if (n.eq.0) return

c   o If one of the leading dimensions takes an illegal value, then the
c     only way to recover, at this point, is if (k.eq.0). dgemm will
c     scale c by beta if k is 0 regardless of lda or ldb. If ldc is
c     also zero, then the programmer really needs to be shot.

      if ((lda*ldb).eq.0) then
         if (k.eq.0) call B_GEMM(transa,transb,m,n,k,
     &                           alpha,a,max(1,lda),
     &                                 b,max(1,ldb),
     &                           beta, c,ldc)
         return
      end if


c ----------------------------------------------------------------------

      call B_GEMM(transa,transb,m,n,k,
     &            alpha,a,lda,
     &                  b,ldb,
     &            beta, c,ldc)

      return
      end

