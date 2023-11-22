cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      subroutine gridgausschebyshev (nr,rmin,rmax,c1,c2,c3,r,w)
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This routine calculates the Gauss-Chebyshev quadrature points and
c weights of for the interval RMIN < R < RMAX.
c
c In:
c      nr     - number of radial points
c      rmin   - interval minium 
c      rmax   - interval maximum
c Scr:
c      c1 (nr), c2 (nr), c3 (nr)
c Out:
c      r      - grid point coordinates
c      w      - grid point weights
c 
      implicit none
c
      integer nr
c
      double precision
     &     rmin, rmax, c1 (nr), c2 (nr), c3 (nr), r (nr), w (nr)
c
c     Local variables
c
      integer i
      double precision slope, intercept, ri
c
c     Construct the grid.
c     
      slope=(rmax-rmin)*.5d0
      intercept=(rmin+rmax)*.5d0
c
      call gridgausslegendre(nr, c1, c2, c3)
      do i=1,nr
         ri = slope * c2 (nr-i+1) + intercept
         r (i) = ri
         w (i) = c3 (nr-i+1)*ri*ri*slope
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      subroutine gridgausslegendre (nn,c,x,a)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This routine calculates the Gauss-Legendre points and weights for the
c Legendre polynomials using subroutines provided in Gaussian Quadrature
c Formulas by A. H. Stroud and D. Secrest. These routines are for Jacobi
c polynomials specialized to the case of alf=bta=0. The routine works
c easily up to order 200.
c
c     nn      [in]   order of the integration
c     c (nn)  [scr] 
c     x (nn)  [out]  integration points
c     a (nn)  [out]  integration weights for the integral int(-1,1) f(x) dx.
c
      implicit none
      integer nn
      double precision c(*),x(*),a(*)
c
      integer n
      double precision eps,csx,csa,tsx,tsa
c
      eps=2.D-15
      do 10 n=1,nn
       c(n)=(n-1.D0)**2/((2.D0*n-1.D0)*(2.D0*n-3.D0))
   10 continue
      call grid_legendre(nn,x,a,c,eps,csx,csa,tsx,tsa)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      subroutine grid_legendre(nn,x,a,c,eps,csx,csa,tsx,tsa)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nn
      double precision c(nn),x(nn),a(nn),csx,csa,tsx,tsa,eps
      double precision fn,beta,cc,xt,dpn,pn1
      integer j,i
C
      fn=nn
      csx=0.D0
      csa=0.D0
      beta=1.D0
      cc=2.D0
      tsx=0.D0
      tsa=cc
      do 1 j=2,nn
    1  cc=cc*c(j)
      do 12 i=1,nn
       if (i-1)  12,2,3
c largest zero
    2  xt=1.D0-2.78D0/(4.D0+fn*fn)
       go to 11
    3  if (i-2)  12,4,5
c second zero
    4  xt=xt-4.1D0*(1.D0+.06D0*(fn-8.D0)/fn)*(1.D0-xt)
       go to 11
    5  if (i-3)  12,6,7
c third zero
    6  xt=xt-1.67D0*(1.D0+.22D0*(fn-8.D0)/fn)*(x(1)-xt)
       go to 11
    7  if (nn-i-1)  10,9,8
c middle zeros
    8  xt=3.D0*x(i-1)-3.D0*x(i-2)+x(i-3)
       go to 11
c second last zero
    9  xt=xt+(xt-x(i-2))/
     &       (0.766D0*(1.D0+.639D0*(fn-4.D0)/(1.D0+.71D0*(fn-4.D0))))
       go to 11
c last zero
   10  xt=xt+(xt-x(i-2))/
     &       (1.67D0*(1.D0+.22D0*(fn-8.D0)/fn))
c
   11  call grid_root(xt,nn,dpn,pn1,c,eps)
       x(i)=xt
       a(i)=cc/(dpn*pn1)
       csx=csx+x(i)
       csa=csa+a(i)
   12 continue
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      subroutine grid_root(x,nn,dpn,pn1,c,eps)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Improves the approximate root "x".  In addition the following quantities
c are obtained:
c
c     dpn = derivative of P(n) at x
c     pn1 = value of P(n-1) at x
c
      implicit none
      integer nn
      double precision c(nn)
      double precision pn,dpn,pn1,eps,x
c
      integer i
      double precision d
c
      do 10 i=1,10
       call grid_recur(pn,dpn,pn1,x,nn,c)
       d=pn/dpn
       x=x-d
       if(abs(d).le.eps) return
   10 continue
      write(6,*)'WARNING:  FAILED TO CONVERGE A ROOT OF THE '//
     &                    'LEGENDRE POLYNOMIAL TO THE DESIRED'
      write(6,*)'          PRECISION (SEE GAUSSLEG.F).'
      write(6,*)'            DESIRED PRECISION:  ',eps
      write(6,*)'            OBTAINED PRECISION:  ',abs(d)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      subroutine grid_recur(pn,dpn,pn1,x,nn,c)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nn
      double precision c(nn),pn,dpn,pn1,x
      double precision p1,p,q,dq,dp1,dp
      integer j
c
      p1=1.D0
      p=x
      dp1=0.D0
      dp=1.D0
      do 1 j=2,nn
       q=x*p-c(j)*p1
       dq=x*dp+p-c(j)*dp1
       p1=p
       p=q
       dp1=dp
       dp=dq
    1 continue
      pn=p
      dpn=dp
      pn1=p1
      end
