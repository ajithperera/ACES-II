      REAL*8 FUNCTION ssssPu(alpha,beta,gamma,delta,u,m)
      REAL*8 alpha,beta,gamma,delta,Pi,zeta,eta,u
      INTEGER m
      REAL*8 Binom,Qk
      EXTERNAL Binom, Qk
      Pi=DACOS(-1.0d0)
      zeta=alpha+beta
      eta=gamma+delta
      Up=(Pi/zeta)**1.5d0
      Uq=(Pi/eta)**1.5d0
      theta=DSQRT(zeta*eta/(zeta+eta))
      T=0.0d0
      SUMQ=0.0d0
      DO k=0,m
         SUMQ=SUMQ+Binom(m,k)*(-2.0d0)**k*(u*theta)**(2*k)*Qk(k)
      ENDDO
      ssssPu=Up*Uq*((2.0d0*u**2*theta*(2*theta**2)**(m+1))/DSQRT(Pi))
     &     *DEXP(-(u*theta)**2)*SUMQ
      RETURN
      END


      REAL*8 FUNCTION Qk(k)
      INTEGER k,fac
      fac=1
      DO i=1,2*k+1,2
         fac=fac*i
      ENDDO
      Qk=1.0d0/fac
      END
         
 
      REAL*8 FUNCTION Binom(n,k)
      INTEGER n,k
      REAL*8 Factorial
      EXTERNAL Factorial
      Binom=Factorial(n)/(Factorial(n-k)*Factorial(k))
      RETURN
      END
 

      REAL*8 FUNCTION Factorial(n)
      INTEGER n
      Factorial=1.0d0
      DO i=1,n
         Factorial=Factorial*DBLE(i)
      ENDDO
      RETURN
      END
      
