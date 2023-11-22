      subroutine eigshf(f,n,nd,ns,alpha,beta)
      implicit double precision (a-h,o-z)
      dimension f(n,n)
      do 10 i=nd+1,nd+ns
      f(i,i) = f(i,i) + alpha
   10 continue
C
      do 20 i=nd+ns+1,n
      f(i,i) = f(i,i) + alpha + beta
   20 continue
      return
      end
