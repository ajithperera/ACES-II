      subroutine mnoz(a,b,c,n1,n2,n3)
      implicit double precision (a-h,o-z)
      dimension a(n1,n3),b(n3,n2),c(n1,n2)
      data zero/0.0d+0/
      call zeroma(c,1,n1*n2)
      do 10 i=1,n1
         do 20 j=1,n2
            xx=zero
            do 30 k=1,n3
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
 30         continue
 20      continue
 10   continue
      return
      end
