      SUBROUTINE CHKSGN(B, N)
C
C  CHECK THE SIGN OF EIGENVECTORS, SUCH THAT THE 'FIRST' LARGEST 
C  COMPONENT IS POSITIVE
C
      implicit none
      integer n, i, j, jmax
      double precision tol, b(n,n), amax,abs
c
      tol=1.d-5
      do 31 i=1,n
       jmax=1
       amax = abs(b(1,i)) + tol
       do 32 j=2,n
        if(abs(b(j,i)).gt.amax) then
          jmax = j
          amax = abs(b(j,i)) + tol
        endif
32     continue
       if(b(jmax,i).lt.0.d0)call vminus(b(1,i),n)
31    continue
c
      return
      end
