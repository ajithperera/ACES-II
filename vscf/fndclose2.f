      subroutine fndclose2(length,v1,test1,v2, test2, ilocmin, 
     $    tol)
c
c Tries to match element in two vectors to be 'closest' to test1 and test2
c error if difference in either is greater than tol: ilocmin = -1
c
      implicit none
c
      integer length, ilocmin, i, ierror, iter
      double precision
     &    v1(length), test1, vlocmin, xclose, tmp,
     $    v2(length), test2, tol
c
c first locate element closest to test1
c
      ierror = 0
      iter = 0
 100  iter = iter + 1
      xclose=1.d+30
      do 10 i=1,length
        tmp=abs(v1(i)-test1) 
        if(tmp.lt.xclose)then
          ilocmin=i
          xclose=tmp
        endif 
 10   continue
c
      if (xclose .gt. 0.05) then
        write(6,*) ' something wrong in fndclose2'
        ierror = 1
      endif
c
c check if element is also close to test2
c
      if (abs(v2(ilocmin) - test2) .gt. tol) then
c
c search again, excluding current value
c
        v1(ilocmin) = 1.0d30
        if (iter .lt. 5) then
          goto 100
        else
          ierror = 1
        endif
      endif
c
      if (ierror .ne. 0) ilocmin = -1
c
      return
      end
c
c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
