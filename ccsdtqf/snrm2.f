      function snrm2(n,a,iskip)

      integer n,iskip
      double precision
     &    snrm2,a(1), sqrt

      integer ntop,i

      ntop=1+(n-1)*iskip
      snrm2=0
      do 10 i=1,ntop,iskip
        snrm2=snrm2+a(i)*a(i)
   10 continue
      snrm2=sqrt(snrm2)
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
