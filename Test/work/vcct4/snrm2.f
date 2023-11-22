      function snrm2(n,a,iskip)

c emulator of cray scilib snrm2 routine.























c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler




































































cYAU - ACES3 stuff . . . we hope - #include <aces.par>















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
