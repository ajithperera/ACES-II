c-----------------------------------------------------------------------
      subroutine x_rks_pbe(r,g,f,dfdra,dfdgaa)
      implicit none
c
c     This subroutine evaluates the PBE exchange functional [1] for
c     the closed shell case. The original PBE expression requires 
c     second derivatives of the density to evaluate the potential. This
c     requirement was lifted by Matt Challacombe.
c
c     The original subroutine was provided by 
c
c        Matt Challacombe
c        Los Alamos National Laboratory
c        The University of California
c
c     and is part of the MondoSCF suite of linear scaling electronic
c     structure codes. It was modified to conform with the design of
c     repository.
c
c     [1] John P. Perdew, Kieron Burke, Matthias Ernzerhof
c         "Generalized gradient approximation made simple"
c         Physical review letters Vol. 77 (1996) 3865-3868.
c
c     Parameters:
c
c     r      the total electron density
c     g      the dot product of the total density gradient with itself
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the alpha
c            electron density
c     dfdgaa On return the derivative of f with respect to the
c            dot product of the alpha density gradient with itself
c
      real*8 r, g
      real*8 f, dfdra, dfdgaa
c
      real*8 a, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11
c
      IF(1.0d-50.lt.R)THEN
         A=g
         R1=R**(1.d0/3.d0)
         R2=R*R1
         R3=R**2
         R4=R1**2
         R5=R3*R4
         R6=7.13182658760049d-3*A
         R7=R5+R6
         R8=R7**2
         R9=1/R8
         R10=A**2
         R11=R3**2
         f=-1.33236001455317d0*R2+(5.93801248171146d-1*R2)/(1.d0+ 
     &       (7.13182658760049d-3*A)/R5)
         dfdra=(-9.03570152478593d-5*R1*R10-8.39954475162682d-3*A*R*R3 
     &       -9.84745021842697d-1*R*R11*R4)*R9
         dfdgaa=-4.23488752945734d-3*R11*R9
         dfdgaa= 0.5d0*dfdgaa
      ELSE
         f      = 0.0d0
         dfdra  = 0.0d0
         dfdgaa = 0.0d0
      ENDIF
      END 
c-----------------------------------------------------------------------
      subroutine x_uks_pbe(ra,rb,gaa,gbb,f,dfdra,dfdrb,dfdgaa,dfdgbb)
      implicit none
c
c     This subroutine evaluates the PBE exchange functional [1] for
c     the open shell case. The original PBE expression requires 
c     second derivatives of the density to evaluate the potential. This
c     requirement was lifted by Matt Challacombe. This subroutine
c     uses x_rks_pbe to evaluate the functional using the spin-scaling
c     relationship of the exact exchange energy (e.g. see [1] Eq. (11)),
c     in analogy to the original PBE code [2].
c
c     The original subroutine was provided by 
c
c        Matt Challacombe
c        Los Alamos National Laboratory
c        The University of California
c
c     and is part of the MondoSCF suite of linear scaling electronic
c     structure codes. It was modified to conform with the design of
c     repository.
c
c     [1] John P. Perdew, Kieron Burke, Matthias Ernzerhof
c         "Generalized gradient approximation made simple"
c         Physical review letters Vol. 77 (1996) 3865-3868.
c
c     [2] John P. Perdew, Kieron Burke, Matthias Ernzerhof
c         "PBE alpha2.1"
c         http://crab.rutgers.edu/~kieron/pubs/PBE.asc
c         September 20, 1996.
c
c     Parameters:
c
c     ra     the alpha electron density
c     rb     the beta  electron density
c     gaa    the dot product of the alpha density gradient with itself
c     gbb    the dot product of the beta  density gradient with itself
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the alpha
c            electron density
c     dfdrb  On return the derivative of f with respect to the beta
c            electron density
c     dfdgaa On return the derivative of f with respect to the
c            dot product of the alpha density gradient with itself
c     dfdgbb On return the derivative of f with respect to the
c            dot product of the beta  density gradient with itself
c
      real*8 ra, rb, gaa, gbb
      real*8 f, dfdra, dfdrb, dfdgaa, dfdgbb
c
      real*8 ra2, rb2, ga2, gb2, fa, fb
c
      if (1.0d-20.lt.ra) then
         ra2 = 2.0d0*ra
         ga2 = 4.0d0*gaa
         call x_rks_pbe(ra2,ga2,fa,dfdra,dfdgaa)
      else
         fa     = 0.0d0
         dfdra  = 0.0d0
         dfdgaa = 0.0d0
      endif
      if (1.0d-20.lt.rb) then
         rb2 = 2.0d0*rb
         gb2 = 4.0d0*gbb
         call x_rks_pbe(rb2,gb2,fb,dfdrb,dfdgbb)
      else
         fb     = 0.0d0
         dfdrb  = 0.0d0
         dfdgbb = 0.0d0
      endif
      f = 0.5d0*(fa+fb)
      END 
c-----------------------------------------------------------------------

