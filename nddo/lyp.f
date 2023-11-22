      subroutine crlyp(r,g,f,dfdra,dfdgaa,dfdgab)
      implicit none
c
c     This subroutine evaluates LYP correlation functional [1] and
c     the ingredients of the corresponding potential for the closed
c     shell case. The implementation is based on the expression by 
c     Miehlich et al. [2]. This expression is simpler than the original 
c     one in the sense that the second derivatives of the density have 
c     been removed through partial integration. 
c
c     The o iginal code was provided by Dr. Phillip Young.
c
c     [1] C. Lee, W. Yang, and R.G. Parr,
c         "Development of the Colle-Salvetti correlation-energy
c          formula into a functional of the electron density"
c         Phys.Rev. Vol. B37 (1988) 785-789.
c
c     [2] B. Miehlich, A. Savin, H. Stoll, and H. Preuss,
c         "Results obtained with the correlation energy density
c          functionals of Becke and Lee, Yang and Parr"
c         Chem.Phys.Lett. Vol. 157 (1989) 200-206.
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
c     dfdgab On return the derivative of f with respect to the dot 
c            product of the alpha density gradient with the beta density
c            gradient
c
      real*8 r, g
      real*8 f, dfdra, dfdgaa, dfdgab
c
      real*8 rhoa13,rho13,rho43
      real*8 rhoa83
      real*8 r_a,rab,raa,rab9,ra,gaa,gab
      real*8 faa,fab,dgaa,dgab
      real*8 d2agaa,d2bgaa
      real*8 d2agab
      real*8 dfaaa,dfaab
      real*8 dfaba
      real*8 del,delp,ome,ro,xtc,xtd,rx,lt1,lt2a,lt2
      real*8 dlt1a,dlt2a
      real*8 zeta, cd
c
      real*8 t1,ta,tb,tc,td,t2,t3,t43,t83,t19,n13
      parameter(n13 = 0.3333333333333333d0)
      parameter(t1  = 0.3646239897876487d+02)   ! p1
      parameter(ta  = 0.04918d0)                ! a
      parameter(tb  = 0.13200d0)                ! b
      parameter(tc  = 0.25330d0)                ! c
      parameter(td  = 0.34900d0)                ! d
      parameter(t2  = ta*tb/td )                ! p2
      parameter(t3  = 4.0d0*ta/td)              ! p3
      parameter(t43 = 4.0d0*n13)
      parameter(t83 = 8.0d0*n13)
      parameter(t19 = 1.0d0/9.0d0)
c
      gaa    = 0.25d0*g
      gab    = gaa
c
      ra     = 0.5d0*r
      cd     = r
      rho13  = cd**n13
      rhoa13 = ra**n13
c
      rho43  = cd*rho13
      rhoa83 = rhoa13*rhoa13*ra*ra
c
      r_a  = 0.5d0
c
      raa = 0.25d0
      rab = 0.25d0
c
      zeta = 0.0d0
c
      xtc = tc/rho13
      xtd = td/rho13
c
      rx = xtd/(1.d0 + xtd)
c
      del = xtc + rx
      ome = t2*rx*exp(-xtc)/rho43
c
      rab9 = t19*rab
c
      faa = rab9*(1.d0 - 3.d0*del - (del - 11.d0)*r_a)
      fab = rab9*(47.d0 - 7.d0*del)
c
      dgaa = -ome*(faa - raa)
      dgab = -ome*(fab - t43)
c
      lt1  = -t3*rx*rab*rho43
      lt2a = -(t1*ome*rab)*rhoa83
      lt2  = 2.0d0*lt2a 
c
      f = lt1 + lt2 + 2.0d0*dgaa*gaa + dgab*gab
      RETURN
c
      delp = n13*(rx**2 - del)/cd
      ro   = n13*(del - 5.d0)/cd
c
      dfaaa = -rab9*( (3.d0 + r_a)*delp + (del - 11.d0)*r_a/cd )
      dfaab = -rab9*( (3.d0 + r_a)*delp - (del - 11.d0)*r_a/cd )
      dfaba = - 7.d0*rab9*delp
c
      dlt1a = (lt1/cd)*(n13*rx + 1.0d0)
      dlt2a = ro*lt2a + (ro + (t83/ra))*lt2a
c
      d2agaa = ro*dgaa - ome * (dfaaa + raa*(2.d0/cd))
      d2bgaa = ro*dgaa - ome * (dfaab - rab*(2.d0/cd))
      d2agab = ro*dgab - ome* dfaba
c
      dfdra = dlt1a + dlt2a +
     &        d2agaa*gaa + d2bgaa*gaa + d2agab*gab
c
      dfdgaa = dgaa
      dfdgab = dgab
c
      end
c-----------------------------------------------------------------------
      subroutine culyp(ra,rb,gaa,gbb,gab,f,dfdra,dfdrb,
     &                     dfdgaa,dfdgbb,dfdgab)
      implicit none
c
c     This subroutine evaluates LYP correlation functional [1] and
c     the ingredient of the corresponding potential. The implementation
c     is based on the expression by Miehlich et al. [2]. This expression
c     is simpler than the original one in the sense that the second
c     derivatives of the density have been removed through partial 
c     integration. 
c
c     The original code was provided by Dr. Phillip Young.
c
c     [1] C. Lee, W. Yang, and R.G. Parr,
c         "Development of the Colle-Salvetti correlation-energy
c          formula into a functional of the electron density"
c         Phys.Rev. Vol. B37 (1988) 785-789.
c
c     [2] B. Miehlich, A. Savin, H. Stoll, and H. Preuss,
c         "Results obtained with the correlation energy density
c          functionals of Becke and Lee, Yang and Parr"
c         Chem.Phys.Lett. Vol. 157 (1989) 200-206.
c
c     Parameters:
c
c     ra     the alpha-electron density
c     rb     the beta-electron density
c     gaa    the dot product of the alpha density gradient with itself
c     gbb    the dot product of the beta density gradient with itself
c     gab    the dot product of the alpha density gradient with 
c            the beta density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to ra
c     dfdrb  On return the derivative of f with respect to rb
c     dfdgaa On return the derivative of f with respect to gaa
c     dfdgbb On return the derivative of f with respect to gbb
c     dfdgab On return the derivative of f with respect to gab
c
      real*8 ra, rb, gaa, gbb, gab
      real*8 f, dfdra, dfdrb, dfdgaa, dfdgbb, dfdgab
c
      real*8 rhoa13,rhob13,rho13,rho43
      real*8 rhoa83,rhob83
      real*8 r_a,r_b,rab,raa,rbb,rab9
      real*8 faa,fbb,fab,dgaa,dgbb,dgab
      real*8 d2agaa,d2bgaa
      real*8 d2agbb,d2bgbb
      real*8 d2agab,d2bgab
      real*8 dfaaa,dfaab
      real*8 dfbba,dfbbb
      real*8 dfaba,dfabb
      real*8 del,delp,ome,ro,xtc,xtd,rx,lt1,lt2a,lt2b,lt2
      real*8 dlt1a,dlt1b,dlt2a,dlt2b
      real*8 zeta, cd
c
      real*8 t1,ta,tb,tc,td,t2,t3,t43,t83,t19,n13
      parameter(n13 = 0.3333333333333333d0)
      parameter(t1  = 0.3646239897876487d+02)   ! p1
      parameter(ta  = 0.04918d0)                ! a
      parameter(tb  = 0.13200d0)                ! b
      parameter(tc  = 0.25330d0)                ! c
      parameter(td  = 0.34900d0)                ! d
      parameter(t2  = ta*tb/td )                ! p2
      parameter(t3  = 4.0d0*ta/td)              ! p3
      parameter(t43 = 4.0d0*n13)
      parameter(t83 = 8.0d0*n13)
      parameter(t19 = 1.0d0/9.0d0)
c
      cd     = ra + rb
      rho13  = cd**n13
      rhoa13 = ra**n13
      rhob13 = rb**n13
c
      rho43  = cd*rho13
      rhoa83 = rhoa13*rhoa13*ra*ra
      rhob83 = rhob13*rhob13*rb*rb
c
      r_a  = ra/cd
      r_b  = rb/cd
c
      rab = r_a*r_b
      raa = r_a*r_a
      rbb = r_b*r_b
c
      zeta = r_a - r_b
c
      xtc = tc/rho13
      xtd = td/rho13
c
      rx = xtd/(1.d0 + xtd)
c
      del = xtc + rx
      ome = t2*rx*exp(-xtc)/rho43
c
      rab9 = t19*rab
c
      faa = rab9*(1.d0 - 3.d0*del - (del - 11.d0)*r_a)
      fbb = rab9*(1.d0 - 3.d0*del - (del - 11.d0)*r_b)
      fab = rab9*(47.d0 - 7.d0*del)
c
      dgaa = -ome*(faa - rbb)
      dgbb = -ome*(fbb - raa)
      dgab = -ome*(fab - t43)
c
      lt1  = -t3*rx*rab*rho43
      lt2a = -(t1*ome*rab)*rhoa83
      lt2b = -(t1*ome*rab)*rhob83
      lt2  = lt2a + lt2b
c
      f = lt1 + lt2 + dgaa*gaa + dgbb*gbb + dgab*gab
c
      delp = n13*(rx**2 - del)/cd
      ro   = n13*(del - 5.d0)/cd
c
      if (ra.le.0.0d0.or.rb.le.0.0d0) then
         dfaaa = 0.0d0
         dfaab = 0.0d0
         dfbba = 0.0d0
         dfbbb = 0.0d0
         dfaba = 0.0d0
         dfabb = 0.0d0
         dlt1a = 0.0d0
         dlt1b = 0.0d0
         dlt2a = 0.0d0
         dlt2b = 0.0d0
      else
         dfaaa = -(zeta/ra)*faa -
     &           rab9*( (3.d0 + r_a)*delp + (del - 11.d0)*r_b/cd )
         dfaab =  (zeta/rb)*faa -
     &           rab9*( (3.d0 + r_a)*delp - (del - 11.d0)*r_a/cd )
         dfbba = -(zeta/ra)*fbb -
     &           rab9*( (3.d0 + r_b)*delp - (del - 11.d0)*r_b/cd )
         dfbbb =  (zeta/rb)*fbb -
     &           rab9*( (3.d0 + r_b)*delp + (del - 11.d0)*r_a/cd )
         dfaba = -(zeta/ra)*fab - 7.d0*rab9*delp
         dfabb =  (zeta/rb)*fab - 7.d0*rab9*delp
c
         dlt1a = (lt1/cd)*(n13*rx + r_b/r_a)
         dlt1b = (lt1/cd)*(n13*rx + r_a/r_b)
         dlt2a = (ro - (zeta/ra))*lt2b +
     &           (ro + (t83 - zeta)/ra)*lt2a
         dlt2b = (ro + (zeta/rb))*lt2a +
     &           (ro + (t83 + zeta)/rb)*lt2b
      endif
c
      d2agaa = ro*dgaa - ome * (dfaaa + rbb*(2.d0/cd))
      d2bgaa = ro*dgaa - ome * (dfaab - rab*(2.d0/cd))
      d2agbb = ro*dgbb - ome * (dfbba - rab*(2.d0/cd))
      d2bgbb = ro*dgbb - ome * (dfbbb + raa*(2.d0/cd))
      d2agab = ro*dgab - ome* dfaba
      d2bgab = ro*dgab - ome* dfabb
c
      dfdra = dlt1a + dlt2a +
     &        d2agaa*gaa + d2agbb*gbb + d2agab*gab
      dfdrb = dlt1b + dlt2b +
     &        d2bgaa*gaa + d2bgbb*gbb + d2bgab*gab
c
      dfdgaa = dgaa
      dfdgbb = dgbb
      dfdgab = dgab
c
      end

