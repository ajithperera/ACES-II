c-----------------------------------------------------------------------
      subroutine xcrb3lyp(r,g,f,dfdra,dfdgaa,dfdgab)
      implicit none
c
c     This subroutine evaluates the B3LYP exchange-correlation 
c     functional [1,2,3] for the closed shell case and its derivative
c     with respect to the alpha density, and the density gradient
c     dot products.
c
c     [1] P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch,
c         "Ab initio calculation of vibrational absorption and circular
c          dichroism spectra using density functional force fields",
c         J.Phys.Chem. Vol 98 (1994) 11623-11627.
c
c     [2] R.H. Hertwig, W. Koch,
c         "On the parameterization of the local correlation functional:
c          What is Becke-3-LYP?",
c         Chem.Phys.Lett. Vol 268 (1997) 345-351.
c
c     [3] "Gaussian 98 User's Reference
c          Density Functional Methods (DFT) Keywords",
c         M.J. Frisch, G.W. Trucks, H.B. Schlegel, G.E. Scuseria, 
c         M.A. Robb, J.R. Cheeseman, V.G. Zakrzewski, J.A. Montgomery, 
c         R.E. Stratmann, J.C. Burant, S. Dapprich, J.M. Millam, 
c         A.D. Daniels, K.N. Kudin, M.C. Strain, O. Farkas, J. Tomasi, 
c         V. Barone, M. Cossi, R. Cammi, B. Mennucci, C. Pomelli, 
c         C. Adamo, S. Clifford, J. Ochterski, G.A. Petersson, 
c         P.Y. Ayala, Q. Cui, K. Morokuma, D.K. Malick, A.D. Rabuck, 
c         K. Raghavachari, J.B. Foresman, J. Cioslowski, J.V. Ortiz, 
c         B.B. Stefanov, G. Liu, A. Liashenko, P. Piskorz, I. Komaromi,
c         R. Gomperts, R.L. Martin, D.J. Fox, T. Keith, M.A. Al-Laham, 
c         C.Y. Peng, A. Nanayakkara, C. Gonzalez, M. Challacombe, 
c         P.M.W. Gill, B.G. Johnson, W. Chen, M.W. Wong, J.L. Andres, 
c         M. Head-Gordon, E.S. Replogle, J.A. Pople,
c         Gaussian, Inc., Pittsburgh PA, 1998.
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
      real*8 f1, dfdra1, dfdgaa1, dfdgab1
c
      real*8 vwnrpa_wt, lyp_wt
      parameter(vwnrpa_wt=0.19d0)
      parameter(lyp_wt   =0.81d0)
c
      f      = 0.0d0
      dfdra  = 0.0d0
      dfdgaa = 0.0d0
      dfdgab = 0.0d0
      call c_rks_vwnrpa(r,f1,dfdra1)
      f      = f      + vwnrpa_wt*f1
      dfdra  = dfdra  + vwnrpa_wt*dfdra1
      call x_rks_b3(r,g,f1,dfdra1,dfdgaa1)
      f      = f      + f1
      dfdra  = dfdra  + dfdra1
      dfdgaa = dfdgaa + dfdgaa1
      call c_rks_lyp(r,g,f1,dfdra1,dfdgaa1,dfdgab1)
      f      = f      + lyp_wt*f1
      dfdra  = dfdra  + lyp_wt*dfdra1
      dfdgaa = dfdgaa + lyp_wt*dfdgaa1
      dfdgab = dfdgab + lyp_wt*dfdgab1
c
      end
c-----------------------------------------------------------------------
      subroutine xcub3lyp(ra,rb,gaa,gbb,gab,f,dfdra,dfdrb,
     +                        dfdgaa,dfdgbb,dfdgab)
      implicit none
c
c     This subroutine evaluates the B3LYP exchange-correlation 
c     functional [1,2,3] for the closed shell case and its derivative
c     with respect to the alpha density, and the density gradient
c     dot products.
c
c     [1] P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch,
c         "Ab initio calculation of vibrational absorption and circular
c          dichroism spectra using density functional force fields",
c         J.Phys.Chem. Vol 98 (1994) 11623-11627.
c
c     [2] R.H. Hertwig, W. Koch,
c         "On the parameterization of the local correlation functional:
c          What is Becke-3-LYP?",
c         Chem.Phys.Lett. Vol 268 (1997) 345-351.
c
c     [3] "Gaussian 98 User's Reference
c          Density Functional Methods (DFT) Keywords",
c         M.J. Frisch, G.W. Trucks, H.B. Schlegel, G.E. Scuseria, 
c         M.A. Robb, J.R. Cheeseman, V.G. Zakrzewski, J.A. Montgomery, 
c         R.E. Stratmann, J.C. Burant, S. Dapprich, J.M. Millam, 
c         A.D. Daniels, K.N. Kudin, M.C. Strain, O. Farkas, J. Tomasi, 
c         V. Barone, M. Cossi, R. Cammi, B. Mennucci, C. Pomelli, 
c         C. Adamo, S. Clifford, J. Ochterski, G.A. Petersson, 
c         P.Y. Ayala, Q. Cui, K. Morokuma, D.K. Malick, A.D. Rabuck, 
c         K. Raghavachari, J.B. Foresman, J. Cioslowski, J.V. Ortiz, 
c         B.B. Stefanov, G. Liu, A. Liashenko, P. Piskorz, I. Komaromi,
c         R. Gomperts, R.L. Martin, D.J. Fox, T. Keith, M.A. Al-Laham, 
c         C.Y. Peng, A. Nanayakkara, C. Gonzalez, M. Challacombe, 
c         P.M.W. Gill, B.G. Johnson, W. Chen, M.W. Wong, J.L. Andres, 
c         M. Head-Gordon, E.S. Replogle, J.A. Pople,
c         Gaussian Inc., Pittsburgh PA, 1998.
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
      real*8 f1, dfdra1, dfdrb1, dfdgaa1, dfdgbb1, dfdgab1
c
      real*8 vwnrpa_wt, lyp_wt
      parameter(vwnrpa_wt=0.19d0)
      parameter(lyp_wt   =0.81d0)
c
      dfdgab = 0.0d0

      call x_uks_b3(ra,rb,gaa,gbb,f,dfdra,dfdrb,dfdgaa,dfdgbb)

      call c_uks_vwnrpa(ra,rb,f1,dfdra1,dfdrb1)
      f      = f      + vwnrpa_wt*f1
      dfdra  = dfdra  + vwnrpa_wt*dfdra1
      dfdrb  = dfdrb  + vwnrpa_wt*dfdrb1
      call c_uks_lyp(ra,rb,gaa,gbb,gab,f1,dfdra1,dfdrb1,
     +               dfdgaa1,dfdgbb1,dfdgab1)
      f      = f      + lyp_wt*f1
      dfdra  = dfdra  + lyp_wt*dfdra1
      dfdrb  = dfdrb  + lyp_wt*dfdrb1
      dfdgaa = dfdgaa + lyp_wt*dfdgaa1
      dfdgab = dfdgab + lyp_wt*dfdgab1
      dfdgbb = dfdgbb + lyp_wt*dfdgbb1
c
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine x_rks_b3(r,g,f,dfdra,dfdgaa)
      implicit none
c
c     This subroutine evaluates the exchange part of the B3xxx
c     functionals excluding the exact exchange term for closed shell
c     cases. In the context of Becke's paper [1] this means that this 
c     subroutine evaluates the part
c
c         (1-a0) E_X^LSDA + aX dE_X^B88
c
c     of equation (2). To complete the exchange part one has to add
c     
c         a0 E_X^exact
c
c     In [1] a0 and aX were given as
c
c         a0 = 0.20
c         aX = 0.72
c
c     In addition to the functional the ingredients of the corresponding
c     potential are calculated as well.
c     
c
c     The original code was provided by Dr. Phillip Young.
c
c     [1] A.D. Becke,
c         "Density-functional thermochemistry. III. The role of 
c          exact exchange",
c         J.Chem.Phys. Vol. 98 (1993) 5648-5652.
c
c     Parameters:
c
c     r      the total electron density
c     g      the dot product of the total density gradient with itself
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the alpha-
c            electron density
c     dfdgaa On return the derivative of f with respect to the
c            dot product of the alpha-density gradient with itself
c
      real*8 a0, aX
      parameter(a0 = 0.20d0)
      parameter(aX = 0.72d0)
c
c     A = -3/4*(6/pi)**(1/3)
      real*8 A
      parameter(A    = -0.930525736349100013861335177933319d0)
c
      real*8 third, third4
      parameter(third  =  0.333333333333333333333333333333333d0)
      parameter(third4 =  1.333333333333333333333333333333333d0)
c
c...  BeckeP : The parameter beta in Becke's paper
c...  BeckeP6: 6 Times the parameter beta in Becke's paper
c
      real*8 BeckeP, BeckeP6
      parameter(BeckeP  = 0.0042d0)
      parameter(BeckeP6 = 0.0252d0)
c
      real*8 r, g
      real*8 f, dfdra, dfdgaa
      real*8 ra, gaa
      real*8 alpha_rho13, alpha_rho43
c
      real*8 alpha_X, alpha_X2, alpha_Xt, alpha_Xhi, alpha_de, alpha_G
      real*8 alpha_Gp, gamma_a
c
      gaa         = 0.25d0*g
      ra          = 0.5d0*r
      alpha_rho13 = ra**third
      alpha_rho43 = ra*alpha_rho13
c
      if (alpha_rho43.gt.0.0d0) then
         gamma_a     = sqrt(gaa)
c...     Xa = the reduced density gradient
         alpha_X     = gamma_a/alpha_rho43
         alpha_X2    = alpha_X * alpha_X
c ...    Sa = sqrt(Xa*Xa + ONE)
         alpha_Xt    = dsqrt(1.0d0+alpha_X2)
c ...    Ha = log(Xa + Sa)
         alpha_Xhi   = log(alpha_X+alpha_Xt)
         alpha_Xt    = alpha_X/alpha_Xt
c
c ...    Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
         alpha_de = 1.0d0/(1.0d0+BeckeP6*alpha_X*alpha_Xhi)
c ...    Ga = P2LDA + Ca
         alpha_G  = (1.0d0-a0)*A-aX*(BeckeP*alpha_X2)*alpha_de
c
c ...    Pa = ( (SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa) )*Ta**2
c
         alpha_Gp =(6.0d0*(BeckeP*alpha_X)**2*(alpha_Xt-alpha_Xhi)-
     &             (2.0d0*BeckeP*alpha_X) ) * (alpha_de*alpha_de)
c        FOURTHIRDS*CD13a*(Ga - Xa*Pa)*fac
         dfdra = third4*alpha_rho13*(alpha_G - aX*alpha_X*alpha_Gp)
c
         if(gamma_a.gt.0.0d0) then
            dfdgaa = (0.5d0/gamma_a)*aX*alpha_Gp
         else
            dfdgaa = 0.0d0
         endif
      else
         dfdra   = 0.0d0
         dfdgaa  = 0.0d0
         alpha_G = 0.0d0
      endif
c
      f = 2.0d0*alpha_rho43*alpha_G
      end
c-----------------------------------------------------------------------
      subroutine x_uks_b3(ra,rb,gaa,gbb,f,dfdra,dfdrb,dfdgaa,dfdgbb)
      implicit none
c
c     This subroutine evaluates the exchange part of the B3xxx
c     functionals excluding the exact exchange term. In the context of
c     Becke's paper [1] this means that this subroutine evaluates the
c     part
c
c         (1-a0) E_X^LSDA + aX dE_X^B88
c
c     of equation (2). To complete the exchange part one has to add
c     
c         a0 E_X^exact
c
c     In [1] a0 and aX were given as
c
c         a0 = 0.20
c         aX = 0.72
c
c     In addition to the functional the ingredients of the corresponding
c     potential are calculated as well.
c     
c
c     The original code was provided by Dr. Phillip Young.
c
c     [1] A.D. Becke,
c         "Density-functional thermochemistry. III. The role of 
c          exact exchange",
c         J.Chem.Phys. Vol. 98 (1993) 5648-5652.
c
c     Parameters:
c
c     ra     the alpha-electron density
c     rb     the beta-electron density
c     gaa    the dot product of the alpha density gradient with itself
c     gbb    the dot product of the beta density gradient with itself
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to ra
c     dfdrb  On return the derivative of f with respect to rb
c     dfdgaa On return the derivative of f with respect to gaa
c     dfdgbb On return the derivative of f with respect to gbb
c
      real*8 a0, aX
      parameter(a0 = 0.20d0)
      parameter(aX = 0.72d0)
c
c     A = -3/4*(6/pi)**(1/3)
      real*8 A
      parameter(A    = -0.930525736349100013861335177933319d0)
c
      real*8 third, third4
      parameter(third  =  0.333333333333333333333333333333333d0)
      parameter(third4 =  1.333333333333333333333333333333333d0)
c
c...  BeckeP : The parameter beta in Becke's paper
c...  BeckeP6: 6 Times the parameter beta in Becke's paper
c
      real*8 BeckeP, BeckeP6
      parameter(BeckeP  = 0.0042d0)
      parameter(BeckeP6 = 0.0252d0)
c
      real*8 ra, rb, gaa, gbb
      real*8 f, dfdra, dfdrb, dfdgaa, dfdgbb
      real*8 alpha_rho13, alpha_rho43, beta_rho13, beta_rho43
c
      real*8 alpha_X, alpha_X2, alpha_Xt, alpha_Xhi, alpha_de, alpha_G
      real*8 alpha_Gp, gamma_a
      real*8 beta_X, beta_X2, beta_Xt, beta_Xhi, beta_de, beta_G
      real*8 beta_Gp, gamma_b
c
      alpha_rho13 = ra**third
      beta_rho13  = rb**third
      alpha_rho43 = ra*alpha_rho13
      beta_rho43  = rb*beta_rho13
c
      if (alpha_rho43.gt.0.0d0) then
         gamma_a     = sqrt(gaa)
c...     Xa = the reduced density gradient
         alpha_X     = gamma_a/alpha_rho43
         alpha_X2    = alpha_X * alpha_X
c ...    Sa = sqrt(Xa*Xa + ONE)
         alpha_Xt    = dsqrt(1.0d0+alpha_X2)
c ...    Ha = log(Xa + Sa)
         alpha_Xhi   = log(alpha_X+alpha_Xt)
         alpha_Xt    = alpha_X/alpha_Xt
c
c ...    Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
         alpha_de = 1.0d0/(1.0d0+BeckeP6*alpha_X*alpha_Xhi)
c ...    Ga = P2LDA + Ca
         alpha_G  = (1.0d0-a0)*A-aX*(BeckeP*alpha_X2)*alpha_de
c
c ...    Pa = ( (SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa) )*Ta**2
c
         alpha_Gp =(6.0d0*(BeckeP*alpha_X)**2*(alpha_Xt-alpha_Xhi)-
     &             (2.0d0*BeckeP*alpha_X) ) * (alpha_de*alpha_de)
c        FOURTHIRDS*CD13a*(Ga - Xa*Pa)*fac
         dfdra = third4*alpha_rho13*(alpha_G - aX*alpha_X*alpha_Gp)
c
         if(gamma_a.gt.0.0d0) then
            dfdgaa = (0.5d0/gamma_a)*aX*alpha_Gp
         else
            dfdgaa = 0.0d0
         endif
      else
         dfdra   = 0.0d0
         dfdgaa  = 0.0d0
         alpha_G = 0.0d0
      endif
c
      if (beta_rho43.gt.0.0d0) then
         gamma_b     = sqrt(gbb)
         beta_X      = gamma_b/beta_rho43
         beta_X2     = beta_X * beta_X
         beta_Xt     = sqrt(beta_X*beta_X+1.0d0)
         beta_Xhi    = log(beta_X+beta_Xt)
         beta_Xt     = beta_X/beta_Xt
c
C        energy
c
c        Tb = 1/(ONE + SIX*(BETA*Xb)*Hb)
         beta_de  = 1.0d0/(1.0d0+BeckeP6*beta_X*beta_Xhi)
c        Gb = P2LDA + Cb
         beta_G   = (1.0d0-a0)*A-aX*(BeckeP*beta_X2)*beta_de
C    
C        potential and gradient potential
c
c        Pb = ( (SIX*(BETA*Xb)**2)*(Xb/Sb - Hb) - TWO*(BETA*Xb) )*Tb**2
         beta_Gp  =(6.0d0*(BeckeP*beta_X)**2*(beta_Xt -beta_Xhi )-
     &             (2.0d0*BeckeP*beta_X ) ) * (beta_de *beta_de )
         dfdrb = third4*beta_rho13 *(beta_G   - aX*beta_X *beta_Gp )
c
c        gaa = Pa/GNa*fac
         if(gamma_b.gt.0.0d0) then
            dfdgbb = (0.5d0/gamma_b)*aX*beta_Gp
         else
            dfdgbb = 0.0d0
         endif
      else
         dfdrb  = 0.0d0
         dfdgbb = 0.0d0
         beta_G = 0.0d0
      endif
c
      f = alpha_rho43*alpha_G+beta_rho43*beta_G
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine c_rks_lyp(r,g,f,dfdra,dfdgaa,dfdgab)
      implicit none
c
c     This subroutine evaluates LYP correlation functional [1] and
c     the ingredients of the corresponding potential for the closed
c     shell case. The implementation is based on the expression by 
c     Miehlich et al. [2]. This expression is simpler than the original 
c     one in the sense that the second derivatives of the density have 
c     been removed through partial integration. 
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
      subroutine c_uks_lyp(ra,rb,gaa,gbb,gab,f,dfdra,dfdrb,
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
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine c_rks_vwnrpa(r,f,dfdra)
      implicit none
c
c     This subroutine evaluates the Vosko, Wilk and Nusair correlation
c     functional for the closed shell case using [4.5] with the RPA 
c     parametrisation given below equation [4.4].
c
c     The original code was obtained from Dr. Phillip Young.
c
c     [1] S.H. Vosko, L. Wilk, and M. Nusair
c         "Accurate spin-dependent electron liquid correlation energies
c          for local spin density calculations: a critical analysis",
c         Can.J.Phys, Vol. 58 (1980) 1200-1211.
c
c     Parameters:
c
c     r      the total electron density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the alpha-
c            electron density
c
      real*8 r
      real*8 f, dfdra
c
      real*8 t4,t5,t6,t7
      real*8 a2,b2,c2,d2
      real*8 P1,P2,P3,P4
      real*8 srho,srho13
      real*8 iv,iv2,inv,i1,i2,i3
      real*8 pp1,pp2
c
      real*8 n13, n43, n16, one, two, four
      parameter(n13 = 0.333333333333333333333333333333d0)
      parameter(n43 = 1.333333333333333333333333333333d0)
      parameter(n16 = 0.166666666666666666666666666666d0)
      parameter(one = 1.0d0)
      parameter(two = 2.0d0)
      parameter(four= 4.0d0)
C 
C VWN (RPA) interpolation parameters
C
C paramagnetic
      a2 =  0.0621814d0 
      b2 =  13.0720d0           ! ch
      c2 =  42.7198d0           ! ch
      d2 = -0.409286d0          ! ch
C
C t4 = (1/(4/3)*pi)**(1/3)
      t4=0.620350490899399531d0
C
C t5 = 0.5/(2**(1/3)-1)
      t5 = 1.92366105093153617d0
C
C t6 = 2/(3*(2**(1/3)-1))
      t6 = 2.56488140124204822d0
C
C t7 = 2.25*(2**(1/3)-1)
      t7 = 0.584822362263464735d0
C
C
C Paramagnetic interpolation constants
C 
      p1 = 0.448998886415767975d-01
      p2 = 0.310907d-01
      p3 = 0.443137364463981956d-02
      p4 = 20.5219723675762680d0
C
C closed shell case
      srho        = r
      srho13      = srho**n13
      iv2         = T4/srho13
      iv          = sqrt(iv2)
C
C paramagnetic
      inv = 1.0d0/(iv2+b2*iv+c2)
      i1  = log(iv2*inv)
      i2  = log((iv-d2)*(iv-d2)*inv)
c corrected b1->b2 ps Apr98
      i3  = atan(P1/(2.0d0*iv+b2))
      pp1 = P2*i1 + P3*i2 + P4*i3
      pp2 = a2*(1.0d0/iv-iv*inv*(1.0d0+b2/(iv-d2))) 
C
      f     = pp1*srho
      dfdra = pp1 - n16*iv*pp2
c
      end
c-----------------------------------------------------------------------
      subroutine c_uks_vwnrpa(ra,rb,f,dfdra,dfdrb)
      implicit none
c
c     This subroutine evaluates the Vosko, Wilk and Nusair correlation
c     functional using [4.5] with the RPA parametrisation given below 
c     equation [4.4].
c
c     The original code was obtained from Dr. Phillip Young.
c
c     [1] S.H. Vosko, L. Wilk, and M. Nusair
c         "Accurate spin-dependent electron liquid correlation energies
c          for local spin density calculations: a critical analysis",
c         Can.J.Phys, Vol. 58 (1980) 1200-1211.
c
c     Parameters:
c
c     ra     the alpha-electron density
c     rb     the beta-electron density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to ra
c     dfdrb  On return the derivative of f with respect to rb
c
      real*8 ra, rb
      real*8 f, dfdra, dfdrb
c
      real*8 t1,t2,t4,t5,t6,t7
      real*8 a1,b1,c1,d1
      real*8 a2,b2,c2,d2
      real*8 a3,b3,c3,d3
      real*8 S1,S2,S3,S4
      real*8 P1,P2,P3,P4
      real*8 F1,F2,F3,F4
      real*8 inter1,inter2
      real*8 alpha_rho13,beta_rho13
      real*8 srho,srho13
      real*8 iv,iv2,inv,i1,i2,i3
      real*8 vwn1,vwn2
      real*8 ss1,ss2,pp1,pp2,ff1,ff2
      real*8 zeta,zeta3,zeta4
      real*8 tau,dtau,v
c
      real*8 n13, n43, n16, one, two, four, tn13, tn43
      parameter(n13 = 0.333333333333333333333333333333d0)
      parameter(n43 = 1.333333333333333333333333333333d0)
      parameter(n16 = 0.166666666666666666666666666666d0)
      parameter(tn13= 1.25992104989487316476721060727823d0)
      parameter(tn43= 2.51984209978974632953442121455646d0)
      parameter(one = 1.0d0)
      parameter(two = 2.0d0)
      parameter(four= 4.0d0)
C 
C VWN (RPA) interpolation parameters
C
C spin stiffness
      a1 = -0.337737278807791058d-01
      b1 =  1.13107d0
      c1 = 13.0045d0
      d1 = -0.00475840d0
C paramagnetic
      a2 =  0.0621814d0 
      b2 =  13.0720d0           ! ch
      c2 =  42.7198d0           ! ch
      d2 = -0.409286d0          ! ch
C ferromagnetic
c try cadpac/nwchem value (.5*a2)
      a3 = 0.0310907d0
      b3 =  20.1231d0           ! ch
      c3 = 101.578d0            ! ch
      d3 = -0.743294d0          ! ch
C
C t4 = (1/(4/3)*pi)**(1/3)
      t4=0.620350490899399531d0
C
C t5 = 0.5/(2**(1/3)-1)
      t5 = 1.92366105093153617d0
C
C t6 = 2/(3*(2**(1/3)-1))
      t6 = 2.56488140124204822d0
C
C t7 = 2.25*(2**(1/3)-1)
      t7 = 0.584822362263464735d0
C
C Spin stiffness interpolation constants
C
      S1 = 7.12310891781811772d0
      S2 = a1 *0.5d0
      S3 = -0.139834647015288626d-04 *0.5d0
      S4 = -0.107301836977671539d-01 *0.5d0
C
C Paramagnetic interpolation constants
C 
c     p1 = sqrt(4.0d0*c2 - b2*b2)
c     p2 = 0.5d0 * a2
c     p3 = -0.5d0 * a2*b2*d2/(d2*d2 + b2*d2 + c2)
c     p4 = 0.5d0 * a2*(two*b2/p1)*((c2-d2*d2)/(d2*d2 + b2*d2 + c2))
      p1 = 0.448998886415767975d-01
      p2 = 0.310907d-01
      p3 = 0.443137364463981956d-02
      p4 = 20.5219723675762680d0
C
C Ferromagnetic interpolation constants
C
c     f1 = sqrt(four*c3 - b3*b3)
c     f2 = 0.5d0 * a3
c     f3 = -0.5d0 * a3*b3*d3/(d3*d3 + b3*d3 + c3)
c     f4 = 0.5d0 * a3*(two*b3/f1)*((c3-d3*d3)/(d3*d3 + b3*d3 + c3))
c
      f1 = 1.17168527770897146d0
      f2 = 0.155453495681285858d-01
      f3 = 0.266730993317173832d-02
      f4 = 0.618818012598993383d0
C
C Interpolation intervals
C 
      inter1 =  1.0d0-1.0d-10
      inter2 = -1.0d0+1.0d-10
C
C open shell case
      alpha_rho13 = ra**n13
      beta_rho13  = rb**n13
      srho        = ra+rb
      srho13      = srho**n13
      iv2         = T4/srho13
      iv          = sqrt(iv2)
C
C spin-stiffness
      inv = 1.0d0/(iv2+b1*iv+c1)
      i1  = log(iv2*inv)
      i2  = log((iv-d1)*(iv-d1)*inv)
      i3  = atan(S1/(2.0d0*iv+b1))
      ss1 = S2*i1 + S3*i2 + S4*i3
      ss2 = a1*(1.0d0/iv-iv*inv*(1.0d0+b1/(iv-d1)))
C
C paramagnetic
      inv = 1.0d0/(iv2+b2*iv+c2)
      i1  = log(iv2*inv)
      i2  = log((iv-d2)*(iv-d2)*inv)
c corrected b1->b2 ps Apr98
      i3  = atan(P1/(2.0d0*iv+b2))
      pp1 = P2*i1 + P3*i2 + P4*i3
      pp2 = a2*(1.0d0/iv-iv*inv*(1.0d0+b2/(iv-d2))) 
C
C ferromagnetic
      inv = 1.0d0/(iv2+b3*iv+c3)
      i1  = log(iv2*inv)
      i2  = log((iv-d3)*(iv-d3)*inv)
      i3  = atan(F1/(2.0d0*iv+b3))
      ff1 = F2*i1 + F3*i2 + F4*i3
      ff2 = a3*(1.0d0/iv-iv*inv*(1.0d0+b3/(iv-d3)))
C
C polarisation function
c     zeta  = (alpha_rho13-beta_rho13)/srho13
c
      zeta  = (ra-rb)/srho
      zeta3 = zeta*zeta*zeta
      zeta4 = zeta3*zeta
      if(zeta.gt.inter1) then
         vwn1 = (tn43-two)*t5
         vwn2 = (tn13)*t6
      elseif(zeta.lt.inter2) then
         vwn1 = (tn43-two)*t5
         vwn2 = -(tn13)*t6
      else
         vwn1 = ((one+zeta)**n43+(one-zeta)**n43-two)*t5
         vwn2 = ((one+zeta)**n13-(one-zeta)**n13)*t6
      endif
      ss1  = ss1*t7
      ss2  = ss2*t7 
c     tau  = ff1-pp1-ss1
c     dtau = ff2-pp2-ss2
      tau  = ff1-pp1
      dtau = ff2-pp2
c
c     v = pp1+vwn1*(ss1+tau*zeta4)
      v = pp1+vwn1*tau
      f = v*srho
c
c     t1 = v - n16*iv*(pp2+vwn1*(ss2+dtau*zeta4))
c     t2 = vwn2*(ss1+tau*zeta4)+vwn1*four*tau*zeta3
      t1 = v - n16*iv*(pp2+vwn1*dtau)
      t2 = vwn2*tau
      dfdra = t1+t2*(one-zeta)
      dfdrb = t1-t2*(one+zeta)
c
      end
c-----------------------------------------------------------------------

