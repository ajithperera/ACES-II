c-----------------------------------------------------------------------
      subroutine xrb3(r,g,f,dfdra,dfdgaa)
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
      subroutine xub3(ra,rb,gaa,gbb,f,dfdra,dfdrb,dfdgaa,dfdgbb)
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

