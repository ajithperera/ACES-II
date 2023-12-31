      subroutine xc_pbe (
     &     iuks, ispotential, iscorrelation,
     &     rho, absgrad, delgrad, laplas, 
     &     ex, ec, vx, vc)
c
c This is a front-end for the PBE routines distributed by
c K. Burke.
c
c In:
c     iuks              - =0 restricted reference =1 unrestricted reference
c     ispotential       - calculate the potential?
c     iscorrelation     - calculate the correlation energy and potential?
c     rho     (iuks+1)  - value of the density
c     absgrad (iuks+1)  - absolute value of the density gradient abs(grad rho)
c     delgrad (iuks+1)  - grad rho . grad (abs(grad rho))
c     laplas  (iuks+1)  - laplasian of the density, grad . grad rho
c Out:
c     ex        - PBE exchange energy
c     ec        - PBE correlation energy      
c     vx        - PBE exchange potential
c     vc        - PBE correlation potential      
c
c Note that unlike this routine modern DFT codes do not calculate the
c value of Vxc at a grid point. Instead, they calculate the the matrix
c elements of Vxc directly using the integrated-by-parts formula (Pople,
c 1992, CPL:199:557). This formula does not require the hessians of the
c density. 
c
c Igor Schweigert, Apr 2005
c $Id: xc_pbe.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     Arguments
c
      logical ispotential, iscorrelation
c
      integer iuks
c
      double precision
     &     rho (iuks+1), absgrad (iuks+1), delgrad (iuks+1),
     &     laplas (iuks+1), ex, ec, vx(iuks+1), vc (iuks+1)
c
c     Local variables
c
      integer flag
c 
      double precision
     &     pi, pi32, third, thrd2, alpha,
     &     kf, ss, uu, vv, rs, zeta, sk, skg, g, t,
     &     e1, e2, v1, v2, v3, v4, dabs
c     
      parameter(alpha=1.91915829267751300662482032624669d0)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(third=1d0/3d0)
      parameter(thrd2=2d0/3d0)
c
      ex = 0d0
      ec = 0d0
      vx (1) = 0d0
      vc (1) = 0d0
      if (iuks.ne.0) then
         vx (2) = 0d0
         vc (2) = 0d0
      endif
      if (dabs (rho(1)+rho(iuks+1)).le.2d-18) goto 999
c      
      if (iuks.eq.0) then
c 
         kf = (3d0*pi**2*rho(1))**third
         ss = absgrad(1) / (2d0*kf*rho(1))
         if (ispotential) then
            uu = delgrad(1) / ((2d0*kf)**3*rho(1)*rho(1))
            vv = laplas(1)  / ((2d0*kf)**2*rho(1))
            flag = 1
         else
            uu = 0d0
            vv = 0d0
            flag = 0
         endif
         call exchpbe (rho(1), ss, uu, vv, 1, flag, ex, vx(1))
c 
         if (.not.iscorrelation) goto 999
         zeta = 0d0
         g  = ((1.d0+zeta)**thrd2+(1.d0-zeta)**thrd2)/2d0
         rs = alpha/kf
         sk = dsqrt(4d0*kf/pi)
         skg= 2d0*sk*g
         t  = absgrad(1) /(skg*rho(1))
         if (ispotential) then
            uu = delgrad(1) /(rho(1)*rho(1)*skg**3)
            vv = laplas(1)  /(rho(1)*skg**2)
            flag = 1
         else
            uu = 0d0
            vv = 0d0
            flag = 0
         endif
         call corpbe (
     &        rs, 0d0, t, uu, vv, 0d0, 1, flag,
     &        e1, v1, v2, e2, v3, v4)
         ec = e1 + e2
         if (ispotential) vc(1) = v1 + v3
c     
      else
         ex = 0d0
         ec = 0d0
         vx (1) = 0d0
         vx (2) = 0d0
         vc (1) = 0d0
         vc (2) = 0d0
      endif
c
 999  continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
      subroutine xc_lsda (
     &     iuks, ispotential, iscorrelation, rho,
     &     ex, ec, vx, vc)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This routine calculates the LSDA exchange-correlation energy and
c potenial using the PBE routines distributed by K.Burke.
c
c
c In:
c     iuks              - =0 restricted reference =1 unrestricted reference
c     ispotential       - calculate the potential?
c     iscorrelation     - calculate the correlation energy and potential?
c     rho     (iuks+1)  - value of the density
c
c Out:
c     ex        - PBE exchange energy
c     ec        - PBE correlation energy      
c     vx        - PBE exchange potential
c     vc        - PBE correlation potential
c
c Igor Schweigert, Apr 2005
c $Id: xc_pbe.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer iuks
c
      logical ispotential, iscorrelation
c
      double precision
     &     rho (iuks+1), ex, ec, vx (iuks+1), vc (iuks+1)
c 
c
c     Local variables
c
      double precision
     &     pi, pi32, third, thrd2, alpha,
     &     rs, zeta, e1, e2, v1, v2, v3, v4, dabs
c     
      parameter(alpha=1.91915829267751300662482032624669d0)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(pi=3.1415926535897932384626433832795d0)
      parameter(third=1d0/3d0)
      parameter(thrd2=2d0/3d0)
c
      ex = 0d0
      ec = 0d0
      vx (1) = 0d0
      vc (1) = 0d0
      if (iuks.ne.0) then
         vx (2) = 0d0
         vc (2) = 0d0
      endif
      if (dabs (rho(1)+rho(iuks+1)).le.2d-18) goto 999
c      
      if (iuks.eq.0) then
c 
         if (ispotential) then
            call exchpbe (rho(1), 0d0, 0d0, 0d0, 0, 1, ex, vx(1))
         else
            call exchpbe (rho(1), 0d0, 0d0, 0d0, 0, 0, ex, v1)
         endif
c 
         if (.not.iscorrelation) goto 999
         rs = alpha / (pi32*rho(1))**third
         if (ispotential) then
            call corpbe (
     &           rs, 0d0, 0d0, 0d0, 0d0, 0d0, 0, 1,
     &           ec, vc(1), v2, e2, v3, v4)
         else
            call corpbe (
     &           rs, 0d0, 0d0, 0d0, 0d0, 0d0, 0, 0,
     &           ec, v1, v2, e2, v3, v4)
         endif
c     
      else
c 
         ex = 0d0
         ec = 0d0
         vx (1) = 0d0
         vx (2) = 0d0
         vc (1) = 0d0
         vc (2) = 0d0
c 
      endif
c
 999  continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                  cc
cc Below are PBE routines distributed by K.Burke: EXCHPBE, CORPBE,  cc
cc GCOR2                                                            cc
cc                                                                  cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=3.14159265358979323846264338327950d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
	ex=exunif
        vx=ex*thrd4
	return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24) 
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn,
     1                  H,DVCUP,DVCDN)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c EU=unpolarized LSD correlation energy
c EURS=dEU/drs
c EP=fully polarized LSD correlation energy
c EPRS=dEP/drs
c ALFM=-spin stiffness, [c](3).
c ALFRSM=-dalpha/drs
c F=spin-scaling factor from [c](9).
c construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,rtrs,EU,EURS)
      CALL gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,rtRS,EP,EPRS)
      CALL gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm-
     1((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c      
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
c                                                                    c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END
