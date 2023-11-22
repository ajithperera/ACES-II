c The following values are used in determining the values of the various
c functionals:
c
c    roa    : The alpha density
c    rob    :
c    ro     : roa+rob
c    roinv  : 1/ro
c    rom    : roa-rob
c    rs     : (3 / 4 pi ro)**(1/3)
c    zeta   : (roa-rob)/ro
c
c    gradx  : d(roa)/dx + d(rob)/dx
c    grady  :
c    gradz  :
c    gradxm : d(roa)/dx - d(rob)/dx
c    gradym :
c    gradzm :
c    gro2   : gradx**2 + grady**2 + gradz**2
c    gro    : sqrt(gro2)
c    gro2a  : [d(roa)/dx]**2 + [d(roa)/dy]**2 + [d(roa)/dz]**2
c    gro2b  :
c
c    hesxx  : d2(roa)/dx2 + d2(rob)/dx2
c    hesyy  :
c    heszz  :
c    hesxy  :
c    hesxz  :
c    hesyz  :
c    hesxxm : d2(roa)/dx2 - d2(rob)/dx2
c    hesyym :
c    heszzm :
c    hesxym :
c    hesxzm :
c    hesyzm :
c    xlap   : hesxx + hesyy + heszz
c    xlapm  : hesxxm + hesyym + heszzm
c    trm1   :
c    trm1m  :
c    trm2   :
c    trm2m  :
c    grdaa  : gradroa.gradroa
c    grdbb  :
c    grdab  :
c    thresh : Parameter which tells what is the cutoff value
c           :  for is for the density, functionals and potetentials.
c           : The flag is cutoff and the default value is 10**(-12).
      M_REAL
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh

c There are several types of exchange and correlation energy functionals
c as well as hybrid schemes.  There are also corresponding potentials
c for use in Kohn-Sham.  Define the number of each type of energy functional
c and parameters to make calling them clearer.
c
c We want to have freedom to sellect potential and functional for
c better flexibility. Usually, the choice of functional determines 
c the poential  to be used in the  KS calculations.
c 
c
c fun_num_exch  : the number of exchange energy functionals (or potentials)
c fun_exch      : the total exchange energy contribution from the current
c                 atom
c tot_exch      : the total exchange energy contribution from all atoms
c coef_exch     : the coefficient of each exchange functional
c ene_exch      : the total energy with this functional
c nam_exch      : the name of the functional
c coef_pot_exch : the coefficient of each exchange component in the KS potential
c coef_pot_corr : the coefficient of each correlation component in the KS potential

      integer fun_num_exch,fun_num_corr,fun_num_hyb
      parameter (fun_num_exch = 5)
      parameter (fun_num_corr = 6)
      parameter (fun_num_hyb  = 1)

      M_REAL
     &    fun_exch(fun_num_exch),fun_corr(fun_num_corr),
     &    fun_hyb(fun_num_hyb),
     &    tot_exch(fun_num_exch),tot_corr(fun_num_corr),
     &    tot_hyb(fun_num_hyb),
     &    ene_exch(fun_num_exch),ene_corr(fun_num_corr),
     &    ene_hyb(fun_num_hyb),
     &    coef_exch(fun_num_exch),coef_corr(fun_num_corr),
     &    coef_pot_exch(fun_num_exch),coef_pot_corr(fun_num_corr)

      M_REAL vxc_ksalpha,vxc_ksbeta

      integer fun_exch_none,fun_corr_none,fun_special,
     &    fun_dft_none,
     &    fun_exch_lda,fun_exch_becke,fun_exch_pbe,fun_exch_pw91,
     &    fun_exch_hf,
     &    fun_corr_vwn,fun_corr_lyp,fun_corr_pbe,
     &    fun_corr_pw91, fun_corr_wl,fun_corr_wi,
     &    fun_hyb_b3lyp

      character*50
     &    nam_exch(fun_num_exch),nam_corr(fun_num_corr),
     &    nam_hyb(fun_num_hyb)
      character*80 nam_func,nam_kspot
      character*10
     &    abb_exch(fun_num_exch),abb_corr(fun_num_corr),
     &    abb_hyb(fun_num_hyb)

      parameter (fun_exch_none = 0)
      parameter (fun_exch_lda  = 1)
      parameter (fun_exch_becke= 2)
      parameter (fun_exch_pbe  = 3)
      parameter (fun_exch_pw91 = 4)
      parameter (fun_exch_hf   = 5)

      parameter (fun_corr_none = 0)
      parameter (fun_corr_vwn  = 1)
      parameter (fun_corr_lyp  = 2)
      parameter (fun_corr_pbe  = 3)
      parameter (fun_corr_pw91 = 4)
      parameter (fun_corr_wl   = 5)
      parameter (fun_corr_wi   = 6)

      parameter (fun_special   = 0)

      parameter (fun_dft_none  = -1)
      parameter (fun_hyb_b3lyp = 1)

      common /dftfunc/ fun_exch,fun_corr,fun_hyb,tot_exch,tot_corr,
     &    tot_hyb,ene_exch,ene_corr,ene_hyb,coef_exch,coef_corr,
     &    roa,rob,ro,rom,rs,zeta,roinv,
     &    gradx,grady,gradz,gradxm,gradym,gradzm,gro,gro2,gro2a,gro2b,
     &    hesxx,hesyy,heszz,hesxy,hesxz,hesyz,
     &    hesxxm,hesyym,heszzm,hesxym,hesxzm,hesyzm,
     &    xlap,xlapm,trm1,trm1m,trm2,trm2m,
     &    grdaa,grdbb,grdab,thresh,vxc_ksalpha,vxc_ksbeta,
     &    coef_pot_exch,coef_pot_corr

      common /dftfuncc/ nam_exch,nam_corr,nam_hyb,nam_func,
     &    abb_exch,abb_corr,abb_hyb,nam_kspot

      save /dftfunc/
      save /dftfuncc/

c Parameters for hybrid methods:
      M_REAL
     &    b3lypa,b3lypb,b3lypc
      parameter (b3lypa=0.20d0)
      parameter (b3lypb=0.72d0)
      parameter (b3lypc=0.81d0)

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End: c
