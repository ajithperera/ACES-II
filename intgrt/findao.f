         subroutine findao(xcoeff,
     &              valao,valgradao,
     &           zvalprim,zgradprim,
     &           iradpt,ncount)

      implicit none





















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler




















































































































































cYAU - ACES3 stuff . . . we hope - #include <aces.par>











c ***NOTE*** This is a genuine (though not serious) limit on what Aces3 can do.
c     12 => s,p,d,f,g,h,i,j,k,l,m,n
      integer maxangshell
      parameter (maxangshell=12)







      integer           iuhf
      common /iuhf_com/ iuhf
      save   /iuhf_com/







c This contains flags that are set in the INTGRT namelist.  See the
c file initintgrt.F for a description of each of them.
c
c The following are exceptions:
c    int_ks           : .true. if we are doing Kohn-Sham
c    int_ks_finaliter : .true. if this is the final iteration of KS
c    int_ks_exch      : which potential to use to calculate exchange
c    int_ks_corr      : which potential to use to calculate correlation
c    int_kspot        : which hybrid functional to use to calculate potential
c                       (if equal to fun_special, use int_ks_exch and
c                       int_ks_corr)
c    int_dft_fun      : which functional to use with any SCF density
c                       Added and modified by Stan Ivanov
c                       (if fun_special the functional is user-defined
c                        if fun_hyb_name then use hybrid functional) 
c    int_printlev     : 0 if we are doing a dft calculation, 1 if we
c                       are doing the final iteration of a KS calculation,
c                       2 if we are doing a KS iteration.
c These are set in the calling routines, NOT in the namelist.
c                     : Additions by S. Ivanov
c     num_acc_ks      : .true.  if numerical accelerator is used for KS 
c                        Default is .true.
c     ks_exact_ex     : .true. if exact LOCAL exchange is used for KS
c                        Deafult is .false.
c     int_tdks        : .true. if time-dependent KS calculation is
c                        requested
c                        Default is .false.
c     int_ks_scf      : .true. if the actual KS SCF energy is being
c                        calculated and printed out. Default is .false.
c
      integer int_numradpts,int_radtyp,int_partpoly,int_radscal,
     &    int_parttyp,int_fuzzyiter,int_defenegrid,int_defenetype,
     &    int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem

      logical int_ks,num_acc_ks,ks_exact_ex,int_tdks,int_ks_scf,
     &        int_ks_finaliter

      double precision
     &    int_radlimit,coef_pot_nonlocal

      common /intgrtflags/  int_numradpts,int_radtyp,int_partpoly,
     &    int_radscal,int_parttyp,int_fuzzyiter,int_defenegrid,
     &    int_defenetype,int_defpotgrid,int_defpottype,int_kspot,
     &    int_ks_exch,int_ks_corr,int_dft_fun,

     &    int_printlev,
     &    int_printscf,int_printint,int_printsize,int_printatom,
     &    int_ks_finaliter,int_printmos,int_printocc,
     &    potradpts, numauxbas,int_ksmem
c

c  prakash added int_ksmem to the common block
      common /intgrtflagsd/ int_radlimit,coef_pot_nonlocal
      common /intgrtflagsl/ int_ks,num_acc_ks,ks_exact_ex,int_tdks,
     &                      int_ks_scf

      save /intgrtflags/
      save /intgrtflagsl/
      save /intgrtflagsd/

c The following are parameters used in the namelist

      integer int_prt_never,int_prt_dft,int_prt_ks,int_prt_always
      parameter (int_prt_never   =1)
      parameter (int_prt_dft     =2)
      parameter (int_prt_ks      =3)
      parameter (int_prt_always  =4)

      integer int_radtyp_handy,int_radtyp_gl
      parameter (int_radtyp_handy=1)
      parameter (int_radtyp_gl   =2)

      integer int_partpoly_equal,int_partpoly_bsrad,
     &    int_partpoly_dynamic
      parameter (int_partpoly_equal  =1)
      parameter (int_partpoly_bsrad  =2)
      parameter (int_partpoly_dynamic=3)

      integer int_radscal_none,int_radscal_slater
      parameter (int_radscal_none  =1)
      parameter (int_radscal_slater=2)

      integer int_parttyp_rigid,int_parttyp_fuzzy
      parameter (int_parttyp_rigid=1)
      parameter (int_parttyp_fuzzy=2)

      integer int_gridtype_leb
      parameter (int_gridtype_leb=1)







c This common block contains the molecule and basis set information
c read in from the JOBARC file.  Since much of this information is
c used in a large number of modules, and since most of the information
c is relatively small compared to the other things held in memory,
c a large percentage of the data stored in the JOBARC file is stored
c here, even though some modules will not use all of it.

c   maxangshell - The maximum number of angular momentum shells.  Since this
c                 is used VERY infrequently, set it high enough to never
c                 cause a problem.
c   spinc(2)    - The characters 'A' and 'B' (useful for alpha/beta labels)
c   natoms      - The number of atoms in the Z-matrix (including X/GH).  After
c                 remove is called, natoms becomes equivalent to nrealatm.
c   natomsx     - The number of atoms in the Z-matrix (including X/GH).  This
c                 does not change.
c   nrealatm    - The number of atoms in the Z-matrix (including GH).
c   naobasfn    - The number of AOs in the basis
c   nbastot     - The number of symmetry adapted orbitals in the basis (the AO
c                 basis may be larger than the SO basis if spherical orbitals
c                 are used since harmonic contaminants are deleted)
c   linear      - 1 if the molecule is linear
c   orientmt    - 3x3 matrix which relates the computational and canonical
c                 orientations
c   nucrep      - Nuclear repulsion energy in a.u.
c   nmproton    - Number of protons in the molecule.
c
c   compptgp    - Point group
c   fullptgp    -
c   compordr    - Order of the point group
c   fullordr    -
c   compnirr    - Number of irreps in the point group
c   fullnirr    -
c   compnorb    - Number of unique atoms (orbits) in the point group
c   fullnorb    -
c   c1symmet    - 1 if the molecule is C1 symmetry
c   nirrep      - The same as compnirr (since nirrep is used so commonly,
c                 this is included for conveniance)
c                 ***NOTE*** nirrep is read in twice and is stored in /sym/
c                            so it is not actually included here
c   totprim     - Total number of primitive functions in the molecule
c   maxshlprim  - Largest number of primitives in a single shell
c   maxshlao    - Largest number of AOs in a single shell
c   maxshlorb   - Largest number of primitive orbitals (primitive functions
c                 times the number of AOs) in a single shell
c   maxangmom   - Largest angular momentum for any atom
c   maxshell    - Larges number of angular momentum shells for any atom
c   noccorb(2)  - The number of alpha and beta occupied orbitals
c   nvrtorb(2)  - The number of alpha and beta virtual orbitals

c The parameter maxorbit is needed because of how dynamic memory is used.
c Two runs of the program are needed.  The first to calculate memory usage,
c the second to use it.  In order to calculate totprim, we have to know the
c orbit population vector (the number of each type of atom).  BUT, this is
c stored in dynamic memory since we do not know how long this vector is.
c In the future, joda or vmol will write this information to JOBARC, and
c this problem will disappear.  In the meantime, we have to introduce a
c genuine limit on the size of the molecule.  It may have no more than
c maxorbit sets of unique atoms.  This limit is ONLY used in the subroutine
c basis, so it probably will disappear when the information in the MOL file
c is put in JOBARC.
c    maxorbit   - the number of symmetry unique atoms

c The following are pointers to real arrays
c
c   zatommass(natoms)  - Atomic mass of all atoms (X=0.0, GH=100.0)
c   zcoord(3,natoms)   - Coordinates of all atoms (computational orientation)
c   zalpha(totprim)    - The alpha for each primitive function
c   zprimcoef(totprim,naobasfn)
c                      - The primitive to AO coefficients
c
c The following are pointers to integer arrays
c
c   patomchrg(natoms)  - Atomic number of all atoms (X=0, GH=110)
c   pfullclss(fullordr)- Class type vector
c   pcompclss(compordr)-
c   pfullpopv(natoms)  - Number of atoms in each orbit
c   pcomppopv(natoms)  -
c   pfullmemb(natoms)  - Atoms sorted by point group orbits
c   pcompmemb(natoms)  -
c   pnprimatom(natoms) - Number of primitive functions for each atom
c   pnshellatom(natoms)- Number of different angular momentum shells for each
c                        atom (takes on values of 1,4,9,16, etc.)
c   pnangatom(natoms)  - The number of different angular momentum for each
c                        atom (takes on values of 1,2,3,4, etc.)
c   pnaoatom(natoms)   - Number of AOs for each atom
c   pnshellprim(maxshell,natoms)
c                      - The number of primitive functions in each shell
c                        of each atom
c   pnshellao(maxshell,natoms)
c                      - The number of AOs in each shell of each atom
c   pprimoff(maxshell,natoms)
c   paooff(maxshell,natoms)
c                      - The primcoef matrix is a block diagonal matrix.
c                        Each shell of each atom has a block.  If you have
c                        a list of all primitive functions, pprimoff(ishell,
c                        iatom) tells the location of the first primitive
c                        function in the block (ishell,iatom) and paooff
c                        contains similar information for the AOs.
c
c ***NOTE***  Because joda stores pfullpopv/pcomppopv as size natoms, we
c             do to, but they should be of size fullnorb/compnorb.  The
c             first ones have real values.  The remaining ones are 0.

      double precision orientmt(3,3),nucrep
      integer natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb(2),
     &    nvrtorb(2),maxshlao,maxangmom,natomsx
      integer patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,
     &    zalpha,zprimcoef,pprimoff,paooff,pnangatom

      common /mol_com/ orientmt,nucrep,
     &    natoms,nrealatm,naobasfn,nbastot,linear,compnirr,
     &    fullnirr,compnorb,fullnorb,compordr,fullordr,nmproton,
     &    c1symmet,totprim,maxshlprim,maxshlorb,maxshell,noccorb,
     &    nvrtorb,maxshlao,maxangmom,natomsx,
     &    patomchrg,zatommass,zcoord,pfullclss,pcompclss,
     &    pfullpopv,pcomppopv,pfullmemb,pcompmemb,pnprimatom,
     &    pnshellatom,pnaoatom,pnshellprim,pnshellao,zalpha,
     &    zprimcoef,pprimoff,paooff,pnangatom
      save   /mol_com/

      character*4 compptgp,fullptgp
      character*1 spinc(2)
      common /molc_com/ compptgp,fullptgp,spinc
      save   /molc_com/







c This contains information about each of the possible grids for
c performing the numerical integration.

c###########################################################################
c MISC
c###########################################################################
c maxgrdatm : The largest atomic number for which the Slater and Bragg-Slater
c             atomic size is known.
c atmrad    : Atomic size using Slater's' rules for the radial integration
c xbsl      : The Bragg-Slater radii (one for each atom)
c             ***NOTE*** This is a genuine constraint.  Only atoms smaller
c             then this (currently 86) may be calculated.
c numangfct : number of angular momentum functions
c minpt     : The number of points used in the integration over
c             interatomic paths to find the minimum density point
c             between two atoms
c pangfct   : a pointer to the array of x,y,z angular momentum for each
c             angular momentum function

      integer maxgrdatm
      parameter (maxgrdatm=86)
      integer minpt
      parameter (minpt=100)

      integer numangfct

      double precision
     &    atmrad(maxgrdatm),xbsl(maxgrdatm),
     &     TA(maxgrdatm),multiEX(maxgrdatm)
      common /grid/  numangfct
      save /grid/

      common /gridd/ atmrad,xbsl,TA,multiEX
      save /gridd/

      integer pangfct
      common /gridp/ pangfct
      save /gridp/

c###########################################################################
c RADGRD file
c###########################################################################
c maxanggrid: The maximum number of different angular grids which can be
c             used in any given calculation.
c             ***NOTE*** This is a genuine constraint, but it must be used
c             since we must be able to keep a record of which angular grids
c             are used (before we have any allocated memory) since we have
c             to know how many grids are used in order to determine how much
c             memory to allocate.  This is set high enough it should never
c             be a problem.
c gridlist  : A list of all grids used (see the comment on maxanggrid).  It
c             is of dimension (maxanggrid,3) to keep track of the type and
c             subtype, and the number of times each grid is used.  The type
c             refers to how the grid is arrived at.
c                1 : Lebedev
c                2 : ?
c             The subtype refers to the degree of the grid of this type.
c numgrid   : The number of different angular grids used in the calculation.
c maxangpts : The maximum number of points in any of the angular grids used.
c maxanggrd : The grid with the maximum number of angular points.
c numradpts : The number of different radial points.
c ntotrad   : The total number of points in all angular grids at all radial
c             points (i.e. the entire integration grid)
c
c iradint   : determines if the Handy method (1) or Gauss-Legendre (2)
c             radial integration is used
c autosiz   : A flag which sets whether the polyhedra are (1) equally
c             sized, (2) sized according to Bragg-Slater radii or
c             (3) automatically sized according to the minimums in
c             density.
c slater    : A flag which determines whether (0) Slater's' rules are
c             used to determine the atomic size and scale the radial
c             integration or (1) no scaling is used.
c rigid     : A flag which determines whether rigid (0) or fuzzy
c             partitioning is used.
c nitr      : The number of iterations of the equations which create the
c             'fuzzy' boundary.

      integer maxanggrid
      parameter (maxanggrid=1000)

      integer gridlist(maxanggrid,3),numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad

      common /radgrd/  gridlist,numgrid,maxangpts,numradpts,
     &    iradint,autosiz,slater,rigid,nitr,maxanggrd,ntotrad
      save /radgrd/

c Memory pointers
c
c pradgrid(numradpts) : The angular grid to use at each radial point.
c pgrdangpts(numgrid) : The number of angular points in each grid.
c zgridxyz(3,maxangpts,numgrid)
c                     : The x,y,z coordinate of each angular point in each grid
c zgridwt(maxangpts,numgrid)
c                     : The weight at each point.
c pintegaxis(natoms,3): Contains information about how much of each axis to
c                       integrate over.  If integaxis(iatom,i) is set to i,
c                       integrate only over the positive half of the i^th
c                       axis.  Otherwise, integrate over the entire axis.

      integer pgrdangpts,pradgrid,zgridxyz,zgridwt,pintegaxis

      common /radgrdp/ pgrdangpts,pradgrid,zgridxyz,zgridwt,
     &    pintegaxis
      save /radgrdp/

c###########################################################################
c Old stuff
c###########################################################################

c polist    : Contains an ordered list of unique atoms
c zatmvc    : The x, y, and z distance between each pair of atoms.
c zrij      : The distance between each pair of atoms.
c zatmpth   : The cartesian coordinates for the path integration between
c              all atom pairs
c zptdis    : The distance from atom i to a point along the path between
c              atoms i and j
c zprsqrd   : The distance squared from each atom to a point along all the
c              paths between all the atoms
c zpthpt    : The cartesian coordinates with respect to each atom for
c              the points along all the paths between all the atoms
c zbslrd    : The Bragg-Slater radii.
c zaij      : Surface shifting parameter dependent on the distance between
c               pairs of atoms.

      integer polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,zptdis,
     &    zprsqrd,zpthpt,zbslrd,zaij
      common /gridold/ polist,ixx,iyy,izz,zatmvc,zrij,zatmpth,
     &    zptdis,zprsqrd,zpthpt,zbslrd,zaij
      save /gridold/







c This commonblock contains values for the functionals involved in the
c numerical integration either for DFT, plotting or fitting

c totele    : The total number of electrons
c etotekin  : The total kinetic energy
c etotenatr : The total nuclear-electron attraction energy
c etottf    : The total Thomas Fermi kinetic energy
c etotw     : The total Weizacker kinetic energy
c xldax     : The total LDA exchange energy
c becke     : The total Becke exchange energy
c lda       : The total LDA correlation energy
c xlyp      : The total LYP correlation energy
c icntr     : The integration center
c idns      : a flag =0 for SCF orbitals and =1 for natural orbitals

      integer icntr,idns
      double precision
     &    totele,etotkin,etotnatr,etottf,etotw,xldax,becke,
     &    lda,xlyp

      common /int_com/  icntr,idns
      save   /int_com/
      common /intr_com/ totele,etotkin,etotnatr,etottf,etotw,xldax,
     &                  becke,lda,xlyp
      save   /intr_com/

c array pointers

c zpcoeff(2): alpha/beta MO to primitive function transformation matrix
c zxocc     : alpha/beta orbital occupation

      integer
     &    zpcoeff(2),zxocc
      common /molecp/
     &    zpcoeff,zxocc
      save /molecp/







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
      double precision
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

      double precision
     &    fun_exch(fun_num_exch),fun_corr(fun_num_corr),
     &    fun_hyb(fun_num_hyb),
     &    tot_exch(fun_num_exch),tot_corr(fun_num_corr),
     &    tot_hyb(fun_num_hyb),
     &    ene_exch(fun_num_exch),ene_corr(fun_num_corr),
     &    ene_hyb(fun_num_hyb),
     &    coef_exch(fun_num_exch),coef_corr(fun_num_corr),
     &    coef_pot_exch(fun_num_exch),coef_pot_corr(fun_num_corr)

      double precision vxc_ksalpha,vxc_ksbeta

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
      double precision
     &    b3lypa,b3lypb,b3lypc
      parameter (b3lypa=0.20d0)
      parameter (b3lypb=0.72d0)
      parameter (b3lypc=0.81d0)





        integer iradpt,ncount,xyz

      double precision
     &    zvalprim(totprim,maxangpts,int_numradpts,
     &    ncount),
     &    zgradprim(totprim,maxangpts,int_numradpts,
     &    ncount,3),
     &    valao(naobasfn,maxangpts,int_numradpts,
     &    ncount),
     & valgradao(naobasfn,maxangpts,int_numradpts,
     &     ncount,3),
     &    xcoeff(naobasfn*totprim)

       call callstack_push('FINDAO')
      call xgemm('t','n',naobasfn,maxangpts,totprim,
     &       1.d0,xcoeff,totprim,
     &        zvalprim(1,1,iradpt,icntr),totprim,
     &      0.d0,valao(1,1,iradpt,icntr),naobasfn)

       do xyz=1,3
      call xgemm('t','n',naobasfn,maxangpts,totprim,
     &                1.d0,xcoeff,totprim,
     &    zgradprim(1,1,iradpt,icntr,xyz),totprim,
     &            0.d0,valgradao(1,1,iradpt,
     &                    icntr,xyz),
     &                             naobasfn)
       end do

      call callstack_pop
      return

      end


 
