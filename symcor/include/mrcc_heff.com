      integer ncases
      parameter (ncases=10)
      character*8 ndimheff(ncases)
      character*8 nameheff(ncases)
      character*8 namenvec(ncases)
      character*8 namehvec(ncases)
      character*8 namehdiab(ncases)
      character*8 nametran(ncases)
      character*8 nametmom(ncases)
      character*8 refheff(ncases)
      character*8 reftmom(ncases)
      character*8 pntheff(ncases)
      character*8 pnttmom(ncases)
c
      data ndimheff /'NEFF_IP2', 'NEFF_EA2', 'NEFF_EE1', 'NEFF_EE3',
     $                'NEFFDIP1', 'NEFFDIP3', 'NEFFDEA1', 'NEFFDEA3',
     $                'NEFFDEE1', 'NEFFDEE3' /
      data nameheff /'HEFF_IP2', 'HEFF_EA2', 'HEFF_EE1', 'HEFF_EE3',
     $                'HEFFDIP1', 'HEFFDIP3', 'HEFFDEA1', 'HEFFDEA3',
     $                'HEFFDEE1', 'HEFFDEE3' /
      data namenvec /'NVEC_IP2', 'NVEC_EA2', 'NVEC_EE1', 'NVEC_EE3',
     $     'NVECDIP1', 'NVECIP3', 'NVECDEA1', 'NVECDEA3',
     $     'NVECDEE1', 'NVECDEE3' /
      data namehvec /'HVEC_IP2', 'HVEC_EA2', 'HVEC_EE1', 'HVEC_EE3',
     $     'HVECDIP1', 'HVECDIP3', 'HVECDEA1', 'HVECDEA3',
     $     'HVECDEE1', 'HVECDEE3' /
c     
      data namehdiab /'HDIA_IP2', 'HDIA_EA2', 'HDIA_EE1', 'HDIA_EE3',
     $     'HDIADIP1', 'HDIADIP3', 'HDIADEA1', 'HDIADEA3',
     $     'HDIADEE1', 'HDIADEE3' /
      data nametran /'TRAN_IP2', 'TRAN_EA2', 'TRAN_EE1', 'TRAN_EE3',
     $     'TRANDIP1', 'TRANDIP3', 'TRANDEA1', 'TRANDEA3',
     $     'TRANDEE1', 'TRANDEE3' /
      data refheff /'HREF_IP2', 'HREF_EA2', 'HREF_EE1', 'HREF_EE3',
     $                'HREFDIP1', 'HREFDIP3', 'HREFDEA1', 'HREFDEA3',
     $                'HREFDEE1', 'HREFDEE3' /
      data pntheff /'HPNT_IP2', 'HPNT_EA2', 'HPNT_EE1', 'HPNT_EE3',
     $                'HPNTDIP1', 'HPNTDIP3', 'HPNTDEA1', 'HPNTDEA3',
     $                'HPNTDEE1', 'HPNTDEE3' /
      data nametmom /'TMOM_IP2', 'TMOM_EA2', 'TMOM_EE1', 'TMOM_EE3',
     $                'TMOMDIP1', 'TMOMDIP3', 'TMOMDEA1', 'TMOMDEA3',
     $                'TMOMDEE1', 'TMOMDEE3' /
      data reftmom /'TREF_IP2', 'TREF_EA2', 'TREF_EE1', 'TREF_EE3',
     $                'TREFDIP1', 'TREFDIP3', 'TREFDEA1', 'TREFDEA3',
     $                'TREFDEE1', 'TREFDEE3' /
      data pnttmom /'TPNT_IP2', 'TPNT_EA2', 'TPNT_EE1', 'TPNT_EE3',
     $                'TPNTDIP1', 'TPNTDIP3', 'TPNTDEA1', 'TPNTDEA3',
     $                'TPNTDEE1', 'TPNTDEE3' /
c
c Effective hamiltonians could be stored on JOBARC for a variety of methods.
c IP / EA / EE / DIP / DEA /DEE indicate the IP/EA-eomcc methods or the STEOM variants. DEE is eomee
c The last digit 1, 2 or 3 indicates the multiplicity of the manifold.
c
c The mrcc program writes the # of states included (ndimheff), and it writes the
c records 'nameheff' (ndim*ndim) dimensional and a transition moment record nametmom,
c of dimension 3*ndim. they contain the effective hamiltonian and the tmoms in 
c the diabatic basis.
c
c The rest of these records is internal to symcor. The reference values (refheff, reftmom)
c need to be
c stored and the finite difference records for all geometry points (pntheff, pnttmom)
c
c finally we will calculate heffgrd, tmomgrd, and possibly hefffcm and tmomfcm later on.
c These will be written to files HGRD_method and TGRD_method etc.
c



