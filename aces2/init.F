
#include "flags.h"

      subroutine init
      implicit none

#include "aces2.h"
      integer linblnk, iTmp
      character arg*20, alaska_inp*11
      integer i, f_iargc, num_args,length
      logical file_exists, gen_quartic

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

#ifdef _DEBUG
      print *, '@ACES2: INITIALIZATION DUMP'
#endif

c ---------------------------------------------------------------------------
c Analyze the IFLAGS arrays first -------------------------------------------
c ---------------------------------------------------------------------------
c                  = 
c .TRUE. if INIT_HESSIAN=MOPAC
      mopac_guess = (iflags2(h_IFLAGS2_init_hessian).eq.2)
#ifdef _DEBUG
      print *, '  mopac_guess = ', mopac_guess
#endif

c .TRUE. if DIRECT=ON and INTEGRALS!=GAMESS
      direct=((iflags2(h_IFLAGS2_direct).eq.1).and.
     &        (iflags(h_IFLAGS_integrals).ne.5))
#ifdef _DEBUG
      print *, '  direct = ', direct
#endif

c .TRUE. if CALC=SCF and EXCITE=NONE
      plain_scf=((iflags(2).eq.0).and.(iflags(87).lt.3))
#ifdef _DEBUG
      print *, '  plain_scf = ', plain_scf
#endif

c .TRUE. if GUESS=NDDO
      nddo_guess = (iflags(h_IFLAGS_guess).eq.2)
#ifdef _DEBUG
      print *, '  nddo_guess = ', nddo_guesS
#endif

c .TRUE. if SCF_TYPE=HF
      hf_scf=(iflags2(153).eq.0)
      ks_scf=(iflags2(153).eq.1)
      hf_dft=(iflags2(153).eq.2)
      sl_oep=(iflags2(153).ge.3)
      dkh   =(iflags2(167).gt.0)
      contract=(iflags2(168).eq.0)
#ifdef _DEBUG
      print *, '  hf_scf = ', hf_scf
      print *, '  ks_scf = ', ks_scf
      print *, '  hf_dft = ', hf_dft
      print *, '  sl_ope = ', sl_oep
      print *, '  dkh    = ', dkh
      print *, 'contract = ', contract 
#endif

c .TRUE. if DIRECT=ON, INTEGRALS=GAMESS, and CALC=MBPT(2)
      dirmp2=((iflags2(h_IFLAGS2_direct).eq.1).and.
     &        (iflags(h_IFLAGS_integrals).eq.5).and.
     &        (iflags(h_IFLAGS_calc).eq.1))
#ifdef _DEBUG
      print *, '  dirmp2 = ', dirmp2
#endif

      fno=(iflags2(144).gt.0)
#ifdef _DEBUG
      print *, '  fno = ', fno
#endif

      geom_opt=(iflags2(h_IFLAGS2_geom_opt).ne.0)
#ifdef _DEBUG
      print *, '  geom_opt = ', geom_opt
#endif

      raman=(iflags2(151).eq.1)
      vib_specs=(iflags(54).ne.0)
#ifdef _DEBUG
      print *, '  raman = ', raman
      print *, '  vib_specs = ', vib_specs
#endif

c Analytical_gradient=.TRUE. if GRAD_CALC=ANALYTICAL
c Numerical_gradient =.TRUE. if GRAD_CALC=NUMERICAL

      analytical_gradient=(iflags2(138).eq.1)
      numerical_gradient =(iflags2(138).eq.2) 
#ifdef _DEBUG
      print *, '  analytical_gradient = ', analytical_gradient
      print *, '   numerical_gradient = ', numerical_gradient 
#endif
c
c .TRUE. if the DERIVATIVE_LEVEL=SECOND
      Analytical_hessian = (iflags(h_IFLAGS_deriv_lev) .eq. 2) .and.
     &                     .NOT. hf_dft
#ifdef _DEBUG
      print *, 'Analytical_hessian   = ', Analytical_hessian
#endif
c
c In general, this is wrong since ACES3 is not just mrcc code.
      mrcc=(iflags2(132).eq.3)
#ifdef _DEBUG
      print *, '  mrcc = ', mrcc
#endif
c
cmn include procedure to do multiple points on a surface
cmn  for the purpose of calculating quartic force constants.
cmn  In MRCC branch a number of such options are available. Here only one.
c
      gen_quartic = iflags2(h_IFLAGS2_gen_quartic).eq.1
     $     .and. .not. mrcc
#ifdef _DEBUG
      print *, '  gen_quartic = ', gen_quartic
#endif
c
c Single Point calculations: Energy;gradients;hessians
c There is no direct keyword that we can acces to choose single point.
c The single point is true if not a geometry optimization, finite difference 
c frequency or finite difference property (energy gradient included) calculation. 
c
      Single_point          = .NOT. (Geom_opt .OR. Vib_specs .OR. 
     &                               Numerical_gradient .OR.
     $                                  gen_quartic) 
      Single_point_energy   = (Single_point .AND. (.NOT. 
     &                         (Analytical_gradient .OR. 
     &                          Numerical_gradient)))
#ifdef _DEBUG
      Print *, 'Single_point = ',          Single_point
      Print *, 'Single_point_energy = ',   Single_point_energy
#endif
c                
c First_order_props; .TRUE. if PROPS = FIRST_ORDER (1)
c 
      First_order_props = (Iflags(h_IFLAGS_props) .EQ. 1)
#ifdef _DEBUG
      print *, 'First_order_props = ', First_order_props
#endif
c
c Higher_order_props; .TRUE. if PROPS = NMR, EOM_NLO, J_FC,....(> 1)
c
      Higher_order_props = (Iflags(h_IFLAGS_props) .GT. 1)
#ifdef _DEBUG_LVL0
      print *, 'Higher_order_props = ', Higher_order_props
#endif
c
      Single_point_gradient = (Analytical_gradient .OR.
     &                         Numerical_gradient) .AND. 
     &                                             .NOT. 
     &                        (First_order_props   .OR. 
     &                         Higher_order_props  .OR.
     &                         Geom_opt            .OR. 
     &                         Vib_specs)
#ifdef _DEBUG 
      Print *, 'Single_point_gradient = ', Single_point_gradient
#endif
c
c NMR_SHIFTS; .TRUE. if PROPS = NMR
c
      NMR_SHIFTS = ((Iflags(h_IFLAGS_props) .EQ. 3) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ.  4) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ.  5) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ.  6) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ. 12))

#ifdef _DEBUG
      print *, 'NMR_SHIFTS = ', NMR_SHIFTS
#endif
c
c NMR_SPNSPN; .TRUE. if PROPS = J_FC,J_SD,J_PSO,J_DS), JSC_ALL
c
      NMR_SPNSPN = ((Iflags(h_IFLAGS_props) .EQ. 8) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ.  9) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ. 10) .OR.
     &             (Iflags(h_IFLAGS_props) .EQ. 13))
#ifdef _DEBUG
      print *, 'NMR_SPNSPN = ', NMR_SPNSPN
#endif
c
c NLO_PROPS; .TRUE. if PROPS = EOM_NLO
c
      NLO_PROPS = (Iflags(h_IFLAGS_props) .EQ. 11) .OR. 
     &            (Higher_order_props)
#ifdef _DEBUG
      print *, 'NLO_PROPS = ', NLO_PROPS
#endif
c
c TDHF; .TRUE. if PROPS = TDHF
c
      TDHF = (Iflags2(h_IFLAGS2_tdhf) .NE. 0)
#ifdef _DEBUG_LVL0
      print *, 'TDHF = ', TDHF
#endif
      Cmpt_props_numrcl = (Analytical_gradient .AND. Raman 
     &                                        .AND. Vib_specs)
     &                .OR. (.NOT. Analytical_gradient   .AND. 
     &                           (Single_point_gradient .OR.
     &                            Higher_order_props))
c
#ifdef _DEBUG
      print *, 'PROPS_NUMRCL = ', Cmpt_props_numrcl 
#endif
c bHyper;bMolden; TRUE if EXTERNAL = HYPERCHEM, MOLDEN
c
      bHyper  = (iflags2(h_IFLAGS2_external) .eq. 1)
      bMolden = (iflags2(h_IFLAGS2_external) .eq. 2)
#ifdef _DEBUG
      print *, 'bHyper  = ', bHyper
      print *, 'bMolden = ', bMolden
#endif
c
c blCCSDT; TRUE if CALCLEVEL = ACCSD(T)
c
      blCCSDT = (iflags(h_IFLAGS_calc).eq.42)
#ifdef _DEBUG
      print *, 'blCCSDT  = ', blCCSDT
#endif
c
c bGExtrap;bEExtrap;bCExtrap TRUE if EXTRAPOLATE = GRADIENT, ENERGY, COMBO
c
      bGExtrap = (iflags2(h_IFLAGS2_extrap) .eq .1)
      bEExtrap = (iflags2(h_IFLAGS2_extrap) .eq .2)
      bCExtrap = (iflags2(h_IFLAGS2_extrap) .eq .3)
#ifdef _DEBUG
      print *, 'bGExtrap  = ', bGExtrap
      print *, 'bEExtrap  = ', bEExtrap
      print *, 'bCExtrap  = ', bCExtrap
#endif
c ---------------------------------------------------------------------------
c Analyze the command-line arguments next -----------------------------------
c ---------------------------------------------------------------------------

c set the defaults and override them with the CL.
      sewswitch=.FALSE.
      integral_package='xvmol'
      if (hf_scf) then
         if ((iflags2(h_IFLAGS2_direct).eq.1).and.
     &       (iflags(h_IFLAGS_integrals).eq.5).and.
     &       (iflags(h_IFLAGS_calc).eq.0)) then
            der_integral_package='xscfgrd'
         else
            der_integral_package='xvdint'
         end if
      else if (ks_scf .or. hf_dft) then
         der_integral_package='xvdint'
      end if

      num_args = f_iargc()

c BEGIN GETARG LOOP
      do i = 1, num_args
         call f_getarg(i,arg)

c SEWARD
         if (arg(1:1).eq.'s') then
            sewswitch=.TRUE.
            integral_package= 'molcas run seward MOLCAS.INP'
         end if

c ALASKA
         if (arg(1:1).eq.'a') then
            inquire(file='ALASKA.INP',exist=file_exists)
            if (file_exists) then
               alaska_inp='ALASKA.INP'
            else
               alaska_inp='MOLCAS.INP'
            endif
            der_integral_package= 'molcas run alaska '//alaska_inp
         end if

c END GETARG LOOP
      end do

c verify consistency
      if ((.not.sewswitch).and.(iflags(56).eq.4)) then
         write(*,*) '@ACES2: VMOL cannot generate MOLCAS files.'
         write(*,*) '        Try running `xaces2 s` instead.'
         call c_exit(1)
      end if
      if (dirmp2.and.analytical_gradient.and.(geom_opt.or.vib_specs))
     &   then
         write(*,*) '@ACES2: There are no gradients for direct MBPT(2).'
         call c_exit(1)
      end if
      if (dirmp2.and.iflags(h_IFLAGS_hfstab).ne.0) then
         write(*,*) '@ACES2: HFSTAB does not work for direct MBPT(2).'
         call c_exit(1)
      end if
      if ( ((.not.plain_scf).and.(.not.hf_scf)).and.
     &     (geom_opt.or.vib_specs.or.raman)             ) then
         write(*,*)
         write(*,*) '@ACES2: WARNING! post-scf energy gradients ',
     &              'using KS-SCF orbitals are not fully tested.'
 1       write(*,*)
      end if

#ifdef _DEBUG
      iTmp = linblnk(integral_package)
      print *, '  integral command is "',integral_package(1:iTmp),'"'
      iTmp = linblnk(der_integral_package)
      print *, '  derivative-integral command is "',
     &         der_integral_package(1:iTmp),'"'
#endif

#ifdef _DEBUG
      print *, '@ACES2: END INITIALIZATION DUMP'
#endif

      return
      end

