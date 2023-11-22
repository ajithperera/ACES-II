
#include "flags.h"

      Subroutine Job_control
c
      Logical Converged, Both_gopt_vfreq, gen_quartic, NORMCO_exsist
      Character*4 Calc_level 
c
#include "aces2.h"
#include "parallel_aces.com"
c
c This is a replacement for the old branch.F. The primary reason for replacing
c branch.F was that it was becoming unmaintainable. Also, the inflexibility
c which hinder adding new options was a handicapp too. 02/2006, Ajith Perera.

      integer icore, is, ius
      character*79 szHyper, szMolden
c
      szHyper   = 'xa2proc hyperchem'
      szMolden  = 'xa2proc molden'
c
c Single point calculations. The possiblities are single point energy (ground & excited), 
c first order properties and higher order properties (the gradients are considered as a
c property). 
c 
      If (Single_point) Then
c
cmn Inserted code to run finite field properties
c
         if (iflags2(h_IFLAGS2_ffield) .gt. 0) then
c
c check that symmetry is off
c
            if (iflags(60) .ne. 1) then
               write(6,*) ' Switch off Symmetry with FField > 0 '
               call ACES_EXIT
            endif
c     
cmn     run a finite_field calculation and terminate
c     
            call finite_field
            go to 999
c
         elseif (iflags2(h_IFLAGS2_ffield) .lt. 0) then
c
cmn     run a finite_field static multipole calculation and terminate
c     
            call finite_field
            go to 999
c
         endif
c
         If (Single_point_energy) Then
c
c Plain_scf is .TRUE. when CALCLEVEL=.SCF. and EXCITE=TDA.
c
            If (Plain_scf) Then
c
c The HF/KS-SCF stability analysis is considerd as a part of SCF and is included.
c If we have not traped the Scf energy only calculations with DROPMO already, we can trap it
c here.
c
               Call Scf_eneg
c
            Else
c
c Interpret call to Scf as a generation of a reference. When interpreted that way,
c the stability analysis should be part of it. The dropmo processing and fno processing 
c or any other procedure that change the size is done in change_orbital_space.F.
c If there is a valid instabilty that we follow, we need to rerun transformation and 
c and processing, handle that and other preparatory stuff in prep4_post_scf.F. 
C Also, For example, ACCSD(T) we need to change the CALCLEVEL=CCSD(T) internaly to get standard
c CCSD(T) calculations done porior to (we need make program think that it is doing both
c ACCSD(T) and CCSD(T) at the same time). That is hard to acomplish
c unless we introduce a new keyword essentialy playing the same role as CALCLEVEL which
c is unpleasent - what we do is to switch to CCSD(T) internaly when we wants it to behave
c as CCSD(T). These are the kind of preparatory things that we shoud try to do in 
c prep4_post_scf.F (eventualy transformation and processing should go in there as well).
c
c The post_scf_eneg.F perhaps need to be more organized to handle some of the developments 
c that are underway now. I consider direct second-order through xdirmp2 (note that
c direct methods dose not need pre transformations or processing) ACCSD(T) as such. 
c There is some DROPMO logic in Scf_eneg which can should be handle by vtran.
c
               If (Iflags(h_IFLAGS_brueckner)  .EQ. 1 .or.
     +             Iflags2(h_IFLAGS2_oo_method) .EQ. 4) Then
C
                  Call Do_brueckner

CSSS                  Converged = .FALSE.
#ifdef _DEBUG
CSSS      print *, "Entering Brueckner block"
#endif
CSSS                  Do while (.NOT. Converged)
CSSS                     Call Scf_eneg
CSSS                     Call Change_orbital_space
CSSS                     Call Prep4_post_scf
CSSS                     Call Post_scf_eneg 
CSSS                     Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
CSSS                     IF (Itest .EQ. 1) Converged = .TRUE.
CCCC                  Enddo
c
c Microiteration for Brueckner orbitals previously handled runvcc.F which was roughly
c equivalent to post_scf_eneg.F in terms of functionalities. Such procedure can be risky.
c For example, in  previous version users can run Brueckner analytical optimizations without
c knowing that Breuckner analytic gradients were unavilable.
c
               Else
#ifdef _DEBUG_LVL0
      print *, "Entering post-SCF single point energy block"
#endif
c The excitation energy code use second moments integrals computed
c by xvporps. The vprops program needs the "IIII" integral file and 
c any for any post SCF calcualtions which require property integral,
c we must ensure that the property integral code run before the transformation
c (vtran) since vtran automatically deletes the "IIII" file (unless
c specificaly asked it not to do it).
c 
                  Call Scf_eneg

                  If (Iflags(h_IFLAGS_excite) .GT. 2) Call Runit
     &                                               ("xvprops")
                  Call Change_orbital_space

                  If (.NOT. Dirmp2) Call Prep4_post_scf

                  If ((iflags(h_IFLAGS_calc) .eq. 45) .or.
     &                (iflags(h_IFLAGS_calc) .eq. 46)) then
                      Call runoombpt
                  Else
                      Call Post_scf_eneg 
                  Endif
             
               Endif
c
c Endif for plain_scf
c
            Endif
c     
         Else if (First_order_props) Then
c
c Properties can be computed both analyticaly and numericaly. By default numerical
c property calculations are not considered as single point calculations.
c
c The call to "xvpoprs" and "xprops" is to compute the first order property integrals and the
c the properties based on the appropriate density. Note that the first-order properties
c does not require anti, bcktrn and the derivative integrals (argument is set to FALSE by
c .NOT. First_order_prps).
c
            If (Plain_scf) Then
               Call Scf_eneg
               If (First_order_props .and. .Not. Hf_dft) 
     &             Call Runit("xvprops")  
               Call Scf_anlytcl_grad(.NOT. First_order_props)
            Else
               If (Iflags(h_IFLAGS_brueckner) .EQ. 1 .AND.
     &             First_order_props) Then
                   Call do_brueckner
               Else
                   Call Scf_eneg
                   If (First_order_props .and. .Not. Hf_dft) 
     &                 Call Runit("xvprops") 
                   Call Change_correl_space
                   Call Prep4_post_scf
                   Call Post_scf_eneg
                   Call Post_scf_anlytcl_grad(.NOT. First_order_props)
               Endif
            Endif
            IF (First_order_props .and. .not. Hf_dft) 
     &             Call Runit("xprops")
c
         Else if (Higher_order_props) Then
#ifdef _DEBUG
      print *, "Entering higher order props block"
#endif
c
c Analytical second-order properties that we currently capable are NMR shifts, polarizabilities
c and NMR coupling constants and TDHF (polarizabilities and hyperpolarizabilites).
c 
            If (Plain_scf) Then
c
c We have only TDHF, polarizabilities, NMR shifts and coupling constants (Juergen's implementation) that can 
c go into this block
c
#ifdef _DEBUG_LVL0
      print *, "Entering scf second derivative  block"
#endif
               Call Scf_eneg
               Call Scf_higher_order_prps(TDHF, NMR_SHIFTS, NMR_SPNSPN,
     &                                    NLO_PROPS)
            Else
#ifdef _DEBUG
      print *, "Entering post-scf second derivative  block"
#endif
               Call Scf_eneg
               Call Prep4_post_scf
               Call Post_scf_eneg
               Call Post_scf_higher_order_prps(NMR_SHIFTS, NMR_SPNSPN, 
     &                                         NLO_PROPS)
            Endif
c
         Else if (single_point_gradient) Then 
c
            If (Plain_scf) Then
               Call Scf_eneg
               Call Scf_anlytcl_grad(.NOT. First_order_props)
            Else
               Call Scf_eneg
               Call Runit("xvprops") 
               Call Change_correl_space
               Call Prep4_post_scf
               Call Post_scf_eneg
               Call Post_scf_anlytcl_grad(.NOT. First_order_props)
            Endif
c
c Endif for single_point_energy
c
         Endif
C
C Else block for anything other than single point calculations.
c      
      Else 
c
c There are at least 3 main type of calculations that falls into this block. 
c They are geometry optimizations, frequency calculations, finite difference
c calculations for properties. 
c
c All the finite difference calculations have a very simple structure. That is Do while (test
c is pass; all the points are done, geometry optimization converged), compute energies or 
c gradients.  So, in principal we can do all of it in one structure with arguments to 
c control various options, or have simple seperate structures for freqeuncies and optimizations.
c The latter is done for frequncies and optimizations here only for the clarity, and the
c numerical properties are handled in a one structure with controling arguments in order
c to illustrate the prior approach.
c 
c Both SCF and Post SCF is handled by opt_anlytcl, opt_numrcl, freq_anlytcl and freq_numrcl,
c and is controled  by the argument string.
c
         If (iflags(h_IFLAGS_calclevel) .Eq. 0) Then
            Calc_level = "SCF "
         Else
            Calc_level = "PSCF" 
         Endif

         Both_gopt_vfreq = Geom_opt .AND. Vib_specs 
         If (Geom_opt .AND. .NOT. Both_gopt_vfreq) then
c
            If (Analytical_gradient) then
#ifdef _DEBUG
      print *, "Entering anlytical optimization block"
#endif

               Call Opt_anlytcl(Calc_level)
            Else
               Call Opt_numrcl(Calc_level)
            End if
c
c Generate HyperChem input file after optimizing the geometry
c
            If (bHyper) Call runit(szHyper)

         Else If (Vib_specs .AND. .NOT. Raman .AND. 
     &            .NOT. Both_gopt_vfreq) Then
#ifdef _DEBUG
      print *, "Entering the freq. block"
#endif
c
c While it is possible to do both IR and Raman inside freq_anlytcl/numrcl, the Raman
c is done in props_numrcl.F (see the comments below) for the purpose of clarity. Let's
c limit freq_xxx routines to handle properties that involve "energy" derivative
c calculations: force constants, cubic, quartic force fields.
c
            If (Analytical_hessian) then
c
C At the moment only availble for HF-SCF
c
#ifdef _DEBUG
      print *, "Entering freq. exact block"
#endif
               Call Freq_anlytcl(Calc_level)
            Else
c
c The first derivative can be numerical or analytical. The second argument
c signals the analytical first derivative or not. 
c
#ifdef _DEBUG
      print *, "Entering freq. finite diff. block"
#endif
               Call Freq_numrcl(Calc_level, Analytical_gradient)
            End if
c        
c Generate HyperChem input file after calculating vibrational frequencies
c
c            If (bHyper) call runit(szHyper)
c            If (bMolden) call runit(szMolden)
c 
         Else If (Geom_opt .AND. Vib_specs .AND. .NOT. Raman) Then
c
c The frequency calculations directly followed by geometry optimizations.
c
#ifdef _DEBUG
      print *, "Entering vib. and freq. block"
#endif
            If (Analytical_gradient) then
               Call Opt_anlytcl(Calc_level)
               Call Prep4_post_opt_freq 
               Call A2getrec(10, 'JOBARC', 'IFLAGS  ', 100, Iflags)
                Print*, "@JOB-CONTRLO", iflags(h_IFLAGS_deriv_lev)
               If (iflags(h_IFLAGS_deriv_lev) .eq. 2) Then
#ifdef _DEBUG
      print *, "Entering analytical hessian block"
#endif
                  Call Freq_anlytcl(Calc_level)
               Else
                  Call Freq_numrcl(Calc_level, Analytical_gradient)
               Endif
            Else
               Call Opt_numrcl(Calc_level)
               Call Prep4_post_opt_freq 
               Call Freq_numrcl(Calc_level, Analytical_gradient)
            End if
c        
c Generate HyperChem input file after calculating vibrational frequencies
c
c            If (bHyper) call runit(szHyper)
c            If (bMolden) call runit(szMolden)
c
cmn
c add block for gen_quartic. This can follow opt_numrcl procedure
c
         Else if (gen_quartic) then
#ifdef _DEBUG
      print *, "Entering gen_quartic block"
#endif
            call Opt_numrcl(Calc_level)
c
         Else if (Cmpt_props_numrcl) Then
C
CSSS(.NOT. Analytical_gradient .AND. (First_order_props
CSSS     &         .OR. Single_point_gradient .OR. Higher_order_props))
CSSS     &         Then
c
c Numerical calculation of first and higher-order properties, We don't do as
c many numerical property calculations. In principal all the analytical 
c property calcualtions can be done numericaly, but the procress can be 
c tedious (extensively used in debugging analytic derivatives). The 
c numerical energy gradients ("new") and Raman Intensities and depolarization
c ratios are computed. In the course of Raman calculations, IR frequencies and 
c intensities are also become available, however, Raman calculations are limited 
c to SCF and CCSD.
c
            If (iflags(h_IFLAGS_calclevel) .Eq. 0) Then
               Calc_level = "SCF "
            Else
               Calc_level = "PSCF" 
            Endif
c
c Single_point_gradient is the identifier for "energy" derivative. 
c
#ifdef _DEBUG
      print *, "Entering numerical single point gradients block" 
#endif
            Call Props_numrcl(Calc_level, Single_point_gradient, Raman)
c
         Endif
      Endif        

c
c Generate HyperChem/MOLDEN input file after calculating vibrational frequencies
c Lets check whether "NORMCO file is present before we try to do this
c 08/2014. 

      INQUIRE(FILE='NORMCO',EXIST=NORMCO_exsist)
      If (bHyper  .and. NORMCO_exsist) call runit(szHyper)
      If (bMolden .and. NORMCO_exsist) call runit(szMolden)
c
 999  continue
C
      return
      End
