#include "flags.h"
      Subroutine Tdee_driver(Work,Maxmem,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)
      Integer AAAA_LENGTH_IJKA
      Integer Isympert(10)
      Integer Nondiag_component
      Character*8 Label(6),Work_label
      Double Precision Mubar_00, M_expect 

      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     +        CC2,CC
      Logical Mu_bar,Mubar_t,T_dot,Tranabci
      Logical Source, Target
      Logical SS,SD,DS,DD
      Logical Restart_present,Restart 
      Logical Done 
      Logical Sing 
      Logical Nondiag

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
#include "parcc.com"
#include "info.com"
#include "tdee_intgrt_vars.com"
#include "lanczos_vars.com"

       Common /METH/MBPT2_DUMMY,MBPT3,M4DQ,M4SDQ,M4SDTQ,
     +              CCD_DUMMY,QCISD,CCSD,UCC
       Common /REFTYPE/MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
       Common /DRVHBAR/ SS, SD, DS, DD
       COMMON/LISTDENS/LDENS

       Character*1 DXYZ(3)
       Character*2 QXYZ(6)

       Data Ione, One /1, 1.0D0/ 
       Data DXYZ /"x", "y", "z"/
       Data QXYZ /"xx", "yy", "zz","xy","xz","yz"/

       Sing   = .False.
       Icontl = Iflags(h_IFLAGS_cc_conv) 
       Iconvg = 1
       Ncycle = 0
       Sing   = (Iflags(h_IFLAGS_calc) .gt. 9)

       Mu_bar  = .True.
       T_dot   = .True.
       Mubar_t = .True.
       Source  = .True.
       Target  = .True.
       Restart = .False.
       Done    = .False.

       Call Tdee_set_vars(Work,Iuhf)

C Check whether the restart file is present. If read the integrator time
C data from it. 

       Inquire(File="Restart_log",Exist=Restart_present)
       If (Restart_present) Then
          Call Tdee_read_from_restart(Work,Maxcor,Irrepx,Done)
          If (Done) Return 

C Set the restart flag and rewinde the Restart_log to record restart
C data for the current run

          Restart = .True. 
          Rewind(10)
          Call Tdee_reset_vars_4restart(Component,Irrepx,Naobfns,
     &                                  Nbfns,Ipert,Work_Label,
     &                                  D_pole,Q_pole)
          Goto 100
       Else
          Call Tdee_read_tintgrt_vars(Work,Maxcor,Nondiag_component)
       Endif 

       If (Nondiag_component .Ne. 0) Nondiag = .True.

       If (LZ) Then 
      
       Write(6,*)
       Write(6,"(2a)")"  Asymmetric Lanczos algorithm is used to",
     &                " to generate a sub-space tridiagonal"
       Write(6,"(a)") "  matrix."
       Write(6,"(2a)")"  Eigenvalues and vectors of the tridigonal",
     &                " matrix is used to generate"
       Write(6,"(a)") "  spectra."
       Else 

       If (Side .EQ. 2) Then

       Write(6,*)
       Write(6,"(4a)") "  The left, id/dt <Mbar(t)| = <Mbar(t)Hbar",
     &                 " multiplication is carrired out for "
       If (D_pole)
     & Write(6,"(2x,2a)")DXYZ(component), " dipole perturbation, M(t)."
       If (Q_pole)
     & Write(6,"(2x,3a)")QXYZ(component), " quad. pole perturbation, ",
     &                     "M(t)."
       If (Nondiag) Then
          If (D_pole)
     &    Write(6,"(1x,4a)")" Non-diagonal contributions are evaluated",
     &                      " and the other pertubation is ",
     &                        DXYZ(Nondiag_component),"."
       If (Q_pole)
     &    Write(6,"(1x,4a)")" Non-diagonal contributions are evaluated",
     &                      " and the other pertubation is ",
     &                        QXYZ(Nondiag_component),"."
       Endif 
       Write(6,*)
       Else
       Write(6,*)
       Write(6,"(2a)") "  The right, id/dt |Mbar(t)> = Hbar|Mbar(t)>",
     &               " multiplication is carrired out for "
       If (D_pole)
     & Write(6,"(2x,2a)")DXYZ(component), " dipole perturbation, M(t)."
       If (Q_pole)
     & Write(6,"(2x,3a)")QXYZ(component), " quad. pole perturbation, ",
     &                      "M(t)."
       If (Nondiag) Then
         If (D_pole)
     &    Write(6,"(1x,4a)")" Non-diagonal contributions are evaluated",
     &                      " and the other pertubation is ",
     &                        DXYZ(Nondiag_component),"."
       If (Q_pole)
     &    Write(6,"(1x,4a)")" Non-diagonal contributions are evaluated",
     &                      " and the other pertubation is ",
     &                        QXYZ(Nondiag_component),"."
       Endif 
       Write(6,*)

       Endif 
       Endif 
C
       If (CCSD) Then
          Coulomb = .False.
          Call Parread(iuhf) 
          write(6,*)
         If (ispar) Then
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
          write(6,"(a,a)") ' Perform a parameterized CC HBAR',
     &                    ' calculations'
          Write(6,*)
          if (coulomb) write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
          write(6,*)
          Fae_scale    = (Paralpha - 1.0D0)
          Fmi_scale    = (Parbeta  - 1.0D0)
          Wmnij_scale  = Pargamma
          Wmbej_scale  = Pardelta
          Gae_scale    = Paralpha  
          Gmi_scale    = Parbeta  
        Else
          write(6,*) ' Perform a regular CC HBAR calculations'
          write(6,*)
          Fae_scale    = 0.0D0
          Fmi_scale    = 0.0D0
          Wmnij_scale  = 1.0D0
          Wmbej_scale  = 1.0D0
          Gae_scale    = 1.0D0
          Gmi_scale    = 1.0D0
         Endif 
      Endif
C
C Obtain the symmetry of dipole moment integrals.
C
      Nbfns   = Nocco(1) + Nvrto(1)
      Call Getrec(20, "JOBARC","NBASTOT ", Ione, Naobfns)

      Call Tdee_get_pert_type(Work,Maxmem,Nbfns,Label,Isympert,
     +                        Iuhf)

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,1x,I4,1x,I4)") "The Nbfns and Naobfns:",
     +                            Nbfns,Naobfns
      Write(6,*)
#endif
      Call Tdee_form_lists(Work,Maxmem,Iuhf)

#ifdef _DEBUG_LVLM
      Write(6,*) "Checking T, L and Hbar lists"
      Call check_t(work,Maxmem,Iuhf)
      Call check_l(work,Maxmem,Iuhf)
#endif

      Call Tdee_built_respn_den(Work,Maxmem,Iuhf)

      Call Tdee_Modf(Work,Maxmem,Iuhf,One)
      If (Iuhf .EQ. 0) Call Tdee_makess(Work,Maxmem,Iuhf)
  
      Tranabci = (CC .AND. Iflags2(123) .EQ. 1 .AND. .NOT. DD)
      Call Modhbar(Work,Maxmem,iuhf,.False.,Tranabci) 

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " Hbar elements"
      Call checkhbar(Work,Maxmem,Iuhf)
#endif

 100  Continue 

C If a non-diagonal component. Save the M_bar_tilde(0) and M_bar(0) to
C built the autocorrelation function.

      If (Nondiag) Call Tdee_form_nondiag_state(Work,Maxmem,
     +                                          Nondiag_component,
     +                                          Iuhf)

C Loop over the 3-components of the diople integrals 

      Do Ixyz = Component, Component 

         If (.NOT. Restart) Then
         Irrepx     = Isympert(Ixyz)
         Ipert      = Ixyz
         Work_label = Label(Ixyz)

#ifdef _DEBUG_LVL0
      write(6,"(a,i2)") " Enter Tdee_reset_lists"
#endif 
         Call Tdee_reset_lists(Work,Maxmem,Iuhf,Irrepx)

         Endif 

         Call Tdee_init_lists(Work,Maxmem,Iuhf,Irrepx,Source,
     +                        Target,Side)

         Length   = Naobfns * (Naobfns+1)/2
         Length21 = Naobfns * Naobfns
         Length22 = Naobfns * Naobfns
         Lenoo    = Irpdpd(Irrepx,21) + Iuhf * Irpdpd(Irrepx,22)
         Lenvv    = Irpdpd(Irrepx,19) + Iuhf * Irpdpd(Irrepx,20)
         Lenvo    = Irpdpd(Irrepx,9)  + Iuhf * Irpdpd(Irrepx,10)
         Lenvo2_aa= Irpdpd(Irrepx,9)  * Irpdpd(Irrepx,9) 
         Lenvo2_bb= Irpdpd(Irrepx,10) * Irpdpd(Irrepx,10)
         Lenvo2_ab= Irpdpd(Irrepx,11) * Irpdpd(Irrepx,11)

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a,a)") "Length,Length21,Length22,Lenoo,Lenvv,Lenv0",
     +                   ",Lenvo2_aa,Lenvo2_bb,Lenvo2_ab: "
        Write(6,"(9(1x,I5))")Length,Length21,Length22,Lenoo,Lenvv,
     +                       Lenvo,Lenvo2_aa, Lenvo2_bb, Lenvo2_ab
        Write(6,*)
#endif 
         Nsize = Irpdpd(Irrepx,9)
         If (Iuhf .EQ. 0) Then
            Nsize = Nsize  + Idsymsz(Irrepx,13,14)
         Else
           Nsize = Nsize + Irpdpd(Irrepx,10)
           Nsize = Nsize + Idsymsz(Irrepx,1,3)
           Nsize = Nsize + Idsymsz(Irrepx,2,4)
           Nsize = Nsize + Idsymsz(Irrepx,13,14)
         Endif

         Ibgn = Ione
         I000 = Ibgn
         I010 = I000 + Lenoo 
         I020 = I010 + Lenvv
         I030 = I020 + Lenvo
         I040 = I030 + Length 
         I050 = I040 + Length21
         I060 = I050 + Length22
         I070 = I060 + Length22
         Iend = I070
         Memleft = Maxmem - Iend

         If (Iend .Gt. Maxmem) Call Insmem("@-Tdee_driver",Iend,
     +                                      Memleft)
         If (D_Pole) Then
         Call Tdee_prep_dipole_ints(Work(I030),Work(I040),Work(I050),
     +                              Work(I060),Work(I000),Work(I010),
     +                              Work(I020),Work(Iend),
     +                              Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                              Naobfns,Iuhf,Ipert,Work_label,
     +                              Irrepx)
         Elseif (Q_pole) Then
         Call Tdee_prep_qdpole_ints(Work(I030),Work(I040),Work(I050),
     +                              Work(I060),Work(I000),Work(I010),
     +                              Work(I020),Work(Iend),
     +                              Memleft,Lenoo,Lenvv,Lenvo,Nbfns,
     +                              Naobfns,Iuhf,Ipert,Work_label,
     +                              Irrepx)
        Endif 

C Built the  Mu_bar = exp(-T)Mexp(T)

         If (Mu_bar) Then

            Ibgn   = I030
            I040   = I030 + Lenoo
            I050   = I040 + Lenvv
            I060   = I050 + Lenvo
            I070   = I060 + Lenvo
            Iend   = I070 
            Iside  = 1

            Memleft = Maxmem - Iend
            Call Dcopy(Lenoo,Work(I000),1,Work(I030),1)
            Call Dcopy(Lenvv,Work(I010),1,Work(I040),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I050),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I060),1)

            Call Tdee_form_mubar_00(Work(I000),Work(I010),Work(I020),
     +                             Work(Iend),Memleft,Iuhf,Irrepx,
     +                             M_expect,Mubar_00)

            Call Tdee_form_mubar_s(Work(I000),Work(I010),Work(I020),
     +                             Work(I030),Work(I040),Work(I050),
     +                             Work(I060),Work(Iend),Memleft,Iuhf,
     +                             Irrepx,Iside,Mubar_00)

            If (Iuhf .Eq. 0) Then

               Call Tdee_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                    Work(I050),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00)
            Else

               Call Tdee_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                    Work(I050),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00)
            Endif 

         Endif 
           
         If (Mubar_t) Then 

            Ibgn   = I070
            I080   = I070 + Lenoo
            I090   = I080 + Lenvv
            I100   = I090 + Lenvo
            I110   = I100 + Lenvo
            Iend   = I110
            Iside  = 2

            Memleft = Maxmem - Iend
            Call Dcopy(Lenoo,Work(I000),1,Work(I070),1)
            Call Dcopy(Lenvv,Work(I010),1,Work(I080),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I090),1)
            Call Dcopy(Lenvo,Work(I020),1,Work(I100),1)

            Call Tdee_form_mubar_s(Work(I030),Work(I040),Work(I050),
     +                             Work(I070),Work(I080),Work(I090),
     +                             Work(I100),Work(Iend),Memleft,Iuhf,
     +                             Irrepx,Iside,Mubar_00-M_expect)

            If (Iuhf .Eq. 0) Then

               Call Tdee_form_mubar_d_rhf(Work(I030),Work(I040),
     +                                    Work(I060),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00-M_expect)
            Else

               Call Tdee_form_mubar_d_uhf(Work(I030),Work(I040),
     +                                    Work(I060),Work(Iend),
     +                                    Memleft,Iuhf,Irrepx,Iside,
     +                                    Mubar_00-M_expect)
            Endif

         Endif 

         Memleft = Maxmem
C
C Open a file to record data for restart (only in the first run)
C
         If (Iflags(h_IOPPAR_restart) .NE. 0) Then
             Inquire(File="Restart_log",Exist=Restart_present)
             If (.Not. Restart_present) Then
                Open(Unit=10,File="Restart_log",Form="formatted",
     +               Status="New")

C The orginal plan was to store in a list but that does not seems
C to work (the ACES IO subsystem can not be closed before aces_fin)
C Now we use the same stratgey used in CC and VEE code to T and 
C R amplitudes. 
CSSS                Call Tdee_form_restart_lists(Work(I000),Memleft,
CSSS     +                                       Iuhf,Irrepx,Nsize)
             Endif 
         Endif

C Save the M_bar_tilde(0) and M_bar(0) to built the autocorrelation function.
C Note that this is a safe measure on since any of the subsequent 
C steps do not touch the list containing M_bar_tilde(0) and M_bar(0).
C When nondiag components are considered reference state is already
C constructed and saved in appropriate place. 

         If (.Not. Nondiag) Then

             Do Iside = 1, 2
                Call Tdee_save_initial_state(Work(I000),Memleft,
     +                                       Irrepx,Iuhf,Iside,Nsize) 
             Enddo 
         Endif 

#ifdef _DEBUG_LVLM
         If (T_dot)Then
            Call Tdee_hbar_mult(Work(I000),Memleft,Iuhf,Irrepx,
     +                          Iside)
         Endif
#endif
#ifdef _DEBUG_LVL0 
         Write(6,*)
         Write(6,"(a,I5)") "The length of the Mubar,Mubar_dot vecs: ",
     +                   Nsize 
#endif 
         Do Iside = Side, Side 

C Iside = 1 The right-hand side AC=<Mbar_tilde(0)| M_bar(t)>
C Iside = 2 The left-hand side  AC=<Mbar_tilde(t)| M_bar(0)>

         If (SG) Then
            Iscrt_bgn = Maxmem    - (42*Nsize+100) + 1
            Imub0_bgn = Iscrt_bgn - 2*Nsize  
            Iac_rl_r  = Imub0_bgn - No_time_steps
            Iac_rl_i  = Iac_rl_r  - No_time_steps
            Memleft   = Iac_rl_i  - 1
            I000      = 1
  
            Call Tdee_init_global_Vars(Irrepx,Iuhf,Iside,Nsize,Memleft)
     +                                 

            If (Iac_rl_i .LE. 0) Then 
                Length = 44*Nsize + 100 + 2*No_time_steps
                Iend   = Maxmem
                Call Insmem("@-Tdee_driver",Length,Iend)
            Endif 

            Call Tdee_time_propagate_SG(Work(Iscrt_bgn),
     +                                  Work(Imub0_bgn),
     +                                  Work(Iac_rl_r),
     +                                  Work(Iac_rl_i), 
     +                                  Work(I000),Memleft,Iuhf,
     +                                  Irrepx,Iside,Nsize,Restart)
         Endif 

         If (RK) Then
            Iscrt_bgn = Maxmem    - (12*Nsize+100) + 1
            Imub0_bgn = Iscrt_bgn - 2*Nsize
            Iac_rl_r  = Imub0_bgn - No_time_steps
            Iac_rl_i  = Iac_rl_r  - No_time_steps
            Memleft   = Iac_rl_i  - 1
            I000      = 1

            Call Tdee_init_global_Vars(Irrepx,Iuhf,Iside,Nsize,Memleft)
     +                                 

            If (Iac_rl_i .LE. 0) Then
                Length = 14*Nsize + 100 + 2*No_time_steps
                Iend   = Maxmem
                Call Insmem("@-Tdee_driver",Length,Iend)
            Endif

            Call Tdee_time_propagate_RK45(Work(Iscrt_bgn),
     +                                    Work(Imub0_bgn),
     +                                    Work(Iac_rl_r),
     +                                    Work(Iac_rl_i), 
     +                                    Work(I000),Memleft,Iuhf,
     +                                    Irrepx,Iside,Nsize,Restart)
         Endif

#ifdef _DEBUG_LVLM
         Call Tdee_ac_debug(Work,Memleft,Irrepx,Iuhf)
#endif 
         Enddo 

C Use the non-symmetric Lanczos algorithm to obtain a sub-space
C that represent the Hbar matirx. This need both right and left
C multiplications. 

         If (LZ) Then

            Call Read_lanczos_vars()
            
            Irl_vec_bgn = Maxmem    - (2*Nsize*(Maxvw+3)) + 1
            Memleft     = Irl_vec_bgn  - 1
            I000        = 1

            If (Irl_vec_bgn .LE. 0) Then
                Length = 2*Nsize*(Maxvw+3)
                Iend   = Maxmem
                Call Insmem("@-Tdee_driver",Length,Iend)
            Endif

            Call Tdee_do_lanczos(Work(Irl_vec_bgn),Work(I000),Memleft,
     +                           Iuhf,Irrepx,Nsize,Maxmem)
         Endif

      Enddo 

      Return
      End
