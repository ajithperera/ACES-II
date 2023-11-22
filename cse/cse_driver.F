      Subroutine Cse_driver(T1old_aa,T1old_bb,W0_aa,W0_bb,W0_ab,
     +                      Fockoo_a,Fockoo_b,Fockvv_a,
     +                      Fockvv_b,Fockov_a,Fockov_b,
     +                      Fockoo_od_a,Fockoo_od_b,
     +                      Fockvv_od_a,Fockvv_od_b,Fmi_a,
     +                      Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                      Dlt_a,Dlt_b,h1bar_aa,h1bar_bb,
     +                      Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +                      Work,Maxcor,Repuls,Ref_energy,Frac_occ)

      Implicit Double Precision(A-H,O-Z)

#include "cse_vars.com"

      Dimension Work(Maxcor)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Wc_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1old_aa(Nvrt_a,Nocc_a)
      Dimension T1old_bb(Nvrt_b,Nocc_b)
      Dimension Dlt_a(Nvrt_a,Nocc_a)
      Dimension Dlt_b(Nvrt_b,Nocc_b)

      Dimension Densp_aa(Nbasis,Nbasis)
      Dimension Densp_bb(Nbasis,Nbasis)
      Dimension Densm_aa(Nbasis,Nbasis)
      Dimension Densm_bb(Nbasis,Nbasis)

      Dimension fockoo_a(Nocc_a,Nocc_a)
      Dimension fockoo_b(Nocc_b,Nocc_b)
      Dimension fockvv_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_b(Nvrt_b,Nvrt_b)
      Dimension fockov_a(Nocc_a,Nvrt_a)
      Dimension fockov_b(Nocc_b,Nvrt_b)

      Dimension fockoo_od_a(Nocc_a,Nocc_a)
      Dimension fockoo_od_b(Nocc_b,Nocc_b)
      Dimension fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_od_b(Nvrt_b,Nvrt_b)

      Dimension fmi_a(Nocc_a,Nocc_a)
      Dimension fmi_b(Nocc_b,Nocc_b)
      Dimension fae_a(Nvrt_a,Nvrt_a)
      Dimension fae_b(Nvrt_b,Nvrt_b)
      Dimension fme_a(Nocc_a,Nvrt_a)
      Dimension fme_b(Nocc_b,Nvrt_b)

      Dimension H1bar_aa(Nocc_a*Nvrt_a,Nocc_a*Nvrt_a) 
      Dimension H1bar_ba(Nocc_a*Nvrt_a,Nocc_a*Nvrt_a) 

      Logical Converged,Fme_on,T_guess_exist 
      Logical Frac_occ,Nr 
      Character*80 Fullname
      
      Call Guesst1(T1old_aa,T1old_bb,Densp_aa,Densp_bb,Nocc_a,Nocc_b,
     +             Nvrt_a,Nvrt_b,Nbasis)

C Check whether the TGUESS file exists. Then read it and use those
C as the starting guess

      Call Gfname("TGUESS",Fullname,Namelen)
      Inquire(File=Fullname(1:Namelen),Exist=T_guess_exist)

      Iunit = 94
      If (T_guess_exist .and. .Not. Brueck) Then
         Open(Unit=Iunit,File=Fullname(1:Namelen),Status="Old",
     +        Form="Unformatted",Access="Sequential")
         Call Read_amps(Iunit,T1old_aa,T1old_bb,Nocc_a,Nocc_b,
     +                  Nvrt_a,Nvrt_b,Work,Maxcor)
         Close(Unit=Iunit,Status="Keep") 
      Endif

      L_aa   = Nocc_a*Nvrt_a
      L_bb   = Nocc_b*Nvrt_b

      It1resid_aa = 1 
      It1resid_bb = It1resid_aa + L_aa 
      Iend        = It1resid_bb + L_bb
      Memleft     = Maxcor - Iend 

      If (Iend .GT. Maxcor) Call Insmem("ccsd_light_driver",Iend,
     +                                   Maxcor)
      Iter      = 0
      Converged = .False. 
      E_old     = 0.0D0
      E_corr(0) = 0.0D0

      Do while ((.NOT. Converged) .AND. Iter .Le. Cc_maxcyc) 
         Iter = Iter + 1

C Fmi(m,i) intermediate 

         Call Form_fmi(T1old_aa,T1old_bb,Fmi_a,Fmi_b,
     +                 Fockoo_od_a,Fockoo_od_b,Fockov_a,
     +                 Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Fae(a,e) intermediate 

         Call Form_fae(T1old_aa,T1old_bb,Fae_a,Fae_b,
     +                 Fockvv_od_a,Fockvv_od_b,Fockov_a,
     +                 Fockov_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

C Fme(m,e) intermediate 

         Call Form_fme(Fme_a,Fme_b,Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b,Fme_on)

C Initialize T1 residual f(a,i)
 
         Call init_2abij(Work(It1resid_aa),Work(It1resid_bb),
     +                   Fockov_a,Fockov_b,Nocc_a,Nocc_b,
     +                   Nvrt_a,Nvrt_b)

C Fae->T1

         Call Fae_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fae_a,Fae_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)

C Fmi->T1 

         Call Fmi_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fmi_a,Fmi_b,Nocc_a,Nocc_b,
     +                 Nvrt_a,Nvrt_b)
         NR = .False.

         If (NR) Then

C foo + fvv -> T1 (diagonals)

            Call fpq_int1(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                 Work(It1resid_bb),Fockoo_a,Fockoo_b,
     +                 Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                 Nvrt_b)

            Call Form_h1bar(Fmi_a,Fmi_b,Fae_a,Fae_b,Fme_a,Fme_b,
     +                      Work(It1resid_aa),Work(It1resid_bb),
     +                      H1bar_aa,H1bar_bb,Fockoo_a,
     +                      Fockoo_b,Fockvv_a,Fockvv_b,Work(Iend),
     +                      Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,Nbasis,
     +                      Memleft)

            Call Newt1_yj(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                    Work(It1resid_bb),H1bar_aa,H1bar_bb,
     +                    Dlt_a,Dlt_b,Dt1,Nocc_a,Nocc_b,Nvrt_a,
     +                    Nvrt_b,Nbasis) 
         Else

C T1old->T1new 

            Call Newt1_nj(T1old_aa,T1old_bb,Work(It1resid_aa),
     +                    Work(It1resid_bb),Fockoo_a,Fockoo_b,
     +                    Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                    Nvrt_b)
         Endif 

         Call Energy_frmt1(T1old_aa,T1old_bb,W0_aa,W0_bb,W0_ab,
     +                     Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,E1ccsd_aa,E1ccsd_bb,E1ccsd_ab)

         Call Dens_frmt1(Densp_aa,Densp_bb,Densm_aa,Densm_bb,
     +                   T1old_aa,T1old_bb,Fockoo_a,Fockoo_b,
     +                   Repuls,Work,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                   Nbasis,Maxcor,Frac_occ,E1) 

         E1_new =  E1ccsd_aa + E1ccsd_bb
         E_new  =  E1_new
         E_diff= Dabs(E_new - E_old)
         E_corr(Iter) = E_new 

#ifdef _DEBUG_LVL0
         Write(6,"(a,3(1x,F15.10))")"E1aa,E1bb            :",
     +                               E1ccsd_aa,
     +                               E1ccsd_bb
         Write(6,"(a,2(1x,F15.10),2x,a)")
     +         "The energy and the energy difference: ", E_new, E_diff,
     +          "at" 
         Write(6,"(a,1x,I3))") "iteration no.: ", Iter
#endif 
CSSS         If (Iter .Eq. 3) Converged = .True.

         If (NR) Then
            IF (E_diff .LT. CC_conv .OR. DT1 .LT. CC_conv) 
     +          Converged = .True.
         Else
            IF (E_diff .LT. CC_conv) Converged = .True.
         Endif 
         E_old = E_new 

      enddo

      Write(6,"(a)") "The FON reference energy is added."
      E_scf = Ref_energy 

      If (Converged) Then  
         Call dump_amps(T1old_aa,T1old_bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
         Write(6,*)
         Write(6,9) 
         Write(6,10)
         Write(6,11)
         Write(6,10)
      Else
         Write(6,*)
         Write(6,"(a,1x,I3,a)") " The CCSD did not converged in", 
     +                          cc_maxcyc, " iterations."
         Write(6,10)
      Endif 

 9    Format(T3, '    Summary of iterative solution of CCSD equations')
 10   Format(6x,59("-")) 
 11   Format(T24,'Correlation',T55,'Total',/,T8,'Iteration',T26,
     +       'Energy',T55,'Energy')
      Do I = 0, Iter
         Write(6,"(5x,I3,10x,F18.12,10x,F19.12)") I, E_corr(i), 
     +                                             E_corr(i)+E_scf
      Enddo
      Write(6,10)

      Return
      End 

