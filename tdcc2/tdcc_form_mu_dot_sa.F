#include "flags.h"
      Subroutine Tdcc_form_mu_dot_sa(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

      Double Precision Mone,Mtwo 
      Integer Hbr_pp_list,Hbr_hh_list,Hbar_phph_list
      Integer Hbar_phph_list_a,Hbar_phph_list_b
 
      Data Zero, One, Mone, Two, Mtwo /0.0D0,  1.0D0, -1.0D0, 
     +                                2.0D0, -2.0D0/
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mu_dot_sa---"
        Write(6,*)
#endif
      Mu0_ph_list = 390
      Hbr_pp_list = 92
      Hbr_hh_list = 91
      Mud_ph_list = 394

      Do Ispin = 1, 1+Iuhf 
         
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Irpdpd(Irrepx,8+Ispin)
         I020 = I010 + Irpdpd(1,18+Ispin) 
         I030 = I020 + Irpdpd(1,20+Ispin) 
         I040 = I030 + Irpdpd(Irrepx,8+Ispin)
         Iend = I040

         If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                       Memleft,Iend)

         Call Getlst(Work(I000),1,1,1,Ispin,Mu0_ph_list)
         Call Getlst(Work(I010),1,1,1,Ispin,Hbr_pp_list)
         Call Getlst(Work(I020),1,1,1,Ispin,Hbr_hh_list)

#ifdef _DEBUG_LVLM
         Call checksum("Mu(ai)",Work(I000), Irpdpd(Irrepx,8+Ispin),S)
         Call checksum("Hbar(ae)",Work(I010), Irpdpd(1,18+Ispin),S)
         Call checksum("Hbar(ij)",Work(I020),Irpdpd(1,20+Ispin),S)
#endif 

         Ioff_m0_bgn = I000 
         Ioff_pp_bgn = I010 
         Ioff_hh_bgn = I020 
         Ioff_zd_bgn = I030 

         Ndim_zd = Irpdpd(Irrepx,8+Ispin)
         Call Dzero(Work(I030),Ndim_zd)

C NDeP Eqn. 25, First term: Mu(e,i) * Hbar(a,e) = Zdot(a,i)
C Stored as Hbar(e,a)

         Do Irrep_i = 1, Nirrep
            Irrep_e = Dirprd(Irrep_i,Irrepx)
            Irrep_a = Irrep_e

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_e = Vrt(Irrep_e,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ioff_pp = Ioff_pp_bgn + (Isymoff(Irrep_e,1,18+Ispin)-1)
            Ioff_m0 = Ioff_m0_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ioff_zd = Ioff_zd_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)

            Nrow = Ndim_a
            Nsum = Ndim_e
            Ncol = Ndim_i

            Call Xgemm("T","N",Nrow,Ncol,Nsum,One,Work(Ioff_pp),
     +                  Nrow,Work(Ioff_m0),Nsum,One,Work(Ioff_zd),
     +                  Nrow)
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(1) :",
     +               Work(Ioff_zd_bgn),Ndim_zd,S)
#endif 
         Enddo

C NDeP Eqn. 25, Second term: Mu(a,m) * Hbar(m,i) = Zdot(a,i)

         Do Irrep_m = 1, Nirrep
            Irrep_a = Dirprd(Irrep_m,Irrepx)
            Irrep_i = Irrep_m

            Ndim_m = Pop(Irrep_m,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)

            Ioff_hh = Ioff_hh_bgn + (Isymoff(Irrep_i,1,20+Ispin)-1)
            Ioff_m0 = Ioff_m0_bgn + (Isymoff(Irrep_m,Irrepx,8+Ispin)-1)
            Ioff_zd = Ioff_zd_bgn + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            
            Nrow = Ndim_a
            Nsum = Ndim_m
            Ncol = Ndim_i

            Call Xgemm("N","N",Nrow,Ncol,Nsum,Mone,Work(Ioff_m0),
     +                 Nrow,Work(Ioff_hh),Nsum,One,Work(Ioff_zd),
     +                 Nrow)
         Enddo 

#ifdef _DEBUG_LVL0
      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(f1) :",
     +               Work(Ioff_zd_bgn),Ndim_zd,S)
#endif 
         Call Putlst(Work(Ioff_zd_bgn),1,1,1,Ispin,Mud_ph_list)

      Enddo 

C NDep Eqn. 25, Third term:

C UHF:  Zdot(A,I) = Mu(E,M)*Hbar(MA,EI) + Mu(e,m)Hbar(mA,eI)
C       Zdot(a,i) = Mu(E,m)*Hbar(ma,Ei) + Mu(E,M)Hbar(Ma,Ei)
C RHF:  Zdot(A,I) = Mu(E,M){(2 Hbar(Ma,Ei) - Hba(Ma,Ie)}

      Do Ispin = 1, 1+Iuhf

         If (IUHF .EQ. 0) Then
            Hbar_phph_list = 56
         Else
            Hbar_phph_list = 58 - Ispin
         Endif 
         
         Ndim_hb_em = Irpdpd(Irrepx,Isytyp(1,Hbar_phph_list))
         Ndim_hb_ai = Irpdpd(Irrepx,Isytyp(2,Hbar_phph_list))
         Ndim_mu_em = Irpdpd(Irrepx,11-Ispin)
         Ndim_zs_ai = Irpdpd(Irrepx,8+Ispin)
         Ndim_zd_ai = Ndim_zs_ai 

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max(Ndim_zd_ai,Ndim_mu_em)
         I020 = I010 + Max(Ndim_zs_ai,Ndim_mu_em)
         I030 = I020 + Ndim_hb_em * Ndim_hb_ai
         Iend = I030
         If (Iend .GT. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                       Memleft,Iend)
         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(em,AI)*mu(e,m) (Ispin=1)
C Form Z(a,i) = Hbar(EM,ai)*mu(E,M) (Ispin=2)

            Call Getlst(Work(I000),1,1,1,Ispin,Mud_ph_list)
            Call Getlst(Work(I010),1,1,1,3-Ispin,Mu0_ph_list)
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                       Hbar_phph_list)
            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

#ifdef _DEBUG_LVLM
      Call checksum("Mu_d(ai)",Work(I000), Irpdpd(Irrepx,8+ispin),S)
      Call checksum("Mu_s(ai)",Work(I010), Irpdpd(Irrepx,8+Ispin),S)
      Call checksum("Hbar(em,ai)",Work(I020),Ndim_hb_em*Ndim_hb_ai,S)
#endif
         Else 

C List 54=Hbar(EM,AI), List 58 = (Em,Ai)

            Call Getlst(Work(I000),1,1,1,1,Mud_ph_list)
            Call Getlst(Work(I010),1,1,1,1,Mu0_ph_list)
            Hbar_phph_list_a = 54
            Hbar_phph_list_b = 58

            Ibgn = I030
            I040 = I030 + Ndim_hb_em * Ndim_hb_ai
            Iend = I040
            If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                          Memleft,Iend)
            
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_b)
            Call Getlst(Work(I030),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_a)
            Call Daxpy(Ndim_hb_em * Ndim_hb_ai,Mtwo,Work(I030),1,
     +                 Work(I020),1)
            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

         Endif

         Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(I010),
     +               Nrow,Work(I020),Nsum,One,Work(I000),
     +               Nrow)

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(3a):",
     +               Work(I000),Ndim_zd_ai,S)
#endif
         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(EM,AI)*mu(E,M) (Ispin=1)
C Form Z(a,i) = Hbar(em,ai)*mu(e,m) (Ispin=2)

            Hbar_phph_list = 53 + Ispin

            Ndim_hb_em = Irpdpd(Irrepx,Isytyp(1,Hbar_phph_list))
            Ndim_hb_ai = Irpdpd(Irrepx,Isytyp(2,Hbar_phph_list))
            Ndim_zs_ai = Irpdpd(Irrepx,8+Ispin)
            Ndim_mu_em = Irpdpd(Irrepx,11-Ispin)

            Ibgn = 1
            I000 = Ibgn
            I010 = I000 + Max(Ndim_zs_ai,Ndim_zs_em)
            I020 = I010 + Max(Ndim_mu_ai,Ndim_mu_em)
            I030 = I020 + Ndim_hb_em * Ndim_hb_ai 
            Iend = I030
            If (Iend .Gt. Memleft) Call Insmem("@-Tdcc_form_mu_0_S",
     +                                         Memleft,Iend)
       
            Call Getlst(Work(I010),1,1,1,Ispin,Mu0_ph_list)
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,Hbar_phph_list)

            Nsum = Ndim_hb_em
            Nrow = 1
            Ncol = Ndim_hb_ai

C Double check the minus sign (Stored as Hbar(ME,AI))
         
            Call Xgemm("N","N",Nrow,Ncol,Nsum,Mone,Work(I010),
     +                  Nrow,Work(I020),Nsum,One,Work(I000),
     +                  Nrow)
         Endif

#ifdef _DEBUG_LVL0
      call Checksum("Tdcc_form_mu_dot_sa,Mu_dot(f2) :",
     +               Work(I000),Ndim_zd_ai,S)
#endif
#ifdef _DEBUG_LVL0

#else
         Call Putlst(Work(I000),1,1,1,Ispin,Mud_ph_list)
#endif 

      Enddo 

      Return
      End
     
           
        
