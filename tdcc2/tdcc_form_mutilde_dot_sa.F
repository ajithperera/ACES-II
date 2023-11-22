#include "flags.h"
      Subroutine Tdcc_form_mutilde_dot_sa(Work,Memleft,Irrepx,Iuhf)

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

      Mut_ph_list  = 392
      Hbr_pp_list  = 92
      Hbr_hh_list  = 91
      Mutd_ph_list = 396

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_dot_sa---"
        Write(6,*)
#endif
      Do Ispin = 1, 1+Iuhf 
         
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Irpdpd(Irrepx,8+Ispin)
         I020 = I010 + Irpdpd(1,18+Ispin) 
         I030 = I020 + Irpdpd(1,20+Ispin) 
         I040 = I030 + Irpdpd(Irrepx,8+Ispin)
         Iend = I040

         If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mu_0_S",
     +                                       Memleft,Iend)

         Call Getlst(Work(I000),1,1,1,Ispin,Mut_ph_list)
         Call Getlst(Work(I010),1,1,1,Ispin,Hbr_pp_list)
         Call Getlst(Work(I020),1,1,1,Ispin,Hbr_hh_list)
         
         Ioff_mt_bgn  = I000 
         Ioff_pp_bgn  = I010 
         Ioff_hh_bgn  = I020
         Ioff_ztd_bgn = I030 

C NDeP Eqn. 25, First term: Mut(i,e) * Hbar(e,a) = Ztdot(a,i)
C Form as Ztdot(a,i)=Mut(e,i)*Hbar(a,e) (stored as)

         Ndim_ztd =  Irpdpd(Irrepx,8+Ispin)
         Call Dzero(Work(I030),Ndim_ztd)

         Do Irrep_i = 1, Nirrep
            Irrep_e = Dirprd(Irrep_i,Irrepx)
            Irrep_a = Irrep_e

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_e = Vrt(Irrep_e,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ioff_pp  = Ioff_pp_bgn  + (Isymoff(Irrep_e,1,
     +                                 18+Ispin)-1)
            Ioff_mt  = Ioff_mt_bgn  + (Isymoff(Irrep_i,Irrepx,
     +                                 8+Ispin)-1)
            Ioff_ztd = Ioff_ztd_bgn + (Isymoff(Irrep_i,Irrepx,
     +                                 8+Ispin)-1)
            Nrow = Ndim_a
            Nsum = Ndim_e
            Ncol = Ndim_i

            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_pp),
     +                  Nrow,Work(Ioff_mt),Nsum,Zero,Work(Ioff_ztd),
     +                  Nrow)
         Enddo

#ifdef _DEBUG_LVLM
      write(6,*)
      call Checksum("Tdcc_form_mutilde_dot_sa,Mutilde_dot(1)  :",
     +               Work(Ioff_ztd_bgn),Ndim_ztd,S)
#endif

C NDeP Eqn. 25, Second term: Mut(m,a) * Hbar(i,m) = Ztdot(a,i)
C Form as Ztdot(a,i) Mut(a,m)*Hbar(i,m)^t

         Do Irrep_m = 1, Nirrep
            Irrep_a = Dirprd(Irrep_m,Irrepx)
            Irrep_i = Irrep_m

            Ndim_m = Pop(Irrep_m,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)

            Ioff_hh  = Ioff_hh_bgn  + (Isymoff(Irrep_i,1,
     +                                 20+Ispin)-1)
            Ioff_mt  = Ioff_mt_bgn  + (Isymoff(Irrep_m,Irrepx,
     +                                 8+Ispin)-1)
            Ioff_ztd = Ioff_ztd_bgn + (Isymoff(Irrep_i,Irrepx,
     +                                 8+Ispin)-1)
            Nrow = Ndim_a
            Nsum = Ndim_m
            Ncol = Ndim_i

            Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_mt),
     +                 Nrow,Work(Ioff_hh),Nsum,One,Work(Ioff_ztd),
     +                 Nrow)
         Enddo 
#ifdef _DEBUG_LVL0
      write(6,*)
      call Checksum("Tdcc_form_mutilde_dot_sa,Mutilde_dot(f):",
     +               Work(Ioff_ztd_bgn),Ndim_ztd,S)
#endif
         Call Putlst(Work(Ioff_ztd_bgn),1,1,1,Ispin,Mutd_ph_list)

      Enddo 

C NDep Eqn. 25, Third term:

C UHF:  Zdot(A,I) = Mu(M,E)*Hbar(EI,MA) + Mu(m,e)Hbar(eI,mA)
C       Zdot(a,i) = Mu(m,e)*Hbar(ei,ma) + Mu(M,E)Hbar(Ei,ma)
C RHF:  Zdot(A,I) = Mu(M,E){(2 Hbar(Ei,Ma) - Hba(Ie,M)}

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
         If (Iend .Gt. Memleft)Call Insmem("@-Tdcc_form_mutilde_0_Sa",
     +                                     Memleft,Iend)
         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(em,AI)*mu(m,e) (Ispin=1)
C Form Z(a,i) = Hbar(EM,ai)*mu(M,E) (Ispin=2)

            Call Getlst(Work(I000),1,1,1,Ispin,Mutd_ph_list)
            Call Getlst(Work(I010),1,1,1,3-Ispin,Mut_ph_list)
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                       Hbar_phph_list)
            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

         Else 

C List 54=Hbar(EM,AI), List 58 = (Em,Ai)

            Call Getlst(Work(I000),1,1,1,1,Mutd_ph_list)
            Call Getlst(Work(I010),1,1,1,1,Mut_ph_list)
            Hbar_phph_list_a = 54
            Hbar_phph_list_b = 58

            Ibgn = I030
            I040 = I030 + Ndim_hb_em * Ndim_hb_ai
            Iend = I040
            If (Iend .Gt. Memleft) Call 
     +                    Insmem("@-Tdcc_form_mutilde_0_Sa",
     +                           Memleft,Iend)
            
            Call Getlst(Work(I020),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_b)
            Call Getlst(Work(I030),1,Ndim_hb_ai,1,Irrepx,
     +                  Hbar_phph_list_a)
            Call Daxpy(Ndim_hb_em*Ndim_hb_ai,Mtwo,Work(I030),1,
     +                 Work(I020),1)

            Nsum = Ndim_mu_em
            Nrow = 1
            Ncol = Ndim_hb_ai

         Endif

         Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(I010),
     +               Nrow,Work(I020),Nsum,One,Work(I000),
     +               Nrow)
#ifdef _DEBUG_LVLM
      write(6,*)
      call Checksum("Tdcc_form_mutilde_dot_sa,Mutilde_dot(3a) :",
     +               Work(I000),Ndim_zd_ai,S)
#endif
         If (Iuhf .NE. 0) Then

C Form Z(A,I) = Hbar(EM,AI)*mu(M,E) (Ispin=1)
C Form Z(a,i) = Hbar(em,ai)*mu(M,E) (Ispin=2)

            Hbar_phph_list = 53 + Ispin

            Ndim_em = Irpdpd(Irrepx,Isytyp(1,Hbar_phph_list))
            Ndim_ai = Irpdpd(Irrepx,Isytyp(2,Hbar_phph_list))
            Ndim_ai = Irpdpd(Irrepx,8+Ispin)
            Ndim_em = Irpdpd(Irrepx,11-Ispin)

            Ibgn = 1
            I000 = Ibgn
            I010 = I000 + Max(Ndim_ai,Ndim_em)
            I020 = I010 + Max(Ndim_ai,Ndim_em)
            I030 = I020 + Ndim_em * Ndim_ai
       
            Call Getlst(Work(I010),1,1,1,Ispin,Mut_ph_list)
            Call Getlst(Work(I020),1,Ndim_ai,1,Irrepx,
     +                  Hbar_phph_list)

            Nsum = Ndim_em
            Nrow = 1
            Ncol = Ndim_ai

C Double check the minus sign (Stored as Hbar(ME,AI))
         
            Call Xgemm("N","T",Nrow,Ncol,Nsum,MOne,Work(I010),
     +                  Nrow,Work(I020),Nsum,One,Work(I000),
     +                  Nrow)
         Endif
#ifdef _DEBUG_LVL0
      call Checksum("Tdcc_form_mutilde_dot_sa,Mutilde_dot(f):",
     +               Work(I000),Ndim_zd_ai,S)
#endif
#ifdef _DEBUG_LVL0

#else
         Call Putlst(Work(I000),1,1,1,Ispin,Mutd_ph_list)
#endif 
      Enddo 

      Return
      End
     
