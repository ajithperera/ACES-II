      Subroutine Tdcc_form_mutilde_dot_dc_rhf(Work,Memleft,Irrepx)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Logical Rhf
      Integer Tmp_list,Hb_ph_list
      Character*4 Spcase(2)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

      Data Spcase /'AAAA','BBBB'/

      Data One /1.0D0/

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_dot_dc_rhf---"
        Write(6,*)
#endif

C The ABAB block.

      Iz2_td_list = 346
      Hb_ph_list  = 93
      Mut_s_list  = 392
      Tmp_list    = 302

      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_td_list),
     +                     Isytyp(2,Iz2_td_list))

C NDeP Eqn. 29, fourth term: Z(Ib,Aj) = Hbar(j,b) * Mut(A,I) and
C Z is built as Z(AI,bj).

      Do Irrep_bj_z = 1, Nirrep
         Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)

         Ndim_ai_z   = Irpdpd(Irrep_ai_z,9)
         Ndim_bj_z   = Irpdpd(Irrep_bj_z,10)

         Ndim_aibj_z = Ndim_bj_z * Ndim_ai_z

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_aibj_z
         I020 = I010 + Irpdpd(1,9)
         Iend = I020 + Irpdpd(1,10)

         Ioff_hb     = I010
         Ioff_mut    = I020
         Ioff_mut_d  = I000

         Call Dzero(Work(I000),Ndim_aibj_z)

C Z(AI,bj) = Mut(b,j) Hbar(A,I)

         If (Irrep_bj_z .EQ. Irrepx) Then

            Call Getlst(Work(I010),1,1,1,1,Hb_ph_list)
            Call Getlst(Work(I020),1,1,1,1,Mut_s_list)

            Nrow = Ndim_ai_z
            NCol = Ndim_bj_z
            Nsum = 1
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ioff_hb),Nrow,
     +                  Work(Ioff_mut),1,One,Work(Ioff_mut_d),Nrow)
         Endif

         Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,
     +               Tmp_list)
      Enddo

      Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                            Isytyp(2,Tmp_list))
      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_td_list),
     +                            Isytyp(2,Iz2_td_list))
      Ibgn = 1
      I000 = Ibgn
      I010 = I000  + Max(Isize_phph,Isize_pphh)
      I020 = I010  + Max(Isize_phph,Isize_pphh)
      I030 = I020  + 100
      Iend = I030

      If (Iend .GT. Memleft) Call Insmem
     +   ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

      Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)
      Call Sstgen(Work(I000),Work(I010),Isize_phph,Vrt(1,1),Pop(1,1),
     +            Vrt(1,2),Pop(1,2),Work(I020),Irrepx,'1324')

#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_dot_dc,mu^t_dot(Ab,Ij):",
     +               Work(I010),Isize_pphh,S)
#endif 
#ifdef _DEBUG_LVL0
      Call Putall(Work(I000),Isize_pphh,Irrepx,Iz2_td_list)
#else
      Call Getall(Work(I000),Isize_phph,Irrepx,Iz2_td_list)
      Call Daxpy(Isize_phph,One,Work(I000),1,Work(I010),1)
      Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_td_list)
#endif 

#ifdef _DEBUG_LVL0
      call checksum("Tdcc_form_mutilde_dot_dc,mu^t_dot(Ab,Ij)(f):",
     +               Work(I010),Isize_pphh,S)
#endif 
      Return
      End
