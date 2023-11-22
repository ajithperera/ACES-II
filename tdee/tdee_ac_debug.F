      Subroutine Tdee_ac_debug(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
      Logical Rhf

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"
     
C Singles left, right dipole functions and derivative lists; 
C Mu : Right
C Mut: Left 

      Mu_s      = 390
      Mu_s_dot  = 394 
      Mut_s     = 392 
      Mut_s_dot = 396

C Corresponding doubles lists; There are no right dipole doubles
C list 

      Mu_d_aa = 314
      Mu_d_bb = 315
      Mu_d_ab = 316

      Mut_d_aa = 324
      Mut_d_bb = 325
      Mut_d_ab = 326
      
      Mu_d_dot_aa = 334
      Mu_d_dot_bb = 335
      Mu_d_dot_ab = 336

      Mut_d_dot_aa = 344
      Mut_d_dot_bb = 345
      Mut_d_dot_ab = 346

      Rhf = (Iuhf .EQ. 0)

      If (Rhf) Then 

         Ndim_ai = Irpdpd(Irrepx,9)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         I030 = I020 + Ndim_ai
         I040 = I030 + Ndim_ai
         I050 = I040 + Ndim_pphh_ab
         I060 = I050 + Ndim_pphh_ab
         I070 = I060 + Ndim_pphh_ab
         I080 = I070 + Ndim_pphh_ab
         Iend = I080
         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,1,Mu_s) 
         Call Getlst(Work(I010),1,1,1,1,Mu_s_dot) 
         Call Getlst(Work(I020),1,1,1,1,Mut_s) 
         Call Getlst(Work(I030),1,1,1,1,Mut_s_dot) 

#ifdef _DEBUG_LVL0
         Write(6,*)
         call checksum("Mu_s(a,i)    :",Work(I000),Ndim_ai,s)
         call checksum("Mu_s_dot(a,i):",Work(I010),Ndim_ai,s)
         call checksum("Mut_s(a,i)   :",Work(I020),Ndim_ai,s)
         call checksum("Mut_s_dot(a,i):",Work(I030),Ndim_ai,s)
#endif 
         Call Getall(Work(I040),Ndim_pphh_ab,Irrepx,Mu_d_ab)
         Call Getall(Work(I050),Ndim_pphh_ab,Irrepx,Mut_d_ab)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_dot_ab)
         Call Getall(Work(I070),Ndim_pphh_ab,Irrepx,Mut_d_dot_ab)

#ifdef _DEBUG_LVL0
         Write(6,*)
         call checksum("Mu_d(Ab,Ij)     :",Work(I040),Ndim_pphh_ab,s)
         call checksum("Mut_d(Ab,Ij)    :",Work(I050),Ndim_pphh_ab,s)
         call checksum("Mu_d_dot(Ab,Ij) :",Work(I060),Ndim_pphh_ab,s)
         call checksum("Mut_d_dot(Ai,Bj):",Work(I070),Ndim_pphh_ab,s)
#endif 
         Return
      Endif 

      Do Ispin = 1, (Iuhf+1)
         Ndim_ai = Irpdpd(Irrepx,8+Ispin)

         Ndim_pphh_aa = Idsymsz(Irrepx,1,3)
         Ndim_pphh_bb = Idsymsz(Irrepx,2,4)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         I030 = I020 + Ndim_ai
         I040 = I030 + Ndim_ai
         I050 = I040 + Ndim_pphh_aa
         I060 = I050 + Ndim_pphh_bb
         I070 = I060 + Ndim_pphh_ab
         I080 = I070 + Ndim_pphh_aa
         I090 = I080 + Ndim_pphh_bb
         I100 = I090 + Ndim_pphh_ab
         Iend = I070
         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,Ispin,Mu_s) 
         Call Getlst(Work(I010),1,1,1,Ispin,Mu_s_dot)
         Call Getlst(Work(I020),1,1,1,Ispin,Mut_s)
         Call Getlst(Work(I030),1,1,1,Ispin,Mut_s_dot)

#ifdef _DEBUG_LVL0
         write(6,*)
         call checksum("Mu_s(a,i)    :",Work(I000),Ndim_ai,s)
         call checksum("Mu_s_dot(a,i):",Work(I010),Ndim_ai,s)
         call checksum("Mut_s(a,i)   :",Work(I020),Ndim_ai,s)
         call checksum("Mut_s_dot(a,i):",Work(I030),Ndim_ai,s)
#endif 
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mu_d_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mu_d_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_ab)

#ifdef _DEBUG_LVL0
         write(6,*)
         call checksum("Mu_d(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mu_d(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mu_d(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
#endif 
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mut_d_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mut_d_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mut_d_ab)

#ifdef _DEBUG_LVL0
         write(6,*)
         call checksum("Mut_d(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mut_d(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mut_d(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
#endif 
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mut_d_dot_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mut_d_dot_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mut_d_dot_ab)

#ifdef _DEBUG_LVL0
         Write(6,*)
         call checksum("Mut_d_dot(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mut_d_dot(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mut_d_dot(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
#endif 
         Call Getall(Work(I040),Ndim_pphh_aa,Irrepx,Mu_d_dot_aa)
         Call Getall(Work(I050),Ndim_pphh_bb,Irrepx,Mu_d_dot_bb)
         Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_dot_ab)

#ifdef _DEBUG_LVL0
         Write(6,*)
         call checksum("Mu_d_dot(A<B,I<J):",Work(I040),Ndim_pphh_aa,s)
         call checksum("Mu_d_dot(a<b,i<j):",Work(I050),Ndim_pphh_bb,s)
         call checksum("Mu_d_dot(Ab,Ij)  :",Work(I060),Ndim_pphh_ab,s)
#endif 
      Enddo

      Return 
      End
