      Subroutine Tdee_built_ac_func(Work,Memleft,Irrepx,Iuhf)

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
      Mut_s     = 392 

      Mu_d_aa_0 = 317
      Mu_d_bb_0 = 318
      Mu_d_ab_0 = 319

      Mu_d_aa_t = 314
      Mu_d_bb_t = 315
      Mu_d_ab_t = 316

      Mut_d_aa_0 = 327
      Mut_d_bb_0 = 328
      Mut_d_ab_0 = 329

      Mut_d_aa_t = 324
      Mut_d_bb_t = 325
      Mut_d_ab_t = 326

      Rhf = (Iuhf .EQ. 0)

      If (Rhf) Then 

         Ndim_ai = Irpdpd(Irrepx,9)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         I030 = I020 + Ndim_ai
         I040 = I030 + Ndim_pphh_ab
         I060 = I040 + Ndim_pphh_ab

         Iend = I060
         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,3,Mu_s_0) 
         Call Getlst(Work(I010),1,1,1,3,Mut_s) 

         ML_t_MR_0 = Ddot(Ndim_ai,Work(I000),1,Work(I010),1)

         Call Getlst(Work(I000),1,1,1,1,Mut_s_0) 
         Call Getlst(Work(I010),1,1,1,3,Mu_s) 

         ML_0_MR_t = Ddot(Ndim_ai,Work(I000),1,Work(I010),1)

         Call Getall(Work(I030),Ndim_pphh_ab,Irrepx,Mu_d_ab_0)
         Call Getall(Work(I040),Ndim_pphh_ab,Irrepx,Mut_d_ab_t)
 
         ML_t_MR_0 = ML_t_MR_0 + Ddot(Ndim_pphh_ab,Work(I030),1,
     +               Work(I040),1)

         Call Getall(Work(I030),Ndim_pphh_ab,Irrepx,Mu_d_ab_t)
         Call Getall(Work(I040),Ndim_pphh_ab,Irrepx,Mut_d_ab_0)
  
         ML_0_MR_t = ML_0_MR_t + Ddot(Ndim_pphh_ab,Work(I030),1,
     +               Work(I040),1)

      Endif 

      ML_t_MR_0 = 0.0D0
      ML_0_MR_t = 0.0D0

      Do Ispin = 1, (Iuhf+1)
         Ndim_ai = Irpdpd(Irrepx,8+Ispin)

         Ndim_pphh_aa = Idsymsz(Irrepx,1,3)
         Ndim_pphh_bb = Idsymsz(Irrepx,2,4)
         Ndim_pphh_ab = Idsymsz(Irrepx,13,14)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_ai
         I020 = I010 + Ndim_ai
         Iend = I020

         If (Iend .GT. Memleft) Then
            Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
         Endif

         Call Getlst(Work(I000),1,1,1,2+Ispin,Mu_s_0)
         Call Getlst(Work(I020),1,1,1,2+Ispin,Mut_s)

         ML_t_MR_0 = ML_t_MR_0 + Ddot(Ndim_ai,Work(I000),1,
     +                                Work(I010),1)
      Enddo

      I030 = I020 + Ndim_pphh_aa
      I040 = I030 + Ndim_pphh_bb
      Iend = I040 

      If (Iend .GT. Memleft) Then
         Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
      Endif

      Call Getall(Work(I020),Ndim_pphh_aa,Irrepx,Mu_d_aa_0)
      Call Getall(Work(I030),Ndim_pphh_aa,Irrepx,Mu_d_aa_t)

      ML_t_MR_0 = ML_t_MR_0 + Ddot(Ndim_ai,Work(I020),1,
     +                                Work(I030),1)
      I030 = I020 + Ndim_pphh_bb
      I040 = I030 + Ndim_pphh_bb

      If (Iend .GT. Memleft) Then
         Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
      Endif

      Call Getall(Work(I020),Ndim_pphh_bb,Irrepx,Mu_d_bb_0)
      Call Getall(Work(I030),Ndim_pphh_bb,Irrepx,Mu_d_bb_t)

      ML_t_MR_0 = ML_t_MR_0 + Ddot(Ndim_ai,Work(I020),1,
     +                                Work(I030),1)
 
      I030 = I020 + Ndim_pphh_ab
      I040 = I030 + Ndim_pphh_ab

      If (Iend .GT. Memleft) Then
          Call Insmem("@-Tdcc_built_ac_func",Iend,Memleft)
      Endif

      Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_ab_0)
      Call Getall(Work(I060),Ndim_pphh_ab,Irrepx,Mu_d_ab_t)

      ML_t_MR_0 = ML_t_MR_0 + Ddot(Ndim_ai,Work(I020),1,
     +                                Work(I030),1)

      Return 
      End
