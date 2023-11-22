#include "flags.h"
      Subroutine Tdcc_form_mutilde_0_d_rhf(Work,Memleft,Doo,Dvv,Dvo,
     +                                     Irrepx,Lenoo,Lenvv,Lenvo,
     +                                     Iuhf,Mu_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Character*4 Spcase(2)

      Double Precision Mu_0,Mone
      Integer Phb_type,Pha_type,Hpb_type,Ppb_type,Hha_type
      Integer Hhb_type,Hpa_type
      Integer Tmp_list 
      Integer Idid(8)
 
      Data Zero, One, Mone, Two, Half /0.0D0, 1.0D0, -1.0D0, 2.0D0,
     +                                 0.50D0/
      Data Spcase /'AAAA','BBBB'/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_0_d_rhf--"
        Write(6,*)
#endif

C The ABAB block.

      Tmp_list   = 302 
      Iz2_t_list = 326
      L2_list    = 137
      L1_list    = 190

      Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                     Isytyp(2,Tmp_list)) 
      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_t_list),
     +                     Isytyp(2,Iz2_t_list))
      
      Phb_type = 10
      Pha_type = 9
      Hpa_type = 16
      Hpb_type = 17
      Ppb_type = 20
      hha_type = 21
      hhb_type = 22

      Ivo_0   = 1 
      Ivv_0   = 1
      Ioo_0   = 1 

C NDeP Eqn. 21, First term: Z(Ib,Aj) = mu(b,j) * L(A,I) and
C Z is built as Z(AI,bj). 

      Do Irrep_bj_z = 1, Nirrep
         Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)
         Irrep_ai_l = Irrep_ai_z
         Irrep_bj_l = Irrep_ai_l

         Ndim_ai_z   = Irpdpd(Irrep_ai_z,Pha_type)
         Ndim_bj_z   = Irpdpd(Irrep_bj_z,Phb_type)
         Ndim_jb_z   = Irpdpd(Irrep_bj_z,Hpb_type)
         Ndim_ai_l   = Irpdpd(Irrep_ai_l,Pha_type)
         Ndim_bj_l   = Irpdpd(Irrep_bj_l,Phb_type)

         Ndim_aibj_l = Ndim_ai_l * Ndim_bj_l 
         Ndim_aibj_z = Ndim_bj_z * Ndim_ai_z

         Max_ph = Max(Ndim_ai_z,Ndim_bj_z,Ndim_jb_z,Ndim_ai_l,
     +                Ndim_bj_l)
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_aibj_z
         I020 = I010 + Irpdpd(1,9)
         Iend = I020

         Ivo_l = I010
         Iz2_t = I000
  
         Call Dzero(Work(I000),Ndim_aibj_z)

         Call Getlst(Work(I010),1,1,1,1,L1_list)
   
         If (Irrep_bj_z .EQ. Irrepx) Then

            Nrow = Ndim_ai_z
            NCol = Ndim_bj_z
            Nsum = 1
            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ivo_l),
     +                  Nrow,Dvo(Ivo_0),1,One,Work(Iz2_t),Nrow)
         Endif 

#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_rhf,mu^t_ab(1) :",
     +               Work(I000),Ndim_aibj_z,S)
#endif
         Ibgn = I010
         I020 = I010 + Ndim_aibj_l
         I030 = I020 + Max_ph
         I040 = I030 + Max_ph
         I050 = I040 + Max_ph
         Iend = I050

         If (Iend .GT. Memleft) Call Insmem
     +      ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

         Call Getlst(Work(I010),1,Ndim_bj_l,1,Irrep_bj_l,L2_list)

C NDeP Eqn. 20, Second term. mu(m,n)delta(m,n)*L2(AI,bj) = Z(AI,bj)

         If (Irrepx .EQ. 1) Then
             Call Daxpy(Ndim_aibj_z,Half*mu_0,Work(Ivo_l),1,
     +                  Work(Iz2_t),1)

         Endif 

#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_rhf,mu^t_ab(2) :",
     +               Work(I000),Ndim_aibj_z,S)
#endif

C NDep Eqn. 20, The fourth term; Construct Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 
C Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 

         Irrep_bm_l = Irrep_bj_l
         Ndim_bm_l  = Ndim_bj_l

         Do irrep_m = 1, Nirrep
            Irrep_j = Dirprd(Irrep_m,Irrepx)
            Irrep_b = Dirprd(Irrep_m,Irrep_bm_l)
    
            Ndim_m = Pop(Irrep_m,2)
            Ndim_j = Pop(Irrep_j,2)
            Ndim_b = Vrt(Irrep_b,2)

            Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,22)-1)
            Ioff_l2 = I010  + Ndim_ai_l*(Isymoff(Irrep_m,Irrep_bm_l,
     +                                Phb_type)-1)
            Iz2_t   = I000  + Ndim_ai_z*(Isymoff(Irrep_j,Irrep_bj_z,
     +                                Phb_type)-1)

            Nrow = Ndim_ai_l * Ndim_b
            Ncol = Ndim_j
            Nsum = Ndim_m
            Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                 Nrow,Doo(Ioff_oo),Ncol,One,Work(Iz2_t),
     +                 Nrow)
         Enddo
#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_rhf,mu^t_ab(3) :",
     +               Work(I000),Ndim_aibj_z,S)
#endif

C Z(AI,bj) -> Z(AI,jb) 

          Call Symtr1(Irrep_bj_z,Vrt(1,2),Pop(1,2),Ndim_ai_z,
     +                Work(I000),Work(I020),Work(I030),
     +                Work(I040))

C L2(AI,ej)-> L2(AI,je)
C Irrep_ej_l = Irrep_bj_l=irrep_bm_l

          Irrep_ej_l = Irrep_bm_l
          Call Symtr1(Irrep_ej_l,Vrt(1,2),Pop(1,2),Ndim_ai_l,
     +                Work(I010),Work(I020),Work(I030),
     +                Work(I040))

C NDep Eqn. 20, third term Z(AI,jb)=L2(AI,je)*Mu(b,e) 

          Do Irrep_e = 1, Nirrep
             Irrep_b = Dirprd(Irrep_e,Irrepx)
             Irrep_j = Dirprd(Irrep_e,Irrep_ej_l)
          
             Ndim_e = Vrt(Irrep_e,2)
             Ndim_b = Vrt(Irrep_b,2)
             Ndim_j = Pop(Irrep_j,2)
  
             Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Ppb_type)-1)
             Ioff_l2 = I010  + Ndim_ai_l*(Isymoff(Irrep_e,Irrep_ej_l,
     +                                         Hpb_type)-1)
             Iz2_t   = I000  + Ndim_ai_z*(Isymoff(Irrep_b,Irrep_bj_z,
     +                                         Hpb_type)-1)
             Nrow = Ndim_ai_l * Ndim_j
             Ncol = Ndim_b
             Nsum = Ndim_e
              
             Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                   Nrow,Dvv(Ioff_vv),Ncol,One,Work(Iz2_t),
     +                   Nrow)
          Enddo 
#ifdef _DEBUG_LVLM
      Write(6,*)
      Ndim_aijb_z = Ndim_ai_z * Ndim_jb_z
      call checksum("Tdcc_form_mutilde_0_d_rhf,mu^t_ab(4) :",
     +               Work(I000),Ndim_aijb_z,S)
#endif

C Z(AI,jb)  -> Z(AI,bj)
C L2(AI,je) -> L2(AI,ej)

          Call Symtr1(Irrep_bj_z,Pop(1,2),Vrt(1,2),Ndim_ai_z,
     +                Work(I000),Work(I020),Work(I030),
     +                Work(I040))
          Call Symtr1(Irrep_ej_l,Pop(1,2),Vrt(1,2),Ndim_ai_l,
     +                Work(I010),Work(I020),Work(I030),
     +                Work(I040))

          Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,Tmp_list)
  
      Enddo 

C Form Z(AI,bj) + Z(bj,AI)

      Call Izero(Idid,8)
      Do Irrep_bj_z = 1, Nirrep

         If (Idid(Irrep_bj_z) .EQ. 0) Then
             Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)
 
             Ndim_ai_z = Irpdpd(Irrep_ai_z,Pha_type)
             Ndim_bj_z = Irpdpd(Irrep_bj_z,Phb_type)
       
             Ibgn = 1
             I000 = Ibgn
             I010 = I000 + Ndim_ai_z * Ndim_bj_z
             I020 = I010 + Ndim_ai_z * Ndim_bj_z
             Iend = I020
             If (Iend .GT. Memleft) Call Insmem
     +          ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)
  
             If (Irrepx .EQ. 1) Then

                Call Getlst(Work(I000),1,Ndim_ai_z,1,Irrep_ai_z,
     +                      Tmp_list)
                Call Transp(Work(I000),Work(I010),Ndim_ai_z,
     +                      Ndim_bj_z)
                Call Daxpy(Ndim_ai_z*Ndim_bj_z,One,Work(I010),1,
     +                     Work(I000),1)
                Call Putlst(Work(I000),1,Ndim_ai_z,1,Irrep_ai_z,
     +                      Tmp_list)
             Else
                Call Getlst(Work(I000),1,Ndim_ai_z,1,Irrep_ai_z,
     +                      Tmp_list)
                Call Transp(Work(I000),Work(I010),Ndim_ai_z,
     +                      Ndim_bj_z)
                Call Getlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,
     +                     Tmp_list)
                Call Daxpy(Ndim_ai_z*Ndim_bj_z,One,Work(I010),1,
     +                     Work(I000),1)

                Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,
     +                      Tmp_list)
                Call Transp(Work(I000),Work(I010),Ndim_bj_z,
     +                      Ndim_ai_z)
                Call Putlst(Work(I000),1,Ndim_ai_z,1,Irrep_ai_z,
     +          Tmp_list)
             Endif 
             Idid (Irrep_bj_z) = 1 
             Idid (Irrep_ai_z) = 1 
         Endif 
      Enddo 

      Ibgn = 1
      I000 = Ibgn
      I010 = I000  + Isize_phph
      I020 = I010  + Isize_phph
      I030 = I020  + Ndim_ai_z
      Iend = I030 
      If (Iend .GT. Memleft) Call Insmem
     +   ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

C Z(AI,bj) -> Z(Ab,Ij)

      Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)
      Call Sstgen(Work(I000),Work(I010),Isize_phph,Vrt(1,1),
     +            Pop(1,1),Vrt(1,1),Pop(1,1),Work(I020),
     +            Irrepx,'1324')
#ifdef _DEBUG_LVL0
      call checksum("Tdcc_form_mutilde_0_d_rhf,mu^t_ab(f) :",
     +               Work(I010),Isize_pphh,S)
#endif
      Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_t_list)

      Return
      End

