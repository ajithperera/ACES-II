#include "flags.h"
      Subroutine Tdcc_form_mutilde_0_d_uhf(Work,Memleft,Doo,Dvv,Dvo,
     +                                     Irrepx,Lenoo,Lenvv,Lenvo,
     +                                     Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Character*4 Spcase(2)

      Double Precision Mu_jj,Mone
      Integer Hp_type,Hh_type,Pp_type,Ph_type
      Integer Tmp_list 
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/
      Data Spcase /'AAAA','BBBB'/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

C  First do the AAAA and BBBB pieces 

      L1_list = 190
      Do Ispin = 1, 2
      
         Ioff_oo = 1 + (Ispin-1)*Irpdpd(Irrepx,20+Ispin)

         Ndim_ph = Irpdpd(Irrepx,8+Ispin)
         Ndim_hp = Irpdpd(Irrepx,15+Ispin)

         Max_hphp = Max(Ndim_ph*Ndim_ph,Ndim_ph*Ndim_hp)

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max_hphp
         I020 = I010 + Irpdpd(1,8+Ispin)
         Iend = I020 
         If (Iend .Gt. Memleft) Call Insmem
     +      ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)
      
         Call Dzero(Work(I000),Max_hphp)
         Call Getlst(Work(I010),1,1,1,Ispin,L1_list)

         Iz2_t   = I000
         Ivo_l   = I010
         Ivo_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,9)
         Ivv_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioo_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,21)

C NDeP Eqn. 21, First term: Z(IB,AJ) = mu(B,J) * L(I,A) and 
C Z(ib,aj) = mu(b,j) * L(i,a)
C Z is built as Z(AI,BJ) or (ai,bj) and L and Dvo are
C stored as A,I or a,i. 

         Ioff_z2  = Iz2_t
         Do Irrep_bj = 1, Nirrep
            Irrep_ai = Dirprd(Irrep_bj,Irrepx)

            If (Irrep_bj .EQ. Irrepx) Then

               Nrow_ai = Irpdpd(Irrep_ai,8+Ispin)
               Ncol_bj = Irpdpd(Irrep_bj,8+Ispin)

               Nrow = Nrow_ai
               Ncol = Ncol_bj
               Nsum = 1
               Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ivo_l),
     +                     Nrow,Dvo(Ivo_0),1,Zero,Work(Ioff_z2),
     +                     Nrow)
               Ioff_z2 = Ioff_z2 + Nrow_ai * Ncol_bj
            Endif 
         Enddo 
         Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Ioff_z2  = Iz2_t
      do irrep_bj=1,Nirrep
         If (irrep_bj .eq.irrepx) then
         irrep_ai=dirprd(irrep_bj,irrepx)
         Ndim_ai=Irpdpd(irrep_ai,8+ispin)
         Ndim_bj=Irpdpd(irrep_bj,8+ispin)
         Ioff_z2=Iz2_t
         call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t_aa(2):",
     +               Work(Ioff_z2),Ndim_ai*Ndim_bj,S)
         Ioff_z2=Ioff_z2 + Ndim_ai*Ndim_bj
         endif
      enddo
#endif
#ifdef _DEBUG_LVLM
         Ioff_z2  = Iz2_t
         Do Irrep_bj = 1, Nirrep
            Irrep_ai = Dirprd(Irrep_bj,Irrepx)
            If (Irrep_bj .EQ. Irrepx) Then
               Nrow_ai = Irpdpd(Irrep_ai,8+Ispin)
               Ncol_bj = Irpdpd(Irrep_bj,8+Ispin)
               Do Irrep_j = 1, Nirrep
                  Irrep_b = Dirprd(Irrep_j,Irrepx)

                  Nrow_b = Vrt(Irrep_b, Ispin)
                  Ncol_j = Pop(Irrep_j, Ispin)

                  Do Irrep_i = 1, Nirrep
                     Irrep_a = Irrep_i

                     Nrow_a = Vrt(Irrep_a, Ispin)
                     Ncol_i = Pop(Irrep_i, Ispin)

                     Ivo_l  = I010  +
     +                        (Isymoff(Irrep_i,1,
     +                                8+Ispin)-1)
                     Ivo_j  = Ivo_0 +
     +                        (Isymoff(Irrep_j,Irrepx,
     +                                 8+Ispin)-1)
                     Call Vecprd_22(Nrow_a,Nrow_b,Ncol_i,
     +                              Ncol_j,Work(Ivo_l),
     +                              Dvo(Ivo_j),Work(Ioff_z2),
     +                              1.0D0)
                  Enddo
                  Ioff_z2 = Ioff_z2 + Nrow_a * Ncol_i
               Enddo
               Ioff_z2 = Ioff_z2 + Nrow_b * Ncol_j

            Endif 
         Enddo
#ifdef _DEBUG_LVL0
      Write(6,*)
      Ioff_z2  = Iz2_t
      do irrep_bj=1,Nirrep
         If (irrep_bj .eq.irrepx) then
         irrep_ai=dirprd(irrep_bj,irrepx)
         Ndim_ai=Irpdpd(irrep_ai,8+ispin)
         Ndim_bj=Irpdpd(irrep_bj,8+ispin)
         Ioff_z2=Iz2_t
         call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t_aa(2):",
     +               Work(Ioff_z2),Ndim_ai*Ndim_bj,S)
         Ioff_z2=Ioff_z2 + Ndim_ai*Ndim_bj
         endif 
      enddo
#endif
#endif

C NDeP Eqn. 20, Second term: mu(M,N)delta(M,N)*L2(IJ,AB) and 
C mu(m,n)delta(m,n)*L2(IJ,AB)
C Z is built as Z(AI,BJ) and Read L2 as (AI,BJ) (and  Z(ai,bj)
C L2(ai,bj)
         Ioff_z2 = Iz2_t 
         If (Irrepx .EQ. 1) Then 
            Do Irrep_bj = 1, Nirrep
               Irrep_ai = Dirprd(Irrep_bj,Irrepx)

               L2_list    = 133 + Ispin 
               Tmp_list   = 299 + Ispin
               Iz2_t_list = 323 + Ispin

               Len_bj = Irpdpd(Irrep_bj,8+Ispin)
               Len_ai = Irpdpd(Irrep_ai,8+Ispin)

               Ndim_aibj = Len_bj * Len_ai
            
               Ibgn = I010
               I020 = I010 + Ndim_aibj
               Iend = I020 
               If (Iend .GT. Memleft) Call Insmem
     +            ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)
   
               Call Getlst(Work(I010),1,Len_bj,1,Irrep_bj,L2_list)
            
               Ivo_l  = I010 + (Isymoff(Irrep_bj,1,Ispin)-1)
           
               Do Irrep_n = 1, Nirrep 
                  Irrep_m = Dirprd(Irrep_n,Irrepx)
                
                  Num_n = Pop(Irrep_n,Ispin)   
                  Num_m = Pop(Irrep_m,Ispin)

                  Ioo_n  = Ioo_0 - 1 + (Isymoff(Irrep_n,Irrepx,
     +                                  20+Ispin)-1)
                  Do N = 1, Num_n
                     mu_jj = Doo(Ioo_n+N+(N-1)*Num_m)
                     Call Daxpy(Ndim_aibj,mu_jj,Work(Ivo_l),1,
     +                          Work(Ioff_z2),1)
                  Enddo 
               Enddo
               Ioff_z2 = Ioff_z2 + Ndim_aibj
            Enddo 
         Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Ioff_z2=Iz2_t
      If (irrepx .eq.1) then
      do irrep_bj=1,Nirrep
         irrep_ai=dirprd(irrep_bj,irrepx)
         Ndim_ai=Irpdpd(irrep_ai,8+ispin)
         Ndim_bj=Irpdpd(irrep_bj,8+ispin)
         call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t_aa(2):",
     +                 Work(Ioff_z2),Ndim_ai*Ndim_bj,S)
         Ioff_z2=Ioff_z2 + Ndim_ai*Ndim_bj
      enddo 
      endif 
#endif

C NDep Eqn. 20, The fourth term -L2(JM,BA)*Mu(I,M) = Z(JI,BA) and
C -L2(jm,ba)*Mu(i,m) = Z(ji,ba)
C Construct Z(AI,BJ)=-L2(AI,BM)*Mu(J,M) and Z(ai,jb)=-L2(ai,bm)*Mu(j,m)
C 
         L2_list = 133 + Ispin
         Ph_type = 8   + Ispin
         Pp_type = 18  + Ispin
         Hh_type = 20  + Ispin
         Hp_type = 15  + Ispin 
         
         Do Irrep_bj = 1, Nirrep
            Irrep_ai = Dirprd(Irrep_bj,Irrepx)
            Irrep_bm = Irrep_ai

C NDeP Eqn. 21, First term: Z(IB,AJ) = mu(B,J) * L(I,A) and
C Z(ib,aj) = mu(b,j) * L(i,a)
C Z is built as Z(AI,BJ) or (ai,bj) and L and Dvo are
C stored as A,I or a,i.

            If (Irrep_bj .EQ. Irrepx) Then

               Nrow_ai = Irpdpd(Irrep_ai,8+Ispin)
               Ncol_bj = Irpdpd(Irrep_bj,8+Ispin)

               Nrow = Nrow_ai
               Ncol = Ncol_bj
               Nsum = 1
               Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ivo_l),
     +                     Nrow,Dvo(Ivo_0),1,Zero,Work(Ioff_z2),
     +                     Nrow)
            Endif 

C NDeP Eqn. 20, Second term: mu(M,N)delta(M,N)*L2(IJ,AB) and
C mu(m,n)delta(m,n)*L2(IJ,AB)
C Z is built as Z(AI,BJ) and Read L2 as (AI,BJ) (and  Z(ai,bj)
C L2(ai,bj)
            If (Irrepx .EQ. 1) Then
 
               L2_list    = 133 + Ispin
               Tmp_list   = 299 + Ispin
               Iz2_t_list = 323 + Ispin

               Len_bj = Irpdpd(Irrep_bj,8+Ispin)
               Len_ai = Irpdpd(Irrep_ai,8+Ispin)

               Ndim_aibj = Len_bj * Len_ai

               Ibgn = I010
               I020 = I010 + Ndim_aibj
               Iend = I020
               If (Iend .GT. Memleft) Call Insmem
     +            ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

               Call Getlst(Work(I010),1,Len_bj,1,Irrep_bj,L2_list)

               Ivo_l  = I010 + (Isymoff(Irrep_bj,1,Ispin)-1)

               Do Irrep_n = 1, Nirrep
                  Irrep_m = Dirprd(Irrep_n,Irrepx)

                  Num_n = Pop(Irrep_n,Ispin)
                  Num_m = Pop(Irrep_m,Ispin)

                  Ioo_n  = Ioo_0 - 1 + (Isymoff(Irrep_n,Irrepx,
     +                                  20+Ispin)-1)
                  Do N = 1, Num_n
                     mu_jj = Doo(Ioo_n+N+(N-1)*Num_m)
                     Call Daxpy(Ndim_aibj,mu_jj,Work(Ivo_l),1,
     +                          Work(Ioff_z2),1)
                  Enddo
               Enddo
               Ioff_z2 = Ioff_z2 + Ndim_aibj
            Enddo
         Endif

C NDep Eqn. 20, The fourth term -L2(JM,BA)*Mu(I,M) = Z(JI,BA) and
C -L2(jm,ba)*Mu(i,m) = Z(ji,ba)
C Construct Z(AI,BJ)=-L2(AI,BM)*Mu(J,M) and Z(ai,jb)=-L2(ai,bm)*Mu(j,m)
C
            Len_bj = Irpdpd(Irrep_bj,Ph_type)
            Len_jb = Irpdpd(Irrep_bj,Hp_type)
            Len_bm = Irpdpd(Irrep_bm,Ph_type)
            Len_ai = Irpdpd(Irrep_ai,Ph_type)
            Max_ph = Max(Len_bj,Len_bm,Len_ai,Len_mb)

            I010 = Ibgn 
            I020 = I010 + Len_ai * Len_bm
            I030 = I020 + Max_ph
            I040 = I030 + Max_ph
            I050 = I040 + Max_ph
            Iend = I050
            If (Iend .GT. Memleft) Call Insmem
     +         ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

            Call Getlst(Work(I010),1,Len_bm,1,Irrep_bm,L2_list)

C Z(AI,BJ)=-L2(AI,BM)*Mu(J,M) and Z(ai,jb)=-L2(ai,bm)*Mu(j,m)

            Do irrep_m = 1, Nirrep
               Irrep_j = Dirprd(Irrep_m,Irrepx)
               Irrep_b = Dirprd(Irrep_m,Irrep_bm)
    
               Ndim_m = Pop(Irrep_m,Ispin)
               Ndim_j = Pop(Irrep_j,Ispin)
               Ndim_b = Vrt(Irrep_b,Ispin)

               Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,20+ispin)-1)
               Ioff_l2 = I010  + Len_ai*(Isymoff(Irrep_m,Irrep_bm,
     +                           Ph_type)-1)
               Iz2_t   = I000  + Len_ai*(Isymoff(Irrep_j,Irrep_bj,
     +                           Ph_type)-1)

               Nrow = Len_ai * Ndim_b
               Ncol = Ndim_j
               Nsum = Ndim_m

               Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                     Nrow,Doo(Ioff_oo),Nsum,One,Work(Iz2_t),
     +                     Nrow)
            Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t_aa(2):",
     +               Work(Iz2_t),Ndim_ai*Ndim_bj,S)
#endif 

C Z(AI,BJ) -> Z(AI,JB) and Z(ai,bj) -> Z(ai,jb)

            Call Symtr1(Irrep_bj,Vrt(1,Ispin),Pop(1,Ispin),Len_ai,
     +                  Work(I000),Work(I020),Work(I030),
     +                  Work(I040))

C L2(AI,EJ)-> L2(AI,JE) and L2(ai,ej)-> L2(ai,je)
C Irrep_je = Irrep_bj

            Call Symtr1(Irrep_bj,Vrt(1,Ispin),Pop(1,Ispin),Len_ai,
     +                  Work(I010),Work(I020),Work(I030),
     +                  Work(I040))

C NDep Eqn. 20, The third term L2(IJ,AE)*Mu(B,E) = Z(IJ,AB) and
C L2(ij,ae)*Mu(b,e) = Z(IJ,AB)
C
C Z(AI,JB)=L2(AI,JE)*Mu(B,E) and Z(ai,jb)=L2(ai,je)*Mu(b,e)

            Do Irrep_e = 1, Nirrep
               Irrep_b = Dirprd(Irrep_e,Irrepx)
               Irrep_j = Dirprd(Irrep_e,Irrep_bj)
        
               Ndim_e = Vrt(Irrep_e,Ispin)
               Ndim_b = Vrt(Irrep_b,Ispin)
               Ndim_j = Pop(Irrep_j,Ispin)

               Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Pp_type)-1)
               Ioff_l2 = I010  + Len_ai*(Isymoff(Irrep_e,Irrep_bj,
     +                           Hp_type)-1)
               Iz2_t   = I000  + Len_ai**(Isymoff(Irrep_b,Irrep_bj,
     +                           Hp_type)-1)
               
               Nrow = Len_ai * Ndim_j
               Ncol = Ndim_b
               Nsum = Ndim_e
               
               Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                     Nrow,Dvv(Ioff_vv),Nsum,One,Work(Iz2_t),
     +                     Nrow)
            Enddo 

C Z(AI,JB) -> Z(AI,BJ) and Z(ai,jb) -> Z(ai,bj)
C Irrep_jb = Irrep_bj

            Call Symtr1(Irrep_bj,Pop(1,Ispin),Vrt(1,Ispin),Len_ai,
     +                  Work(I000),Work(I020),Work(I030),Work(I040))

            Call Putlst(Work(I000),1,Len_bj,1,Irrep_bj,Tmp_list) 

        Enddo

        Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                       Isytyp(2,Tmp_list))
        Ibgn = 1
        I000 = Ibgn
        I010 = I000 + Isize_phph
        I020 = I010 + Isize_phph
        I030 = I020 + 100
        Iend = I030 
        If (Iend .GT. Memleft) Call Insmem
     +     ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

        Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)

        Call Dsst03i(Work(I000),Work(I010),Isize_phph, 
     +               Isize_Phph,Work(I020),Spcase(Ispin),Irrepx)

        Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_t_list),
     +                       Isytyp(2,Iz2_t_list))

        Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_t_list)
        
      Enddo 

C The ABAB block.

      Tmp_list   = 302 
      Iz2_t_list = 326

      Ndim_ph_a = Irpdpd(Irrepx,9)
      Ndim_ph_b = Irpdpd(Irrepx,10)
      Ndim_hp_a = Irpdpd(Irrepx,16)
      Ndim_hp_b = Irpdpd(Irrepx,17)

      Max_hphp = Max(Ndim_ph_a*Ndim_ph_b,Ndim_ph_a*Ndim_hp_b,
     +               Ndim_hp_a*Ndim_hp_b)

      Ibgn = 1
      I000 = Ibgn
      I010 = I000 + Max_hphp
      I020 = I010 + Irpdpd(1,9)
      Iend = I020

      Call Getlst(Work(I010),1,1,1,1,L1_list)

      Ivo_l   = I010
      Iz2_t   = I000
      Ivo_0   = 1 +  Irpdpd(Irrepx,9)
      Ivv_0   = 1 +  Irpdpd(Irrepx,19)
      Ioo_0   = 1 +  Irpdpd(Irrepx,21)

C NDeP Eqn. 21, First term: Z(Ib,Aj) = mu(b,j) * L(I,A) and
C Z is built as Z(AI,bj) or (ai,bj) and L and Dvo are
C stored as A,I and b,j.

      Do Irrep_bj = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_bj,Irrepx)

         Nrow_ai = Irpdpd(Irrep_ai,9)
         Ncol_bj = Irpdpd(Irrep_bj,10)

         Do Irrep_j = 1, Nirrep
            Irrep_b = Dirprd(Irrep_j,Irrepx)

            Nrow_b = Vrt(Irrep_a, 2)
            Ncol_j = Pop(Irrep_i, 2)

            Do Irrep_i = 1, Nirrep
               Irrep_a = Irrep_i
  
               Nrow_a = Vrt(Irrep_a, 1)
               Ncol_i = Pop(Irrep_i, 1)

               Ivo_l  = I000  +
     +                  (Isymoff(Irrep_i,1,9)-1)
               Ivo_j  = Ivo_0 +
     +                  (Isymoff(Irrep_j,Irrepx,
     +                   10)-1)

               Call Vecprd_22(Nrow_a,Nrow_b,Ncol_i,Ncol_j,
     +                        Work(Ivo_l),Dvo(Ivo_j),
     +                        Work(Iz2_t),1.0D0)
            Enddo
         Enddo
         Iz2_t = Iz2_t + Nrow_ai * Ncol_bj
      Enddo

      L2_list = 137

      Ioff_l2 = I010 
      Iz2_t   = I000

      Do Irrep_bj = 1, Nirrep
         Irrep_ai = Irrep_bj
          
         Ndim_bj   = Irpdpd(Irrep_bj,10)
         Ndim_ai   = Irpdpd(Irrep_ai,9)
         Ndim_aibj = Ndim_bj * Ndim_ai

         Ibgn = I010
         I020 = I010 + Ndim_aibj 
         Iend = I020 

         If (Iend .GT. Memleft) Call Insmem
     +      ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

            Call Getlst(Work(I010),1,Ndim_bj,1,Irrep_bj,L2_list)

C NDeP Eqn. 20, Second term. mu(m,n)delta(m,n)*L2(AI,bj) = Z(AI,bj)

            Ivo_l  = I010 + (Isymoff(Irrep_bj,1,10)-1)

            Do Irrep_n = 1, Nirrep
               Irrep_m = Dirprd(Irrep_n,Irrepx)

               Num_n = Pop(Irrep_n,Ispin)
               Num_m = Pop(Irrep_m,Ispin)

               Ndim_nm = Num_n * Num_m

               Ioo_n  = Ioo_0 + (Isymoff(Irrep_n,Irrepx,
     +                                   22)-1)
               Do N = 1, Num_n
                 mu_jj = Doo(Ioo_n+N+(N-1)*Num_m)
                 Call Daxpy(Ndim_nm,mu_jj,Work(Ivo_l),1,Work(Iz2_t),1)
               Enddo
               Iz2_t = Iz2_t + Ndim_aibj
            Enddo
         Enddo

C NDep Eqn. 20, The fourth term; Construct Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 

      L2_list = 137
      Pha_type = 9
      Phb_type = 10
      Hpb_type = 17
      Ppb_type = 20
      hha_type = 21
      hhb_type = 22
         
      Do Irrep_bj = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_bj,Irrepx)
         Irrep_bm = Irrep_ai
  
         Len_bj = Irpdpd(Irrep_bj,Phb_type)
         Len_jb = Irpdpd(Irrep_bj,Hpb_type)
         Len_ai = Irpdpd(Irrep_ai,Pha_type)
         Len_ia = Irpdpd(Irrep_ai,Hpa_type)
         Len_bm = Irpdpd(Irrep_bm,Phb_type)
         Max_ph = Max(Len_bj,Len_bm,Len_ai,Len_jb) 

         I010 = Ibgn 
         I020 = I010 + Len_ai * Len_bm
         I030 = I020 + Max_ph
         I040 = I030 + Max_ph
         I050 = I040 + Max_ph
         Iend = I050

         If (Iend .GT. Memleft) Call Insmem
     +      ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

         Call Getlst(Work(I010),1,Len_bj,1,Irrep_bj,L2_list)

C Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 

         Do irrep_m = 1, Nirrep
            Irrep_j = Dirprd(Irrep_m,Irrepx)
            Irrep_b = Dirprd(Irrep_m,Irrep_bm)
    
            Ndim_m = Pop(Irrep_m,2)
            Ndim_j = Pop(Irrep_j,2)
            Ndim_b = Vrt(Irrep_b,2)

            Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,22)-1)
            Ioff_l2 = I010  + Len_ai*(Isymoff(Irrep_m,Irrep_bm,
     +                                Phb_type)-1)
            Iz2_t   = I000  + Len_ai*(Isymoff(Irrep_j,Irrep_bj,
     +                                Phb_type)-1)

            Nrow = Len_ai * Ndim_b
            Ncol = Ndim_j
            Nsum = Ndim_m

            Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                 Nrow,Doo(Ioff_oo),Nsum,One,Work(Iz2_t),
     +                 Nrow)
         Enddo

C Z(AI,bj) -> Z(AI,jb) 

          Call Symtr1(Irrep_bj,Vrt(1,2),Pop(1,2),Len_ai,
     +                Work(I000),Work(I020),Work(I030),
     +                Work(I040))

C L2(AI,ej)-> L2(AI,je)
C Irrep_ej = Irrep_bj

          Call Symtr1(Irrep_bj,Vrt(1,2),Pop(1,2),Len_ai,
     +                Work(I010),Work(I020),Work(I030),
     +                Work(I040))

C NDep Eqn. 20, third term Z(AI,jb)=L2(AI,ej)*Mu(b,e) 

          Do Irrep_e = 1, Nirrep
             Irrep_b = Dirprd(Irrep_e,Irrepx)
             Irrep_j = Dirprd(Irrep_e,Irrep_bj)
          
             Ndim_e = Vrt(Irrep_e,2)
             Ndim_b = Vrt(Irrep_b,2)
             Ndim_j = Pop(Irrep_j,2)
  
             Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Ppb_type)-1)
             Ioff_l2 = I010  + Len_ai*(Isymoff(Irrep_e,Irrep_bj,
     +                                         Phb_type)-1)
             Iz2_t   = I000  + Len_ai**(Isymoff(Irrep_b,Irrep_bj,
     +                                         Hpb_type)-1)
            
             Nrow = Len_ai * Ndim_j
             Ncol = Ndim_b
             Nsum = Ndim_e
              
             Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                   Nrow,Dvv(Ioff_vv),Nsum,One,Work(Iz2_t),
     +                   Nrow)
          Enddo 

          Ibgn = I020
          I030 = I020 + Len_ai * Len_jb
          Iend = I030
          If (Iend .GT. Memleft) Call Insmem
     +       ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

C Z(AI,jb) -> Z(jb,AI)
        
          Call Transp(Work(I000),Work(I020),Len_bj,Len_ai)
          Call Dcopy(Len_ai*Len_jb,Work(I020),1,Work(I000),1)

C NDep Eqn. 20, The fourth term; Construct Z(bj,AI)=-L2(bj,AM)*Mu(I,M)
C Read L2 as L2(AM,bj)

          Call Getlst(Work(I010),1,Len_bj,1,Irrep_bj,L2_list)

C L(AM,bj) ->  L(AM,jb)
C Irrep_jb = Irrep_bj

          Len_am = Len_ai

          Ibgn = I030
          I040 = I030 + Max_ph
          I050 = I040 + Max_ph
          I060 = I050 + Max_ph
          Iend = I060
          If (Iend .GT. Memleft) Call Insmem
     +       ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)
   
          Call Symtr1(Irrep_bj,Vrt(1,2),Pop(1,2),Len_am,
     +                Work(I010),Work(I030),Work(I040),
     +                Work(I050))

C L(AM,jb) -> L(jb,AM)
        
          Call Transp(Work(I010),Work(I020),Len_jb,Len_am)
          Call Dcopy(Len_am*Len_jb,Work(I020),1,Work(I010),1)

C Construct Z(jb,AI)=-L2(jb,AM)*Mu(I,M)
C Irrep_am = Irrep_ai

           Do irrep_m = 1, Nirrep
              Irrep_i = Dirprd(Irrep_m,Irrepx)
              Irrep_a = Dirprd(Irrep_m,Irrep_ai)

              Ndim_m = Pop(Irrep_m,1)
              Ndim_i = Pop(Irrep_i,1)
              Ndim_a = Vrt(Irrep_a,1)

              Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,21)-1)
              Ioff_l2 = I010  + Len_jb*(Isymoff(Irrep_m,Irrep_ai,
     +                                          Pha_type)-1)
              Iz2_t   = I000  + Len_jb*(Isymoff(Irrep_i,Irrep_ai,
     +                                          Pha_type)-1)
              Nrow = Len_jb * Ndim_a
              Ncol = Ndim_i
              Nsum = Ndim_m

              Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                    Nrow,Doo(Ioff_oo),Nsum,One,Work(Iz2_t),
     +                    Nrow)
           Enddo

C Z(jb,AI) -> Z(jb,IA)
        
           Call Symtr1(Irrep_ai,Vrt(1,1),Pop(1,1),Len_jb,
     +                 Work(I000),Work(I020),Work(I030),
     +                 Work(I040))

C L(jb,EI) -> L(jb,IE)
C Irrep_ei = Irrep_ai,Len_ei=Len_ai

           Call Symtr1(Irrep_ai,Vrt(1,1),Pop(1,1),Len_jb,
     +                 Work(I010),Work(I020),Work(I030),
     +                 Work(I040))

C NDep Eqn. 20, third term Z(jb,IA)=L2(jb,IE)*Mu(A,E)

           Do Irrep_e = 1, Nirrep
              Irrep_a = Dirprd(Irrep_e,Irrepx)
              Irrep_i = Dirprd(Irrep_e,Irrep_ai)

              Ndim_e = Vrt(Irrep_e,1)
              Ndim_a = Vrt(Irrep_a,1)
              Ndim_i = Pop(Irrep_j,1)

              Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Ppa_type)-1)
              Ioff_l2 = I010  + Len_jb*(Isymoff(Irrep_i,Irrep_ai,
     +                                          Hpa_type)-1)
              Iz2_t   = I000  + Len_jb*(Isymoff(Irrep_a,Irrep_ai,
     +                                  Hpa_type)-1)

              Nrow = Len_bj * Ndim_a
              Ncol = Ndim_a
              Nsum = Ndim_e

              Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                    Nrow,Dvv(Ioff_vv),Nsum,One,Work(Iz2_t),
     +                    Nrow)
           Enddo

C Z(jb,IA) -> Z(IA,jb) -> Z(IA,bj) -> Z(AI,bj)

           Call Transp(Work(I000),Work(I010),Len_ia,Len_jb)
           Call Symtr1(Irrep_bj,Pop(1,2),Vrt(1,2),Len_ia,Work(I010),
     +                 Work(I020),Work(I030),Work(I040))

           Call Symtr3(Irrep_ai,Pop(1,1),Vrt(1,1),Irrep_ai,Irrep_bj,
     +                 Work(I010),Work(I020),Work(I030),Work(I040))

           Call Putlst(Work(I010),1,Len_bj,1,Irrep_bj,Tmp_list)

      Enddo 

      Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                            Isytyp(1,Tmp_list))

      Ibgn = 1
      I000 = Ibgn
      I010 = I000  + Isize_phph
      I020 = I010  + Isize_phph
      I030 = I020  + 100 
      Iend = I030 

      If (Iend .GT. Memleft) Call Insmem
     +   ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

      Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)
      Call Sstgen(Work(I000),Work(I010),Isize_phph,Vrt(1,1),Pop(1,1),
     +            Vrt(1,2),Pop(1,2),Work(1020),Irrepx,'1324')

      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_t_list),
     +                            Isytyp(2,Iz2_t_list))

      Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_t_list)

      Return
      End
        
