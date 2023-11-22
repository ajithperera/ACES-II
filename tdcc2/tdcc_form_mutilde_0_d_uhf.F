#include "flags.h"
      Subroutine Tdcc_form_mutilde_0_d_uhf(Work,Memleft,Doo,Dvv,Dvo,
     +                                     Irrepx,Lenoo,Lenvv,Lenvo,
     +                                     Iuhf,Mu_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Character*4 Spcase(2)

      Double Precision Mu_0,Mone
      Integer Hp_type,Hh_type,Pp_type,Ph_type
      Integer Pha_type,Phb_type,Hpb_type,Hpa_type,Ppb_type,Ppa_type
      Integer Hha_type,Hhb_type
      Integer Tmp_list 

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

      Data Zero, One, Mone, Two, Half /0.0D0, 1.0D0, -1.0D0, 2.0D0,
     +                                 0.50D0/
      Data Spcase /'AAAA','BBBB'/

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_0_d_uhf---"
        Write(6,*)
#endif
C  First do the AAAA and BBBB pieces 
 
      L1_list = 190
      Do Ispin = 1, 2
      
         Ivo_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,9)
         Ivv_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioo_0   = 1 + (Ispin-1) * Irpdpd(Irrepx,21)

         L2_list    = 133 + Ispin
         Ph_type    = 8   + Ispin
         Pp_type    = 18  + Ispin
         Hh_type    = 20  + Ispin
         Hp_type    = 15  + Ispin 
         Tmp_list   = 299 + Ispin
         Iz2_t_list = 323 + Ispin 

C Target, Z(AI,Bj): Ispin =1
C Target, Z(ai,bj): Ispin =2
         
         Do Irrep_bj_z = 1, Nirrep
            Irrep_ai_z = Dirprd(Irrep_bj_z,Irrepx)
            Irrep_ai_l = Irrep_ai_z
            Irrep_bj_l = Irrep_ai_l

            Ndim_ai_z = Irpdpd(Irrep_ai_z,Ph_type)
            Ndim_bj_z = Irpdpd(Irrep_bj_z,Ph_type)
            Ndim_jb_z = Irpdpd(Irrep_bj_z,Hp_type)
            Ndim_ai_l = Irpdpd(Irrep_ai_l,Ph_type)
            Ndim_bj_l = Irpdpd(Irrep_bj_l,Ph_type)
            Max_pq  = Max(Ndim_ai_z,Ndim_bj_z,Ndim_jb_z,Ndim_ai_l,
     +                    Ndim_bj_l)

            Max_aibj    = Max_pq  * Max_pq
            Ndim_aibj_z = Ndim_ai_z * Ndim_bj_z
            Ndim_aibj_l = Ndim_ai_l * Ndim_bj_l

            Ibgn = 1
            I000 = Ibgn
            I010 = I000 + Max_aibj
            I020 = I010 + Irpdpd(1,8+Ispin)
            Iend = I020
            If (Iend .Gt. Memleft) Call Insmem
     +         ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)

            Call Dzero(Work(I000),Max_aibj)
            Call Getlst(Work(I010),1,1,1,Ispin,L1_list)

            Iz2_t   = I000
            Ivo_l   = I010

C NDeP Eqn. 21, First term: Z(IB,AJ) = mu(B,J) * L(I,A) and
C Z(ib,aj) = mu(b,j) * L(i,a)
C Z is built as Z(AI,BJ) or (ai,bj) and L and Dvo are
C stored as A,I or a,i.

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,a,6(1x,i4))")"Ph,Pp,Hp_type,Irrep_bj_z,Irrep_ai_z",
     +        "Irrep_ai_l:",Ph_type,Pp_type,Hp_type,Irrep_bj_z,
     +         Irrep_ai_z,Irrep_ai_l
      write(6,"(a,a,4(1x,i4))")"Irrep_bj_l,Ndim_ai_z,Ndim_bj_z,",
     +        "Ndim_jb_z:",
     +         Irrep_bj_l,Ndim_ai_z,Ndim_bj_z,Ndim_jb_z
      write(6,"(a,5(1x,i4))")"Ndim_bj_l,Max_pq,Ivo_0,Ioo_0,Ivv_0:",
     +        Ndim_bj_l,Max_pq,Ivo_0,Ioo_0,Ivv_0 
      write(6,"(a,a,5(1x,i4))")"Max_aibj,Ndim_aibj_z,Ndim_aibj_l,",
     +        "Iz2_t,Ivo_l:",
     +         Max_aibj,Ndim_aibj_z,Ndim_aibj_l,Iz2_t,Ivo_l
#endif 
C The term Mut(B,J) L(I,A) performed as
C Z(AI,BJ) = Mut(B,J) * L(A,I) ISPIN=1
C Z(ai,bj) = Mut(b,j) * L(a,i) ISPIN=2

            If (Irrep_bj_z .EQ. Irrepx) Then

               Nrow = Ndim_ai_z
               Ncol = Ndim_bj_z
               Nsum = 1
               Call Xgemm("N","N",Nrow,Ncol,Nsum,-Two,Work(Ivo_l),
     +                     Nrow,Dvo(Ivo_0),1,One,Work(Iz2_t),
     +                     Nrow)
            Endif 

#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(AA/bb)(1):",
     +               Work(I000),Ndim_aibj_z,S)
#endif 
C L2(AI,BM), Ispin=1, L2(ai,bm) Ispin=2
C Irrep_bm_l=Irrep_bj_l,Ndim_bm_l=Ndim_bj_l

            Ibgn = I010
            I020 = I010 + Ndim_aibj_l
            I030 = I020 + Max_pq
            I040 = I030 + Max_pq
            I050 = I040 + Max_pq

            If (Iend .GT. Memleft) Call Insmem
     +         ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

            Irrep_bm_l = Irrep_bj_l
            Ndim_bm_l  = Ndim_bj_l
            Call Getlst(Work(I010),1,Ndim_bm_l,1,Irrep_bm_l,L2_list)

C NDeP Eqn. 20, Second term: mu(M,N)delta(M,N)*L2(IJ,AB) and
C mu(m,n)delta(m,n)*L2(IJ,AB)
C Z is built as Z(AI,BJ) and Read L2 as (AI,BJ) (and  Z(ai,bj)
C L2(ai,bj)
            If (Irrepx .EQ. 1) Then
               
               Call Daxpy(Ndim_aibj_l,Half*mu_0,Work(I010),1,
     +                    Work(Iz2_t),1)

            Endif
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,Mu^t(AA,bb)(2):",
     +               Work(I000),Ndim_aibj_z,S)
#endif 

C NDep Eqn. 20, The third and fouth terms -L2(JM,BA)*Mu(I,M) = Z(JI,BA) and
C -L2(jm,ba)*Mu(i,m) = Z(ji,ba)
C Construct Z(AI,BJ)=-L2(AI,BM)*Mu(J,M) and Z(ai,jb)=-L2(ai,bm)*Mu(j,m)
C
C Z(AI,BJ)=-L2(AI,BM)*Mu(J,M) and Z(ai,bj)=-L2(ai,bm)*Mu(j,m)

            Do irrep_m = 1, Nirrep
               Irrep_j = Dirprd(Irrep_m,Irrepx)
               Irrep_b = Dirprd(Irrep_m,Irrep_bm_l)
    
               Ndim_m = Pop(Irrep_m,Ispin)
               Ndim_j = Pop(Irrep_j,Ispin)
               Ndim_b = Vrt(Irrep_b,Ispin)

               Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,Hh_type)-1)
               Ioff_l2 = I010  + Ndim_ai_l*(Isymoff(Irrep_m,Irrep_bm_l,
     +                                      Ph_type)-1)
               Iz2_t   = I000  + Ndim_ai_z*(Isymoff(Irrep_j,Irrep_bj_z,
     +                                      Ph_type)-1)

               Nrow = Ndim_ai_l * Ndim_b
               Ncol = Ndim_j
               Nsum = Ndim_m

               Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                     Nrow,Doo(Ioff_oo),Ncol,One,Work(Iz2_t),
     +                     Nrow)
            Enddo
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,Mu^t(AA,bb)(3):",
     +               Work(I000),Ndim_aibj_z,S)
#endif 

C Z(AI,BJ) -> Z(AI,JB) and Z(ai,bj) -> Z(ai,jb)

            Call Symtr1(Irrep_bj_z,Vrt(1,Ispin),Pop(1,Ispin),
     +                  Ndim_ai_z,Work(I000),Work(I020),Work(I030),
     +                  Work(I040))

C L2(AI,EJ)-> L2(AI,JE) and L2(ai,ej)-> L2(ai,je)
C Irrep_ej_l=Irrep_bj_l=irrep_ai_z

            Irrep_ej_l = Irrep_bj_l
            Call Symtr1(Irrep_ej_l,Vrt(1,Ispin),Pop(1,Ispin),
     +                  Ndim_ai_l,Work(I010),Work(I020),Work(I030),
     +                  Work(I040))

C NDep Eqn. 20, The third term L2(IJ,AE)*Mu(B,E) = Z(IJ,AB) and
C L2(ij,ae)*Mu(b,e) = Z(IJ,AB)
C
C Z(AI,JB)=L2(AI,JE)*Mu(B,E) and Z(ai,jb)=L2(ai,je)*Mu(b,e)

            Do Irrep_e = 1, Nirrep
               Irrep_b = Dirprd(Irrep_e,Irrepx)
               Irrep_j = Dirprd(Irrep_e,Irrep_ej_l)
        
               Ndim_e = Vrt(Irrep_e,Ispin)
               Ndim_b = Vrt(Irrep_b,Ispin)
               Ndim_j = Pop(Irrep_j,Ispin)

               Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Pp_type)-1)
               Ioff_l2 = I010  + Ndim_ai_l*(Isymoff(Irrep_e,Irrep_ej_l,
     +                                      Hp_type)-1)
               Iz2_t   = I000  + Ndim_ai_z*(Isymoff(Irrep_b,Irrep_bj_z,
     +                                      Hp_type)-1)
               
               Nrow = Ndim_ai_l * Ndim_j
               Ncol = Ndim_b
               Nsum = Ndim_e
               
               Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                     Nrow,Dvv(Ioff_vv),Ncol,One,Work(Iz2_t),
     +                     Nrow)
            Enddo 
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,Mu^t(AA,bb)(4):",
     +               Work(I000),Ndim_aibj_z,S)
#endif 

C Z(AI,JB) -> Z(AI,BJ) and Z(ai,jb) -> Z(ai,bj)
C Irrep_jb = Irrep_bj

            Call Symtr1(Irrep_bj_z,Pop(1,Ispin),Vrt(1,Ispin),Ndim_ai_z,
     +                  Work(I000),Work(I020),Work(I030),Work(I040))

            Call Putlst(Work(I000),1,Ndim_bj_z,1,Irrep_bj_z,Tmp_list) 
        Enddo

        Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                       Isytyp(2,Tmp_list))
        Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_t_list),
     +                       Isytyp(2,Iz2_t_list))
        Ibgn = 1
        I000 = Ibgn
        I010 = I000 + Max(Isize_phph,Isize_pphh)
        I020 = I010 + Max(Isize_phph,Isize_pphh)
        I030 = I020 + 100
        Iend = I030 
        If (Iend .GT. Memleft) Call Insmem
     +     ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

        Call Getall(Work(I000),Isize_phph,Irrepx,Tmp_list)
        Call Dsst03i(Work(I000),Work(I010),Isize_Phph, 
     +               Isize_Phph,Work(I020),Spcase(Ispin),Irrepx)
        Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_t_list)

#ifdef _DEBUG_LVL0
      Call Getall(Work(I000),Isize_pphh,Irrepx,Iz2_t_list)
      call checksum("Tdcc_form_mutilde_0_d_uhf,Mu^t(AA,bb)(f):",
     +               Work(I000),Isize_pphh,S)
#endif
      Enddo 

C The ABAB block.

      Tmp_list   = 302 
      Iz2_t_list = 326
      L2_list    = 137

      Pha_type = 9
      Phb_type = 10
      Hpa_type = 16
      Hpb_type = 17
      Ppb_type = 20
      Hha_type = 21
      Hhb_type = 22
      Ppa_type = 19

      Ivo_0   = 1 +  Irpdpd(Irrepx,9)
      Ivv_0   = 1 +  Irpdpd(Irrepx,19)
      Ioo_0   = 1 +  Irpdpd(Irrepx,21)

C Target, Z(AI,bj)

      Do Irrep_bj_z1 = 1, Nirrep
         Irrep_ai_z1 = Dirprd(Irrep_bj_z1,Irrepx)

         Ndim_ai_z1  = Irpdpd(Irrep_ai_z1,Pha_type)
         Ndim_ia_z1  = Irpdpd(Irrep_ai_z1,Hpa_type)
         Ndim_bj_z1  = Irpdpd(Irrep_bj_z1,Phb_type) 
         Ndim_jb_z1  = Irpdpd(Irrep_bj_z1,Hpb_type)

         Ndim_aibj_z = Ndim_ai_z1 * Ndim_bj_z1

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Ndim_aibj_z
         Iend = I010 
         If (Iend .Gt. Memleft) Call Insmem
     +     ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)

         Call Dzero(Work(I000),Ndim_aibj_z)
   
         Iz2_t   = I000

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,7(1x,i2))")"Pha,Phb,Hpa,Hpb,Ppb,Hha,Hhb:",
     +Pha_type,Phb_type,Hpa_type,Hpb_type,Ppb_type,Hha_type,hhb_type
      write(6,"(a,a,7(1x,i4))")"Irrep_bj_z1,Irrep_ai_z1,Ndim_ai_z1,",
     +        "Ndim_ia_z1,Ndim_bj_z1,Ndim_jb_z1,Ndim_aibj_z:",
     +         Irrep_bj_z1,Irrep_ai_z1,Ndim_ai_z1,Ndim_ia_z1,
     +         Ndim_bj_z1,Ndim_jb_z1,Ndim_aibj_z
      write(6,"(a,4(1x,i2))")"Iz2_t,Ivo_0,Ivv_0,Ioo_0:",Iz2_t,Ivo_0,
     +                        Ivv_0,Ioo_0
#endif
C NDeP Eqn. 21, First term: Z(Ib,Aj) = mu(b,j) * L(I,A) and
C Z is built as Z(AI,bj) or (ai,BJ) and L and Dvo are
C stored as A,I and b,j.

        If (Irrep_bj_z1 .EQ. Irrepx) Then

C Z(AI,bj)=L(I,A)*Mu(b,j)

           Ibgn = I010 
           I020 = I010 + Irpdpd(1,9)
           Iend = I020
           If (Iend .Gt. Memleft) Call Insmem
     +        ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)

           Call Getlst(Work(I010),1,1,1,1,L1_list)

           Ivo_l   = I010
           Ivo_0   = 1 +  Irpdpd(Irrepx,9)

           Nrow = Ndim_ai_z1
           Ncol = Ndim_bj_z1
           Nsum = 1
           Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(Ivo_l),
     +                 Nrow,Dvo(Ivo_0),1,One,Work(Iz2_t),
     +                 Nrow)
        Endif 

        If (Irrep_ai_z1 .EQ. Irrepx) Then

C Z(AI,BJ)=L(j,b)*Mu(A,I)
        
           Ibgn = I010
           I020 = I010 + Irpdpd(1,10)
           Iend = I020
           If (Iend .Gt. Memleft) Call Insmem
     +        ("Tdcc_form_mutilde_0_d_uhf",Iend,Memleft)

           Call Getlst(Work(I010),1,1,1,2,L1_list)

           Ivo_l   = I010
           Ivo_0   = 1 

           Nrow = Ndim_ai_z1
           Ncol = Ndim_bj_z1
           Nsum = 1
           Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Dvo(Ivo_0),
     +                 Nrow,Work(Ivo_l),1,One,Work(Iz2_t),
     +                 Nrow)
        Endif

#ifdef _DEBUG_LVLM
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(1) :",
     +               Work(I000),Ndim_aibj_z,S)
#endif
C Read L2 as L2(AI,bm)
 
        Irrep_ai_l1 = Irrep_ai_z1
        Irrep_bm_l1 = Irrep_ai_z1
        Irrep_ai_l2 = Irrep_bj_z1
        Irrep_bm_l2 = Irrep_bj_z1

        Ndim_ai_l1 = Irpdpd(Irrep_ai_l1,Pha_type)
        Ndim_bm_l1 = Irpdpd(Irrep_bm_l1,Phb_type)
        Ndim_mb_l1 = Irpdpd(Irrep_bm_l1,Hpb_type)
        Ndim_ai_l2 = Irpdpd(Irrep_ai_l2,Pha_type)
        Ndim_bm_l2 = Irpdpd(Irrep_bm_l2,Phb_type)
        Ndim_mb_l2 = Irpdpd(Irrep_bm_l2,Hpb_type)

        Max_pq   = Max(Ndim_ai_z1,Ndim_bj_z1,Ndim_ai_l1,Ndim_bm_l1,
     +                 Ndim_ai_l2,Ndim_bm_l2,Ndim_mb_l2)
        Max_aibm = Max(Ndim_ai_l1*Ndim_bm_l1,Ndim_ai_l2*Ndim_bm_l2,
     +                 Ndim_ai_l1*Ndim_mb_l1)

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,a,7(1x,i4))")"Irrep_ai_l1,Irrep_bm_l1,Irrep_ai_l2,",
     +"Irrep_bm_l2:",Irrep_ai_l1,Irrep_bm_l1,Irrep_ai_l2,Irrep_bm_l2
      write(6,"(a,a,8(1x,i4))")"Ndim_ai_l1,Ndim_bm_l1,Ndim_mb_l1",
     +"Ndim_ai_l2,Ndim_bm_l2,Ndim_mb_l2,Max_pq,Max_aibm:",Ndim_ai_l1,
     + Ndim_bm_l1,Ndim_mb_l1,Ndim_ai_l2,Ndim_bm_l2,Ndim_mb_l2,
     + Max_pq,Max_aibm
#endif
        Ibgn = I010
        I020 = I010 + Max_aibm
        I030 = I020 + Max_pq
        I040 = I030 + Max_pq
        I050 = I040 + Max_pq
        Iend = I050

        If (Iend .GT. Memleft) Call Insmem
     +     ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

        Ndim_aibm_l1 = Ndim_ai_l1 * Ndim_bm_l1

        Call Getlst(Work(I010),1,Ndim_bm_l1,1,Irrep_bm_l1,L2_list)

        If (Irrepx .EQ. 1) Then

            Call Daxpy(Ndim_aibm_l1,Mu_0,Work(I010),1,Work(Iz2_t),1)

        Endif

#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(2) :",
     +               Work(I000),Ndim_aibj_z,S)
#endif

C NDep Eqn. 20, The fourth term; Construct Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 
C Z(AI,bj)=-L2(AI,bm)*Mu(j,m) 

        Do irrep_m = 1, Nirrep
           Irrep_j = Dirprd(Irrep_m,Irrepx)
           Irrep_b = Dirprd(Irrep_m,Irrep_bm_l1)
   
           Ndim_m = Pop(Irrep_m,2)
           Ndim_j = Pop(Irrep_j,2)
           Ndim_b = Vrt(Irrep_b,2)

           Ioff_oo = Ioo_0 + (Isymoff(Irrep_m,Irrepx,HHb_type)-1)
           Ioff_l2 = I010  + Ndim_ai_l1*(Isymoff(Irrep_m,Irrep_bm_l1,
     +                                   Phb_type)-1)
           Iz2_t   = I000  + Ndim_ai_z1*(Isymoff(Irrep_j,Irrep_bj_z1,
     +                                  Phb_type)-1)

           Nrow = Ndim_ai_l1 * Ndim_b
           Ncol = Ndim_j
           Nsum = Ndim_m

           Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                 Nrow,Doo(Ioff_oo),Ncol,One,Work(Iz2_t),
     +                 Nrow)
        Enddo

#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(31):",
     +               Work(I000),Ndim_aibj_z,S)
#endif

C Z(AI,bj) -> Z(AI,jb) 

          Call Symtr1(Irrep_bj_z1,Vrt(1,2),Pop(1,2),Ndim_ai_z1,
     +                Work(I000),Work(I020),Work(I030),
     +                Work(I040))

C L2(AI,ej)-> L2(AI,je)
C Irrep_ej_l1=Irrep_bm_l1

          Irrep_ej_l1 = Irrep_bm_l1
          Call Symtr1(Irrep_ej_l1,Vrt(1,2),Pop(1,2),Ndim_ai_l1,
     +                Work(I010),Work(I020),Work(I030),
     +                Work(I040))

C NDep Eqn. 20, third term Z(AI,jb)=L2(AI,je)*Mu(b,e) 

          Do Irrep_e = 1, Nirrep
             Irrep_b = Dirprd(Irrep_e,Irrepx)
             Irrep_j = Dirprd(Irrep_e,Irrep_ej_l1)
          
             Ndim_e = Vrt(Irrep_e,2)
             Ndim_b = Vrt(Irrep_b,2)
             Ndim_j = Pop(Irrep_j,2)
  
             Ioff_vv = Ivv_0 + (Isymoff(Irrep_e,Irrepx,Ppb_type)-1)
             Ioff_l2 = I010  + Ndim_ai_l1*(Isymoff(Irrep_e,Irrep_ej_l1,
     +                                         Hpb_type)-1)
             Iz2_t   = I000  + Ndim_ai_z1*(Isymoff(Irrep_b,Irrep_bj_z1,
     +                                         Hpb_type)-1)
             Nrow = Ndim_ai_l1 * Ndim_j
             Ncol = Ndim_b
             Nsum = Ndim_e

             Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                   Nrow,Dvv(Ioff_vv),Ncol,One,Work(Iz2_t),
     +                   Nrow)
          Enddo 
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(41):",
     +               Work(I000),Ndim_aibj_z,S)
#endif
          Ibgn = I020
          I030 = I020 + Ndim_ai_z1 * Ndim_jb_z1
          Iend = I030
          If (Iend .GT. Memleft) Call Insmem
     +       ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

C Z(AI,jb) -> Z(jb,AI)
        
          Call Transp(Work(I000),Work(I020),Ndim_jb_z1,Ndim_ai_z1)
          Call Dcopy(Ndim_ai_z1*Ndim_jb_z1,Work(I020),1,Work(I000),1)

C NDep Eqn. 20, The fourth term; Construct Z(bj,AI)=-L2(bj,AM)*Mu(I,M)
C Read L2 as L2(AM,bj)
C Irrep_am_l2 = Irrep_ai_l2
C Irrep_bj_l2 = Irrep_bm_l2

          Irrep_am_l2 = Irrep_ai_l2
          Irrep_bj_l2 = Irrep_bm_l2 
          Ndim_bj_l2  = Ndim_bm_l2 
          Ndim_am_l2  = Ndim_ai_l2 
          Ndim_jb_l2  = Ndim_mb_l2

          Call Getlst(Work(I010),1,Ndim_bj_l2,1,Irrep_bj_l2,L2_list)

C L(AM,bj) ->  L(AM,jb)

          Ibgn = I020
          I030 = I020 + Ndim_am_l2 * Ndim_jb_l2
          I040 = I030 + Max_pq
          I050 = I040 + Max_pq
          I060 = I050 + Max_pq
          Iend = I060
          If (Iend .GT. Memleft) Call Insmem
     +       ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)
   
          Call Symtr1(Irrep_bj_l2,Vrt(1,2),Pop(1,2),Ndim_am_l2,
     +                Work(I010),Work(I030),Work(I040),
     +                Work(I050))

C L(AM,jb) -> L(jb,AM)
        
          Call Transp(Work(I010),Work(I020),Ndim_jb_l2,Ndim_am_l2)
          Call Dcopy(Ndim_am_l2*Ndim_jb_l2,Work(I020),1,Work(I010),
     +               1)

C Construct Z(jb,AI)=-L2(jb,AM)*Mu(I,M)

           Do irrep_m = 1, Nirrep
              Irrep_i = Dirprd(Irrep_m,Irrepx)
              Irrep_a = Dirprd(Irrep_m,Irrep_am_l2)

              Ndim_m = Pop(Irrep_m,1)
              Ndim_i = Pop(Irrep_i,1)
              Ndim_a = Vrt(Irrep_a,1)

              Ioff_oo = 1     + (Isymoff(Irrep_m,Irrepx,Hha_type)-1)
              Ioff_l2 = I010  + Ndim_jb_l2*(Isymoff(Irrep_m,
     +                                      Irrep_am_l2,Pha_type)-1)
              Iz2_t   = I000  + Ndim_jb_z1*(Isymoff(Irrep_i,
     +                                      Irrep_ai_z1,Pha_type)-1)
              Nrow = Ndim_jb_l2 * Ndim_a
              Ncol = Ndim_i
              Nsum = Ndim_m

              Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ioff_l2),
     +                    Nrow,Doo(Ioff_oo),Ncol,One,Work(Iz2_t),
     +                    Nrow)
           Enddo
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(32):",
     +               Work(I000),ndim_aibj_z,S)
#endif
C Z(jb,AI) -> Z(jb,IA)
        
           Call Symtr1(Irrep_ai_z1,Vrt(1,1),Pop(1,1),Ndim_jb_z1,
     +                 Work(I000),Work(I030),Work(I040),
     +                 Work(I050))

C L(jb,EI) -> L(jb,IE)

           Irrep_ei_l2 = Irrep_am_l2

           Call Symtr1(Irrep_ei_l2,Vrt(1,1),Pop(1,1),Ndim_jb_l2,
     +                 Work(I010),Work(I030),Work(I040),
     +                 Work(I050))

C NDep Eqn. 20, third term Z(jb,IA)=L2(jb,IE)*Mu(A,E)

           Do Irrep_e = 1, Nirrep
              Irrep_a = Dirprd(Irrep_e,Irrepx)
              Irrep_i = Dirprd(Irrep_e,Irrep_ei_l2)

              Ndim_e = Vrt(Irrep_e,1)
              Ndim_a = Vrt(Irrep_a,1)
              Ndim_i = Pop(Irrep_i,1)

              Ioff_vv = 1     + (Isymoff(Irrep_e,Irrepx,Ppa_type)-1)
              Ioff_l2 = I010  + Ndim_jb_l2*(Isymoff(Irrep_e,
     +                                      Irrep_ei_l2,Hpa_type)-1)
              Iz2_t   = I000  + Ndim_jb_z1*(Isymoff(Irrep_a,
     +                                      Irrep_ai_z1,Hpa_type)-1)

              Nrow = Ndim_jb_l2 * Ndim_i
              Ncol = Ndim_a
              Nsum = Ndim_e

              Call Xgemm("N","T",Nrow,Ncol,Nsum,One,Work(Ioff_l2),
     +                    Nrow,Dvv(Ioff_vv),Ncol,One,Work(Iz2_t),
     +                    Nrow)
           Enddo
#ifdef _DEBUG_LVLM
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(42):",
     +               Work(I000),Ndim_aibj_z,S)
#endif
           Ibgn = I010
           I020 = I010 + Ndim_ai_z1 * Ndim_jb_z1
           I030 = I020 + Max_pq
           I040 = I030 + Max_pq
           I050 = I040 + Max_pq
           Iend = I050
           If (Iend .GT. Memleft) Call Insmem
     +        ("@-Tdvee_form_mutilde_0_d_uhf",Iend,Memleft)

C Z(jb,IA) -> Z(IA,jb) -> Z(IA,bj) -> Z(AI,bj)

           Call Transp(Work(I000),Work(I010),Ndim_ia_z1,Ndim_jb_z1)
           Call Symtr1(Irrep_bj_z1,Pop(1,2),Vrt(1,2),Ndim_ia_z1,
     +                 Work(I010),Work(I020),Work(I030),Work(I040))

           Call Symtr3(Irrep_ai_z1,Pop(1,1),Vrt(1,1),Ndim_ai_z1,
     +                 Ndim_bj_z1,Work(I010),Work(I020),Work(I030),
     +                 Work(I040))
           Call Putlst(Work(I010),1,Ndim_bj_z1,1,Irrep_bj_z1,Tmp_list)

      Enddo 

      Isize_phph = Idsymsz(Irrepx,Isytyp(1,Tmp_list),
     +                            Isytyp(2,Tmp_list))

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
     +            Vrt(1,2),Pop(1,2),Work(I020),Irrepx,'1324')

      Isize_pphh = Idsymsz(Irrepx,Isytyp(1,Iz2_t_list),
     +                            Isytyp(2,Iz2_t_list))
#ifdef _DEBUG_LVL0
      Write(6,*)
      call checksum("Tdcc_form_mutilde_0_d_uhf,mu^t(Ab)(f)   :",
     +               Work(I010),Isize_pphh,S)
#endif
      Call Putall(Work(I010),Isize_pphh,Irrepx,Iz2_t_list)

      Return
      End
        
