#include "flags.h"
      Subroutine Tdcc_form_mutilde_0_s(Work,Memleft,Doo,Dvv,Dvo,
     +                                 Mu_0_t_S,Irrepx,Lenoo,Lenvv,
     +                                 Lenvo,Iuhf,Mu_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)

      Double Precision Mu_0, Mu_0_t_s(Lenvo)
      Double Precision Mone,Zero,One,Two
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

C Make mu_tilde(i,a) (stored as mu_tilde(a,i))

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_0_s---"
        Write(6,*)
#endif 
      L1_list      = 190 
      Mut_ph_list  = 392
      Ivo_0         = 1

      Do Ispin = 1, 1+Iuhf 
         
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Irpdpd(1,8+Ispin)
         I020 = I010 + Irpdpd(1,18+Ispin) * Irpdpd(1,20+Ispin)
         Iend = I020

         If (Iend .Gt. Memleft) Call Insmem("-Tdvee_form_mu_0_S",
     +                                       Iend,Memleft)
         Call Getlst(Work(I000),1,1,1,Ispin,L1_list)

         Ndim_vo = Irpdpd(Irrepx,8+Ispin)

         Ivo_l0 = I000
         Ivo_z0 = 1 + (Ispin-1) * Irpdpd(Irrepx,9)
         Ivv_d0 = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioo_d0 = 1 + (Ispin-1) * Irpdpd(Irrepx,21)

#ifdef _DEBUG_LVL0
      call Checksum("Tdcc_form_mutilde_0_s,Dvo      :",
     +               Dvo(Ivo_z0),Ndim_vo,S)
#endif
         Call Dzero(Mu_0_t_s(Ivo_z0),Irpdpd(Irrepx,8+Ispin))

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)

            Ndim_i = Pop(Irrep_i,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
               
            Ivo_z = Ivo_z0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)
            Ivo_d = Ivo_z

C NDeP Eqn. 20, First term: z(a,i) = Mu(a,i) (and Z(A,I)=Mu(A,I))

            Ndim_ai = Ndim_i * Ndim_a
            Call Dcopy(Ndim_ai,Dvo(Ivo_d),1,Mu_0_t_s(Ivo_z),1)

C NDep Eqn. 20 Second term: z(A,I) = Mu(E,A) * L(I,E) (and 
C z(a,i) = Mu(e,a) * L(i,e)). The L is stored as L(E,I).
            
            Irrep_e = Irrep_i

            Ndim_e = Vrt(Irrep_e,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)
               
            Nrow = Ndim_a
            Nsum = Ndim_e
            Ncol = Ndim_i

            Ivv_d = Ivv_d0 + (Isymoff(Irrep_a,Irrepx,18+Ispin)-1)
            Ivo_l = Ivo_l0 + (Isymoff(Irrep_i,1,8+Ispin)-1)

            Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Dvv(Ivv_d),Nrow,
     +                  Work(Ivo_l),Nsum,One,Mu_0_t_s(Ivo_z),Nrow)
         Enddo

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_0_s,Mu^t(1-2):",
     +               Mu_0_t_s(Ivo_z0),Ndim_vo,S)
#endif 
#ifdef _DEBUG_LVLM
      do ispin1=1, 1+iuhf
      ioff = 1 + (Ispin-1)*Irpdpd(Irrepx,21)
      do irrep_i=1,Nirrep
         irrep_a=dirprd(irrep_i,Irrepx)
         ndim_i = Pop(irrep_i,ispin1)
         Ndim_a = Vrt(irrep_a,ispin1)
         Write(6,*)
         write(6,"(a,a,5(1x,i2))") "irrep_i,Irrep_a,Ndim_i,",
     +                     "Ndim_a,Ioff: ",irrep_i,Irrep_a,
     +                      Ndim_i,Ndim_a,Ioff
         call output(Mu_0_t_S(ioff),1,ndim_a,1,Ndim_i,Ndim_a,
     +                Ndim_i,1)
         Ioff = Ioff + ndim_i*Ndim_a
      enddo
      enddo
#endif 

C NDeP Eqn. 20, Two third terms: 
C a. z(a,i) = -L(m,a)*Mu(i,m)  (and z(A,I) = -L(M,A)*Mu(I,M))

         Ivo_l0 = I000
         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)
            Irrep_m = Irrep_a

            Ndim_m = Pop(Irrep_m,Ispin)
            Ndim_a = Vrt(Irrep_a,Ispin)
            Ndim_i = Pop(Irrep_i,Ispin)
            
            Ioo_d = Ioo_d0 + (Isymoff(Irrep_m,Irrepx,20+Ispin)-1)
            Ivo_l = Ivo_l0 + (Isymoff(Irrep_m,1,8+Ispin)-1)
            Ivo_z = Ivo_z0 + (Isymoff(Irrep_i,Irrepx,8+Ispin)-1)

            Nrow = Ndim_a
            Nsum = Ndim_m
            Ncol = Ndim_i

#ifdef _DEBUG_LVLM
      write(6,*)
      Write(6,"(a,1x,I4)") "Work(Ivo_l),Ivo_l: ",Ivo_l
      call output(Work(Ivo_l),1,Ndim_a,1,Ndim_i,Ndim_a,Ndim_i,1)
      Write(6,"(a,1x,I4)") "Doo(Ioo_d),Ioo_d: ",Ioo_d
      call output(Doo(Ioo_d),1,Ndim_i,1,Ndim_m,Ndim_i,Ndim_m,1)
      Write(6,"(a,1x,I4)") "Mu_0_t_s(Ivo_z),Ivo_z: ",Ivo_z
      call output(Mu_0_t_s(Ivo_z),1,Ndim_a,1,Ndim_i,Ndim_a,
     +            Ndim_m,1)
#endif     
            Call Xgemm("N","T",Nrow,Ncol,Nsum,Mone,Work(Ivo_l),
     +                 Nrow,Doo(Ioo_d),Ncol,One,Mu_0_t_s(Ivo_z),
     +                 Nrow)
C           call checksum("check1",Mu_0_t_s(Ivo_z),Nrow*Ncol,s)
         Enddo 

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_0_s,Mu^t(3a) :",
     +               Mu_0_t_s(Ivo_0),Ndim_vo,S)
#endif
#ifdef _DEBUG_LVLM
      do ispin1=1, 1+iuhf
      ioff = 1 + (Ispin1-1)*Irpdpd(Irrepx,21)
      do irrep_i=1,Nirrep
         irrep_a=dirprd(irrep_i,Irrepx)
         ndim_i = Pop(irrep_i,ispin1)
         Ndim_a = Vrt(irrep_a,ispin1)
         Write(6,*)
         write(6,"(a,a,5(1x,i2))") "irrep_i,Irrep_a,Ndim_i,",
     +                     "Ndim_a,Ioff: ",irrep_i,Irrep_a,
     +                      Ndim_i,Ndim_a,Ioff
         call output(Mu_0_t_S(ioff),1,ndim_a,1,Ndim_i,Ndim_a,
     +                Ndim_i,1)
         Ioff = Ioff + ndim_i*Ndim_a
      enddo
      enddo
#endif

C z(a,i) = L(m,a)*Mu(i,m)*delta(i,m) (and z(A,I) = L(M,A) * 
C                                         Mu(I,M)Delta(I,M)
         If (Irrepx .EQ. 1) Then
             Call Daxpy(Lenvo,Mu_0,Work(Ivo_l0),1,Mu_0_t_s(Ivo_z0),1)
         Endif 

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_0_s,Mu^t(3b) :",
     +               Mu_0_t_s(Ivo_0),Ndim_vo,S)
#endif 
      Enddo

C NDep Eqn. 20, The third term Z(A,I) = L(MI,EA)*mu(E,M) 
C UHF  Z(A,I) = L(MI,EA) * mu(E,M) and Z(a,i) = L(mi,ea) * mu(e,m)
c Lists 137,136: L(AI,em) and L(ai,EM)
C RHF  Z(A,I) = {2*L(AI,em) -L(Am,eI)} * mu(e,m)

      Do Ispin = 1, 1+Iuhf 

C UHF Ispin=1, L(AI,em)*Mu(e,m); Ispin=2, L(ai,EM)*Mu(E,M)
C     Ispin=2, L(AI,EM)*Mu(E,M);          L(ai,em)*Mu(e,m)
C
C RHF Ispin=1  {2*L(AI,em)-L(Am,eI)}Mu(E,M)

         Ivo_z  = 1 + (Ispin-1) * Irpdpd(Irrepx,8+Ispin)

         L2_list_ab0 = 138 - Ispin
         L2_list_ab1 = 139
         L2_list_pq0 = 133 + Ispin

C Irrep_em=Irrep_ai=Irrepx

         Irrep_em = Irrepx
         Irrep_ai = Irrepx
         Irrep_ei = Irrepx

         Len_em_pq = Irpdpd(Irrepx,11-Ispin)
         Len_ai_pq = Irpdpd(Irrepx,8+Ispin)
         Len_em_pp = Irpdpd(Irrepx,8+Ispin)
         Len_ai_pp = Irpdpd(Irrepx,8+Ispin)
         Len_am_ab = Irpdpd(Irrepx,11)
         Len_ei_ab = Irpdpd(Irrepx,12)

         Max_pq = Max(Len_em_pq,Len_ai_pq,Len_ai_pp,Len_em_pp,
     +                Len_am_ab,Len_ei_ab)
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max_pq*Max_pq 
         I020 = I010 + Max_pq*Max_pq
         Iend = I020 

         If (Iend .GT. Memleft) Call Insmem("-Tdvee_form_mutilde_0_s",
     +                                       Iend,Memleft)

C UHF Ispin=1, Z(A,I) = L(AI,em)*Mu(e,m)
C     Ispin=2  Z(a,i) = L(ai,EM)*Mu(E,M)
C RHF Ispin=1  Z(A,I) = {2*L(AI,em) -L(Am,eI)} Mu(e,m)

         Call Getlst(Work(I000),1,Len_em_pq,1,Irrep_em,L2_list_ab0)

         If (Iuhf .EQ. 0) Then
             Len_aiem_ab = Len_am_ab * Len_ei_ab
             Call Getlst(Work(I010),1,Len_ei_ab,1,Irrep_ei,
     +                    L2_list_ab1)
             Call Dscal(Len_aiem_ab,Two,Work(I000),1)
             Call Daxpy(Len_aiem_ab,Mone,Work(I010),1,Work(I000),1)
         Endif

         If (Ispin .EQ. 1) Then
            Ivo_d0 = 1
            Ivo_d1 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_z0 = 1
            Ivo_z1 = 1 
            Nrow0  = Len_ai_pp
            Ncol0  = 1
            Nsum0  = Len_em_pp
            If (Iuhf .Eq. 0) Then 
                Nrow1 = Len_am_ab
                Ncol1 = 1
                Nsum1 = Len_ei_ab 
            Else
                Nrow1 = Len_ai_pq
                Ncol1 = 1
                Nsum1 = Len_em_pq
            Endif 
              
         Elseif (Ispin .EQ. 2) Then
            Ivo_d0 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_d1 = 1 
            Ivo_z0 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Ivo_z1 = 1 + Iuhf * Irpdpd(Irrepx,9)
            Nrow0  = Len_ai_pp
            Ncol0  = 1
            Nsum0  = Len_em_pp
            Nrow1  = Len_ai_pq
            Ncol1  = 1
            Nsum1  = Len_em_pq
         Endif 

#ifdef _DEBUG_LVLM
      write(6,*)
      Write(6,"(a,3(1x,I4))") "Nrow1,Ncol1,Nsum1: ", 
     +                         Nrow1,Ncol1,Nsum1
      Write(6,"(a,1x,I4)") "Work(I000),Ivo_l: ", I000
      call checksum("L2_abab",Work(I000),Nrow1*Nsum1,s)
      Write(6,"(a,1x,I4)") "Dvo(Ivo_d1),Ivo_d1: ",Ivo_d1
      call checksum("Dvo",Dvo(Ivo_d1),Nsum1*Ncol1,s)
      Write(6,"(a,1x,I4)") "Mu_0_t_s(Ivo_z1),Ivo_z1: ",Ivo_z1
      call checksum("Mu_0_t_s",Mu_0_t_s(Ivo_z1),Nrow1*Ncol1,s)
#endif
            Call Xgemm("N","N",Nrow1,Ncol1,Nsum1,One,Work(I000),Nrow1,
     +                 Dvo(Ivo_d1),Nrow1,One,Mu_0_t_s(Ivo_z1),Nrow1)

            If (Iuhf .NE. 0) Then

C UHF Ispin=1, L(AI,EM)*Mu(E,M)
C     Ispin=2  L(ai,em)*Mu(e,m)

                Call Getlst(Work(I000),1,Len_em_pp,1,Irrep_em,
     +                      L2_list_pq0)
                Call Xgemm("N","N",Nrow0,Ncol0,Nsum0,Mone,Work(I000),
     +                     Nrow0,Dvo(Ivo_d0),Nrow0,One,
     +                     Mu_0_t_s(Ivo_z0),Nrow0)
            Endif 

         Call Putlst(Mu_0_t_s(Ivo_z0),1,1,1,Ispin,Mut_ph_list)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Ndim_vo = Irpdpd(Irrepx,8+Ispin)
      call Checksum("Tdcc_form_mutilde_0_s,Mu^t(f)  :",
     +               Mu_0_t_s(Ivo_z0),Ndim_vo,S)
#endif 

      Enddo 
  
      Return
      End
