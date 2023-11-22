      Subroutine Tdcc_form_mutilde_dot_da(Work,Memleft,Irrepx,Iuhf,
     +                                    Ioff_vo)
   
      Implicit Double Precision (A-H, O-Z)
      Integer Hbar_hhhp_list,z2_pphh_list
      Logical Rhf 
C
      Dimension Work(Memleft),Ioff_vo(8,2),Ioff_h(8),Ioff_z(8)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

      Data One,Onem,Zero/1.0D0,-1.0D0,0.0D0/

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_dot_da---"
        Write(6,*)
#endif
      Rhf = .False.
      If (Iuhf .Eq. 0) Rhf = .True.

C All the expression has to be transposed 
C Z(Ab,Ij) = P(Ab) {Mu(A,M) * Hbar(Ij,Mb)}

      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)
         Irrep_mb = Irrep_ij
         
         Hbar_hhhp_list = 10
         z2_pphh_list   = 346
     
         Ndim_ij_h = Irpdpd(Irrep_ij,Isytyp(1,Hbar_hhhp_list))
         Ndim_mb_h = Irpdpd(Irrep_mb,Isytyp(2,Hbar_hhhp_list))
         Ndim_ab_z = Irpdpd(Irrep_ab,Isytyp(1,z2_pphh_list))
         Ndim_ij_z = Irpdpd(Irrep_ij,Isytyp(2,z2_pphh_list))

         Max_ijmb_h = Max(Ndim_ij_h,Ndim_mb_h)
         Max_ijab_z = Max(Ndim_ij_z,Ndim_ab_z)
         Max_both   = Max(Max_ijmb_h,Max_ijab_z)

         Ndim_ijmb_h = Ndim_ij_h * Ndim_mb_h
         Ndim_ijab_z = Ndim_ij_z * Ndim_ab_z

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max(3*Max_both,Ndim_ijab_z)
         I020 = I010 + Max(Ndim_ijmb_h,Ndim_ijab_z,3*Max_both)
         Iend = I020

         Call Dzero(Work(I000),Ndim_ijab_z)

         If (Iend .LT. Memleft) Then
            Call Getlst(Work(I010),1,Ndim_mb_h,1,Irrep_mb,
     +                  Hbar_hhhp_list)

C List 110 = Hbar(Ij,Mb); A permutation of Ket indices are needed. 
C Hbar(Ij,Mb) -> Hbar(Ij,bM)

            Ibgn  = I000
            IPtr1 = Ibgn
            Iptr2 = Iptr1 + Max_ijmb_h
            Iptr3 = Iptr2 + Max_ijmb_h
            Iptr4 = Iptr3 + Max_ijmb_h
            Iend  = Iptr4

            Call Symtr1(Irrep_mb,Pop(1,1),Vrt(1,2),Ndim_ij_h,
     +                  Work(I010),Work(Iptr1),Work(Iptr2),
     +                  Work(Iptr3))

C Z(Ij,bA) = SUM W(Ij,bM) * T(A,M)^t
            
            Loc_h = 0
            Loc_z = 0
            Do Irrep_hz = 1, Nirrep
               Irrep_bh = Dirprd(Irrep_hz,Irrep_mb)
               Irrep_bz = Dirprd(Irrep_hz,Irrep_ab)
               
               Ioff_h(Irrep_hz) = I010 + Loc_h
               Ioff_z(Irrep_hz) = I000 + Loc_z
          
               Loc_h = Loc_h+Ndim_ij_h*Vrt(Irrep_bh,2)*Pop(Irrep_hz,1)
               Loc_z = Loc_z+Ndim_ij_z*Vrt(Irrep_bz,2)*Vrt(Irrep_hz,1)
            Enddo 

C Z(Ij,bA) = SUM W(Ijb,M) * T(A,M)^t

            Do Irrep_m = 1, Nirrep
               Irrep_a = Dirprd(Irrep_m,Irrepx)
               Irrep_b = Dirprd(Irrep_m,Irrep_mb)
               
                Nrow = Ndim_ij_h*Vrt(Irrep_b,2)
                Ncol = Vrt(Irrep_a,1) 
                Nsum = Pop(Irrep_m,1)

                Loc_h = Ioff_h(Irrep_m)
                Loc_z = Ioff_z(Irrep_a)
                Loc_m = Ioff_vo(Irrep_m,1)
 
                Call Xgemm("N","T",Nrow,Ncol,Nsum,Onem,Work(loc_h),
     +                     Nrow,Work(Loc_m),Ncol,Zero,Work(loc_z),
     +                     Nrow)
            Enddo

         Else 

           Write(6,"(a,a)") "Insufficient amount of memory to perform",
     +                      " the contraction in memory."
           Call Errex

         Endif 
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_dot_da Z(Ab,Ij)(1):",
     +               Work(I000),Ndim_ijab_z,S)
#endif

C Z(Ij,bA) --> Z(Ij,Ab)

       Ibgn  = I010 
       Iptr1 = Ibgn
       Iptr2 = Iptr1 + Max_ijab_z
       Iptr3 = Iptr2 + Max_ijab_z
       Iptr4 = Iptr3 + Max_ijab_z
       Iend  = Iptr4

       Call Symtr1(Irrep_ab,Vrt(1,2),Vrt(1,1),Ndim_ij_z,Work(I000),
     +             Work(Iptr1),Work(Iptr2),Work(Iptr3))

       If (Rhf) Then

C Z(Ij,Ab) -> Z(Ab,Ij)

           Call Transp(Work(I000),Work(I010),Ndim_ab_z,Ndim_ij_z)
           Call Dcopy(Ndim_ab_z*Ndim_ij_z,Work(I010),1,Work(I000),1)
           Call Symrhf3(Irrep_ab,Irrep_ij,Vrt(1,1),Pop(1,1),
     +                  Ndim_ab_z,Work(I000),Work(Iptr1),Work(Iptr2),  
     +                  Work(Iptr3))
#ifdef _DEBUG_LVL0

#else
           Call Getlst(Work(I010),1,Ndim_ij_z,1,Irrep_ij,z2_pphh_list)
           Call Daxpy(Ndim_ab_z*Ndim_ij_z,One,Work(I010),1,
     +                Work(I000),1)
#endif 
C The target, Z(Ab,Ij) 

           Call Putlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,
     +                 z2_pphh_list)
       Else

C Z(Ij,bA) = P(bA) {Mu(b,m) * Hbar(Ij,Am)}
C List 109 = Hbar(Ij,Am) no permutation of ket indices needed. 

           Hbar_hhhp_list = 107+2
           Irrep_am = Irrep_ij

           Ndim_ij_h = Irpdpd(Irrep_ij,Isytyp(1,Hbar_hhhp_list))
           Ndim_am_h = Irpdpd(Irrep_am,Isytyp(2,Hbar_hhhp_list))

           Ndim_ijam_h = Ndim_ij_h * Ndim_am_h
            
           Ibgn = I010 
           I020 = I010 + Max(Ndim_ijam_h,Ndim_ijab_z,3*Max_ijab_z)
           Iend = I020

           If (Iend .LT. Memleft) Then
              Call Getlst(Work(I010),1,Ndim_am_h,1,Irrep_am,
     +                    Hbar_hhhp_list)

C Z(Ij,Ab) = Hbar(Ij,Am) * Mu(b,m)^t

            Loc_h = 0
            Loc_z = 0
            Do Irrep_hz = 1, Nirrep
               Irrep_ah = Dirprd(Irrep_hz,Irrep_am)
               Irrep_az = Dirprd(Irrep_hz,Irrep_ab)

               Ioff_h(Irrep_hz) = I010 + Loc_h
               Ioff_z(Irrep_hz) = I000 + Loc_z

               Loc_h = Loc_h+Ndim_ij_h*Vrt(Irrep_ah,1)*Pop(Irrep_hz,2)
               Loc_z = Loc_z+Ndim_ij_z*Vrt(Irrep_az,1)*Vrt(Irrep_hz,2)
            Enddo

            Do Irrep_m = 1, Nirrep
               Irrep_b = Dirprd(Irrep_m,Irrepx)
               Irrep_a = Dirprd(Irrep_m,Irrep_am)

                Nrow = Ndim_ij_z*Vrt(Irrep_a,1)
                Ncol = Vrt(Irrep_b,1)
                Nsum = Pop(Irrep_m,2)

                Loc_h = Ioff_h(Irrep_m)
                Loc_z = Ioff_z(Irrep_b)
                Loc_m = Ioff_vo(Irrep_m,1)

                Call Xgemm("N","T",Nrow,Ncol,Nsum,Onem,Work(loc_h),
     +                     Nrow,Work(Loc_m),Ncol,One,Work(loc_z),
     +                     Nrow)
            Enddo

C Z(Ij,Ab) -> Z(Ab,Ij)

            Call Transp(Work(I000),Work(I010),Ndim_ab_z,Ndim_ij_z)
            Call Dcopy(Ndim_ab_z*Ndim_ij_z,Work(I010),1,Work(I000),1)

            Call Putlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,z2_pphh_list)

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mu_dot_da Z(Ab,Ij)(2)     :",
     +               Work(I000),Ndim_ijab_z,S)
#endif
         Else

           Write(6,"(a,a)") "Insufficient amount of memory to perform",
     +                      " the contraction in memory."
           Call Errex

         Endif

       Endif

      Enddo 
#ifdef _DEBUG_LVL0
      I000 = 1
      Isize_abij_z=Idsymsz(Irrepx,Isytyp(1,z2_pphh_list),
     +           Isytyp(2,z2_pphh_list))
      call getall(Work(I000),Isize_abij_z,Irrepx,z2_pphh_list)
      call checksum("Tdcc_form_mu_dot_da Z(Ab,Ij)(f)  :",
     +               Work(i000),Isize_abij_z,S)
#endif 
      If (Rhf) Return

      Do Ispin = 1, 1+Iuhf 

C Ispin=1: Z(AB,IJ) = P(AB) {Mu(B,M) * Hbar(IJ,MA)}
C Ispin=2: Z(ab,ij) = P(ab) {Mu(b,m) * Hbar(ij,ma)}

         Do Irrep_ij = 1, Nirrep
            Irrep_ab = Dirprd(Irrep_ij,Irrepx)
            Irrep_ma = Irrep_ij

            Hbar_hhhp_list = 6   + Ispin
            z2_pphh_list   = 343 + Ispin

            Ndim_ij_c_h = Irpdpd(Irrep_ij,Isytyp(1,Hbar_hhhp_list))
            Ndim_ma_f_h = Irpdpd(Irrep_ma,Isytyp(2,Hbar_hhhp_list))
            Ndim_ab_c_z = Irpdpd(Irrep_ab,Isytyp(1,z2_pphh_list))
            Ndim_ij_c_z = Irpdpd(Irrep_ij,Isytyp(2,z2_pphh_list))
            Ndim_ab_f_z = Irpdpd(Irrep_ab,18+Ispin)

            Max_ijma_h = Max(Ndim_ij_c_h,Ndim_ma_f_h)
            Max_ijab_z = Max(Ndim_ij_c_z,Ndim_ab_c_z,Ndim_ab_f_z)

            Ndim_ijma_cf_h = Ndim_ij_c_h * Ndim_ma_f_h
            Ndim_ijab_cf_z = Ndim_ij_c_z * Ndim_ab_f_z
            Ndim_ijab_cc_z = Ndim_ij_c_z * Ndim_ab_c_z
    
            Ibgn = 1
            I000 = Ibgn 
            I010 = I000 + Max(Ndim_ijab_cf_z,Ndim_ijab_cc_z,
     +                        3*Max_ijma_h)
            I020 = I010 + Max(Ndim_ijma_cf_h,Ndim_ijab_cf_z)
            Iend = I020

            Call Dzero(Work(I000),Ndim_ijab_cf_z)
               
            If (Iend .LT. Memleft) Then

C Hbar(i<j,ma) or Hbar(I<J,MA) 

               Call Getlst(Work(I010),1,Ndim_ma_f_h,1,Irrep_ma,
     +                     Hbar_hhhp_list)

C Hbar(i<j,am) or Hbar(I<J,MA)

             Ibgn  = I000
             IPtr1 = Ibgn
             Iptr2 = Iptr1 + Max_ijma_h
             Iptr3 = Iptr2 + Max_ijma_h
             Iptr4 = Iptr3 + Max_ijma_h
             Iend  = Iptr4

             Call Symtr1(Irrep_ma,Pop(1,Ispin),Vrt(1,Ispin),
     +                   Ndim_ij_c_h, Work(I010),Work(Iptr1),
     +                   Work(Iptr2),Work(Iptr3))

C Z(I<J,AB) = Hbar(I<J,AM) * Mu(B,M)^t or Z(i<j,ab)=Hbar(i<j,am) * Mu(b,m)^t

             Loc_h = 0
             Loc_z = 0
             Do Irrep_hz = 1, Nirrep
                Irrep_ah = Dirprd(Irrep_hz,Irrep_ma)
                Irrep_az = Dirprd(Irrep_hz,Irrep_ab)
 
                Ioff_h(Irrep_hz) = I010 + Loc_h
                Ioff_z(Irrep_hz) = I000 + Loc_z
 
                Loc_h = Loc_h+Ndim_ij_c_h*Vrt(Irrep_ah,Ispin)*
     +                                Pop(Irrep_hz,Ispin)
                Loc_z = Loc_z+Ndim_ij_c_z*Vrt(Irrep_az,Ispin)*
     +                                Vrt(Irrep_hz,Ispin)
             Enddo

             Do Irrep_m = 1, Nirrep
                 Irrep_b = Dirprd(Irrep_m,Irrepx)
                 Irrep_a = Dirprd(Irrep_m,Irrep_ma)
 
                 Nrow = Ndim_ij_c_h*Vrt(Irrep_a,Ispin)
                 Ncol = Vrt(Irrep_b,Ispin)
                 Nsum = Pop(Irrep_m,Ispin)
 
                 Loc_h = Ioff_h(Irrep_m)
                 Loc_z = Ioff_z(Irrep_b)
                 Loc_m = Ioff_vo(Irrep_m,Ispin)

                 Call Xgemm("N","T",Nrow,Ncol,Nsum,Onem,Work(loc_h),
     +                      Nrow,Work(Loc_m),Ncol,Zero,Work(loc_z),
     +                      Nrow)
                 Call Dscal(Nrow*Ncol,Onem,Work(loc_z),1)
             Enddo

C Z(I<J,AB) -> Z(I<j,A<B)

             Call Assym2(Irrep_ab,Vrt(1,Ispin),Ndim_ij_c_z,Work(I000))

C Z(I<J,A<B) -> Z(A<B,I<J)

             Call Transp(Work(I000),Work(I010),Ndim_ab_c_z,
     +                   Ndim_ij_c_z)
             Call Dcopy(Ndim_ij_c_z*Ndim_ab_c_z,Work(I010),1,
     +                  Work(I000),1)
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_dot_da Z(A<B,I<J)(3):",
     +               Work(I000),Ndim_ijab_cc_z,S)
#endif
             Call Putlst(Work(I000),1,Ndim_ij_c_z,1,Irrep_ij,
     +                   z2_pphh_list)
         Else 

            Write(6,"(a,a)") "Insufficient amount of memory to perform",
     +                       " the contraction in memory."
            Call Errex
C
         Endif 
C
        Enddo
#ifdef _DEBUG_LVL0
      I000 = 1
      Isize_abij_z=Idsymsz(Irrepx,Isytyp(1,z2_pphh_list),
     +           Isytyp(2,z2_pphh_list))
      call getall(Work(I000),Isize_abij_z,Irrepx,z2_pphh_list)
      call checksum("Tdcc_form_mu_dot_da Z(A<B,I<J)(f):",Work(i000),
     +               Isize_abij_z,S)
#endif
      Enddo
C
      Return
      End
