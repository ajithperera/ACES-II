      Subroutine Tdcc_form_mutilde_dot_db(Work,Memleft,Irrepx,Iuhf,
     +                                    Ioff_vo)
    
      Implicit Double Precision (A-H, O-Z)
      Logical Rhf 
C
      Dimension Work(Memleft),Ioff_vo(8,2)
      Integer E, Dis_bgn
      Integer Hbar_ppph_list,z2_pphh_list

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"
C
      Data One, Onem, Zero/1.0D0, -1.0D0, 0.0D0/

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_form_mutilde_dot_db---"
        Write(6,*)
#endif
      Rhf = .False.
      If (Iuhf .Eq. 0) Rhf = .True.

C All the expression must be transposed.
C Contaction in spin-orbitals Z(AB,IJ) = P(IJ) {Mu(E,I) * Hbar(AB,JE)}
C ABAB spin block
C Z(Ab,Ij) =  Mu(E,I) * Hbar(Ab,Ej) + Mu(e,j) Hbar(Ab,Ie)

      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)
         Irrep_ej = Irrep_ab
         
         Hbar_ppph_list = 30
         z2_pphh_list   = 346
     
         Ndim_ab_h = Irpdpd(Irrep_ab,Isytyp(1,Hbar_ppph_list))
         Ndim_ej_h = Irpdpd(Irrep_ej,Isytyp(2,Hbar_ppph_list))
         Ndim_ab_z = Irpdpd(Irrep_ab,Isytyp(1,z2_pphh_list))
         Ndim_ij_z = Irpdpd(Irrep_ij,Isytyp(2,z2_pphh_list))

         Max_abej_h = Max(Ndim_ab_h,Ndim_ej_h)
         Max_abij_z = Max(Ndim_ab_z,Ndim_ij_z)
         Max_both   = Max(Max_abej_h,Max_abij_z)

         Ndim_abej_h = Ndim_ab_h * Ndim_ej_h
         Ndim_abij_z = Ndim_ij_z * Ndim_ab_z
     
         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max(3*Max_both,Ndim_abij_z)
         I020 = I010 + Max(Ndim_abej_h,Ndim_abij_z,3*Max_both)
         Iend = I020
         
         If (Iend .LT. Memleft) Then
            Call Getlst(Work(I010),1,Ndim_ej_h,1,Irrep_ej,
     +                  Hbar_ppph_list)

C List 130 = Hbar(Ab,Ej); Hbar(Ab,Ej) -> Hbar(Ab,jE)

            Ibgn  = I000
            IPtr1 = Ibgn
            Iptr2 = Iptr1 + Max_both
            Iptr3 = Iptr2 + Max_both
            Iptr4 = Iptr3 + Max_both
            Iend  = Iptr4

            Call Symtr1(Irrep_ej,Vrt(1,1),Pop(1,2),Ndim_ab_h,
     +                  Work(I010),Work(Iptr1),Work(Iptr2),
     +                  Work(Iptr3))
            Call Dzero(Work(I000),Ndim_abij_z)

C Z(Ab,jI) =  Hbar(Ab,jE) * Mu(E,I)
            
            Do Irrep_i = 1, Nirrep
               Irrep_e = Dirprd(Irrep_i,Irrepx)
               Irrep_j = Dirprd(Irrep_e,Irrep_ej)
               
                Nrow = Ndim_ab_h*Pop(Irrep_j,2)
                Ncol = Pop(Irrep_i,1) 
                Nsum = Vrt(Irrep_e,1)

                Loc_h = I010 + Ndim_ab_h*(Isymoff(Irrep_e,
     +                                    Irrep_ej,25)-1)
                Loc_z = I000 + Ndim_ab_z*(Isymoff(Irrep_i,
     +                                    Irrep_ij,24)-1)
                Loc_m = Ioff_vo(Irrep_i,1)

                If (Nsum .NE. 0) Then
                Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(loc_h),
     +                     Nrow,Work(Loc_m),Nsum,Zero,Work(loc_z),
     +                     Nrow)
                Endif 
            Enddo
  
C Z(Ab,jI) -> Z(Ab,Ij)

            Ibgn  = I010
            Iptr1 = Ibgn
            Iptr2 = Iptr1 + Max_both
            Iptr3 = Iptr2 + Max_both
            Iptr4 = Iptr3 + Max_both
            Iend  = Iptr4

            Call Symtr1(Irrep_ij,Pop(1,2),Pop(1,1),Ndim_ab_z,
     +                  Work(I000),Work(IPtr1),Work(IPtr2),
     +                  Work(Iptr3))

         Else 
#ifdef _DEBUG_LVLM
      write(6,"(a)") "Enter the out-of-core block-1"
#endif
C  Hbar(Ab,Ej) * Mu(E,I) per given j (stored as Ab,Ej).

            Call Dzero(Work(I000),Ndim_ab_z*Ndim_ij_z)
            Dis_bgn = 1
            Loc_z = I000

            Do Irrep_j = 1, Nirrep 
               Irrep_e = Dirprd(Irrep_j,Irrep_ej)
               Irrep_i = Dirprd(Irrep_j,Irrep_ij)
                
               Ndim_i = Pop(Irrep_i,1)
               Ndim_j = Pop(Irrep_j,2)
               Ndim_e = Vrt(Irrep_e,1)
               
               Iend = I010 + Max(Ndim_abij_z,Ndim_ab_h*Ndim_e)

               If (Iend .GT. Memleft) Then
                  Call Insmem("@-Tdvee_form_mu_dot_db",Iend,
     +                         Memleft)
               Endif 
  
               Do j=1, Ndim_j
                  Call Getlst(Work(I010),Dis_bgn,Ndim_e,1,Irrep_ej,
     +                        Hbar_ppph_list)
                  Nrow = Ndim_ab_h
                  Ncol = Ndim_i
                  Nsum = Ndim_e

                  Loc_m = Ioff_vo(Irrep_i,1)
                  If (Nsum .NE. 0) Then
                  Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(I010),
     +                        Nrow,Work(Loc_m),Nsum,Zero,Work(loc_z),
     +                        Nrow)
                  Endif 
                  Dis_bgn = Dis_bgn + Ndim_e
                  Loc_z   = Loc_z   + Ndim_ab_z*Ndim_i
               Enddo 
            Enddo    

         Endif 
#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_dot_db Z(Ab,Ij)(1):",
     +               Work(I000),Ndim_abij_z,S)
#endif
         If (Rhf) Then

C Z(Ab,Ij)
            Ibgn  = I010
            Iptr1 = Ibgn
            Iptr2 = Iptr1 + Max_both
            Iptr3 = Iptr2 + Max_both
            Iptr4 = Iptr3 + Max_both
            Iend  = Iptr4

           Call Symrhf3(Irrep_ab,Irrep_ij,Vrt(1,1),Pop(1,1),
     +                  Ndim_ab_z,Work(I000),Work(Iptr1),Work(IPtr2),
     +                  Work(Iptr3))

C The target, Z(Ab,Ij) 

#ifdef _DEBUG_LVL0

#else 
           Call Getlst(Work(I010),1,Ndim_ij_z,1,Irrep_ij,z2_pphh_list)
           Call Daxpy(Ndim_ab_z*Ndim_ij_z,One,Work(I010),1,Work(I000),1)
#endif 
           Call Putlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,z2_pphh_list)


         Else

C Z(Ab,Ij) = Mu(e,j) Hbar(Ab,Ie)
C List 129 = Hbar(Ab,Ie) no permutation of ket indices needed. 

           Hbar_ppph_list = 127+2

           Irrep_ie = Irrep_ab
           Ndim_ab_h = Irpdpd(Irrep_ab,Isytyp(1,Hbar_ppph_list))
           Ndim_ie_h = Irpdpd(Irrep_ie,Isytyp(2,Hbar_ppph_list))

           Ndim_abie_h = Ndim_ab_h * Ndim_ie_h
           Max_pqrs    = Max(Ndim_abie_h,Ndim_abij_z,3*Max_both)

           Ibgn = I010 
           I020 = I010 + Max_pqrs
           Iend = I020

         If (Iend .LT. Memleft) Then
            Call Getlst(Work(I010),1,Ndim_ie_h,1,Irrep_ie,
     +                  Hbar_ppph_list)

C Z(Ab,Ij) = Hbar(Ab,Ie) * Mu(e,j)

            Do Irrep_j = 1, Nirrep
               Irrep_e = Dirprd(Irrep_j,Irrepx)
               Irrep_i = Dirprd(Irrep_e,Irrep_ie)

                Nrow = Ndim_ab_h*Pop(Irrep_i,1)
                Ncol = Pop(Irrep_j,2)
                Nsum = Vrt(Irrep_e,2)

                Loc_h = I010 + Ndim_ab_h*(Isymoff(Irrep_e,
     +                                    Irrep_ie,18)-1)
                Loc_z = I000 + Ndim_ij_z*(Isymoff(Irrep_i,
     +                                    Irrep_ij,14)-1)
                Loc_m = Ioff_vo(Irrep_j,2)
                If (Nsum .NE. 0) Then
                Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(loc_h),
     +                     Nrow,Work(Loc_m),Nsum,One,Work(loc_z),
     +                     Nrow)
                Endif 
            Enddo 

         Else
#ifdef _DEBUG_LVLM
      write(6,"(a)") "Enter the out-of-core block-2"
#endif

C Z(Ab,Ij) = Hbar(Ab,Ie) * Mu(e,j) per given e 

            Ibgn    = 1
            I000    = Ibgn
            Dis_bgn = 1

            Do Irrep_e = 1, Nirrep
               Irrep_i = Dirprd(Irrep_e,Irrep_ie)
               Irrep_j = Dirprd(Irrep_e,Irrepx)

               Ndim_i = Pop(Irrep_i,1)
               Ndim_j = Pop(Irrep_j,2)
               Ndim_e = Vrt(Irrep_e,2)

               Loc_z = I000 + Ndim_ab_z*(Isymoff(Irrep_j,Irrep_ij,14)-1)
               Iend  = I010 + Max(Ndim_abij_z,Ndim_ab_h*Ndim_i)

 
               If (Iend .Gt. Memleft) Then
                  Call Insmem("@-Tdvee_form_mu_dot_db",Iend,
     +                         Memleft)
               Endif
               
               Do E=1, Ndim_e
                  Call Getlst(Work(I010),Dis_bgn,Ndim_i,1,Irrep_ie,
     +                        Hbar_ppph_list)
                  Nrow = Ndim_ab_h*Ndim_i
                  Ncol = Ndim_j
                  Nsum = 1

                  Loc_m = Ioff_vo(Irrep_j,2) + (E-1)
                  If (Nsum .NE. 0) Then
                  Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(I010),
     +                        Nrow,Work(Loc_m),Ndim_e,One,
     +                        Work(loc_z),Nrow)
                  Endif 
                  Dis_bgn = Dis_bgn + Ndim_i
               Enddo
            Enddo 

         Endif

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_dot_db Z(Ab,Ij)(2):",
     +               Work(I000),Ndim_abij_z,S)
#endif
#ifdef _DEBUG_LVL0
            Call Putlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,
     +                       z2_pphh_list)
#else
            Call Sumsym2(Work(I000),Work(I010),Ndim_ab_z*Ndim_ij_z,1,
     +                   Irrep_ij,z2_pphh_list)
#endif 
        Endif
      
      Enddo 

#ifdef _DEBUG_LVL0
      I000 = 1
      Isize_abij_z=Idsymsz(Irrepx,Isytyp(1,z2_pphh_list),
     +           Isytyp(2,z2_pphh_list))
      call getall(Work(I000),Isize_abij_z,Irrepx,z2_pphh_list)
      call checksum("Tdcc_form_mutilde_dot_db Z(Ab,Ij)(f)  :",
     +               Work(i000),Isize_abij_z,S)
#endif

      If (Rhf) Return
C
      Do Ispin = 1, 1+Iuhf 

C AAAA and BBBB spin blocks
C Z(AB,IJ) =  Mu(E,I) * Hbar(AB,EJ) (and Mu(e,i) * Hbar(ab,ej))

      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)
         Irrep_ej = Irrep_ab

         Hbar_ppph_list = 26  + Ispin 
         z2_pphh_list   = 343 + Ispin

         Ndim_ab_c_h = Irpdpd(Irrep_ab,Isytyp(1,Hbar_ppph_list))
         Ndim_ej_f_h = Irpdpd(Irrep_ej,Isytyp(2,Hbar_ppph_list))
         Ndim_ab_c_z = Irpdpd(Irrep_ab,Isytyp(1,z2_pphh_list))
         Ndim_ij_c_z = Irpdpd(Irrep_ij,Isytyp(2,z2_pphh_list))
         Ndim_ij_f_z = Irpdpd(Irrep_ij,20+Ispin)

         Max_pq = Max(Ndim_ab_c_h,Ndim_ej_f_h,Ndim_ij_c_z,Ndim_ij_f,
     +                Ndim_ab_c_z)

         Ndim_abej_cf_h = Ndim_ab_c_h * Ndim_ej_f_h
         Ndim_abij_cf_z = Ndim_ij_f_z * Ndim_ab_c_z
         Ndim_abij_cc_z = Ndim_ij_c_z * Ndim_ab_c_z

         Ibgn = 1
         I000 = Ibgn
         I010 = I000 + Max(3*Max_pq,Ndim_abij_cf_z)
         I020 = I010 + Max(Ndim_abej_cf_h,Ndim_abij_cf_z)
         Iend = I020

C List 127,128 = Hbar(A<B,EJ) (and Hbar(a<b,ej))

         If (Iend .LT. Memleft) Then
            Call Getlst(Work(I010),1,Ndim_ej_f_h,1,Irrep_ej,
     +                  Hbar_ppph_list)

C Hbar(A<B,EJ) (and Hbar(a<b,ej)) -> Hbar(A<B,JE) (and Hbar(a<b,je))

            Ibgn = I000
            I020 = Ibgn
            I030 = I020 + Max_pq
            I040 = I030 + Max_pq
            I050 = I040 + Max_pq
            Iend = I050

            Call Symtr1(Irrep_ej,Vrt(1,Ispin),Pop(1,Ispin),
     +                  Ndim_ab_c_h,Work(I010),Work(I020),Work(I030),
     +                  Work(I040))
            Call Dzero(Work(I000),Ndim_abij_cf_z)

C Hbar(A<B,JE) Mu(E,I) (and Hbar(a<b,je) Mu(e,i))

            Do Irrep_i = 1, Nirrep 
               Irrep_e = Dirprd(Irrep_i,Irrepx)
               Irrep_j = Dirprd(Irrep_e,Irrep_ej)
               
               Nrow = Ndim_ab_c_h*Pop(Irrep_j,Ispin)   
               Ncol = Pop(Irrep_i,Ispin)
               Nsum = Vrt(Irrep_e,Ispin)

               Loc_h = I010 + Ndim_ab_c_h*(Isymoff(Irrep_e,Irrep_ej,
     +                 15+Ispin)-1)
               Loc_z = I000 + Ndim_ab_c_z*(Isymoff(Irrep_i,Irrep_ij,
     +                 20+Ispin)-1)
               Loc_m = Ioff_vo(Irrep_i,Ispin)

               If (Nsum .NE. 0) Then
               Call Xgemm("N","N",Nrow,Ncol,Nsum,Onem,Work(loc_h),
     +                     Nrow,Work(Loc_m),Nsum,Zero,Work(loc_z),
     +                     Nrow)
               Endif 
            Enddo 

         Else

#ifdef _DEBUG_LVLM
      write(6,"(a)") "Enter the out-of-core block"
#endif

C Hbar(A<B,EJ) Mu(E,I) (and Hbar(a<b,ej) Mu(e,i)) per J
            
            Call Dzero(Work(I000),Ndim_abij_cf_z)
            Ibgn    = 1
            I000    = Ibgn
            Dis_bgn = 1
            Loc_z   = I000

            Do Irrep_j = 1, Nirrep
               Irrep_e = Dirprd(Irrep_j,Irrep_ej)
               Irrep_i = Dirprd(Irrep_e,Irrepx)

               Ndim_i = Pop(Irrep_i,Ispin)
               Ndim_j = Pop(Irrep_j,Ispin)
               Ndim_e = Vrt(Irrep_e,Ispin)

               Iend = I010 + Max(Ndim_abij_cf_z,Ndim_ab_cf_h*Ndim_e)

               If (Iend .Gt. Memleft) Then
                  Call Insmem("@-Tdvee_form_mu_dot_db",Iend,
     +                         Memleft)
               Endif

               Do j=1, Ndim_j
                  Call Getlst(Work(I010),Dis_bgn,Ndim_e,1,Irrep_ej,
     +                        Hbar_ppph_list)
                  Nrow = Ndim_ab_c_h
                  Ncol = Ndim_i
                  Nsum = Ndim_e

                  Loc_m = Ioff_vo(Irrep_i,Ispin)
                  If (Nsum .NE. 0) Then
                  Call Xgemm("N","N",Nrow,Ncol,Nsum,One,Work(I010),
     +                        Nrow,Work(Loc_m),Nsum,Zero,Work(loc_z),
     +                        Nrow)
                  Endif  
                  Dis_bgn = Dis_bgn + Ndim_e
                  Loc_z   = Loc_z   + Ndim_ab_c_z*Ndim_i
               Enddo
            Enddo 
          
         Endif 

C Antisymmetrize and pack:Z(A<B,IJ) -> Z(A<B,I<J) (Z(a<b,ij) -> Z(a<b,j<j) 

         Call Assym2(Irrep_ij,Pop(1,Ispin),Ndim_ab_c_z,Work(I000))

#ifdef _DEBUG_LVLM
      call Checksum("Tdcc_form_mutilde_dot_db Z(AB,IJ)(3):",
     +               Work(I000),Ndim_abij_cc_z,S)
#endif

#ifdef _DEBUG_LVL0
         Call Putlst(Work(I000),1,Ndim_ij_c_z,1,Irrep_ij,
     +               z2_pphh_list)
#else
         If (Ndim_abij_cc_z .GT. 0) Then
            Call Sumsym2(Work(I000),Work(I010),Ndim_abij_cc_z,1,
     +                   Irrep_ij,z2_pphh_list)
         Endif 
#endif 
       Enddo 

#ifdef _DEBUG_LVL0
      I000 = 1
      Isize_abij_z=Idsymsz(Irrepx,Isytyp(1,z2_pphh_list),
     +           Isytyp(2,z2_pphh_list))
      call getall(Work(I000),Isize_abij_z,Irrepx,z2_pphh_list)
      call checksum("Tdcc_form_mutilde_dot_db Z(A<B,I<J)(f):",
     +               Work(i000),Isize_abij_z,S)
#endif

      Enddo
C
      Return
      End
