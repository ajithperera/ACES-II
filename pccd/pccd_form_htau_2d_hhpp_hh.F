      Subroutine Pccd_form_htau_2d_hhpp_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                                      Nocc,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/

C V(ce,kj)*Gamma(ce,ki) (Gamma stored as ce,ki)
   
      Irrepx = Ione

      Do Irrep_ki = 1, Nirrep
         Irrep_ce = Dirprd(Irrep_ki,Irrepx)
         Irrep_kj = Dirprd(Irrep_ce,Irrepx)

         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_kj = Irpdpd(Irrep_kj,21)
         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_ki = Irpdpd(Irrep_ki,21)

         I000 = Ione
         I010 = I000 + Nrow_ce*Ncol_kj
         I020 = I010 + Max(Ncol_kj,Nrow_ce)
         Iend = I020 + Max(Ncol_kj,Nrow_ce)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_kj,1,Irrep_kj,List_v)
         Call Spinad3(Irrep_ce,Vrt(1,1),Nrow_ce,Ncol_kj,Work(I000),
     +                Work(I010),Work(I020))
         Iend = I010 + Nrow_ce*Ncol_ki
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ij",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_ki,2,Irrep_ki,List_g)

C V(ce,kj)(t)*G(ce,ki) (stored as (ce,ki)) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_i = 1, Nirrep
            Irrep_k  = Dirprd(Irrep_i,Irrep_ki)
            Irrep_j  = Dirprd(Irrep_k,Irrep_kj)

            Ni = Pop(irrep_i,1)
            Nj = Pop(irrep_j,1)
            Nk = Pop(irrep_k,1)

            Nsum = Nrow_ce*Nk
            Nrow = Nj
            Ncol = Ni

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nsum*Nj
            Koff = Koff + Nsum*Ni
         Enddo 
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_oo :",Hoo_qp,Nocc,"OO","D")
#endif

#ifdef _EXPLICITLY_PERMUTE
C V(ki,ce)*Gamma(kj,ce) (V is stored as ce,ki and Gamma is stored as ce,kj)

      Do Irrep_kj = 1, Nirrep
         Irrep_ce = Dirprd(Irrep_kj,Irrepx)
         Irrep_ki = Dirprd(Irrep_ce,Irrepx)

         Ncol_kj = Irpdpd(Irrep_kj,21)
         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_ki = Irpdpd(Irrep_ki,21)

         I000 = Ione
         I010 = I000 + Ncol_ki*Nrow_ce
         I030 = I020 + Ncol_ki
         Iend = I030 + Ncol_ki
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pphh_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ki,1,Irrep_ki,List_v)
         Call Spinad3(Irrep_ce,Vrt(1,1),Nrow_ce,Ncol_ki,Work(I000),
     +                Work(I020),Work(I030))

         Iend = I010 + Ncol_kj*Nrow_ce
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_kj,2,Irrep_kj,List_g)

C  V(ki,ce)*Gamma(ce,kj) -> S(i,j)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_j = 1, Nirrep
            Irrep_k  = Dirprd(Irrep_j,Irrep_kj)
            Irrep_i  = Dirprd(Irrep_k,Irrep_ki)

            Ni = Pop(irrep_i,1)
            Nj=  Pop(irrep_j,1)
            Nk = Pop(irrep_k,1)

            Nsum = Nrow_ce*Nk
            Nrow = Ni
            Ncol = Nk
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hoo_pq(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nsum*Ni
            Koff = Koff + Nsum*Nj
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_oo :",Hoo_pq,Nocc,"OO","F")
#endif
#endif 
      Return
      End

