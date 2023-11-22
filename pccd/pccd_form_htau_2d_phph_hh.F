      Subroutine Pccd_form_htau_2d_phph_hh(Hoo_pq,Hoo_qp,Work,Maxcor,
     +                                     Nocc,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(ek,cj)*Gamma(ek,ci) 

      Irrepx = Ione
  
      Do Irrep_ek = 1, Nirrep
         Irrep_ci = Dirprd(Irrep_ek,Irrepx)
         Irrep_cj = Dirprd(Irrep_ek,Irrepx)

         Nrow_ek = Irpdpd(Irrep_ek,9)
         Ncol_cj = Irpdpd(Irrep_cj,9)
         Nrow_ci = Irpdpd(Irrep_ci,9)
         Ncol_ek = Irpdpd(Irrep_ek,9)

         I000 = Ione
         I010 = I000 + Nrow_ek*Ncol_cj 
         Iend = I010 + Nrow_ci*Ncol_ek
  
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_cj,1,Irrep_cj,List_v)
         Call Getlst(Work(I010),1,Ncol_ek,2,Irrep_ek,List_g)

C V(ek,cj)^(t)*G(ek,ci) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_i = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_i,Irrep_ci)
            Irrep_j  = Dirprd(Irrep_c,Irrep_cj)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nc = Vrt(irrep_c,1)
            
            Nsum = Nrow_ek*Nc
            Nrow = Nj
            Ncol = Ni
            Icheck = Min(Nsum,Nrow,Ncol)
            
            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),Nrow)
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
C V(ci,ek)*Gamma(cj,ek) (ignore labelling)

      Do Irrep_cj = 1, Nirrep
         Irrep_ek = Dirprd(Irrep_cj,Irrepx)
         Irrep_ci = Dirprd(Irrep_ek,Irrepx)

         Nrow_ci = Irpdpd(Irrep_ci,9)
         Ncol_ek = Irpdpd(Irrep_ek,9)
         Nrow_ek = Irpdpd(Irrep_ek,9)
         Ncol_cj = Irpdpd(Irrep_cj,9)

         I000 = Ione
         I010 = I000 + Nrow_ci*Ncol_ek
         Iend = I010 + Nrow_ek*Ncol_cj

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_hh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_cj,1,Irrep_cj,List_v)
         Call Getlst(Work(I010),1,Ncol_ek,2,Irrep_ek,List_g)

C  V(ek,ci)*Gamma(ek,cj) -> S(i,j)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_j = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_j,Irrep_cj)
            Irrep_i  = Dirprd(Irrep_c,Irrep_ci)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nc = Vrt(irrep_c,1)
            
            Nsum = Ncol_ek*Nc
            Nrow = Ni
            Ncol = Nj
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

