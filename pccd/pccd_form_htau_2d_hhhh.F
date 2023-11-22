      Subroutine Pccd_form_htau_2d_hhhh(Hoo_pq,Hoo_qp,Work,Maxcor,Nocc,
     +                                  Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Two,Onem,Dnull,Ione /1.0D0,2.0D0,-1.0D0,0.0D0,1/

C V(mk,nj)*Gamma(mk,ni)
C Note that this is commented because of the OO gradient term added in 
C Pccd_form_htau_pq (Pccd_form_htau_oo). 

      Call Dzero(Hoo_qp,Nocc*Nocc)
      Call Dzero(Hoo_pq,Nocc*Nocc)

      Irrepx = Ione
      Do Irrep_ni = 1, Nirrep
         Irrep_mk = Dirprd(Irrep_ni,Irrepx)
         Irrep_nj = Dirprd(Irrep_mk,Irrepx)

         Nrow_mk = Irpdpd(Irrep_mk,21)
         Ncol_nj = Irpdpd(Irrep_nj,21)
         Nrow_ni = Irpdpd(Irrep_ni,21)
         Ncol_mk = Irpdpd(Irrep_mk,21)

         I000 = Ione
         I010 = I000 + Nrow_mk*Ncol_nj
         I020 = I010 + Max(Nrow_mk,Ncol_nj)
         Iend = I020 + Max(Nrow_mk,Ncol_nj)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_nj,1,Irrep_nj,List_v)
         Call Spinad1(Irrep_nj,Pop(1,1),Nrow_mk,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_ni*Ncol_mk
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)
         
         Call Getlst(Work(I010),1,Ncol_mk,2,Irrep_mk,List_g)

C V(mkn,j)^(t)*G(mkn,i) -> S(j,i)

         Ioff = Ione
         Joff = I000
         Koff = I010

         Do irrep_i = 1, Nirrep
            Irrep_n  = Dirprd(Irrep_i,Irrep_ni)
            Irrep_j  = Dirprd(Irrep_n,Irrep_nj)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nn = Pop(irrep_n,1)

            Nsum = Nrow_mk*Nn
            Ncol = Ni
            Nrow = Nj

            Icheck = Min(Nsum,Ni,Nj)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                   Nsum,Work(Koff),Nsum,One,Hoo_qp(Ioff),
     +                   Nrow)
            Endif 

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nrow_mk*Nn*Nj
            Koff = Koff + Ncol_mk*Nn*Ni
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      call pccd_check_htau("Htau_oo :",Hoo_qp,Nocc,"OO","D")
   
CSSS      call output(Hoo_qp,1,Nocc,1,Nocc,Nocc,Nocc,1)
#endif 
#ifdef _EXPLICITLY_PERMUTE
      Call Dzero(Hoo_pq,Nocc*Nocc)

C V(mi,nk)*Gamma(mj,nk)->(nk,mi)*Gamma(nk,mj)

      Do Irrep_mj = 1, Nirrep
         Irrep_nk = Dirprd(Irrep_mj,Irrepx)
         Irrep_mi = Dirprd(Irrep_nk,Irrepx)

         Nrow_mi = Irpdpd(Irrep_mi,21)
         Ncol_nk = Irpdpd(Irrep_nk,21)
         Nrow_nk = Irpdpd(Irrep_nk,21)
         Ncol_mj = Irpdpd(Irrep_mj,21)

         I000 = Ione
         I010 = I000 + Nrow_mi*Ncol_nk
         I020 = I010 + Ncol_nk
         Iend = I020 + Ncol_nk
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Nrow_mi,1,Irrep_mi,List_v)
         Call Spinad1(Irrep_mi,Pop(1,1),Ncol_nk,Work(I000),
     +                Work(I010),Work(I020))
        
         Iend = I010 + Nrow_nk*Ncol_mj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_mj,2,Irrep_mj,List_g)

C V(nkm,i)(t)*G(nkm,j) -> S(i,j)

         Ioff = Ione 
         Joff = I000
         Koff = I010
         Do irrep_j = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_j,Irrep_mj)
            Irrep_i  = Dirprd(Irrep_m,Irrep_mi)

            Nj = Pop(irrep_j,1)
            Ni = Pop(irrep_i,1)
            Nm = Pop(irrep_m,1)

            Nsum = Ncol_nk*Nm
            Nrow = Ni
            Ncol = Nj
            Icheck = Min(Nsum,Ni,Nj)

            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                      Nsum,Work(Koff),Nsum,One,Hoo_pq(Ioff),
     +                      Nrow)
            Endif

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Ncol_nk*Ni*Nm
            Koff = Koff + Nrow_nk*Nj*Nm 
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_oo :",Hoo_pq,Nocc,"OO","D")
CSSS      call output(Hoo_pq,1,Nocc,1,Nocc,Nocc,Nocc,1)
#endif 
#endif 

      Return
      End
