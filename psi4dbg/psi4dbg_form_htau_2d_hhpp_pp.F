      Subroutine Psi4dbg_form_htau_2d_hhpp_pp(Hvv_pq,Hvv_qp,Work,
     +                                        Maxcor,Nvrt,Nbas,
     +                                        List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/

C V(cb,km)*Gamma(ca,km) (Gamma stored as ca,km and V stored as cb,km)
   
      Irrepx = Ione
      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hvv_pq,Nvrt*Nvrt)

      Do Irrep_km = 1, Nirrep
         Irrep_ca = Dirprd(Irrep_km,Irrepx)
         Irrep_cb = Dirprd(Irrep_km,Irrepx)

         Nrow_cb = Irpdpd(Irrep_cb,13)
         Ncol_km = Irpdpd(Irrep_km,14)
         Nrow_ca = Irpdpd(Irrep_ca,13)
         Ncol_km = Irpdpd(Irrep_km,14)

         I000 = Ione
         I010 = I000 + Nrow_cb*Ncol_km
         I020 = I010 + Nrow_cb*Ncol_km
         I030 = I020 + Max(Ncol_km,Nrow_cb)
         Iend = I030 + Max(Ncol_km,Nrow_cb)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_pp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_km,1,Irrep_km,List_v)
         Call Spinad3(Irrep_cb,Vrt(1,1),Nrow_cb,Ncol_km,Work(I010),
     +                Work(I020),Work(I030))
C V(cb,km) -> V(km,cb)
         Call Transp(Work(I010),Work(I000),Ncol_km,Nrow_cb)

         I020 = I010 + Nrow_ca*Ncol_km
         Iend = I020 + Nrow_ca*Ncol_km
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ab",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_km,2,Irrep_km,List_g)
C G(cb,km) -> G(km,cb)
         Call Transp(Work(I020),Work(I010),Ncol_km,Nrow_ca)

C V(km,cb)(t)*G(km,ca) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_a,Irrep_ca)
            Irrep_b  = Dirprd(Irrep_c,Irrep_cb)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nc = Vrt(irrep_c,1)

            Nsum = Ncol_km*Nc
            Nrow = Nb
            Ncol = Na

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
            Endif
            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Nb
            Koff = Koff + Nsum*Na
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_vv :",Hvv_qp,Nvrt,"VV","D")
CSSS      call checksum("Htau_vv :",Hvv_qp,Nvrt*Nvrt)
#endif

#ifdef _EXPLICITLY_PERMUTE

C V(km,ca)*G(km,cb) (V(km,ca) stored as ca,km and Gamma stored as cb,km)

      Do Irrep_km = 1, Nirrep
         Irrep_cb = Dirprd(Irrep_km,Irrepx)
         Irrep_ca = Dirprd(Irrep_km,Irrepx)

         Nrow_km = Irpdpd(Irrep_km,21)
         Ncol_ca = Irpdpd(Irrep_ca,19)
         Nrow_cb = Irpdpd(Irrep_cb,19)
         Ncol_km = Irpdpd(Irrep_km,21)

         I000 = Ione
         I010 = I000 + Nrow_km*Ncol_ca
         I020 = I010 + Nrow_km
         Iend = I020 + Nrow_km
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_pp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Nrow_km,1,Irrep_km,List_v)
         Call Spinad3(Irrep_ca,Vrt(1,1),Ncol_ca,Nrow_km,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_cb*Ncol_km
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ab",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_km,2,Irrep_km,List_g)

C V(km,ca)(t)*G(km,cb)-> S(a,b)

         Ioff = Ione
         Joff = I000
         Koff = I010

         Do irrep_b = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_b,Irrep_cb)
            Irrep_a  = Dirprd(Irrep_c,Irrep_ca)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nc = Vrt(irrep_c,1)

            Nsum = Nrow_km*Nc
            Nrow = Na
            Ncol = Nb

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvv_pq(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Nb
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_vv :",Hvv_pq,Nvrt,"VV","F")
#endif
#endif 

      Return
      End

