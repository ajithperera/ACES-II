      Subroutine Pccd_form_htau_2d_pppp(Hvv_pq,Hvv_qp,Work,Maxcor,Nvrt,
     +                                  Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Two,Onem,Dnull,Ione /1.0D0,2.0D0,-1.0D0,0.0D0,1/

C V(ce,da)*Gamma(ce,db) (ignore the labels).

      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hvv_pq,Nvrt*Nvrt)

      Irrepx = Ione
      Do Irrep_ce = 1, Nirrep
         Irrep_db = Dirprd(Irrep_ce,Irrepx)
         Irrep_da = Dirprd(Irrep_ce,Irrepx)

         Nrow_ce = Irpdpd(Irrep_ce,19)
         Ncol_db = Irpdpd(Irrep_db,19)
         Nrow_da = Irpdpd(Irrep_da,19)
         Ncol_ce = Irpdpd(Irrep_ce,19)

         I000 = Ione
         I010 = I000 + Nrow_ce*Ncol_db
         I020 = I010 + Max(Ncol_db,Nrow_ce)
         I030 = I020 + Max(Ncol_db,Nrow_ce)
         Iend = I030 + Max(Ncol_db,Nrow_ce)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pppp",
     +                                      Iend,Maxcor)
 
         Call Getlst(Work(I000),1,Ncol_db,1,Irrep_db,List_v)

         Call Spinad3(Irrep_ce,Vrt(1,1),Nrow_ce,Ncol_db,Work(I000),
     +                Work(I010),Work(I020),Work(I030))

         Iend = I010 + Nrow_da*Ncol_ce
         Call Getlst(Work(I010),1,Ncol_ce,2,Irrep_ce,List_g)

C V(ce,da)(t)*G(ce,db) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_d  = Dirprd(Irrep_a,Irrep_da)
            Irrep_b  = Dirprd(Irrep_d,Irrep_db)

            Nb = Vrt(Irrep_b,1)
            Na = Vrt(Irrep_a,1)
            Nd = Vrt(Irrep_d,1)

            Nsum = Nrow_ce*Nd
            Nrow = Nb
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt .0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                      Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
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
      Call Dzero(Hvv_pq,Nvrt*Nvrt)

C V(ce,da)*Gamma(ce,db) (ignore the labels)

      Do Irrep_cb = 1, Nirrep
         Irrep_de = Dirprd(Irrep_cb,Irrepx)
         Irrep_ca = Dirprd(Irrep_de,Irrepx)

         Nrow_ca = Irpdpd(Irrep_ca,19)
         Ncol_de = Irpdpd(Irrep_de,19)
         Nrow_de = Irpdpd(Irrep_de,19)
         Ncol_cb = Irpdpd(Irrep_cb,19)

         I000 = Ione
         I010 = I000 + Nrow_ca*Ncol_de
         I020 = I010 + Ncol_de
         Iend = I020 + Ncol_de 
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pppp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_de,1,Irrep_de,List_v)
         Call Spinad3(Irrep_ca,Vrt(1,1),Nrow_ca,Ncol_de,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_de*Ncol_cb
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pppp",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_cb,2,Irrep_cb,List_g)

C V(ce,da)*G(ce,db) -> S(a,b)

         Ioff = Ione
         Joff = I000
         Koff = I010
          
         Do irrep_b = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_b,Irrep_cb)
            Irrep_a  = Dirprd(Irrep_c,Irrep_ca)

            Na = Vrt(Irrep_a,1)
            Nb = Vrt(Irrep_b,1)
            Nc = Vrt(Irrep_c,1)

            Nsum = Ncol_de*Nc
            Nrow = Na
            Ncol = Nb
            Icheck = Min(Nsum,Na,Nb)

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
