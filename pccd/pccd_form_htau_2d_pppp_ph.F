      Subroutine Pccd_form_htau_2d_pppp_ph(Hov,Work,Maxcor,Nocc,
     +                                     Nvrt,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hov(Nvrt*Nocc)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(ci,db)*G(ca,db) = V(db,ci)*G(db,ca)

      Irrepx = Ione

      Do Irrep_ca = 1, Nirrep
         Irrep_db = Dirprd(Irrep_ca,Irrepx)
         Irrep_ci = Dirprd(Irrep_db,Irrepx)

         Nrow_db = Irpdpd(Irrep_db,15)
         Ncol_ci = Irpdpd(Irrep_ci,11)
         Ncol_ca = Irpdpd(Irrep_ca,15)

         I000 = Ione
         I010 = I000 + Nrow_db*Ncol_ci
         I020 = I010 + Max(Ncol_ci,Nrow_db)
         Iend = I020 + Max(Ncol_ci,Nrow_db)

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pppp_ph",
     +                                      Iend,Maxcor)
 
         Call Getlst(Work(I000),1,Ncol_ci,1,Irrep_ci,List_v)

         Call Spinad3(Irrep_db,Vrt(1,1),Nrow_db,Ncol_ci,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_db*Ncol_ca
         Call Getlst(Work(I010),1,Ncol_ca,2,Irrep_ca,List_g)

C V(db,ci)(t)*G(db,ca) -> S(i,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_c  = Dirprd(Irrep_a,Irrep_ca)
            Irrep_i  = Dirprd(Irrep_c,Irrep_ci)

            Nc = Vrt(Irrep_c,1)
            Na = Vrt(Irrep_a,1)
            Ni = Pop(Irrep_i,1)

            Nsum = Nrow_db*Nc
            Nrow = Ni
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt .0) Then

                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                      Nsum,Work(Koff),Nsum,One,Hov(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Na*Ni
            Joff = Joff + Nsum*Ni
            Koff = Koff + Nsum*Na
 
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htau_ov :",Hov,Nocc*Nvrt)
      ndim=(Nocc*Nvrt)**(1/2)
      call pccd_check_htau("Htau_ov :",Hov,Ndim,"OV","D")
#endif
      Return
      End
