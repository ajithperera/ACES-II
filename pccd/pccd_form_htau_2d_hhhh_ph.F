      Subroutine Pccd_form_htau_2d_hhhh_ph(Hvo,Work,Maxcor,Nocc,
     +                                     Nvrt,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hvo(Nocc*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(kj,la)*Gamma(kj,li)

      Irrepx = Ione

      Do Irrep_li = 1, Nirrep
         Irrep_kj = Dirprd(Irrep_li,Irrepx)
         Irrep_la = Dirprd(Irrep_kj,Irrepx)

         Nrow_kj = Irpdpd(Irrep_kj,14)
         Ncol_la = Irpdpd(Irrep_la,18)
         Ncol_li = Irpdpd(Irrep_li,14)

         I000 = Ione
         I010 = I000 + Nrow_kj*Ncol_la
         I020 = I010 + Max(Ncol_la,Nrow_kj)
         Iend = I020 + Max(Ncol_la,Nrow_kj)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh_ph",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_la,1,Irrep_la,List_v)
         Call Spinad3(Irrep_kj,Pop(1,1),Nrow_kj,Ncol_la,Work(I000),
     +                Work(I010),Work(I020))

         Iend = I010 + Nrow_kj*Ncol_li
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhhh_ph",
     +                                      Iend,Maxcor)
         
         Call Getlst(Work(I010),1,Ncol_li,2,Irrep_li,List_g)

C V(kjl,a)^(t)*G(kjl,i) -> S(a,i)

         Ioff = Ione
         Joff = I000
         Koff = I010

         Do irrep_i = 1, Nirrep
            Irrep_l  = Dirprd(Irrep_i,Irrep_li)
            Irrep_a  = Dirprd(Irrep_l,Irrep_la)

            Na = Vrt(irrep_a,1)
            Ni = Pop(irrep_i,1)
            Nl = Pop(irrep_l,1)

            Nsum = Nrow_kj*Nl
            Nrow = Na
            Ncol = Ni

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then

               Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                   Nsum,Work(Koff),Nsum,One,Hvo(Ioff),
     +                   Nrow)
            Endif 

            Ioff = Ioff + Ni*Na
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Ni
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htau_vo :",Hvo,Nocc*Nvrt)
      ndim=(Nocc*Nvrt)**(1/2)
      call pccd_check_htau("Htau_vo :",Hvo,Ndim,"VO","D")
#endif

      Return
      End
