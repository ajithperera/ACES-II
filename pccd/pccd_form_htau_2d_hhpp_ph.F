      Subroutine Pccd_form_htau_2d_hhpp_ph(Hov,Hvo,Work,Maxcor,
     +                                     Nocc,Nvrt,Nbas,List_v1,
     +                                     List_v2,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Hov(Nocc*Nvrt)
      Dimension Hvo(Nocc*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D00,1/

C V(fe,na)*Gamma(fe,ni) (Gamma stored as fe,ni and V is stored as fe,an)
C V(fe,an)*Gamma(fe,in)

      Call Dzero(Hvo,Nocc*Nvrt)
      Irrepx = Ione
      Ispin  = Ione

      Do Irrep_ni = 1, Nirrep
         Irrep_fe = Dirprd(Irrep_ni,Irrepx)
         Irrep_na = Dirprd(Irrep_fe,Irrepx)

         Nrow_fe = Irpdpd(Irrep_fe,15)
         Ncol_na = Irpdpd(Irrep_na,18)
         Ncol_ni = Irpdpd(Irrep_ni,14)

         I000 = Ione
         I010 = I000 + Nrow_fe*Ncol_na
         I020 = I010 + Max(Nrow_fe,Ncol_na)
         I030 = I020 + Max(Nrow_fe,Ncol_na)
         Iend = I030 + Max(Nrow_fe,Ncol_na)

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ph",
     +                                      Iend,Maxcor)

C V(fe,an) -> V(fe,na)
         
         Call Getlst(Work(I000),1,Ncol_na,1,Irrep_na,List_v1)
         Call Symtr1(Irrep_na,Vrt(1,1),Pop(1,1),Nrow_fe,Work(I000),
     +               Work(I010),Work(I020),Work(I030))

         I020 = I010 + Max(Ncol_na,Nrow_fe)
         Iend = I020 + Max(Ncol_na,Nrow_fe)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ph",
     +                                      Iend,Maxcor)
         Call Spinad3(Irrep_fe,Vrt(1,1),Nrow_fe,Ncol_na,Work(I000),
     +                Work(I010),Work(I020))

         I020 = I010 + Nrow_fe*Ncol_ni
         I030 = I020 + Max(Nrow_fe,Ncol_ni)
         I040 = I030 + Max(Nrow_fe,Ncol_ni)
         Iend = I030 + Max(Nrow_fe,Ncol_ni)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ph",
     +                                      Iend,Maxcor)

C G(fe,in) -> G(fe,ni)
         Call Getlst(Work(I010),1,Ncol_ni,2,Irrep_ni,List_g)
         Call Symtr1(Irrep_ni,Pop(1,Ispin),Pop(1,Ispin),Nrow_fe,
     +               Work(I010),Work(I020),Work(I030),Work(I040))

C V(fe,na)(t)*G(fe,ni) (stored as (fe,ni)) -> S(a,i))

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_i = 1, Nirrep
            Irrep_n  = Dirprd(Irrep_i,Irrep_ni)
            Irrep_a  = Dirprd(Irrep_n,Irrep_na)

            Na = Vrt(irrep_a,1)
            Ni = Pop(irrep_i,1)
            Nn = Pop(irrep_n,1)

            Nsum = Nrow_fe*Nn
            Nrow = Na
            Ncol = Ni

            Icheck =  Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvo(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Na*Ni
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Ni
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htau_vo :",Hvo,Nocc*Nvrt)
      ndim=(Nocc*Nvrt)**(1/2)
      call pccd_check_htau("Htau_vo :",Hvo,Ndim,"VO","D")
#endif

      Call Dzero(Hov,Nocc*Nvrt)

C V(fi,nm)*Gamma(fa,nm)[V(nm,fi)*Gamma(nm,fa)]  (V is stored as nm,if) 
C V(if,nm)*G(af,nm)

      Do Irrep_fa = 1, Nirrep
         Irrep_nm = Dirprd(Irrep_fa,Irrepx)
         Irrep_fi = Dirprd(Irrep_nm,Irrepx)

         Nrow_nm = Irpdpd(Irrep_nm,14)
         Ncol_fa = Irpdpd(Irrep_fa,15)
         Ncol_fi = Irpdpd(Irrep_fi,11)

         I000 = Ione
         I010 = I000 + Nrow_nm*Ncol_fi
         I020 = I010 + Max(Ncol_fi,Nrow_nm)
         I030 = I020 + Max(Ncol_fi,Nrow_nm)
         Iend = I030 + Max(Ncol_fi,Nrow_nm)
           
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_pphh_ph",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_fi,1,Irrep_fi,List_v2)
         Call Spinad3(Irrep_nm,Pop(1,1),Nrow_nm,Ncol_fi,Work(I000),
     +                Work(I010),Work(I020))

C V(nm,if) -> (nm,fi)
         Call Symtr1(Irrep_fi,Pop(1,Ispin),Vrt(1,Ispin),Nrow_nm,
     +               Work(I000),Work(I010),Work(I020),Work(I030))

         I020 = I010 + Nrow_nm*Ncol_fa
         I030 = I020 + Nrow_nm*Ncol_fa
         I040 = I030 + Max(Nrow_nm,Ncol_fa)
         I050 = I040 + Max(Nrow_nm,Ncol_fa)
         Iend = I050 + Max(Nrow_nm,Ncol_fa)
      
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ph",
     +                                      Iend,Maxcor)
C G stored as (af,nm)

         Call Getlst(Work(I020),1,Nrow_nm,2,Irrep_nm,List_g)

C G(af,nm) -> G(nm,fa)

         Call Transp(Work(I020),Work(I010),Nrow_nm,Ncol_fa)
         Call Symtr1(Irrep_fa,Vrt(1,Ispin),Vrt(1,Ispin),Nrow_nm,
     +               Work(I010),Work(I030),Work(I040),Work(I050))

C  V(nm,fi)*Gamma(nm,fa) -> S(i,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do irrep_a = 1, Nirrep
            Irrep_f  = Dirprd(Irrep_a,Irrep_fa)
            Irrep_i  = Dirprd(Irrep_f,Irrep_fi)

            Ni = Pop(irrep_i,1)
            Na = Vrt(irrep_a,1)
            Nf = Vrt(irrep_f,1)

            Nsum = Nrow_nm*Nf
            Nrow = Ni
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hov(Ioff),Nrow)
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

