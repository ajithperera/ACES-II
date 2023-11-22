      Subroutine Pccd_form_htau_2d_phph_ph_1(Hov,Hvo,Work,Maxcor,
     +                                       Nocc,Nvrt,Nbas,List_v1,
     +                                       List_v2,List_g,Fact)
      Implicit Double Precision(A-H,O-Z)

      Dimension Htau_pq(Nbas,Nbas)
      Dimension Work(Maxcor)
      Dimension Hov(Nocc*Nvrt)
      Dimension Hvo(Nocc*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione,Two,Half /1.0D0,-1.0D0,0.0D0,1,
     +                              2.0D0,0.50D0/

C V(ek,fa)*Gamma(ek,fi)=V(ek,fa)Gamma(ek,fi) (V is stored as fa,ek and G stored as ek,fi)

      Irrepx = Ione
      Ispin  = Ione
  
      Do Irrep_fi = 1, Nirrep
         Irrep_ek = Dirprd(Irrep_fi,Irrepx)
         Irrep_fa = Dirprd(Irrep_ek,Irrepx)

         Nrow_ek = Irpdpd(Irrep_ek,11)
         Ncol_fa = Irpdpd(Irrep_fa,13)
         Ncol_fi = Irpdpd(Irrep_fi,11)

         I000 = Ione
         I010 = I000 + Nrow_ek*Ncol_fa
         I020 = I010 + Nrow_ek*Ncol_fa
         I030 = I020 + Max(Ncol_fa,Nrow_ek)
         I040 = I030 + Max(Ncol_fa,Nrow_ek)
         Iend = I040 + Max(Ncol_fa,Nrow_ek)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hhpp_ph",
     +                                      Iend,Maxcor)

C Read in as V(fa,ek)

         Call Getlst(Work(I010),1,Nrow_ek,1,Irrep_ek,List_v1)
CSSS         Call Spinad3(Irrep_fa,Vrt(1,1),Ncol_fa,Nrow_ek,Work(I010),
CSSS     +                Work(I020),Work(I030))
C V(fa,ek) -> V(ek,fa)

         Call Transp(Work(I010),Work(I000),Nrow_ek,Ncol_fa)
#ifdef _DEBUG_LVLM
         Print*,List_v1
         call output(Work(I010),1,Ncol_fa,1,Nrow_ek,Ncol_fa,Nrow_ek,1)
#endif 

C Gamma(ek,fi) read in as (ek,fi)

         Iend = I010 + Nrow_ek*Ncol_fi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_fi,2,Irrep_fi,List_g)
#ifdef _DEBUG_LVLM
         Print*,List_g
         call output(Work(I010),1,Nrow_ek,1,Ncol_fa,Nrow_ek,Ncol_fa,1)
#endif 
C V(ek,fa)(t)*G(ek,fi) -> S(a,i)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_i = 1, Nirrep
            Irrep_f  = Dirprd(Irrep_i,Irrep_fi)
            Irrep_a  = Dirprd(Irrep_f,Irrep_fa)

            Na = Vrt(irrep_a,1)
            Ni = Pop(irrep_i,1)
            Nf = Vrt(irrep_f,1)
            
            Nsum = Nrow_ek*Nf
            Nrow = Na
            Ncol = Ni
            Icheck = Min(Nsum,Nrow,Ncol)
            
            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hvo(Ioff),Nrow)
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

C V(ni,mf)*Gamma(na,mf)=V(mf,ni)[stored as V(ni,mf)]*Gamma(mf,na)
C V(ni,fm)*Gamma(fm,an)->V(fm,ni)*Gamma(fm,na)

      Do Irrep_na = 1, Nirrep
         Irrep_fm = Dirprd(Irrep_na,Irrepx)
         Irrep_ni = Dirprd(Irrep_fm,Irrepx)

         Nrow_ni = Irpdpd(Irrep_ni,14)
         Nrow_na = Irpdpd(Irrep_na,18)
         Ncol_fm = Irpdpd(Irrep_fm,11)

         I000 = Ione
         I010 = I000 + Nrow_ni*Ncol_fm
         I020 = I010 + Nrow_ni*Ncol_fm
         I030 = I020 + Max(Ncol_fm,Nrow_ni)
         I040 = I030 + Max(Ncol_fm,Nrow_ni)
         Iend = I040 + Max(Ncol_fm,Nrow_ni)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_ph",
     +                                      Iend,Maxcor)
       
         Call Getlst(Work(I010),1,Ncol_fm,1,Irrep_fm,List_v2)
CSSS         Call Spinad3(Irrep_ni,Pop(1,1),Nrow_ni,Ncol_fm,Work(I010),
CSSS     +                Work(I020),Work(I030))

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_ph",
     +                                      Iend,Maxcor)
C V(ni,mf)-> V(ni,fm) -> V(fm,ni)
          
          Call Symtr1(Irrep_fm,Pop(1,Ispin),Vrt(1,Ispin),Nrow_ni,
     +                Work(I010),Work(I020),Work(I030),Work(I040))
          
          Call Transp(Work(I010),Work(I000),Ncol_fm,Nrow_ni)

 	 I020 = I010 + Nrow_na*Ncol_fm
         I030 = I020 + Max(Ncol_fm,Nrow_na)
         I040 = I030 + Max(Ncol_fm,Nrow_na)
         Iend = I040 + Max(Ncol_fm,Nrow_na)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Nrow_na,2,Irrep_na,List_g)

C Gamma(fm,an) -> Gamma(fm,na)

          Call Symtr1(Irrep_na,Vrt(1,Ispin),Pop(1,Ispin),Ncol_fm,
     +                Work(I010),Work(I020),Work(I030),Work(I040))

C V(fm,ni)*Gamma(fm,na) -> S(i,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_a = 1, Nirrep
            Irrep_n  = Dirprd(Irrep_a,Irrep_na)
            Irrep_i  = Dirprd(Irrep_n,Irrep_ni)

            Ni = Pop(irrep_i,1)
            Na = Vrt(irrep_a,1)
            Nn = Pop(irrep_n,1)
            
            Nsum = Ncol_fm*Nn
            Nrow = Ni
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hov(Ioff),Nrow)
            Endif

            Ioff = Ioff + Ni*Na
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

