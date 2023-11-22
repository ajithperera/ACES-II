      Subroutine Pccd_form_htau_pqrs_hhhh_0(Work,Maxcor,List_v1,List_v2,
     +                                      List_g1,List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Four,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                      2.0D0,4.0D0,-2.0D0/

      Irrepx = Ione
      Ispin  = Ione 

      Nsize_t = Idsymsz(Irrepx,14,14)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C <jl|mn>*G(ik,mn)->H(jl,ik)
C <mn|ik>*G(mn,jl)->H(ik,jl)

      Ioff = I000 
      Do Irrep_mn = 1, Nirrep
         Irrep_jl = Dirprd(Irrep_mn,Irrepx) 
         Irrep_ik = Dirprd(Irrep_mn,Irrepx)

         Nrow_jl = Irpdpd(Irrep_jl,14)
         Nrow_ik = Irpdpd(Irrep_ik,14)
         Ncol_mn = Irpdpd(Irrep_mn,14)

         I020 = I010 + Nrow_jl*Ncol_mn
         I030 = I020 + Max(Ncol_mn,Nrow_jl)
         Iend = I030 + Max(Ncol_mn,Nrow_jl)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_mn,1,Irrep_mn,List_v1)

         Call Spinad3(Irrep_jl,Pop(1,1),Nrow_jl,Ncol_mn,Work(I000),
     +                Work(I010),Work(I020))

         I030 = I020 + Nrow_ik*Ncol_mn
         Iend = I030 + Nrow_jl*Nrow_ik
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_mn,2,Irrep_mn,List_g1)

         Icheck = Min(Nrow_jl,Nrow_ik,Ncol_mn)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_jl,Nrow_ik,Ncol_mn,One,Work(I010),
     +                  Nrow_jl,Work(I020),Nrow_ik,Dnull,Work(Ioff),
     +                  Nrow_jl)
         Endif 

C H(jl,ik) + H(jl,ik)^t

         Call Transp(Work(Ioff),Work(I030),Nrow_ik,Nrow_jl)
         Call Daxpy(Nrow_ik*Nrow_jl,One,Work(I030),1,Work(Ioff),1)
         Ioff = Ioff + Nrow_ik*Nrow_jl
      Enddo 

C <ef|jk>*G(ef,il) + <il|ef>*G(jk,ef)

      Ioff = I000
      Do Irrep_il = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_il,Irrepx)
         Irrep_jk = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,15)
         Ncol_il = Irpdpd(Irrep_il,14)
         Ncol_jk = Irpdpd(Irrep_jk,14)

         I020 = I010 + Nrow_ef*Ncol_jk
         I030 = I020 + Max(Nrow_ef,Ncol_jk)
         Iend = I030 + Max(Nrow_ef,Ncol_jk)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_jk,1,Irrep_jk,List_v2)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_jl,Work(I010),
     +                Work(I020),Work(I030))

         I030 = I020 + Nrow_ef*Ncol_il
         Iend = I030 + Ncol_il*Ncol_jk
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I020),1,Ncol_il,2,Irrep_il,List_g2)

         Icheck = Min(Nrow_ef,Ncol_il,Ncol_jk)

C <ef|jk>*G(ef,il)->H(jk,il)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_jk,Ncol_il,Nrow_ef,One,Work(I010),
     +                  Nrow_ef,Work(I020),Nrow_ef,One,Work(Ioff),
     +                  Ncol_jk)
         Endif

C H(jk,il) + H(jk,il)^(t) 

         Call Transp(Work(Ioff),Work(I030),Ncol_il,Ncol_jk)
         Call Daxpy(Ncol_il*Ncol_jk,One,Work(I030),1,Work(Ioff),1)
         Ioff = Ioff + Ncol_il*Ncol_jk
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("OOOO-C  :",Work(I000),Nsize_t)
#endif
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 

