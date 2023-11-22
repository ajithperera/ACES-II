      Subroutine Pccd_form_htau_pqrs_hhhh_1(Work,Maxcor,List_v,
     +                                      List_g,List_h)
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

      Nsize = Idsymsz(Irrepx,14,14)

C <lm|jn>*G(km,in) + <ml|jn>*G(mk,in) -> H(lj,ki)
C <in|km>*G(jn,lm) + <in|mk>*G(jn,ml) -> H(ki,lj)

C First spin-adapt <lm|jn>, then <lm|jn> -> <lj|mn> 
      I020  = I010 + Nsize
      I030  = I020 + Nsize
      Iend  = I030 + Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh-1",
     +                                   Iend,Maxcor)
      Joff = I020 
      Do Irrep_r = 1, Nirrep
         Irrep_l = Dirprd(Irrep_r,Irrepx)
         Nrow = Irpdpd(Irrep_l,14)
         Ncol = Irpdpd(Irrep_r,14)
         I040 = I030 + Max(Nrow,Ncol)
         Iend = I040 + Max(Nrow,Ncol)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh-1",
     +                                   Iend,Maxcor)
         Call Getlst(Work(Joff),1,Ncol,1,Irrep_r,List_v)
         Call Spinad3(Irrep_l,Pop(1,Ispin),Nrow,Ncol,Work(Joff),
     +                Work(I030),Work(I040))
         Joff = Joff + Nrow*Ncol
      Enddo
         
C G(km,in) -> G(ki,mn)

      Call SStgen(Work(I020),Work(I010),Nsize,Pop(1,Ispin),
     +                Pop(1,Ispin),Pop(1,Ispin),Pop(1,Ispin),
     +                Work(I020),Irrepx,"1324")

      I030  = I020 + Nsize
      I040  = I030 + Nsize
      Iend  = I040 +  Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)
      Call Getall(Work(I030),Nsize,Irrepx,List_g)
      Call SStgen(Work(I030),Work(I020),Nsize,Pop(1,Ispin),
     +                Pop(1,Ispin),Pop(1,Ispin),Pop(1,Ispin),
     +                Work(I040),Irrepx,"1324")

C <lj|mn>*G(ki,mn) -> H(lj,ki)

      Ioff = I000
      Joff = I010
      Koff = I020 
      Do Irrep_mn = 1, Nirrep
         Irrep_lj = Dirprd(Irrep_mn,Irrepx)
         Irrep_ki = Dirprd(Irrep_mn,Irrepx)

          Ncol_mn = Irpdpd(Irrep_mn,14)
          Nrow_lj = Irpdpd(Irrep_lj,14)
          Nrow_ki = Irpdpd(Irrep_ki,14)

          Iend  = I030 + Nrow_lj*Nrow_ki

          Icheck = Min(Ncol_mn,Nrow_lj,Nrow_ki)
          If (Icheck .Gt. 0) Then
             Call Dgemm("N","T",Nrow_lj,Nrow_ki,Ncol_mn,Twom,Work(Joff),
     +                   Nrow_lj,Work(Koff),Nrow_ki,One,Work(Ioff),
     +                   Nrow_lj)
          Endif

          Call Transp(Work(Ioff),Work(I030),Nrow_ki,Nrow_lj)
          Call Daxpy(Nrow_ki*Nrow_lj,One,Work(I030),1,Work(Ioff),1)

#ifdef _DEBUG_LVLM
         Call checksum("OOOO    :",Work(Ioff),Nrow_ki*Nrow_lj)
#endif
         Ioff = Ioff + Nrow_ki*Nrow_lj
         Joff = Joff + Nrow_lj*Ncol_mn
         Koff = Koff + Nrow_ki*Ncol_mn
      Enddo
#ifdef _DEBUG_LVL0
      call checksum("OOOO-C  :",Work(I000),Nsize_t)
#endif
      I020 = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 

