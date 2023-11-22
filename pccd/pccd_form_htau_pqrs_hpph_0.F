      Subroutine Pccd_form_htau_pqrs_hpph_0(Work,Maxcor,List_v1,
     +                                      List_v2,List_g1,List_g2,
     +                                      List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two /1.0D0,-1.0D0,0.0D0,1,2.0D0/

C List_v1=18,List_g=125
C <em|ia>*G(em,bj) -> H(ia,bj)
C <em|ia> is stored as (ei,am) and G(em,bj) is stored as (em,bj) 
C Here factor two account for spin-adapation 

      Irrepx = Ione
      Ispin  = Ione
      Ioff   = Ione

      Nsize_t = Idsymsz(Irrepx,11,11)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C I(ei,am) -> I(em,ai)

      Nsize = Idsymsz(Irrepx,11,11)
      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hpph_0",
     +                                   Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_v1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(em,ai)*G(em,bj) -> H(bi,aj)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_bj = Dirprd(Irrep_em,Irrepx)
         Irrep_ai = Dirprd(Irrep_em,Irrepx)

         Ncol_ai = Irpdpd(Irrep_ai,11)
         Ncol_bj = Irpdpd(Irrep_bj,11)
         Nrow_em = Irpdpd(Irrep_em,11)

         Iend = I020 + Ncol_em*Nrow_bj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hpph_0",
     +                                   Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_bj,1,Irrep_bj,List_g1)
         Icheck = Min(Nrow_em,Ncol_ai,Ncol_bj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ai,Ncol_bj,Nrow_em,Two,Work(Joff),
     +                  Ncol_ai,Work(I020),Nrow_em,One,Work(Ioff),
     +                  Ncol_ai)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("OVVO    :",Work(Ioff),Ncol_bj*Ncol_ai)
#endif
         Ioff = Ioff + Ncol_bj*Ncol_ai
         Joff = Joff + Ncol_ai*Nrow_em
      Enddo

C List_v2=25,List_g2=118
C <me|ia>*G(me,bj) -> H(ia,bj)
C <me|ia> is stored as (am,ie) and G(me,bj) is stored as bj,em

C I(am,ei) -> I(ai,em)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hpph_0",
     +                                   Iend,Maxcor)
      Call Getall(Work(I020),Nsize,Irrepx,List_v2)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(ai,em)*G(bj,em)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_bj = Dirprd(Irrep_em,Irrepx)
         Irrep_ai = Dirprd(Irrep_em,Irrepx)

         Nrow_ai = Irpdpd(Irrep_ai,11)
         Nrow_bj = Irpdpd(Irrep_bj,11)
         Ncol_em = Irpdpd(Irrep_em,11)

         Iend = I020 + Ncol_em*Nrow_bj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hpph_0",
     +                                   Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_g2)
         Icheck = Min(Nrow_ai,Nrow_bj,Ncol_em)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_ai,Nrow_bj,Ncol_em,Two,Work(Joff),
     +                  Ncol_em,Work(I020),Nrow_bj,One,Work(Ioff),
     +                  Nrow_ai)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("OVVO    :",Work(Ioff),Nrow_bj*Nrow_ai)
#endif
         Ioff = Ioff + Nrow_bj*Nrow_ai
         Joff = Joff + Nrow_ai*Ncol_em
      Enddo

      Iend = I000 + Nsize_t
      Maxcor = Maxcor - Iend

C H(ia,bj) + H(ia,bj)^(t)

      Call Pccd_symmetrize_phph(Work(I000),Nsize_t,Work(Iend),Maxcor,
     +                          .False.)
#ifdef _DEBUG_LVL0
      call checksum("OVVO-F  :",Work(I000),Nsize_t)
#endif

      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)
      Return
      End 


