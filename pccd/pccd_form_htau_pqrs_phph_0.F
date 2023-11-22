      Subroutine Pccd_form_htau_pqrs_phph_0(Work,Maxcor,List_v1,List_v2,
     +                                      List_v3,List_g1,List_g2,
     +                                      List_g3,List_h)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two /1.0D0,-1.0D0,0.0D0,1,2.0D0/

C List_v1=23,List_g1=123
C <em|ai>*G(em,bj) -> H(ai,bj)
C <em|ai> is stored as (am,ei) and G(em,bj) is stored as (bj,em)

      Irrepx = Ione
      Ispin  = Ione

      Nsize_t = Idsymsz(Irrepx,11,11)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

      Nsize = Idsymsz(Irrepx,11,11)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                   Iend,Maxcor)

C I(am,ei) -> (ai,em)

      Call Getall(Work(I020),Nsize,Irrepx,List_v1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

      Ioff = I000
      Joff = I010

C I(ai,em)*G(bj,em) -> H(ai,bj)

      Do Irrep_em = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_em,Irrepx) 
         Irrep_bj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_ai = Irpdpd(Irrep_ai,11)
         Nrow_bj = Irpdpd(Irrep_bj,11)

         Iend = I020 + Ncol_em*Nrow_bj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_0",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_g1)

         Icheck = Min(Nrow_ai,Nrow_bj,Ncol_em)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_ai,Nrow_bj,Ncol_em,One,Work(Joff),
     +                  Nrow_ai,Work(I020),Nrow_bj,One,Work(Ioff),
     +                  Nrow_ai)
         Endif 
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_ai*Nrow_bj)
#endif 
         Ioff = Ioff + Nrow_bj*Nrow_ai
         Joff = Joff + Nrow_ai*Ncol_em
      Enddo 
 
C List_v2=25,List_g2=123
C <em|ai>*G(em,bj) -> H(ai,bj)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                   Iend,Maxcor)
C I(am,ei) -> (ai,em)

      Call Getall(Work(I020),Nsize,Irrepx,List_v2)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(ai,em)*G(bj,em) -> H(ai,bj)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_em,Irrepx)
         Irrep_bj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_ai = Irpdpd(Irrep_ai,11)
         Nrow_bj = Irpdpd(Irrep_bj,11)

         Iend = I020 + Ncol_em*Nrow_bj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_0",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_g1)

         Icheck = Min(Nrow_ai,Nrow_bj,Ncol_em)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_ai,Nrow_bj,Ncol_em,One,Work(Joff),
     +                  Nrow_ai,Work(I020),Nrow_bj,One,Work(Ioff),
     +                  Nrow_ai)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_ai*Nrow_bj)
#endif
         Ioff = Ioff + Nrow_bj*Nrow_ai
         Joff = Joff + Nrow_ai*Ncol_em
      Enddo

C List_v3=18,List_g3=118
C <me|ai>*G(me,bj) + <bj|me>*G(ai,me) 
C <me|ai> is ordered as (am,ei) and G(me,bj) is oredered as (bm,ej)

C I(am,ei) -> I(ai,em)

       I020 = I010 + Nsize
       I030 = I020 + Nsize
       Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                   Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_v3)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C G(bm,ej) -> G(bj,em)

      I030 = I020 + Nsize
      I040 = I030 + Nsize
      Iend = I040 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                 Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                    Iend,Maxcor)

      Call Getall(Work(I030),Nsize,Irrepx,List_v1)
      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I040),Irrepx,"1432")

C  I(ai,em)*G(bj,em)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_em = 1, Nirrep
         Irrep_ai = Dirprd(Irrep_em,Irrepx)
         Irrep_bj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_ai = Irpdpd(Irrep_ai,11)
         Nrow_bj = Irpdpd(Irrep_bj,11)

         Icheck = Min(Ncol_em,Nrow_ai,Nrow_bj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_ai,Nrow_bj,Ncol_em,One,Work(Joff),
     +                  Nrow_ai,Work(Koff),Nrow_bj,One,Work(Ioff),
     +                  Nrow_ai)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_ai*Nrow_bj)
#endif
         Ioff = Ioff + Nrow_bj*Nrow_ai
         Joff = Joff + Nrow_ai*Ncol_em
         Koff = Koff + Nrow_bj*Ncol_em
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("VOVO-C  :",Work(I000),Nsize_t)
#endif

      Iend   = I000 + Nsize_t
      Maxcor = Maxcor - Iend 

C H(ai,bj) + H(ai,bj)^t

      Call Pccd_symmetrize_phph(Work(I000),Nsize_t,Work(Iend),Maxcor,
     +                          .False.)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h) 

#ifdef _DEBUG_LVL0
      call checksum("VOVO-F  :",Work(I000),Nsize_t)
#endif
      Return
      End 


