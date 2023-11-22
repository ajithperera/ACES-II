      Subroutine Pccd_form_htau_pqrs_phph_1(Work,Maxcor,List_v1,List_v2,
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
C <ej|am>*G(ei,bm)=<je|ma>G(ie,mb) -> H(ai,bj)
C <je|ma> is ordered (aj,em) and G(ie,mb) ordered as bm,ei

      Irrepx = Ione
      Ispin  = Ione
      Nsize_t = Idsymsz(Irrepx,11,11)

      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C G(bm,ei) -> G(bi,em)

       Nsize = Idsymsz(Irrepx,11,11)
       I020 = I010 + Nsize
       I030 = I020 + Nsize
       Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                    Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_g1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(aj,em)*G(bi,em)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx) 
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)
         Nrow_aj = Irpdpd(Irrep_aj,11)

         Iend  = I020 + Ncol_em*Nrow_aj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                    Iend,Maxcor)
         Icheck = Min(Ncol_em,Nrow_bi,Nrow_aj)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_v1)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_aj,Nrow_bi,Ncol_em,Onem,Work(I020),
     +                  Nrow_aj,Work(Joff),Nrow_bi,One,Work(Ioff),
     +                  Nrow_aj)
         Endif 
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_aj*Nrow_bi)
#endif 
         Ioff = Ioff + Nrow_aj*Nrow_bi
         Joff = Joff + Nrow_aj*Ncol_em
      Enddo 
 
C List_v2=25,List_g2=123
C <ej|am>*G(ei,bm)=<je|am>G(ie,mb) -> H(ai,bj)
C <je|am> is ordered (am,ej) and G(ie,mb) ordered as bm,ei

C G(bm,ei) -> G(bi,em)

       I020 = I010 + Nsize
       I030 = I020 + Nsize
       Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                   Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_g2)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(aj,em)*G(bi,em)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx)
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)
         Nrow_aj = Irpdpd(Irrep_aj,11)

         Iend = I020 +  Ncol_em*Nrow_aj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                    Iend,Maxcor)

         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_v2)
         Icheck = Min(Ncol_em,Nrow_bi,Nrow_aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_aj,Nrow_bi,Ncol_em,Onem,Work(I020),
     +                  Nrow_aj,Work(Joff),Nrow_bi,One,Work(Ioff),
     +                  Nrow_aj)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_aj*Nrow_bi)
#endif
         Ioff = Ioff + Nrow_aj*Nrow_bi
         Joff = Joff + Nrow_aj*Ncol_em
      Enddo

C List_v3=18, List_g3=118
C <je|am>*G(ie,bm) <je|am> is stored as (aj,em) and G(ie,bm) is stored as
C (bi,em)
C I(aj,em)*G(bi,em) -> H(aj,bi)

      Ioff = I000
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx)
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)
         Nrow_aj = Irpdpd(Irrep_aj,11)

         I020 = I010 + Ncol_em*Nrow_aj
         Iend = I020 + Ncol_em*Nrow_bi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_1",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_em,1,Irrep_em,List_v3)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_g3)

         Icheck = Min(Ncol_em,Nrow_bi,Nrow_aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_aj,Nrow_bi,Ncol_em,Onem,Work(I010),
     +                  Nrow_aj,Work(I020),Nrow_bi,One,Work(Ioff),
     +                  Nrow_aj)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_aj*Nrow_bi)
#endif
         Ioff = Ioff + Nrow_aj*Nrow_bi
      Enddo
#ifdef _DEBUG_LVL0
      call checksum("VOVO-C  :",Work(I000),Nsize_t)
#endif

C H(aj,bi)-> H(ai,bj) -> H(ai,bj) + H(ai,bj)^t

      Iend   = I010 + Nsize_t
      Maxcor = Maxcor - Iend 
      Call Pccd_symmetrize_phph(Work(I000),Nsize_t,Work(Iend),Maxcor,
     +                           .True.)

      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

#ifdef _DEBUG_LVL0
      call checksum("VOVO-F  :",Work(I000),Nsize_t)
#endif

      Return
      End 

