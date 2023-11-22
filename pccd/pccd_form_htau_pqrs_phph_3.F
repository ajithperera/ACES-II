      Subroutine Pccd_form_htau_pqrs_phph_3(Work,Maxcor,List_v1,List_v2,
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
C <bm|ei>*G(am,ej)=<mb|ie>*G(ma,je)-> H(ai,bj)
C <mb|ie> is ordered (em,bi) and G(ma,je) ordered as (ej,am)

      Irrepx = Ione
      Ispin  = Ione
      Nsize_t = Idsymsz(Irrepx,11,11)

      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C G(ej,am) -> G(em,aj)

       Nsize = Idsymsz(Irrepx,11,11)
       I020 = I010 + Nsize
       I030 = I020 + Nsize
       Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_3",
     +                                   Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_g1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(em,bi)*G(em,aj)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx) 
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Nrow_em = Irpdpd(Irrep_em,11)
         Ncol_bi = Irpdpd(Irrep_bi,11)
         Ncol_aj = Irpdpd(Irrep_aj,11)

         Iend = I020 + Ncol_bi*Nrow_em
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_3",
     +                                   Iend,Maxcor)

         Call Getlst(Work(I020),1,Ncol_bi,1,Irrep_bi,List_v1)

         Icheck = Min(Nrow_em,Ncol_bi,Ncol_aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_bi,Ncol_aj,Nrow_em,Onem,Work(I020),
     +                  Nrow_em,Work(Joff),Nrow_em,One,Work(Ioff),
     +                  Ncol_aj)
         Endif 
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Ncol_aj*Ncol_bi)
#endif 
         Ioff = Ioff + Ncol_aj*Ncol_bi
         Joff = Joff + Ncol_aj*Nrow_em
      Enddo 
 
C List_v2=25,List_g2=125
C <bm|ei>*G(am,ej)=<mb|ie>*G(ma,je)-> H(ai,bj)
C <mb|ie> is ordered (em,bi) and G(ma,je) ordered as (ej,am)

C G(ej,am) -> G(em,aj)

       I020 = I010 + Nsize
       I030 = I020 + Nsize
       Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
       If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_3",
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

         Iend = I020 + Ncol_bi*Nrow_em
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_3",
     +                                   Iend,Maxcor)

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
         Joff = Joff + Nrow_bi*Ncol_em
      Enddo

C List_v3=18,List_g3=118
C <bm|ie>*G(am,je) <bm|ie> is stored as (bi,em) and G(am,je) is 
C stored as (aj,em) 
C I(bi,em)*G(aj,em) -> H(bi,aj)

      Ioff = I000
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx)
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)
         Nrow_aj = Irpdpd(Irrep_aj,11)

         I020 = I010 + Ncol_em*Nrow_aj
         Iend = I020 + Ncol_em*Nrow_bi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_3",
     +                                   Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_em,1,Irrep_em,List_v3)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_g3)

         Icheck = Min(Ncol_em,Nrow_bi,Nrow_aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bi,Nrow_aj,Ncol_em,Onem,Work(I010),
     +                  Nrow_bi,Work(I020),Nrow_aj,One,Work(Ioff),
     +                  Nrow_bi)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_aj*Nrow_bi)
#endif
         Ioff = Ioff + Nrow_aj*Nrow_bi
      Enddo
#ifdef _DEBUG_LVL0
      call checksum("VOVO-C  :",Work(I000),Nsize_t)
#endif

      Iend = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)

#ifdef _DEBUG_LVL0
      call checksum("VOVO-F  :",Work(I000),Nsize_t)
#endif
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h) 

      Return
      End 


