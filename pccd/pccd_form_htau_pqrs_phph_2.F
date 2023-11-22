      Subroutine Pccd_form_htau_pqrs_phph_2(Work,Maxcor,List_v1,List_g1,
     +                                      List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Twom /1.0D0,-1.0D0,0.0D0,1,
     +                              -2.0D0/

C List_v1=18,List_g1=116
C <mj|ae>*G(mi,be)-><ae|mj>*G(be,mi)->H(aj,bi)
C <ae|mj> stored as (am,ej) and G(be,mi) is stored as (be,mi)
C The factor of two accounts for <jm|ae>*G(im,be) 

      Irrepx = 1
      Ispin  = 1 
      Nsize_t = Idsymsz(Irrepx,11,11)

      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

      Nsize = Idsymsz(Irrepx,11,11)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_0",
     +                                   Iend,Maxcor)

C I(am,ej)-> I(aj,em)

      Call Getall(Work(I020),Nsize,Irrepx,List_v1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C G(be,mi) spin-adapt and then reorder to (bi,em)

      I030 = I020 + Nsize
      I040 = I030 + Nsize 
      Iend = I040 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),
     +                  Nocco(1)*Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_0",
     +                                   Iend,Maxcor)

      Joff = I030
      Do Irrep_mi = 1, Nirrep
         Irrep_be = Dirprd(Irrep_mi,Irrepx)
         
         Nrow_be = Irpdpd(Irrep_be,15)
         Ncol_mi = Irpdpd(Irrep_mi,14)

         I050 = I040 + Max(Nrow_be,Ncol_mi)
         Iend = I050 + Max(Nrow_be,Ncol_mi)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_0",
     +                                      Iend,Maxcor)
         Call Getlst(Work(Joff),1,Ncol_mi,1,Irrep_mi,List_g1)
         Call Spinad3(Irrep_be,Vrt(1,1),Nrow_be,Ncol_mi,Work(Joff),
     +                Work(I040),Work(I050))
         Joff = Joff + Nrow_be*Ncol_mi
      Enddo 

C G(be,mi) -> G(bi,em)

      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Vrt(1,Ispin),
     +            Pop(1,Ispin),Pop(1,Ispin),Work(I040),Irrepx,"1423")

C I(aj,em)*G(bi,em)
    
      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_em = 1, Nirrep
         Irrep_bi = Dirprd(Irrep_em,Irrepx)
         Irrep_aj = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)
         Nrow_aj = Irpdpd(Irrep_aj,11)

         Icheck = Min(Ncol_em,Nrow_bi,Nrow_aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_aj,Nrow_bi,Ncol_em,Twom,Work(Joff),
     +                  Nrow_aj,Work(Koff),Nrow_bi,One,Work(Ioff),
     +                  Nrow_aj)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_aj*Nrow_bi)
#endif
         Ioff = Ioff + Nrow_aj*Nrow_bi
         Joff = Joff + Nrow_aj*Ncol_em
         Koff = Koff + Nrow_bi*Ncol_em
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("VOVO-C  :",Work(I000),Nsize_t)
#endif

C H(aj,bi)-> H(ai,bj) -> H(ai,bj) + H(ai,bj)^t

      Iend   = I010 + Nsize_t
      Maxcor = Maxcor - Iend
      Call  Pccd_symmetrize_phph(Work(I000),Nsize_t,Work(Iend),Maxcor,
     +                           .True.)

      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)

#ifdef _DEBUG_LVL0
      call checksum("VOVO-F  :",Work(I000),Nsize_t)
#endif
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 


