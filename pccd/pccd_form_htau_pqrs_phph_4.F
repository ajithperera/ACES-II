      Subroutine Pccd_form_htau_pqrs_phph_4(Work,Maxcor,List_v1,List_g1,
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
C <be|mi>*G(ae,mj)->H(aj,bi)
C <be|mi> stored as (bm,ei) and G(ae,mj) is stored as (ae,mj)
C The factor of two accounts for <be|im>*G(ae,jm)

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
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_4",
     +                                   Iend,Maxcor)

C I(bm,ei) -> I(bi,em)

      Call Getall(Work(I020),Nsize,Irrepx,List_v1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C G(ae,mj) spin-adapt and then reorder to (aj,em)

      Nsize = Idsymsz(Irrepx,11,11)

      I030 = I020 + Nsize
      I040 = I030 + Nsize 
      Iend = I040 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),
     +                  Nocco(1)*Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_4",
     +                                   Iend,Maxcor)

      Joff = I030
      Do Irrep_mj = 1, Nirrep
         Irrep_ae = Dirprd(Irrep_mj,Irrepx)
         
         Nrow_ae = Irpdpd(Irrep_ae,15)
         Ncol_mj = Irpdpd(Irrep_mj,14)

         I050 = I040 + Max(Nrow_ae,Ncol_mj)
         Iend = I050 + Max(Nrow_ae,Ncol_mj)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph_4",
     +                                   Iend,Maxcor)
     
         Call Getlst(Work(Joff),1,Ncol_mj,1,Irrep_mj,List_g1)
         Call Spinad3(Irrep_ae,Vrt(1,1),Nrow_ae,Ncol_mj,Work(Joff),
     +                Work(I040),Work(I050))
         Joff = Joff + Nrow_ae*Ncol_mj
      Enddo 

C G(be,mj) -> G(bj,em)

      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Vrt(1,Ispin),
     +            Pop(1,Ispin),Pop(1,Ispin),Work(I040),Irrepx,"1423")

C I(ai,em)*G(bj,em)
    
      Ioff = I000
      Joff = I010 
      Koff = I020
      Do Irrep_em = 1, Nirrep
         Irrep_bj = Dirprd(Irrep_em,Irrepx)
         Irrep_ai = Dirprd(Irrep_em,Irrepx)

         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bj = Irpdpd(Irrep_bj,11)
         Nrow_ai = Irpdpd(Irrep_ai,11)

         Icheck = Min(Ncol_em,Nrow_bj,Nrow_ai)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_ai,Nrow_bj,Ncol_em,Twom,Work(Joff),
     +                  Nrow_ai,Work(Koff),Nrow_bj,One,Work(Ioff),
     +                  Nrow_ai)
         Endif
#ifdef _DEBUG_LVLM
         call checksum("VOVO    :",Work(Ioff),Nrow_ai*Nrow_bj)
#endif
         Ioff = Ioff + Nrow_ai*Nrow_bj
         Joff = Joff + Nrow_ai*Ncol_em
         Koff = Koff + Nrow_bj*Ncol_em
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("VOVO-C  :",Work(I000),Nsize_t)
#endif

      Iend = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I010),Nsize_t,Irrepx,List_h)
#ifdef _DEBUG_LVL0
      call checksum("VOVO-F  :",Work(I000),Nsize_t)
#endif

      Return
      End 


