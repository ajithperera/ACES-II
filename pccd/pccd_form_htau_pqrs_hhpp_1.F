      Subroutine Pccd_form_htau_pqrs_hhpp_1(Work,Maxcor,List_v1,
     +                                      List_v2,List_g1,List_h)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                  2.0D0,-2.0D0/

C <bm|ie>*G(jm,ae)+<mb|ie>*G(mj.ae)->H(bi,aj) <bm|ie>=<be,im> is stored as (bi,em) and 
C G(jm,ae) is stored as (ae,jm).Reorder to aj,em

      Irrepx = Ione
      Ispin  = Ione
      Nsize_t  = Idsymsz(Irrepx,11,11)

      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

      Nsize = Idsymsz(Irrepx,15,14)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_1",
     +                                   Iend,Maxcor)

C G(jm,ae) stored as ae,jm and reorder to aj,em

      Call Getall(Work(I020),Nsize,Irrepx,List_g1)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Vrt(1,Ispin),
     +            Pop(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1324")

C I(bi,em)*G(aj,em) -> H(bi,aj)

      Joff = I010
      Ioff = Ione
      Do Irrep_em = 1, Nirrep
         Irrep_aj = Dirprd(Irrep_em,Irrepx) 
         Irrep_bi = Dirprd(Irrep_em,Irrepx)

         Nrow_aj = Irpdpd(Irrep_aj,11)
         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)

         Iend = I020 + Nrow_bi*Ncol_em
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh_1",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_v1)

         Icheck = Min(Nrow_aj,Ncol_em,Nrow_bi)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bi,Nrow_aj,Ncol_em,Two,Work(I020),
     +                  Nrow_bi,Work(Joff),Nrow_aj,One,Work(Ioff),
     +                  Nrow_bi)
         Endif 
         Ioff = Ioff + Nrow_bi*Nrow_aj
         Joff = Joff + Nrow_aj*Ncol_em
      Enddo 

C <mb|ie>*G(mj,ae) (<mb|ie> ordered as (em,bi) and G(mj,ae) is 
C ordered as (ae,mj). This is the direct contribution: spin-adapt the
C integral or G.

      Nsize = Idsymsz(Irrepx,15,14)

      I020 = I010 + Nsize 
      I030 = I020 + Nsize

      Koff = I020 
      Do Irrep_mj = 1, Nirrep
         Irrep_ae = Dirprd(irrep_mj,Irrepx)

         Nrow_ae = Irpdpd(Irrep_ae,15)
         Ncol_mj = Irpdpd(Irrep_mj,14)

         I040 = I030 + Max(Nrow_ae,Ncol_mj)
         Iend = I040 + Max(Nrow_ae,Ncol_mj)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh_1",
     +                                      Iend,Maxcor)

	 Call Getlst(Work(Koff),1,Ncol_mj,2,Irrep_mj,List_g1)
         Call Spinad3(Irrep_ae,Vrt(1,1),Nrow_ae,Ncol_mj,Work(Koff),
     +                Work(I030),Work(I040))
         Koff = Koff + Nrow_ae*Ncol_mj
      Enddo 

C G(ae,mj) (spin-adpated) -> G(aj,em)

      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),
     +                  Nocco(1)*Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh_1",
     +                                   Iend,Maxcor)

C G(ae,mj) stored as ae,mj and reorder to aj,em

      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Vrt(1,Ispin),
     +            Pop(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1423")

C <mb|ie>*G(mj,ae), Reordered to I(bi,em)*G(aj,em)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
         Irrep_aj = Dirprd(Irrep_em,Irrepx)
         Irrep_bi = Dirprd(Irrep_em,Irrepx)

         Nrow_aj = Irpdpd(Irrep_aj,11)
         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)

         I030 = I020 + Nrow_bi*Ncol_em
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh_1",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_em,1,Irrep_em,List_v2)

C I(bi,em)*G(aj,em) -> H(bi,aj)

         Icheck = Min(Nrow_aj,Ncol_em,Nrow_bi)
         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bi,Nrow_aj,Ncol_em,Twom,Work(I020),
     +                  Nrow_bi,Work(Joff),Nrow_aj,One,Work(Ioff),
     +                  Nrow_bi)
         Endif
         Ioff = Ioff + Nrow_bi*Nrow_aj
         Joff = Joff + Nrow_aj*Ncol_em
      Enddo

C H(bi,aj) -> H(ab,ij)

      I020 = I010 + Nsize_t
      Iend = I020 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh_1",
     +                                   Iend,Maxcor)
      Call Sstgen(Work(I000),Work(I010),Nsize_t,Vrt(1,Ispin),
     +            Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),Work(I020),
     +            Irrepx,"1324")

      Call Pccd_transp_pphh(Work(I010),Work(I000),Nsize_t)
C <bm|ie>*G(jm,ae) + <mb|ie>*G(mj,ae)
C <am|je>*G(ie,bm) + <am|ej>*G(ie,mb) Note that the second-contribution is
C identical to the first and therefore scale the firsr contribution by
C factor two to account for the full term.

      Call Dscal(Nsize_t,Two,Work(I000),1) 
#ifdef _DEBUG_LVL0
         Call checksum("OOVV-C  :",Work(I000),Nsize_t)
#endif
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 


