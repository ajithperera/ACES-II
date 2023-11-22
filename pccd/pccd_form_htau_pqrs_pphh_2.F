      Subroutine Pccd_form_htau_pqrs_pphh_2(Work,Maxcor,List_v1,
     +                                      List_v2,List_v3,List_g1,
     +                                      List_g2,List_g3,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                  2.0D0,-2.0D0/

C <im|be>*G(am,je)+<im|eb>*G(am,ej)  (This is the exchange term,List_v1=18,List_g1=118)
C <im|be>=<be|im> stored as (bi,em) G(am,je)=G(je,am) is stored as 
C (aj,em) 

      Irrepx = Ione
      Ispin  = Ione
      Nsize_t  = Idsymsz(Irrepx,15,14)

      I000 = Ione
      I010 = I000 + Nsize_t
      Iend = I010 + Nsize_t
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_2",
     +                                   Iend,Maxcor)
      Call Dzero(Work(I000),Nsize_t)

C <bi|em>*G(aj,em) -> H(bi,aj)

      Ioff   = I000
      Do Irrep_em = 1, Nirrep
         Irrep_aj = Dirprd(Irrep_em,Irrepx) 
         Irrep_bi = Dirprd(Irrep_em,Irrepx)

         Nrow_aj = Irpdpd(Irrep_aj,11)
         Ncol_em = Irpdpd(Irrep_em,11)
         Nrow_bi = Irpdpd(Irrep_bi,11)

         I020 = I010 + Nrow_aj*Ncol_em
         Iend = I020 + Nrow_bi*Ncol_em
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_2",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_em,1,Irrep_em,List_v1)
         Call Getlst(Work(I020),1,Ncol_em,2,Irrep_em,List_g1)

         Icheck = Min(Nrow_aj,Ncol_em,Nrow_bi)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bi,Nrow_aj,Ncol_em,Two,Work(I010),
     +                  Nrow_bi,Work(I020),Nrow_aj,ONe,Work(Ioff),
     +                  Nrow_bi)
         Endif 
         Ioff = Ioff + Nrow_bi*Nrow_aj
      Enddo 
 
C <im|eb>*G(am,ej) (This is the direct contribution <im|eb>=<eb|im>,list_v2=19
C  is stored as (ei,bm) and G(am,ej)=(ma,je),list_g2=123 is stored ad ej,am)

      Nsize = Idsymsz(Irrepx,11,11)
      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nvrto(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_2",
     +                                   Iend,Maxcor)
C (ei,bm) -> (em,bi)

      Call Getall(Work(I020),Nsize,Irrepx,List_v2)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C G(ej,am) reorder to (em,aj))

      I030 = I020 + Nsize
      I040 = I030 + Nsize
      Iend = I040 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nvrto(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_0",
     +                                   Iend,Maxcor)
      Call Getall(Work(I030),Nsize,Irrepx,List_v2)
      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(em,bi)*G(em,aj) -> H(bi,aj)

      Ioff = I000
      Joff = I010 
      Koff = I020 
      Do Irrep_em = 1, Nirrep
        Irrep_bi = Dirprd(Irrep_em,Irrepx)
        Irrep_aj = Dirprd(Irrep_em,Irrepx)

        Nrow_bi = Irpdpd(Irrep_bi,11)
        Ncol_em = Irpdpd(Irrep_em,11)
        Nrow_aj = Irpdpd(Irrep_aj,11)

        Icheck = Min(Nrow_aj,Ncol_em,Nrow_bi)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Nrow_bi,Nrow_aj,Ncol_em,Twom,Work(Joff),
     +                  Ncol_em,Work(Koff),Ncol_em,One,Work(Ioff),
     +                  Nrow_bi)
         Endif
         Ioff = Ioff + Nrow_bi*Nrow_aj
         Joff = Joff + Nrow_aj*Ncol_em
         Koff = Koff + Nrow_bi*Ncol_em
      Enddo

C <im|eb>*G(am,ej) (This is the direct contribution <im|eb>=<eb|im>,list_v2=21
C  is stored as (em,bi) and G(am,ej)=(ma,je),list_g2=125 is stored ad ej,am)

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nvrto(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_2",
     +                                   Iend,Maxcor)

C G(ej,am) reorder to (em,aj))

      Call Getall(Work(I020),Nsize,Irrepx,List_g3)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1432")

C I(em,bi)*G(em,aj) -> H(bi,aj)

      Ioff = I000
      Joff = I010
      Do Irrep_em = 1, Nirrep
        Irrep_bi = Dirprd(Irrep_em,Irrepx)
        Irrep_aj = Dirprd(Irrep_em,Irrepx)

        Nrow_bi = Irpdpd(Irrep_bi,11)
        Ncol_em = Irpdpd(Irrep_em,11)
        Nrow_aj = Irpdpd(Irrep_aj,11)

        Iend = I020 + Nrow_aj*Ncol_em
        If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_2",
     +                                     Iend,Maxcor)
        Call Getlst(Work(I020),1,Ncol_em,2,Irrep_em,List_v3)

        Icheck = Min(Nrow_aj,Ncol_em,Nrow_bi)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Nrow_bi,Nrow_aj,Ncol_em,Twom,Work(I020),
     +                  Ncol_em,Work(Joff),Ncol_em,One,Work(Ioff),
     +                  Nrow_bi)
         Endif
         Ioff = Ioff + Nrow_bi*Nrow_aj
         Joff = Joff + Nrow_bi*Ncol_em
      Enddo

C H(aj,bi) -> H(ab,ij) 

      I020 = I010 + Nsize_t
      Iend = I020 + Max(Nocco(1)*Nvrto(1),Nvrto(1)*Nvrto(1),
     +                  Nocco(1)*Nvrto(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp_0",
     +                                   Iend,Maxcor)
      Call Sstgen(Work(I000),Work(I010),Nsize_t,Vrt(1,Ispin),
     +            Pop(1,Ispin),Vrt(1,Ispin),Pop(1,Ispin),Work(I020),
     +            Irrepx,"1324")

C <im|be>*G(am,je) + <im|eb>*G(am,ej) 
C <jm|ae>*G(bm,ie) + <mj|ae>*G(mb,ie). Note that the second-contribution is 
C identical to the first and therefore scale the firsr contribution by
C factor two to account for the full term.

      Call Dscal(Nsize_t,Two,Work(I010),1)
#ifdef _DEBUG_LVL0
         Call checksum("VVOO-C  :",Work(I010),Nsize_t)
#endif
      Call Getall(Work(I000),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I000),1,Work(I010),1)
      Call Putall(Work(I010),Nsize_t,Irrepx,List_h)

C Antisymmetrize H(ab,ij)

      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nrow_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)

         I000 = Ione
         I010 = I000 + Nrow_ab*Ncol_ij
         Iend = I010 + Nrow_ab*Ncol_ij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_h)
         Call Assym2a(Irrep_ab,Vrt(1,Ispin),Nrow_ab,Ncol_ij,
     +                Work(I000),Work(Iend),Work(Iend))
         Call Transp(Work(I000),Work(I010),Ncol_ij,Nrow_ab)
         Call Assym2a(Irrep_ij,Pop(1,Ispin),Ncol_ij,Nrow_ab,
     +                Work(I010),Work(Iend),Work(Iend))
         Call Transp(Work(I010),Work(I000),Nrow_ab,Ncol_ij)
#ifdef _DEBUG_LVL0
         Call checksum("VVOO-F  :",Work(I000),Nrow_ab*Ncol_ij)
#endif
         Call Putlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_h)
      Enddo

      Return
      End 


