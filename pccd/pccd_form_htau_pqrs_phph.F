      Subroutine Pccd_form_htau_pqrs_phph(Work,Maxcor,List_v1,
     +                                    List_v2,List_g1,
     +                                    List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C <ke|ib>*G(ke,aj)=> <ek|bi>*G(ek,ja) [<ke|ib> stored as ek,bi and G(ek,ja)
C is stored as G(ek,aj) 

      Irrepx = Ione
      Ispin  = Ione
      Do Irrep_aj = 1, Nirrep
         Irrep_ek = Dirprd(Irrep_aj,Irrepx) 
         Irrep_bi = Dirprd(Irrep_ek,Irrepx)

         Nrow_ek = Irpdpd(Irrep_ek,11)
         Ncol_bi = Irpdpd(Irrep_bi,11)
         Ncol_aj = Irpdpd(Irrep_aj,11)

         I000 = Ione
         I010 = I000 + Nrow_ek*Ncol_bi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_bi,1,Irrep_bi,List_v1)

         I020 = I010 + Nrow_ek*Ncol_aj
         I030 = I020 + Nrow_ek
         I040 = I030 + Nrow_ek
         Iend = I040 + Nrow_ek
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_aj,2,Irrep_aj,List_g1)

C G(ek,aj) -> (ek,ja)

         Call Symtr1(Irrep_aj,Vrt(1,Ispin),Pop(1,Ispin),Nrow_ek,
     +               Work(I010),Work(I020),Work(I030),Work(I040))

         Iend = I020 + Ncol_ib*Ncol_aj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_phph",
     +                                      Iend,Maxcor)
         Icheck = Min(Ncol_ij,Ncol_ab,Nrow_kl)

C <ek|bi>*G(ek,ja)-> H(bi,aj)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_bi,Ncol_aj,Nrow_ek,One,Work(I000),
     +                  Nrow_ek,Work(I010),Nrow_ek,Dnull,Work(I020),
     +                  Ncol_bi)
         Endif 
#ifdef _DEBUG_LVL0
         call checksum("VOVO    :",Work(I020),Ncol_ib*Ncol_aj)
#endif 
         Call Putlst(Work(I020),1,Ncol_aj,1,Irrep_aj,List_h)
      Enddo 

      Return
      End 


