      Subroutine Pccd_form_htau_pqrs_pphh_0(Work,Maxcor,List_v1,List_v2,
     +                                      List_g1,List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione,Two /1.0D0,-1.0D0,0.0D0,1,
     +                              2.0D0/

      Irrepx = Ione
      Ispin  = Ione

C <ef||ab>*G(ef,ij) -> H(ab,ij)

      Do Irrep_ij = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ij,Irrepx)
         Irrep_ab = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,15)
         Ncol_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)

         I000 = Ione
         I010 = I000 + Nrow_ef*Ncol_ab
         I020 = I010 + Max(Ncol_ab,Nrow_ef)
         Iend = I020 + Max(Ncol_ab,Nrow_ef)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ab,1,Irrep_ab,List_v2)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_ab,Work(I000),
     +                Work(I010),Work(I020))

         I020 = I010 + Nrow_ef*Ncol_ij
         Iend = I020 + Ncol_ab*Ncol_ij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_ij,2,Irrep_ij,List_g2)

         Icheck = Min(Ncol_ab,Ncol_ij,Nrow_ef)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ab,Ncol_ij,Nrow_ef,One,Work(I000),
     +                  Nrow_ef,Work(I010),Nrow_ef,Dnull,Work(I020),
     +                  Ncol_ab)
         Endif
         Call Putlst(Work(I020),1,Ncol_ij,1,Irrep_ij,List_h)
      Enddo

C <mn||ab>*G(mn,ij)->H(ab,ij) (<mn|ab> is stored as ab,mn)

      Do Irrep_ij = 1, Nirrep
         Irrep_mn = Dirprd(Irrep_ij,Irrepx) 
         Irrep_ab = Dirprd(Irrep_mn,Irrepx)

         Nrow_mn = Irpdpd(Irrep_mn,14)
         Ncol_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)

         I000 = Ione
         I010 = I000 + Nrow_mn*Ncol_ab
         I020 = I010 + Nrow_mn*Ncol_ab
         I030 = I020 + Max(Nrow_mn,Ncol_ab)
         Iend = I030 + Max(Nrow_mn,Ncol_ab)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Nrow_mn,1,Irrep_mn,List_v1)

         Call Spinad3(Irrep_ab,Vrt(1,1),Ncol_ab,Nrow_mn,Work(I010),
     +                Work(I020),Work(I030))

C <ab|mn> -> <mn|ab>

         Call Transp(Work(I010),Work(I000),Nrow_mn,Ncol_ab)

         I020 = I010 + Nrow_mn*Ncol_ij
         Iend = I020 + Ncol_ab*Ncol_ij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_ij,2,Irrep_ij,List_g1)

C <mn|ab>*G(mn,ij) -> H(ab,ij)
  
         Icheck = Min(Ncol_ab,Ncol_ij,Nrow_mn)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ab,Ncol_ij,Nrow_mn,One,Work(I000),
     +                  Nrow_mn,Work(I010),Nrow_mn,Dnull,Work(I020),
     +                  Ncol_ab)
         Endif 
         Iend = I000 + Ncol_ab*Ncol_ij
         If (Iend .Gt. I020) Then
            Write(6,"(2a)") " Breach of array bound integrity in",
     +                      " pccd_form_htau_pqrs_pphh.F"
            Call Errex
         Endif
  
         Call Getlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_h)
         Call Daxpy(Ncol_ab*Ncol_ij,One,Work(I000),1,Work(I020),1)
         Call Putlst(Work(I020),1,Ncol_ij,1,Irrep_ij,List_h)
      Enddo 

#ifdef _DEBUG_LVL0
      Nsize_t = Idsymsz(Irrepx,15,14)
      Call Getall(work(I000),Nsize_t,Irrepx,List_h)
      Call checksum("VVOO-C  :",Work(I000),Nsize_t)
#endif 
      Return
      End 


