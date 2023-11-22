      Subroutine Pccd_form_htau_pqrs_hhpp(Work,Maxcor,List_v1,
     +                                    List_v2,List_g1,
     +                                    List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C <kl||ij>*G(kl,ab)  G is stores as ab,kl

      Irrepx = 1
      Ispin  = 1
      Do Irrep_ab = 1, Nirrep
         Irrep_kl = Dirprd(Irrep_ab,Irrepx) 
         Irrep_ij = Dirprd(Irrep_kl,Irrepx)

         Nrow_kl = Irpdpd(Irrep_kl,14)
         Ncol_ij = Irpdpd(Irrep_ij,14)
         Ncol_ab = Irpdpd(Irrep_ab,15)

         I000 = Ione
         I010 = I000 + Nrow_kl*Ncol_ij
         I020 = I010 + Max(Ncol_ij,Nrow_kl)
         Iend = I020 + Max(Ncol_ij,Nrow_kl)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_v1)

         Call Spinad3(Irrep_kl,Pop(1,1),Nrow_kl,Ncol_ij,Work(I000),
     +                Work(I010),Work(I020))
         I020 = I010 + Nrow_kl*Ncol_ab
         Iend = I020 + Ncol_ij*Ncol_ab
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Nrow_kl,2,Irrep_kl,List_g1)

         Icheck = Min(Ncol_ij,Ncol_ab,Nrow_kl)

C <ij|kl>*G(ab,kl)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Ncol_ij,Ncol_ab,Nrow_kl,One,Work(I000),
     +                  Nrow_kl,Work(I010),Ncol_ab,Dnull,Work(I020),
     +                  Ncol_ij)
         Endif 
         Call Putlst(Work(I020),1,Ncol_ab,1,Irrep_ab,List_h)
      Enddo 

C <ef||ij>*G(ef,ab)

      Irrepx = 1
      Do Irrep_ab = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ab,Irrepx)
         Irrep_ij = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)
         Ncol_ab = Irpdpd(Irrep_ab,15)

         I000 = Ione
         I010 = I000 + Nrow_ef*Ncol_ij
         I020 = I010 + Ncol_ij
         Iend = I020 + Ncol_ij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_v2)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_ij,Work(I000),
     +                Work(I010),Work(I020))

         I020 = I010 + Nrow_ef*Ncol_ab
         Iend = I020 + Ncol_ij*Ncol_ab
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhpp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_ab,2,Irrep_ab,List_g2)

         Icheck = Min(Ncol_ij,Ncol_ab,Nrow_ef)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ij,Ncol_ab,Nrow_ef,One,Work(I000),
     +                  Nrow_ef,Work(I010),Nrow_ef,Dnull,Work(I020),
     +                  Ncol_ij)
         Endif
         Iend = I000 + Ncol_ab*Ncol_ij
         If (Iend .Gt. I020) Then
            Write(6,"(2a)") " Breach of array bound integrity in",
     +                      " pccd_form_htau_pqrs_pphh.F"
            Call Errex
         Endif

         Call Getlst(Work(I000),1,Ncol_ab,1,Irrep_ab,List_h)
         Call Daxpy(Ncol_ab*Ncol_ij,One,Work(I000),1,Work(I020),1)
         Call Putlst(Work(I020),1,Ncol_ab,1,Irrep_ab,List_h)
      Enddo

C Antisymmetrize H(ij,ab)

      Do Irrep_ab = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_ab,Irrepx)

         Nrow_ij = Irpdpd(Irrep_ij,14)
         Ncol_ab = Irpdpd(Irrep_ab,13)

         I000 = Ione
         I010 = I000 + Nrow_ij*Ncol_ab
         Iend = I010 + Nrow_ij*Ncol_ab
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pphh",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I000),1,Ncol_ab,1,Irrep_ab,List_h)
         Call Assym2a(Irrep_ij,Pop(1,Ispin),Nrow_ij,Ncol_ab,
     +                Work(I000),Work(Iend),Work(Iend))
         Call Transp(Work(I000),Work(I010),Ncol_ab,Nrow_ij)
         Call Assym2a(Irrep_ab,Vrt(1,Ispin),Ncol_ab,Nrow_ij,
     +                Work(I010),Work(Iend),Work(Iend))
         Call Transp(Work(I010),Work(I000),Nrow_ij,Ncol_ab)

#ifdef _DEBUG_LVL0
         Call checksum("OOVV    :",Work(I000),Nrow_ij*Ncol_ab)
#endif
         Call Putlst(Work(I000),1,Ncol_ab,1,Irrep_ab,List_h)
      Enddo

      Return
      End 


