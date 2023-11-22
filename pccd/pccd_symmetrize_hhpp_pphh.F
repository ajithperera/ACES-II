      Subroutine Pccd_symmetrize_hhpp_pphh(Work,Maxcor,List_h1,
     +                                    List_h2)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C List_h1 = H(ab,ij), List_h2=H(ij,ab)

      Irrepx = 1

      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nrow_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)
       
         I000 = Ione
         I010 = I000 + Nrow_ab*Ncol_ij
         I020 = I010 + Nrow_ab*Ncol_ij
         I030 = I020 + Nrow_ab*Ncol_ij
         Iend = I030 + Nrow_ab*Ncol_ij
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_symmetriz_hhpp_pphh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_h1)
         Call Getlst(Work(I010),1,Nrow_ab,1,Irrep_ab,List_h2)
         Call Dcopy(Nrow_ab*Ncol_ij,Work(I010),1,Work(I020),1)
         
         Call Transp(Work(I000),Work(I030),Ncol_ij,Nrow_ab)
         Call Daxpy(Ncol_ij*Nrow_ab,One,Work(I030),1,Work(I010),1)
         Call Putlst(Work(I010),1,Nrow_ab,1,Irrep_ab,List_h2)

         Call Transp(Work(I020),Work(I030),Nrow_ab,Ncol_ij)
         Call Daxpy(Ncol_ij*Nrow_ab,One,Work(I030),1,Work(I000),1)
         Call Putlst(Work(I000),1,Ncol_ij,1,Irrep_ij,List_h1)
        
      Enddo 

#ifdef _DEBUG_LVL0
         Nsize_t = Idsymsz(Irrepx,15,14)
         Call Getall(work(I000),Nsize_t,Irrepx,List_h1)
         Call checksum("VVOO-C  :",Work(I000),Nsize_t)
         Nsize_t = Idsymsz(Irrepx,14,15)
         Call Getall(work(I000),Nsize_t,Irrepx,List_h1)
         Call checksum("OOVV-C  :",Work(I000),Nsize_t)
#endif

      Return
      End
