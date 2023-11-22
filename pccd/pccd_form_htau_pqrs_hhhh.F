      Subroutine Pccd_form_htau_pqrs_hhhh(Work,Maxcor,List_v,
     +                                    List_g,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Four/1.0D0,-1.0D0,0.0D0,1,
     +                                  2.0D0,4.0D0/

      Irrepx = Ione
      Ispin  = Ione 

      Nsize_s = Idsymsz(Irrepx,14,14)
      Nsize_t = Idsymsz(Irrepx,14,14)

C Reorder HHHH lists from 12|12 11|22

      I000  = Ione 
      I010  = I000 + Nsize_t
      I020  = I010 + Nsize_s
      Iend  = I020 + Nocco(Ispin)*Nocco(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                   Iend,Maxcor)
      Call Getall(Work(I010),Nsize_s,Irrepx,List_v)
      Call SStgen(Work(I010),Work(I000),Nsize_s,Pop(1,Ispin),
     +                Pop(1,Ispin),Pop(1,Ispin),Pop(1,Ispin),
     +                Work(I020),Irrepx,"1324")
      I020  = I010 + Nsize_t
      I030  = I020 + Nsize_s
      Iend  = I030 + Nocco(Ispin)*Nocco(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)
      Call Getall(Work(I020),Nsize_s,Irrepx,List_g)
      Call SStgen(Work(I020),Work(I010),Nsize_s,Pop(1,Ispin),
     +                Pop(1,Ispin),Pop(1,Ispin),Pop(1,Ispin),
     +                Work(I030),Irrepx,"1324")
      Ioff = I000
      Joff = I010

C <mn|jk>*G(mn,il) (four identical contributions)

      Do Irrep_il = 1, Nirrep
         Irrep_mn = Dirprd(Irrep_il,Irrepx)
         Irrep_jk = Dirprd(Irrep_mn,Irrepx)

          Nrow_mn = Irpdpd(Irrep_mn,14)
          Ncol_jk = Irpdpd(Irrep_jk,14)
          Ncol_il = Irpdpd(Irrep_il,14)

          I030  = I020 + Max(Nrow_mn,Ncol_jk)
          Iend  = I030 + Max(Nrow_mn,Ncol_jk)
          If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                       Iend,Maxcor)
          Call Spinad3(Irrep_mn,Pop(1,1),Nrow_mn,Ncol_jk,Work(Ioff),
     +                 Work(I020),Work(I030))

          Iend  = I020 + Ncol_jk*Ncol_il
          If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)
          Icheck = Min(Nrow_mn,Ncol_jk,Ncol_il)

          If (Icheck .Gt. 0) Then
             Call Dgemm("T","N",Ncol_jk,Ncol_il,Nrow_mn,Four,Work(Ioff),
     +                   Nrow_mn,Work(Joff),Nrow_mn,Dnull,Work(I020),
     +                   Ncol_jk)
          Endif

#ifdef _DEBUG_LVL0
         Call checksum("OOOO    :",Work(I020),Ncol_jk*Ncol_il)
#endif
         Call Putlst(Work(I020),1,Ncol_il,1,Irrep_il,List_h)
         Ioff = Ioff + Nrow_mn*Ncol_jk
         Joff = Joff + Nrow_mn*Ncol_il
      Enddo

#ifdef _DEBUG_LVLM
      write(6,*)
#endif

C <jl||mn>*G(ik,mn) (two identical contributions)

      Do Irrep_mn = 1, Nirrep
         Irrep_jl = Dirprd(Irrep_mn,Irrepx) 
         Irrep_ik = Dirprd(Irrep_mn,Irrepx)

         Nrow_jl = Irpdpd(Irrep_jl,14)
         Nrow_ik = Irpdpd(Irrep_ik,14)
         Ncol_mn = Irpdpd(Irrep_mn,14)

         I000 = Ione
         I010 = I000 + Nrow_jl*Ncol_mn
         I020 = I010 + Max(Ncol_mn,Nrow_jl)
         Iend = I020 + Max(Ncol_mn,Nrow_jl)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_mn,1,Irrep_mn,List_v)

         Call Spinad3(Irrep_jl,Pop(1,1),Nrow_jl,Ncol_mn,Work(I000),
     +                Work(I010),Work(I020))

         I020 = I010 + Nrow_ik*Ncol_mn
         Iend = I020 + Nrow_jl*Nrow_ik
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_mn,2,Irrep_mn,List_g)

         Icheck = Min(Nrow_jl,Nrow_ik,Ncol_mn)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_jl,Nrow_ik,Ncol_mn,Two,Work(I000),
     +                  Nrow_jl,Work(I010),Nrow_ik,Dnull,Work(I020),
     +                  Nrow_jl)
         Endif 

         Iend = I000 + Nrow_jl*Nrow_ik
         If (Iend .Gt. I020) Then
            Write(6,"(2a)") " Breach of array bound integrity in",
     +                      " pccd_form_htau_pqrs_hhhh.F"
            Call Errex
         Endif

#ifdef _DEBUG_LVL0
         Call checksum("OOOO    :",Work(I020),Nrow_jl*Nrow_ik)
#endif
         Call Getlst(Work(I000),1,Nrow_ik,1,Irrep_ik,List_h)
         Call Daxpy(Nrow_jl*Nrow_ik,One,Work(I000),1,Work(I020),1)
         Call Putlst(Work(I020),1,Nrow_ik,1,Irrep_ik,List_h)
      Enddo 

C Read the H(ij,kl)

#ifdef _DEBUG_LVLM
      write(6,*)
#endif
      Do Irrep_kl = 1, Nirrep
         Irrep_ij = Dirprd(Irrep_kl,Irrepx)

         Nrow_ij = Irpdpd(Irrep_ij,14)
         Ncol_kl = Irpdpd(Irrep_kl,14)

         I000 = Ione
         I010 = I000 + Nrow_ij*Ncol_kl 
         Iend = I010 + Nrow_ij*Ncol_kl
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I000),1,Ncol_kl,1,Irrep_kl,List_h)
#ifdef _DEBUG_LVL0
         Call checksum("OOOO    :",Work(I000),Nrow_ij*Ncol_kl)
#endif
         Call Assym2a(Irrep_ij,Pop(1,Ispin),Nrow_ij,Ncol_kl,
     +                Work(I000),Work(Iend),Work(Iend))
         Call Transp(Work(I000),Work(I010),Ncol_kl,Nrow_ij)
         Call Assym2a(Irrep_kl,Pop(1,Ispin),Ncol_kl,Nrow_ij,
     +                Work(I010),Work(Iend),Work(Iend))
         Call Transp(Work(I010),Work(I000),Nrow_ij,Ncol_kl)
#ifdef _DEBUG_LVL0
         Call checksum("OOOO    :",Work(I000),Nrow_ij*Ncol_kl)
#endif
         Call Putlst(Work(I000),1,Ncol_kl,1,Irrep_kl,List_h)
      Enddo

      Return
      End 

