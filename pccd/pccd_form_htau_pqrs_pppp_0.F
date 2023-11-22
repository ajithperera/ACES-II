      Subroutine Pccd_form_htau_pqrs_pppp_0(Work,Maxcor,List_v1,List_v2,
     +                                      List_g1,List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Four,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                       2.0D0,4.0D0,-2.0D0/
      write(6,*) 

      Irrepx = Ione
      Ispin  = Ione
    
      Nsize_t = Idsymsz(Irrepx,15,15)

      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C <ef|ac>*G(ef,bd)->H(ac,bd)
C <bd|ef>*G(ac,ef)->H(bd,ac)

      Ioff = Ione
      Do Irrep_bd = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_bd,Irrepx)
         Irrep_ac = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,13)
         Ncol_ac = Irpdpd(Irrep_ac,13)
         Ncol_bd = Irpdpd(Irrep_bd,13)

         I020 = I010 + Nrow_ef*Ncol_ac
         I030 = I020 + Max(Ncol_ac,Nrow_ef)
         Iend = I030 + Max(Ncol_ac,Nrow_ef)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_ac,1,Irrep_ac,List_v1)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_ac,Work(I010),
     +                Work(I020),Work(I030))

         I030 = I020 + Nrow_ef*Ncol_bd
         Iend = I030 + Ncol_ac*Ncol_bd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I020),1,Ncol_bd,2,Irrep_bd,List_g1)

         Icheck = Min(Nrow_ef,Ncol_ac,Ncol_bd)

C <ef|ac>*G(ef,bd)->H(ac,bd)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ac,Ncol_bd,Nrow_ef,One,Work(I010),
     +                  Nrow_ef,Work(I020),Nrow_ef,One,Work(Ioff),
     +                  Ncol_ac)
         Endif

C H(ac,bd)+H(ac,bd)^t

         Call Transp(Work(Ioff),Work(I030),Ncol_bd,Ncol_ac)
         Call Daxpy(Ncol_bd*Ncol_ac,One,Work(I030),1,Work(Ioff),1)
         Ioff = Ioff + Ncol_bd*Ncol_ac
      Enddo

C <bd|mn>*G(ac,mn) + <mn|ac>*G(mn,bd) 

      Ioff = Ione
      Do Irrep_mn = 1, Nirrep
         Irrep_bd = Dirprd(Irrep_mn,Irrepx)
         Irrep_ac = Dirprd(Irrep_mn,Irrepx)

         Ncol_mn = Irpdpd(Irrep_mn,14)
         Nrow_ac = Irpdpd(Irrep_ac,13)
         Nrow_bd = Irpdpd(Irrep_bd,13)

         I020 = I010 + Nrow_bd*Ncol_mn
         I030 = I020 + Max(Nrow_bd,Ncol_mn)
         Iend = I030 + Max(Nrow_bd,Ncol_mn)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_mn,1,Irrep_mn,List_v2)

         Call Spinad3(Irrep_bd,Vrt(1,1),Nrow_bd,Ncol_mn,Work(I010),
     +                Work(I020),Work(I030))

         I030 = I020 + Ncol_mn*Nrow_ac
         Iend = I030 + Nrow_ac*Nrow_bd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_0",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I020),1,Ncol_mn,2,Irrep_mn,List_g2)

         Icheck = Min(Ncol_mn,Nrow_ac,Nrow_bd)

C <bd|mn>*G(ac,mn)->H(bd,ac)

         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bd,Nrow_ac,Ncol_mn,One,Work(I010),
     +                  Nrow_bd,Work(I020),Nrow_ac,One,Work(Ioff),
     +                  Nrow_bd)
         Endif

C H(bd,ac)+H(bd,ac)^t

         Call Transp(Work(Ioff),Work(I030),Nrow_ac,Nrow_bd)
         Call Daxpy(Nrow_bd*Nrow_ac,One,Work(I030),1,Work(Ioff),1)
         Ioff = Ioff + Nrow_bd*Nrow_ac
      Enddo 
#ifdef _DEBUG_LVL0
      call checksum("VVVV-C  :",Work(I000),Nsize_t)
#endif
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

      Return
      End 

