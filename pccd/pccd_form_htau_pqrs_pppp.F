      Subroutine Pccd_form_htau_pqrs_pppp(Work,Maxcor,List_v,
     +                                    List_g,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Four/1.0D0,-1.0D0,0.0D0,1,2.0D0,
     +                              4.0D0/
      write(6,*) 
      Irrepx = Ione
      Ispin  = Ione

C <eb|fc>*G(ea,fd)+<be|cf>*G(ae,cd)+<eb|cf>*G(ea,df)+<be|fc>*G(ae,dc)

      Nsize_s = Idsymsz(Irrepx,15,15)
      Nsize_t = Idsymsz(Irrepx,15,15)

C Reorder PPPP lists from 12|12 to  11|22

      I000  = Ione
      I010  = I000 + Nsize_t
      I020  = I010 + Nsize_s
      Iend  = I020 + Nvrto(Ispin)*Nvrto(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                   Iend,Maxcor)
      Call Getall(Work(I010),Nsize_s,Irrepx,List_v)
      Call Sstgen(Work(I010),Work(I000),Nsize_s,Vrt(1,Ispin),
     +            Vrt(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),Work(I020),
     +            Irrepx,"1324")

      I020  = I010 + Nsize_t
      I030  = I020 + Nsize_s
      Iend  = I030 + Nvrto(Ispin)*Nvrto(Ispin)
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)
      Call Getall(Work(I030),Nsize_s,Irrepx,List_g)
      Call Sstgen(Work(I030),Work(I010),Nsize_s,Vrt(1,Ispin),
     +            Vrt(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),Work(I030),
     +            Irrepx,"1324")

      Ioff = I000
      Joff = I010
      Do Irrep_ad = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ad,Irrepx)
         Irrep_bc = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,13)
         Ncol_ad = Irpdpd(Irrep_ad,13)
         Ncol_bc = Irpdpd(Irrep_bc,13)

         I030  = I020 + Nrow_ef
         Iend  = I030 + Nrow_ef
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_bc,Work(Ioff),
     +                Work(I020),Work(I030))

         I030 = I020 + Ncol_ad*Ncol_bc
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)

C <ef|bc>*G(ef,ad) -> H(bc,ad)

         Icheck = Min(Nrow_ef,Ncol_bc,Ncol_ad)
         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_bc,Ncol_ad,Nrow_ef,Two,Work(Ioff),
     +                  Nrow_ef,Work(Joff),Nrow_ef,Dnull,Work(I020),
     +                  Ncol_bc)
         Endif 

#ifdef _DEBUG_LVLM
         Call checksum("VVVV    :",Work(I020),Ncol_ad*Ncol_bc)
#endif
         Call Putlst(Work(I020),1,Ncol_ad,1,Irrep_ad,List_h)
         Ioff = Ioff + Nrow_ef*Ncol_bc
         Joff = Joff + Nrow_ef*Ncol_ad
      Enddo 

#ifdef _DEBUG_LVLM
      write(6,*) 
#endif 
C <ef|ac>*G(ef,bd)+<bd|ef>*G(ac|ef)

      Do Irrep_bd = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_bd,Irrepx)
         Irrep_ac = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,13)
         Ncol_ac = Irpdpd(Irrep_ac,13)
         Ncol_bd = Irpdpd(Irrep_bd,13)

         I000 = Ione
         I010 = I000 + Nrow_ef*Ncol_ac
         I020 = I010 + Ncol_ac
         Iend = I020 + Ncol_ac
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ac,1,Irrep_ac,List_v)

         Call Spinad3(Irrep_ef,Vrt(1,1),Nrow_ef,Ncol_ac,Work(I000),
     +                Work(I010),Work(I020))

         I020 = I010 + Nrow_ef*Ncol_bd
         Iend = I020 + Ncol_ac*Ncol_bd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_bd,2,Irrep_bd,List_g)

         Icheck = Min(Nrow_ef,Ncol_ac,Ncol_bd)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ac,Ncol_bd,Nrow_ef,Two,Work(I000),
     +                  Nrow_ef,Work(I010),Nrow_ef,Dnull,Work(I020),
     +                  Ncol_ac)
         Endif

         Iend = I000 + Ncol_ac*Ncol_bd
         If (Iend .Gt. I020) Then
            Write(6,"(2a)") " Breach of array bound integrity in",
     +                      " pccd_form_htau_pqrs_pppp.F"
            Call Errex
         Endif
         
#ifdef _DEBUG_LVLM
         Call checksum("VVVV    :",Work(I020),Ncol_bd*Ncol_ac)
#endif
         Call Getlst(Work(I000),1,Ncol_bd,1,Irrep_bd,List_h)
         Call Daxpy(Ncol_bd*Ncol_ac,One,Work(I000),1,Work(I020),1)
         Call Putlst(Work(I020),1,Ncol_bd,1,Irrep_bd,List_h)
      Enddo

C Read the H(ab,cd)  

#ifdef _DEBUG_LVLM
      write(6,*) 
#endif 
      Do Irrep_cd = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_cd,Irrepx)

         Nrow_ab = Irpdpd(Irrep_ab,13)
         Ncol_cd = Irpdpd(Irrep_cd,13)

         I000 = Ione
         I010 = I000 + Nrow_ab*Ncol_cd
         Iend = I010 + Nrow_ab*Ncol_cd
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I000),1,Ncol_cd,1,Irrep_cd,List_h)
         Call Assym2a(Irrep_ab,Vrt(1,Ispin),Nrow_ab,Ncol_cd,
     +                Work(I000),Work(Iend),Work(Iend))
         Call Transp(Work(I000),Work(I010),Ncol_cd,Nrow_ab)
         Call Assym2a(Irrep_cd,Vrt(1,Ispin),Ncol_cd,Nrow_ab,
     +                Work(I010),Work(Iend),Work(Iend))
         Call Transp(Work(I010),Work(I000),Nrow_ab,Ncol_cd)
#ifdef _DEBUG_LVL0
         Call checksum("VVVV    :",Work(I000),Ncol_ab*Ncol_cd)
#endif 
         Call Putlst(Work(I000),1,Ncol_cd,1,Irrep_cd,List_h)
      Enddo 

      Return
      End 


