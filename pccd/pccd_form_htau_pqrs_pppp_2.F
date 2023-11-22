      Subroutine Pccd_form_htau_pqrs_pppp_2(Work,Maxcor,List_v1,List_v2,
     +                                      List_v3,List_g1,List_g2,
     +                                      List_g3,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Four,Twom/1.0D0,-1.0D0,0.0D0,1,
     +                                       2.0D0,4.0D0,-2.0D0/
      Irrepx = Ione
      Ispin  = Ione

      Nsize_t = Idsymsz(Irrepx,15,15)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C <bn|cm>*G(an,dm) + <dm|an>*G(cm,bn) -> H(bc,ad)
C <bn|mc>*G(an,md) + <md|am>*G(mc,bn) -> H(bc,ad). The second contributions are identical 
C to the firsti, hence the factor two scalling. 

      Do Icomp = 1, 2

C <bn|cm> is stored as cn,bm -> <cb|nm>

      Ioff = I000
      If (Icomp .eq.1) Then
      Nsize = Idsymsz(Irrepx,9,9)
      List_v = List_v1
      List_g = List_g1
      ElseIf (Icomp .eq.2) Then
      Nsize = Idsymsz(Irrepx,11,11)
      List_v = List_v2
      List_g = List_g2
      Endif 

      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                    Iend,Maxcor)

      Call Getall(Work(I020),Nsize,Irrepx,List_v)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1324")

C G(an,dm) stored as (dm,an) G(dm,an)^t -> G(an,dm) -> G(ad,nm)

      I030  = I020 + Nsize
      I040  = I030 + Nsize

      Joff = I030
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         If (Icomp .Eq. 1) Then
            Nrow = Irpdpd(Irrepl,9)
            Ncol = Irpdpd(Irrepr,9)
         Elseif (Icomp .Eq. 2) Then
            Nrow = Irpdpd(Irrepl,11)
            Ncol = Irpdpd(Irrepr,11)
         Endif 
         Iend = I040 + Nrow*Ncol
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I040),1,Ncol,1,Irrepr,List_g)
         Call Transp(Work(I040),Work(Joff),Ncol,Nrow)
         Joff = Joff + Nrow*Ncol
      Enddo 
      Iend  = I040 + Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                   Iend,Maxcor)
      
      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I040),Irrepx,"1324")

C <bc|nm>*G(ad,nm) -> H(bc,ad)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_nm = 1, Nirrep
         Irrep_ad = Dirprd(Irrep_nm,Irrepx)
         Irrep_bc = Dirprd(Irrep_nm,Irrepx)
         
         Nrow_ad = Irpdpd(Irrep_ad,15)
         Nrow_bc = Irpdpd(Irrep_bc,15)
         Ncol_nm = Irpdpd(Irrep_nm,14)

         Icheck = Min(Nrow_ad,Nrow_bc,Ncol_nm)
         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bc,Nrow_ad,Ncol_nm,Twom,Work(Joff),
     +                  Nrow_bc,Work(Koff),Nrow_ad,One,Work(Ioff),
     +                  Nrow_bc)
         Endif 
         Call Transp(Work(Ioff),Work(I030),Nrow_ad,Nrow_bc)
         Call Daxpy(Nrow_ad*Nrow_bc,One,Work(I030),1,Work(Ioff),1)

#ifdef _DEBUG_LVLM
         Call checksum("VVVV    :",Work(Ioff),Nrow_ad*Nrow_bc)
#endif
         Ioff = Ioff + Nrow_ad*Nrow_bc
         Joff = Joff + Nrow_bc*Ncol_nm
         Koff = Koff + Nrow_ad*Ncol_nm
      Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      call checksum("VVVV-C  :",Work(I000),Nsize_t)
#endif

C <bn|mc>*G(an,md) + <md|an>*G(mc,bn)
C <bn|mc> is ordered as (bm,cn), reorder to (bc,mn)
C G(an,md)=(na,dm)is ordered as (dn,am) reorder to (da,mn)

C I(bm,cn)->I(bc,mn)

      Nsize = Idsymsz(Irrepx,11,11)
      I020 = I010 + Nsize
      I030 = I020 + Nsize
      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                    Iend,Maxcor)
      Call Getall(Work(I020),Nsize,Irrepx,List_v3)
      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1324")

C G(dn,am) -> G(da,nm) -> G(ad,mn)

      I030  = I020 + Nsize
      I040  = I030 + Nsize
      Iend  = I040 + Max(Nvrto(Ispin)*Nvrto(Ispin),
     +                   Nocco(Ispin)*Nocco(Ispin),
     +                   Nocco(Ispin)*Nvrto(Ispin))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                   Iend,Maxcor)
      Call Getall(Work(I030),Nsize,Irrepx,List_g3)
      Call Sstgen(Work(I030),Work(I020),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I040),Irrepx,"1324")

C <bc|nm>*G(da,mn) -> <bc|mn>*G(ad,mn) -> H(bc,ad)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_mn = 1, Nirrep
         Irrep_ad = Dirprd(Irrep_mn,Irrepx)
         Irrep_bc = Dirprd(Irrep_mn,Irrepx)

         Nrow_ad = Irpdpd(Irrep_ad,15)
         Nrow_bc = Irpdpd(Irrep_bc,15)
         Ncol_mn = Irpdpd(Irrep_mn,14)

C <bc|mn>*G(ad,mn) 

         Icheck = Min(Nrow_ad,Nrow_bc,Ncol_mn)
         If (Icheck .Gt. 0) Then
            Call Dgemm("N","T",Nrow_bc,Nrow_ad,Ncol_mn,Twom,Work(Joff),
     +                  Nrow_bc,Work(Koff),Nrow_ad,One,Work(Ioff),
     +                  Nrow_bc)
         Endif

         Call Transp(Work(Ioff),Work(I030),Nrow_ad,Nrow_bc)
         Call Daxpy(Nrow_ad*Nrow_bc,One,Work(I030),1,Work(Ioff),1)
#ifdef _DEBUG_LVLM
         Call checksum("VVVV    :",Work(Ioff),Nrow_ad*Nrow_bc)
#endif
         Ioff = Ioff + Nrow_ad*Nrow_bc
         Joff = Joff + Nrow_bc*Ncol_mn
         Koff = Koff + Nrow_ad*Ncol_mn
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("VVVV-C  :",Work(I000),Nsize_t)
#endif

      I020 = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)
#ifdef _DEBUG_LVL0
      call checksum("VVVV-C  :",Work(I000),Nsize_t)
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
         Call checksum("VVVV-F  :",Work(I000),Nrow_ab*Ncol_cd)
#endif 
         Call Putlst(Work(I000),1,Ncol_cd,1,Irrep_cd,List_h)
      Enddo 

      Return
      End 


