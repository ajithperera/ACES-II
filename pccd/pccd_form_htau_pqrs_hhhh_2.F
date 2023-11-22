      Subroutine Pccd_form_htau_pqrs_hhhh_2(Work,Maxcor,List_v1,List_v2,
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

      Nsize_t = Idsymsz(Irrepx,14,14)
      I000 = Ione
      I010 = I000 + Nsize_t
      Call Dzero(Work(I000),Nsize_t)

C <le|jf>*G(ke,if) + <if|ke>*G(jf,le) -> H(lj,ki)+H(ki,lj)
C <el|jf>*G(ek,if) + <if|ek>*G(jf,el) -> H(lj,ki)+H(ki,lj). The second contributions
C are identical to the first, hence fator two scalling. 

      Do Icomp = 1, 2

C <le|jf> is stored as fl,ej (fl,ej)^t -> ej,fl -> <ef|jl>

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

      Joff = I020
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         If (Icomp .Eq. 1) Then
            Nrow = Irpdpd(Irrepl,9)
            Ncol = Irpdpd(Irrepr,9)
         Elseif (Icomp .Eq. 2) Then
            Nrow = Irpdpd(Irrepl,11)
            Ncol = Irpdpd(Irrepr,11)
         Endif
         Iend = I030 + Nrow*Ncol
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hhhh_2",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I030),1,Ncol,1,Irrepr,List_v)
         Call Transp(Work(I030),Work(Joff),Ncol,Nrow)
         Joff = Joff + Nrow*Ncol
      Enddo

      Iend = I030 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
      If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_pppp_2",
     +                                    Iend,Maxcor)

      Call Sstgen(Work(I020),Work(I010),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +            Vrt(1,Ispin),Pop(1,Ispin),Work(I030),Irrepx,"1324")

C G(ke,if) stored as (fi,ek) G(fi,ek)^t -> G(ek,fi) -> G(ef,ki)

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

C <ef|jl>*G(ef,ki) -> H(jl,ki)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_ki = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ki,Irrepx)
         Irrep_jl = Dirprd(Irrep_ef,Irrepx)
         
         Ncol_ki = Irpdpd(Irrep_ki,14)
         Ncol_jl = Irpdpd(Irrep_jl,14)
         Nrow_ef = Irpdpd(Irrep_ef,15)

         Icheck = Min(Ncol_ki,Ncol_jl,Nrow_ef)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_jl,Ncol_ki,Nrow_ef,Twom,Work(Joff),
     +                  Nrow_ef,Work(Koff),Nrow_ef,One,Work(Ioff),
     +                  Ncol_jl)
         Endif 
         Call Transp(Work(Ioff),Work(I030),Ncol_ki,Ncol_jl)
         Call Daxpy(Ncol_ki*Ncol_jl,One,Work(I030),1,Work(Ioff),1)

#ifdef _DEBUG_LVLM
         Call checksum("OOOO    :",Work(Ioff),Ncol_ki*Ncol_jl)
#endif
         Ioff = Ioff + Ncol_ki*Ncol_jl
         Joff = Joff + Nrow_ef*Ncol_jl
         Koff = Koff + Nrow_ef*Ncol_ki
      Enddo 
      Enddo 

#ifdef _DEBUG_LVL0
      call checksum("OOOO-C  :",Work(I000),Nsize_t)
#endif

C <el|jf>*G(ek,if) + <if|ek>*G(if,el) -> H(lj,ki)+H(ki,lj)

C I(ej,fl)->I(ef,jl)

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

C G(ek,if) -> G(fk,ei) -> G(fe,ki) -> G(ef,ki)

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

C I(ef,jl)*G(ef,ki) -> H(jl,ki)

      Ioff = I000
      Joff = I010
      Koff = I020
      Do Irrep_ki = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_ki,Irrepx)
         Irrep_jl = Dirprd(Irrep_ef,Irrepx)

         Ncol_ki = Irpdpd(Irrep_ki,14)
         Ncol_jl = Irpdpd(Irrep_jl,14)
         Nrow_ef = Irpdpd(Irrep_ef,15)

         Icheck = Min(Ncol_ki,Ncol_jl,Nrow_ef)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_jl,Ncol_ki,Nrow_ef,Twom,Work(Joff),
     +                  Nrow_ef,Work(Koff),Nrow_ef,One,Work(Ioff),
     +                  Ncol_jl)
         Endif
         Call Transp(Work(Ioff),Work(I030),Ncol_ki,Ncol_jl)
         Call Daxpy(Ncol_ki*Ncol_jl,One,Work(I030),1,Work(Ioff),1)

#ifdef _DEBUG_LVLM
         Call checksum("OOOO    :",Work(Ioff),Ncol_ki*Ncol_jl)
#endif
         Ioff = Ioff + Ncol_ki*Ncol_jl
         Joff = Joff + Nrow_ef*Ncol_jl
         Koff = Koff + Nrow_ef*Ncol_ki
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("OOOO-C  :",Work(I000),Nsize_t)
#endif

      I020 = I010 + Nsize_t
      Call Getall(Work(I010),Nsize_t,Irrepx,List_h)
      Call Daxpy(Nsize_t,One,Work(I010),1,Work(I000),1)
      Call Putall(Work(I000),Nsize_t,Irrepx,List_h)

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
#ifdef _DEBUG_LVLM
         Call checksum("OOOO    :",Work(I000),Nrow_ij*Ncol_kl)
#endif
         Call Assym2a(Irrep_ij,Pop(1,Ispin),Nrow_ij,Ncol_kl,
     +                Work(I000),Work(Iend),Work(Iend))
         Call Transp(Work(I000),Work(I010),Ncol_kl,Nrow_ij)
         Call Assym2a(Irrep_kl,Pop(1,Ispin),Ncol_kl,Nrow_ij,
     +                Work(I010),Work(Iend),Work(Iend))
         Call Transp(Work(I010),Work(I000),Nrow_ij,Ncol_kl)
#ifdef _DEBUG_LVL0
         Call checksum("OOOO-F  :",Work(I000),Nrow_ij*Ncol_kl)
#endif
         Call Putlst(Work(I000),1,Ncol_kl,1,Irrep_kl,List_h)
      Enddo

      Return
      End 


