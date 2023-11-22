      Subroutine Pccd_form_htau_2d_phph_pp(Hvv_pq,Hvv_qp,Work,Maxcor,
     +                                     Nvrt,Nbas,List_v,List_g) 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C V(mb,ne)*Gamma(ma,ne)
C V(bm,en)*G(am,en)=V(en,bm)G(en,am)

      Irrepx = Ione
      Ispin  = Ione
  
      Do Irrep_am = 1, Nirrep
         Irrep_en = Dirprd(Irrep_am,Irrepx)
         Irrep_bm = Dirprd(Irrep_en,Irrepx)

         Nrow_en = Irpdpd(Irrep_en,11)
         Ncol_bm = Irpdpd(Irrep_bm,11)
         Ncol_am = Irpdpd(Irrep_am,11)

         I000 = Ione
         I010 = I000 + Nrow_en*Ncol_bm
         I020 = I010 + Nrow_en*Ncol_am
         I030 = I020 + Nrow_en
         I040 = I030 + Nrow_en
         Iend = I040 + Nrow_en

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_vv",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_bm,1,Irrep_bm,List_v)
         Call Getlst(Work(I010),1,Ncol_am,2,Irrep_am,List_g)

         Call Symtr1(Irrep_bm,Vrt(1,Ispin),Pop(1,Ispin),Nrow_en,
     +               Work(I000),Work(I020),Work(I030),Work(I040))
         Call Symtr1(Irrep_am,Vrt(1,Ispin),Pop(1,Ispin),Nrow_en,
     +               Work(I010),Work(I020),Work(I030),Work(I040))

C V(en,mb)^(t)*G(en,ma) -> S(b,a)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_a = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_a,Irrep_am)
            Irrep_b  = Dirprd(Irrep_m,Irrep_bm)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nm = Pop(irrep_m,1)
            
            Nsum = Nrow_en*Nm
            Nrow = Nb
            Ncol = Na
            Icheck = Min(Nsum,Nrow,Ncol)
            
            If (Icheck .Ne. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hvv_qp(Ioff),Nrow)
            Endif 

            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Nb
            Koff = Koff + Nsum*Na
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_vv :",Hvv_qp,Nvrt,"VV","D")
CSSS      call checksum("Htau_vv :",Hvv_qp,Nvrt*Nvrt)
#endif

#ifdef _EXPLICITLY_PERMUTE
C V(me,fa)*Gamma(me,fb)

      Do Irrep_mb = 1, Nirrep
         Irrep_ne = Dirprd(Irrep_mb,Irrepx)
         Irrep_ma = Dirprd(Irrep_ne,Irrepx)

         Nrow_ma = Irpdpd(Irrep_ma,16)
         Ncol_ne = Irpdpd(Irrep_ne,16)
         Nrow_ne = Irpdpd(Irrep_ne,16)
         Ncol_mb = Irpdpd(Irrep_mb,16)

         I000 = Ione
         I010 = I000 + Nrow_ma*Ncol_ne
         Iend = I010 + Nrow_ne*Ncol_mb

         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hphp_vv",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ne,1,Irrep_ne,List_v)
         Call Getlst(Work(I010),1,Ncol_mb,2,Irrep_mb,List_g)

C  V(ma,ne)*Gamma(ne,mb) -> S(a,b)

         Ioff = Ione
         Joff = I000
         Koff = I010
         Do Irrep_b = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_b,Irrep_mb)
            Irrep_a  = Dirprd(Irrep_m,Irrep_ma)

            Na = Vrt(irrep_a,1)
            Nb = Vrt(irrep_b,1)
            Nm = Pop(irrep_m,1)
            
            Nsum = Ncol_ne*Nm
            Nrow = Na
            Ncol = Nb
            Icheck = Min(Nsum,Nrow,Ncol)

            If (Icheck .Gt. 0) Then
                Call Dgemm("T","N",Nrow,Ncol,Nsum,One,Work(Joff),
     +               Nsum,Work(Koff),Nsum,One,Hvv_pq(Ioff),Nrow)
            Endif

            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Nb
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      call pccd_check_htau("Htau_vv :",Hvv_pq,Nvrt,"VV","F")
#endif
#endif

      Return
      End

