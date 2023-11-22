      Subroutine pccd_form_htau_2d_hf_pp_2(Hvv,Dpq,Dhf,Dcc,Work,Maxcor,
     +                                     Nocc,Nvrt,Nbas,List_v,Fact,
     +                                     Spin)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,F,Bn,An,Fn
      Character*4 Spin

      Dimension Work(Maxcor)
      Dimension Hvv(Nvrt*Nvrt)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
      Dimension Dcc(Nbas,Nbas)
      Dimension Ioffo(8),Ioffv(8)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione,Two,Quart,Half/1.0D0,-1.0D0,0.0D0,
     +                                        1,2.0D0,0.250D0,
     +                                        0.50D0/
      Ioffo(1) = 0
      Ioffv(1) = Nocco(1)
      Irrepx   = Ione
      Ispin    = Ione

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo

#ifdef _NO_NRML_ORD
C W(an,fm)(t)*D(fm,nb) -> Hvv(a,b)
C Dhf(n,m)*Dco(f,b) is D(fm,nb)

      Ioff = Ione
      Do Irrep_nb = 1, Nirrep
         Irrep_fm = Dirprd(Irrep_nb,Irrepx)

         Nnb = Irpdpd(Irrep_nb,18)
         Nfm = Irpdpd(Irrep_fm,11)
         Iend = Ioff + Nnb*Nfm
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_pp_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nnb*Nfm)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_n = Dirprd(Irrep_b,Irrep_nb)
            Do B = 1, Vrt(Irrep_b,Ispin)
               Do N = 1, Pop(Irrep_n,Ispin)
                  Do Irrep_m = 1, Nirrep
                     Irrep_f = Dirprd(Irrep_m,Irrep_fm)
                     Do M = 1, Pop(Irrep_M,Ispin)
                        Do F  = 1, Vrt(Irrep_f,Ispin)
                           Bn = B + Ioffv(Irrep_b)
                           Nn = N + Ioffo(Irrep_n)
                           Mn = M + Ioffo(Irrep_m)
                           Fn = F + Ioffv(Irrep_f)
C                           If (Spin .Eq. "ABAB") Then
                           Work(Ioff+Icount)=Dhf(Nn,Mn)*Dpq(Fn,Bn) 
C For pCCD this is zero
C                           Elseif (Spin .Eq. "AAAA") Then
C                           Work(Ioff+Icount)=Dhf(Nn,Mn)*Dpq(Fn,Bn) 
C                           Endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built  Dhf(n,m)*Dco(f,b) = D(fm,nb)"
          Call output(Work(Ioff),1,Nfm,1,Nnb,Nfm,Nnb,1)
#endif
          Ioff = Ioff + Nnb*Nfm 
      Enddo

C Get W(fm,an)
      Kbgn = Ione
      Do Irrep_an = 1, Nirrep
         Irrep_fm = Dirprd(Irrep_an,Irrepx)
         Irrep_nb = Dirprd(Irrep_fm,Irrepx)

         Nrow_fm = Irpdpd(Irrep_fm,11)
         Ncol_an = Irpdpd(Irrep_an,11)
         Ncol_nb = Irpdpd(Irrep_nb,18)

         I010 = Iend 
         I020 = I010 + Nrow_fm*Ncol_an
         I030 = I020 + Max(Nrow_fm,Ncol_an)
         I040 = I030 + Max(Nrow_fm,Ncol_an)
         Iend = I040 + Max(Nrow_fm,Ncol_an)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph_2",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_an,1,Irrep_an,List_v)

C W(fm,an) -> W(fm,na)

         Call Symtr1(Irrep_an,Vrt(1,Ispin),Pop(1,Ispin),Nrow_fm,
     +                Work(I010),Work(I020),Work(I030),Work(I040))

         Ioff = Ione
         Joff = I010
         Koff = Kbgn

#ifdef _DEBUG_LVLM
	  Write(6,"(a)") "Request  Dhf(n,m)*Dco(f,b) = D(fm,nb)"
          Call output(Work(koff),1,Nrow_fm,1,Ncol_nb,Nrow_fm,
     +                Ncol_nb,1)
#endif

C W(fm,na)*D(fm,nb) -> Hvv(a,b)
         Do Irrep_b = 1, Nirrep 
            Irrep_n  = Dirprd(Irrep_b,Irrep_nb)
            Irrep_a  = Dirprd(Irrep_n,Irrep_an)

            Nn = Pop(Irrep_n,Ispin)
            Na = Vrt(Irrep_a,Ispin)
            Nb = Vrt(Irrep_b,Ispin)

            Nsum = Nrow_fm*Nn
            Nrow = Na
            Ncol = Nb

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvv(Ioff),
     +                     Nrow)
            Endif

            Ioff = Ioff + Na*Nb
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Nb
         Enddo
         Kbgn = Kbgn + Nrow_fm*Ncol_nb
      Enddo
#endif 
#ifdef _DEBUG_LVL0
      call checksum("Htau_vv :",Hvv,Nvrt*Nvrt)
      call output(Hvv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
#endif
      Return
      End

 

