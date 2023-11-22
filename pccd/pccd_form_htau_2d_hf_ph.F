      Subroutine pccd_form_htau_2d_hf_ph(Hvo,Hov,Dpq,Dhf,Work,Maxcor,
     +                                   Nocc,Nvrt,Nbas,List_v1,
     +                                   List_v2,Fact)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,Bn,An

      Dimension Work(Maxcor)
      Dimension Hvo(Nocc*Nvrt)
      Dimension Hov(Nocc*Nvrt)
      Dimension Dpq(Nbas,Nbas)
      Dimension Dhf(Nbas,Nbas)
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

      Do Irrep =2, Nirrep
         Ioffo(Irrep)=Ioffo(Irrep-1)+Pop(irrep-1,1)
         Ioffv(Irrep)=Ioffv(Irrep-1)+Vrt(irrep-1,1)
      Enddo

      Call Dzero(Dhf,Nbas*Nbas)
      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         Do I = 1, N
            Dhf(Ioffo(Irrep)+I,Ioffo(Irrep)+I) =  Two
         Enddo
      Enddo

      Ispin  = Ione
      Irrepx = Ione

C W(ck,da)*D(ck,di) -> Hvo(a,i)
C Dco(d,c)*Dhf(k,i) is D(ck,di) 

      Ioff = Ione
      Do Irrep_di = 1, Nirrep
         Irrep_ck = Dirprd(Irrep_di,Irrepx)

         Ndi = Irpdpd(Irrep_di,11)
         Nck = Irpdpd(Irrep_ck,11)
         Iend = Ioff + Ndi*Nck
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Ndi*Nck)
         Icount = Inull
         Do Irrep_i = 1, Nirrep
            Irrep_d = Dirprd(Irrep_i,Irrep_di)
            Do I = 1, Pop(Irrep_i,Ispin)
               Do D = 1, Vrt(Irrep_d,Ispin)
                  Do Irrep_k = 1, Nirrep
                     Irrep_c = Dirprd(Irrep_k,Irrep_ck)
                     Do K = 1, Pop(Irrep_k,Ispin)
                        Do C = 1, Vrt(Irrep_c,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Dn = D + Ioffv(Irrep_d)
                           Kn = K + Ioffo(Irrep_k)
                           Cn = C + Ioffv(Irrep_c)
                           Work(Ioff+Icount)=Dpq(Dn,Cn)*Dhf(Kn,In)
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dco(c,d)*Dhf(k,i) = D(ck,di)"
          Call output(Work(Ioff),1,Nck,1,Ndi,Nck,Ndi,1)
#endif 
          Ioff = Ioff + Ndi*Nck
      Enddo 

C Get W(da,ck) and transpose to form W(ck,da)

      Kbgn = Ione
      Do Irrep_ck = 1, Nirrep
         Irrep_da = Dirprd(Irrep_ck,Irrepx)

         Nrow_da = Irpdpd(Irrep_da,13)
         Ncol_ck = Irpdpd(Irrep_ck,11)

         I010 = Iend
         I020 = I010 + Nrow_da*Ncol_ck
         I030 = I020 + Nrow_da*Ncol_ck
         I040 = I030 + Max(Nrow_da,Ncol_ck)
         Iend = I040 + Max(Nrow_da,Ncol_ck)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_ck,1,Irrep_ck,List_v1)

         Call Spinad3(Irrep_da,Vrt(1,1),Nrow_da,Ncol_ck,Work(I020),
     +                Work(I030),Work(I040))
         Call Transp(Work(I020),Work(I010),Ncol_ck,Nrow_da)

C W(ck,da)*D(ck,di) -> Hvo(a,i)

         Ioff = Ione
         Joff = I010
         Koff = Kbgn
         Irrep_di = Dirprd(Irrep_ck,Irrepx)
         Ncol_di  = Irpdpd(Irrep_di,11)

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Request Dco(c,d)*Dhf(k,i) = D(ck,di)"
         Call output(Work(Koff),1,Ncol_ck,1,Ncol_di,Ncol_ck,
     +               Ncol_di,1)
#endif 
         Do irrep_i = 1, Nirrep
            Irrep_d  = Dirprd(Irrep_i,Irrep_di)
            Irrep_a  = Dirprd(Irrep_d,Irrep_da)

            Nd = Vrt(irrep_d,Ispin)
            Na = Vrt(irrep_a,Ispin)
            Ni = Pop(irrep_i,Ispin)

            Nsum = Ncol_ck*Nd
            Nrow = Na
            Ncol = Ni

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvo(Ioff),
     +                     Nrow)
            Endif

            Ioff = Ioff + Ni*Na
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Ni
          Enddo
          Kbgn = Kbgn + Ncol_di*Ncol_ck 
      Enddo 
#ifdef _DEBUG_LVL0
      call checksum("Htau_vo :",Hvo,Nocc*Nvrt)
#endif 

C W(nk,ma)*D(nk,mi) -> Hvo(a,i)
C Dco(m,n)*Dhf(k,i) is D(nk,mi)

      Ioff = Ione
      Do Irrep_mi = 1, Nirrep
         Irrep_nk = Dirprd(Irrep_mi,Irrepx)

         Nnk = Irpdpd(Irrep_nk,14)
         Nmi = Irpdpd(Irrep_mi,14)
         Iend = Ioff + Nnk*Nmi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nnk*Nmi)
         Icount = Inull
         Do Irrep_i = 1, Nirrep
            Irrep_m = Dirprd(Irrep_i,Irrep_mi)
            Do I = 1, Pop(Irrep_i,Ispin)
               Do M = 1, Pop(Irrep_m,Ispin)
                  Do Irrep_k = 1, Nirrep
                     Irrep_n = Dirprd(Irrep_k,Irrep_nk)
                     Do K = 1, Pop(Irrep_k,Ispin)
                        Do N  = 1, Pop(Irrep_n,Ispin)
                           In = I + Ioffo(Irrep_i)
                           Mn = M + Ioffo(Irrep_m)
                           Kn = K + Ioffo(Irrep_k)
                           Nn = N + Ioffo(Irrep_n)
                           Work(Ioff+Icount)=Dpq(Mn,Nn)*Dhf(Kn,In)-
     +                                       Dhf(Mn,Nn)*Dhf(Kn,In)
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dco(m,n)*Dhf(k,i) = D(nk,mi)"
          Call output(Work(Ioff),1,Nnk,1,Nmi,Nnk,Nmi,1)
#endif
          Ioff = Ioff + Nnk*Nmi
      Enddo

C Get W(nk,ma) 

      Kbgn = Ione
      Do Irrep_ma = 1, Nirrep
         Irrep_nk = Dirprd(Irrep_ma,Irrepx)

         Nrow_nk = Irpdpd(Irrep_nk,14)
         Ncol_ma = Irpdpd(Irrep_ma,18)

         I010 = Iend
         I020 = I010 + Nrow_nk*Ncol_ma
         I030 = I020 + Max(Nrow_nk,Ncol_ma)
         Iend = I030 + Max(Nrow_nk,Ncol_ma)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_ma,1,Irrep_ma,List_v2)

         Call Spinad3(Irrep_nk,Pop(1,1),Nrow_nk,Ncol_ma,Work(I010),
     +                Work(I020),Work(I030))

C W(nk,ma)*D(nk,mi) -> Hvo(a,i)

         Ioff = Ione
         Joff = I010
         Koff = Kbgn
         Irrep_mi = Dirprd(Irrep_nk,Irrepx)
         Ncol_mi  = Irpdpd(Irrep_mi,14)

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Request Dco(c,d)*Dhf(k,i) = D(ck,di)"
         Call output(Work(Koff),1,Nrow_nk,1,Ncol_mi,Nrow_nk,
     +               Ncol_mi,1)
#endif
         Do irrep_i = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_i,Irrep_mi)
            Irrep_a  = Dirprd(Irrep_m,Irrep_ma)

            Nm = Pop(irrep_m,Ispin)
            Na = Vrt(irrep_a,Ispin)
            Ni = Pop(irrep_i,Ispin)

            Nsum = Nrow_nk*Nm
            Nrow = Na
            Ncol = Ni

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hvo(Ioff),
     +                     Nrow)
            Endif

            Ioff = Ioff + Ni*Na
            Joff = Joff + Nsum*Na
            Koff = Koff + Nsum*Ni
          Enddo
          Kbgn = Kbgn + Nrow_nk*Ncol_mi
      Enddo
#ifdef _DEBUG_LVL0
      call checksum("Htau_vo :",Hvo,Nocc*Nvrt)
#endif 

C W(nb,mi)*D(nb,ma) -> Hov(i,a)
C Dco(b,a)*Dhf(m,n) is D(nb,ma)

      Ioff = Ione
      Do Irrep_ma = 1, Nirrep
         Irrep_nb = Dirprd(Irrep_ma,Irrepx)

         Nma = Irpdpd(Irrep_ma,18)
         Nnb = Irpdpd(Irrep_nb,19)
         Iend = Ioff + Nma*Nnb
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nma*Nnb)
         Icount = Inull
         Do Irrep_a = 1, Nirrep
            Irrep_m = Dirprd(Irrep_a,Irrep_ma)
            Do A = 1, Vrt(Irrep_a,Ispin)
               Do M = 1, Pop(Irrep_m,Ispin)
                  Do Irrep_b = 1, Nirrep
                     Irrep_n = Dirprd(Irrep_b,Irrep_nb)
                     Do B = 1, Vrt(Irrep_b,Ispin)
                        Do N = 1, Pop(Irrep_n,Ispin)
                           An = A + Ioffv(Irrep_a)
                           Mn = M + Ioffo(Irrep_m)
                           Bn = B + Ioffv(Irrep_b)
                           Nn = N + Ioffo(Irrep_n)
                           Work(Ioff+Icount)=Dpq(Bn,An)*Dhf(Mn,Nn)
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dco(b,a)*Dhf(m,n) = D(nb,ma)"
          Call output(Work(Ioff),1,Nnb,1,Nma,Nnb,Nma,1)
#endif
          Ioff = Ioff +  Nma*Nnb
      Enddo

C Get W(mi,nb)

      Kbgn = Ione
      Do Irrep_nb = 1, Nirrep
         Irrep_mi = Dirprd(Irrep_nb,Irrepx)

         Nrow_mi = Irpdpd(Irrep_mi,14)
         Ncol_nb = Irpdpd(Irrep_nb,18)

         I010 = Iend
         I020 = I010 + Nrow_mi*Ncol_nb
         I030 = I020 + Nrow_mi*Ncol_nb
         I040 = I030 + Max(Nrow_mi,Ncol_nb)
         Iend = I040 + Max(Nrow_mi,Ncol_nb)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_nb,1,Irrep_nb,List_v2)

         Call Spinad3(Irrep_mi,Pop(1,1),Nrow_mi,Ncol_nb,Work(I020),
     +                Work(I030),Work(I040))
         Call Transp(Work(I020),Work(I010),Ncol_nb,Nrow_mi)

C W(nb,mi)*D(nb,ma) -> Hov(i,a)

         Ioff = Ione
         Joff = I010
         Koff = Kbgn
         Irrep_ma = Dirprd(Irrep_nb,Irrepx)
         Ncol_ma  = Irpdpd(Irrep_ma,18)

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "REquest Dco(b,a)*Dhf(m,n) = D(nb,ma)"
          Call output(Work(Koff),1,Ncol_nb,1,Ncol_ma,Ncol_nb,
     +               Ncol_ma,1)
#endif
         Do irrep_a = 1, Nirrep
            Irrep_m  = Dirprd(Irrep_a,Irrep_ma)
            Irrep_i  = Dirprd(Irrep_m,Irrep_mi)

            Nm = Pop(irrep_m,Ispin)
            Na = Vrt(irrep_a,Ispin)
            Ni = Pop(irrep_i,Ispin)

            Nsum = Ncol_nb*Nm
            Nrow = Ni
            Ncol = Na

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hov(Ioff),
     +                     Nrow)
            Endif

            Ioff = Ioff + Ni*Na
            Joff = Joff + Nsum*Ni
            Koff = Koff + Nsum*Na
          Enddo
          Kbgn = Kbgn + Ncol_nb*Ncol_ma
      Enddo

#ifdef _DEBUG_LVL0
      call checksum("Htau_ov :",Hov,Nocc*Nvrt)
      Print*, "HF-OV"
      call output(Hov,1,nocc,1,nvrt,nocc,nvrt,1)
#endif
      Return
      End

 

