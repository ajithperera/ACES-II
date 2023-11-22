      Subroutine pccd_form_htau_2d_hf_hh_2(Hoo,Dpq,Dhf,Dcc,Work,Maxcor,
     +                                     Nocc,Nvrt,Nbas,List_v,Fact,
     +                                     Spin)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,Bn,An
      Character*4 Spin

      Dimension Work(Maxcor)
      Dimension Hoo(Nocc*Nocc)
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
C W(an,bi)(t)*D(an,bj) -> Hoo(i,j)
C Dco(a,b)*Dhf(n,j) is D(an,bj)
#ifdef _DEBUG_LVLM
      write(6,"(a)") "pccd_form_htau_2d_hf_hh_2:density matrices"
      call output(Dpq,1,nbas,1,nbas,nbas,nbas,1)
      call output(Dhf,1,nbas,1,nbas,nbas,nbas,1)
      call output(Dcc,1,nbas,1,nbas,nbas,nbas,1)
#endif 
     
      Ioff = Ione
      Do Irrep_bj = 1, Nirrep
         Irrep_an = Dirprd(Irrep_bj,Irrepx)

         Nbj = Irpdpd(Irrep_bj,11)
         Nan = Irpdpd(Irrep_an,11)
         Iend = Ioff + Nbj*Nan
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph_2",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nan*Nbj)
         Icount = Inull
         Do Irrep_j = 1, Nirrep
            Irrep_b = Dirprd(Irrep_j,Irrep_bj)
            Do J = 1, Pop(Irrep_j,Ispin)
               Do B = 1, Vrt(Irrep_b,Ispin)
                  Do Irrep_n = 1, Nirrep
                     Irrep_a = Dirprd(Irrep_n,Irrep_an)
                     Do N = 1, Pop(Irrep_n,Ispin)
                        Do A  = 1, Vrt(Irrep_a,Ispin)
                           Jn = J + Ioffo(Irrep_j)
                           Nn = N + Ioffo(Irrep_n)
                           Bn = B + Ioffv(Irrep_b)
                           An = A + Ioffv(Irrep_a)
C                           If (Spin .Eq. "ABAB") Then
                           Work(Ioff+Icount)=Dpq(An,Bn)*Dhf(Nn,Jn) 
C For pCCD this is zero
C                           Elseif (Spin .Eq. "AAAA") Then
C                           Work(Ioff+Icount)=Dpq(An,Bn)*Dhf(Nn,Jn) -
C                           Endif

                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dco(a,b)*Dhf(k,j) = D(an,bj)"
          Call output(Work(Ioff),1,Nan,1,Nbj,Nan,Nbj,1)
#endif
          Ioff = Ioff + Nan*Nbj
      Enddo

C Get W(an,bi)
      Kbgn = Ione
      Do Irrep_bi = 1, Nirrep
         Irrep_an = Dirprd(Irrep_bi,Irrepx)
         Irrep_bj = Dirprd(Irrep_an,Irrepx)

         Nrow_an = Irpdpd(Irrep_an,11)
         Ncol_bi = Irpdpd(Irrep_bi,11)
         Ncol_bj = Irpdpd(Irrep_bj,11)

         I010 = Iend 
         Iend = I010 + Nrow_an*Ncol_bi
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph_2",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_bi,1,Irrep_bi,List_v)

         Ioff = Ione
         Joff = I010
         Koff = Kbgn

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Reqest Dco(a,b)*Dhf(k,j) = D(an,bj)"
          Call output(Work(koff),1,Nrow_an,1,Ncol_bj,Nrow_an,
     +                Ncol_bj,1)
#endif

C W(an,bi)*D(an,bj) -> Hoo(i,j)
         Do Irrep_j = 1, Nirrep 
            Irrep_b  = Dirprd(Irrep_j,Irrep_bj)
            Irrep_i  = Dirprd(Irrep_b,Irrep_bi)

            Nb = Vrt(Irrep_b,Ispin)
            Ni = Pop(Irrep_i,Ispin)
            Nj = Pop(Irrep_j,Ispin)

            Nsum = Nrow_an*Nb
            Nrow = Ni
            Ncol = Nj

            Icheck = Min(Nsum,Nrow,Ncol)
            If (Icheck .Ne. 0)  then
               Call Dgemm("T","N",Nrow,Ncol,Nsum,Fact,Work(Joff),
     +                     Nsum,Work(Koff),Nsum,One,Hoo(Ioff),
     +                     Nrow)
            Endif

            Ioff = Ioff + Ni*Nj
            Joff = Joff + Nsum*Ni
            Koff = Koff + Nsum*Nj
         Enddo
         Kbgn = Kbgn + Nrow_an*Ncol_bj
      Enddo
#endif 
#ifdef _DEBUG_LVL0
      call checksum("Htau_oo :",Hoo,Nocc*Nocc)
      call output(Hoo,1,Nocc,1,Nocc,Nocc,Nocc,1)
#endif
      Return
      End

 

