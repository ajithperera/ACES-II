      Subroutine pccd_form_htau_2d_hf_hh_1(Hoo,Dpq,Dhf,Dcc,Work,Maxcor,
     +                                     Nocc,Nvrt,Nbas,List_v)

      Implicit Double Precision(A-H,O-Z)
      Integer D,C,Dn,Cn
      Integer B,A,Bn,An

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

      Call Dcopy(Nbas*Nbas,Dpq,1,Dcc,1)
      Do Irrep = 1, Nirrep
         N = Pop(Irrep,1)
         Do I = 1, N
            Dcc(Ioffo(Irrep)+I,Ioffo(Irrep)+I) =  
     +         Dpq(Ioffo(Irrep)+I,Ioffo(Irrep)+I) - Two
         Enddo
      Enddo

#ifdef _NO_NRML_ORD
C W(ni,mk)(t)*D(mk,nj) -> Hoo(i,j)
C Dhf(n,m)*Dhf(k,j) is D(mk,nj)
      Fact = Quart 
      write(6,"(2a)") "@-pccd_form_htau_2d_hf_hh_1: The correlated",
     +                " density matrix"
CSSS      call output(dcc,1,nbas,1,nbas,nbas,nbas,1)
CSSS      call output(dpq,1,nbas,1,nbas,nbas,nbas,1)
#ifdef _HF_ONLY
      Call dzero(Hoo,Nocc*Nocc)
      Call Dcopy(Nbas*Nbas,Dpq,1,Dhf,1)
      Fact = Quart 
      write(6,"(2a)") "@-pccd_form_htau_2d_hf_hh_1: The HF",
     +                " density matrix"
CSSS      call output(dpq,1,nbas,1,nbas,nbas,nbas,1)
#endif 
      Ioff = Ione
      Do Irrep_nj = 1, Nirrep
         Irrep_mk = Dirprd(Irrep_nj,Irrepx)

         Nmk = Irpdpd(Irrep_mk,14)
         Nnj = Irpdpd(Irrep_nj,14)
         Iend = Ioff + Nmk*Nnj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nmk*Nnj)
         Icount = Inull
         Do Irrep_j = 1, Nirrep
            Irrep_n = Dirprd(Irrep_j,Irrep_nj)
            Do J = 1, Pop(Irrep_j,Ispin)
               Do N = 1, Pop(Irrep_n,Ispin)
                  Do Irrep_k = 1, Nirrep
                     Irrep_m = Dirprd(Irrep_k,Irrep_mk)
                     Do K = 1, Pop(Irrep_k,Ispin)
                        Do M  = 1, Pop(Irrep_m,Ispin)
                           Jn = J + Ioffo(Irrep_j)
                           Nn = N + Ioffo(Irrep_n)
                           Kn = K + Ioffo(Irrep_k)
                           Mn = M + Ioffo(Irrep_m)
#ifdef _HF_ONLY
                           Work(Ioff+Icount)=Dhf(Nn,Mn)*Dhf(Kn,Jn) 
#else
                           Work(Ioff+Icount)=Dpq(Nn,Mn)*Dpq(Kn,Jn) - 
     +                                       Dcc(Nn,Mn)*Dcc(Kn,Jn)
#endif 
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dhf(m,n)*Dhf(k,j) = D(mk,nj)"
          Call output(Work(Ioff),1,Nmk,1,Nnj,Nmk,Nnj,1)
#endif
          Ioff = Ioff + Nmk*Nnj
      Enddo

C Get W(ni,mk)
      Kbgn = Ione
      Do Irrep_mk = 1, Nirrep
         Irrep_nj = Dirprd(Irrep_mk,Irrepx)
         Irrep_ni = Dirprd(Irrep_mk,Irrepx)

         Nrow_ni = Irpdpd(Irrep_ni,14)
         Ncol_mk = Irpdpd(Irrep_mk,14)
         Ncol_nj = Irpdpd(Irrep_nj,14)

         I010 = Iend 
         I020 = I010 + Nrow_ni*Ncol_mk
         I030 = I020 + Nrow_ni*Ncol_mk
         I040 = I030 + Max(Nrow_ni,Ncol_mk)
         Iend = I040 + Max(Nrow_ni,Ncol_mk)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_ph",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I020),1,Ncol_mk,1,Irrep_mk,List_v)
         Call Spinad1(Irrep_mk,Pop(1,1),Nrow_ni,Work(I020),
     +                Work(I030),Work(I040))
         Call Transp(Work(I020),Work(I010),Ncol_mk,Nrow_ni)

         Ioff = Ione
         Joff = I010
         Koff = Kbgn

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Request Dhf(m,n)*Dhf(k,i) = D(mk,nj)"
          Call output(Work(koff),1,Ncol_mk,1,Ncol_nj,Ncol_mk,
     +                Ncol_nj,1)
#endif

C W(mk,ni)*D(mk,nj) -> Hoo(i,j)
         Do Irrep_j = 1, Nirrep 
            Irrep_n  = Dirprd(Irrep_j,Irrep_nj)
            Irrep_i  = Dirprd(Irrep_n,Irrep_ni)

            Nn = Pop(irrep_n,Ispin)
            Ni = Pop(irrep_i,Ispin)
            Nj = Pop(irrep_j,Ispin)

            Nsum = Ncol_mk*Nn
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
         Kbgn = Kbgn + Ncol_nj*Ncol_mk
      Enddo
#endif 

#ifdef _DEBUG_LVL0
      call checksum("Htau_oo :",Hoo,Nocc*Nocc)
#ifdef _HF_ONLY
       Write(6,"(a)") "The SCF contribution to OO block"
       Call output(Hoo,1,Nocc,1,Nocc,Nocc,Nocc,1)
       Call Dzero(Work,Nocc*Nocc)
       Lenoo=Nocc*Nocc
       I000=Ione
       I010=I000+Lenoo
       I020=I010+Lenoo
       Iend=I020+Lenoo
       Call Putblk(Work(I000),Hoo,"OO",Work(I010),Work(I020),
     +             Nocc,Irpdpd(1,20+Ispin),Lenoo,Ispin)
      call Pccd_sum(Work(I000),Nocc,E,"D")
      write(6,*)
      write(6,"(a,(2x,F15.10))") "The reference electronic energy: ", e
      write(6,*)
#endif 
#endif
      Return
      End

 

