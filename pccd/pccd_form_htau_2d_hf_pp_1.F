      Subroutine pccd_form_htau_2d_hf_pp_1(Hvv,Dpq,Dhf,Dcc,Work,Maxcor,
     +                                     Nocc,Nvrt,Nbas,List_v,Fact)

      Implicit Double Precision(A-H,O-Z)
      Integer B,D,E,F,Bn,Dn,En,Fn

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
C W(ef,da)^(t)*D(ef,db) -> Hvv(a,b)
C Dco(d,e)*Dco(f,b) is D(ef,db)

      Ioff = Ione
      Do Irrep_db = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_db,Irrepx)

         Nef = Irpdpd(Irrep_ef,13)
         Ndb = Irpdpd(Irrep_db,13)
         Iend = Ioff + Nef*Ndb
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_pp",
     +                                     Iend,Maxcor)
         Call Dzero(Work(Ioff),Nef*Ndb)
         Icount = Inull
         Do Irrep_b = 1, Nirrep
            Irrep_d = Dirprd(Irrep_b,Irrep_db)
            Do B = 1, Vrt(Irrep_b,Ispin)
               Do D = 1, Vrt(Irrep_d,Ispin)
                  Do Irrep_f = 1, Nirrep
                     Irrep_e = Dirprd(Irrep_f,Irrep_ef)
                     Do F = 1, Vrt(Irrep_f,Ispin)
                        Do E  = 1, Vrt(Irrep_e,Ispin)
                           Bn = B + Ioffv(Irrep_b)
                           Dn = D + Ioffv(Irrep_d)
                           Fn = F + Ioffv(Irrep_f)
                           En = E + Ioffv(Irrep_e)
                           Work(Ioff+Icount)=Dpq(Dn,En)*Dpq(Fn,Bn)- 
     +                                       Dcc(Dn,En)*Dcc(Fn,Bn)
                           Icount = Icount + 1
                        Enddo
                      Enddo
                  Enddo
               Enddo
            Enddo
          Enddo
#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Built Dco(d,e)*Dco(f,b) = D(ef,db)"
          Call output(Work(Ioff),1,Nef,1,Ndb,Nef,Ndb,1)
#endif
          Ioff = Ioff + Nef*Ndb
      Enddo

C Get W(ef,da)
      Kbgn = Ione
      Do Irrep_da = 1, Nirrep
         Irrep_ef = Dirprd(Irrep_da,Irrepx)
         Irrep_db = Dirprd(Irrep_ef,Irrepx)

         Nrow_ef = Irpdpd(Irrep_ef,13)
         Ncol_da = Irpdpd(Irrep_da,13)
         Ncol_db = Irpdpd(Irrep_db,13)

         I010 = Iend 
         I020 = I010 + Nrow_ef*Ncol_da
         I030 = I020 + Max(Nrow_ef,Ncol_da)
         Iend = I030 + Max(Nrow_ef,Ncol_da)
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_2d_hf_pp",
     +                                      Iend,Maxcor)
         Call Getlst(Work(I010),1,Ncol_da,1,Irrep_da,List_v)
         Call Spinad1(Irrep_da,Vrt(1,1),Nrow_ef,Work(I010),
     +                Work(I020),Work(I030))

         Ioff = Ione
         Joff = I010
         Koff = Kbgn

#ifdef _DEBUG_LVLM
          Write(6,"(a)") "Request Dco(d,e)*Dco(f,b) = D(ef,db)"
          Call output(Work(koff),1,Nrow_ef,1,Ncol_db,Nrow_ef,
     +                Ncol_db,1)
#endif

C W(ef,da)*D(ef,db) -> Hvv(a,b)
         Do Irrep_b = 1, Nirrep 
            Irrep_d  = Dirprd(Irrep_b,Irrep_db)
            Irrep_a  = Dirprd(Irrep_a,Irrep_da)

            Nd = Vrt(irrep_d,Ispin)
            Na = Vrt(irrep_a,Ispin)
            Nb = Vrt(irrep_b,Ispin)

            Nsum = Ncol_ef*Nd
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
         Kbgn = Kbgn + Nrow_ef*Ncol_db
      Enddo
#endif 

#ifdef _DEBUG_LVL0
      call checksum("Htau_vv :",Hvv,Nvrt*Nvrt)
#endif
      Return
      End

 

