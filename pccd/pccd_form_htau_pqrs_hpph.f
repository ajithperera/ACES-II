










      Subroutine Pccd_form_htau_pqrs_hpph(Work,Maxcor,List_v1,
     +                                    List_v2,List_g1,
     +                                    List_g2,List_h)
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C <ke||ib>*G(ke,aj)

      Irrepx = 1
      Do Irrep_aj = 1, Nirrep
         Irrep_ke = Dirprd(Irrep_aj,Irrepx) 
         Irrep_ib = Dirprd(Irrep_ke,Irrepx)

         Nrow_ke = Irpdpd(Irrep_ke,16)
         Ncol_ib = Irpdpd(Irrep_ib,16)
         Ncol_aj = Irpdpd(Irrep_aj,9)

         I000 = Ione
         I010 = I000 + Nrow_ke*Ncol_ib
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hphp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I000),1,Ncol_ib,1,Irrep_ib,List_v1)

         I020 = I010 + Nrow_ke*Ncol_aj
         Iend = I020 + Ncol_ib*Ncol_aj
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_form_htau_pqrs_hphp",
     +                                      Iend,Maxcor)

         Call Getlst(Work(I010),1,Ncol_aj,2,Irrep_aj,List_g1)

         Icheck = Min(Ncol_ij,Ncol_ab,Nrow_kl)

         If (Icheck .Gt. 0) Then
            Call Dgemm("T","N",Ncol_ib,Ncol_aj,Nrow_ke,One,Work(I000),
     +                  Nrow_ke,Work(I010),Nrow_ke,Dnull,Work(I020),
     +                  Ncol_ib)
         Endif 
         
         Call Putlst(Work(I020),1,Ncol_aj,1,Irrep_aj,List_h)
      Enddo 

      Return
      End 


