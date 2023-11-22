      Subroutine Pccd_check_ijka(Work,Maxcor,List_v1,LIst_v2)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione,Two /1.0D0,-1.0D0,0.0D0,1,
     +                              2.0D0/

C <Ij,Ka> 

       Ispin  = 1
       Irrepx = 1
       Do Irrep_ka = 1, Nirrep
          Irrep_ij = Dirprd(Irrep_ka,Irrepx)

          Nrow_ij  = Irpdpd(Irrep_ij,14)
          Ncol_ka  = Irpdpd(Irrep_ka,18)
          Nrow_ijc = Irpdpd(Irrep_ij,3)

          I000 = Ione
          I010 = I000 + Nrow_ij*Ncol_ka
          I020 = I010 + Nrow_ij*Ncol_ka

          Call Getlst(Work(I000),1,Ncol_ka,2,Irrep_ka,List_v2)
          call checksum("1",Work(I000),Ncol_ka*Nrow_ij)
          Call Transp(Work(I000),Work(I010),Ncol_ka,Nrow_ij)
          Call Assym2(Irrep_ij,Pop(1,Ispin),Ncol_ka,Work(I010))
          Call Symexp(Irrep_ij,Pop(1,Ispin),Ncol_ka,Work(I010))
          Write(6,"(a)") "<IJ||kA> constructed from <Ij|ka>"
          call output(Work(I010),1,Ncol_ka,1,Nrow_ij,Ncol_ka,
     +                Nrow_ij,1)
      Enddo

      Do Irrep_ka = 1, Nirrep
          Irrep_ij = Dirprd(Irrep_ka,Irrepx)

          Nrow_ij  = Irpdpd(Irrep_ij,14)
          Ncol_ka  = Irpdpd(Irrep_ka,18)
          Nrow_ijc = Irpdpd(Irrep_ij,3)

          I000 = Ionw
          I020 = I010 + Nrow_ijc*Ncol_ka
          I030 = I020 + Nrow_ij*Ncol_ka
          I040 = I030 + Nrow_ij*Ncol_ka

          Call Getlst(Work(I010),1,Ncol_ka,2,Irrep_ka,List_v1)
          Call Transp(Work(I010),Work(I020),Ncol_ka,Nrow_ijc)
          Call Symexp(Irrep_ij,Pop(1,Ispin),Ncol_ka,Work(I020))
          Write(6,"(a)") "<IJ||kA> read from list 7 and expanded"
          call output(Work(I020),1,Ncol_ka,1,Nrow_ij,Ncol_ka,
     +                Nrow_ij,1)

      Enddo
      Return
      End  
