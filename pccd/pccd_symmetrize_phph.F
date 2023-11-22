      Subroutine Pccd_symmetrize_phph(H_in,Nsize,Work,Maxcor,Ex)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Dimension H_in(Nsize)
      Logical Ex

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "info.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

C Construct H(ai,bj) = H(ai,bj) + H(bi,aj)

C First H(ai,bj) -> H(aj,bi), then transpose and sum.
  
      Irrepx = Ione 
      Ispin  = Ione
      
      If (Ex) Then
         I000 = Ione
         I010 = I000 + Nsize
         Iend = I010 + Max(Nvrto(1)*Nocco(1),Nvrto(1)*Nvrto(1),Nocco(1)*
     +                  Nocco(1))
         If (Iend .Gt. Maxcor) Call Insmem("Pccd_Pccd_symmetrize_phph",
     +                                   Iend,Maxcor)
         Call Sstgen(H_in,Work(I000),Nsize,Vrt(1,Ispin),Pop(1,Ispin),
     +               Vrt(1,Ispin),Pop(1,Ispin),Work(I010),Irrepx,
     +               "1432") 
      
         Ioff   = I000 
         Do Irrep_rph= 1, Nirrep
            Irrep_lph = Dirprd(Irrep_rph,Irrepx)

            Nrow_ph = Irpdpd(Irrep_lph,11)
            Ncol_ph = Irpdpd(Irrep_rph,11)
            Iend    = I010 + Ncol_ph*Nrow_ph
            If (Iend .Gt. Maxcor) 
     +          Call Insmem("Pccd_Pccd_symmetrize_phph",Iend,Maxcor)
            Call Transp(Work(Ioff),Work(I010),Ncol_ph,Nrow_ph)
            Call Daxpy(Nrow_ph*Ncol_ph,One,Work(I010),1,H_in(Ioff),1)
            Ioff = Ioff + Nrow_ph*Ncol_ph
         Enddo 
      Else
         I000 = Ione
         Ioff = Ione
         Do Irrep_rph= 1, Nirrep
            Irrep_lph = Dirprd(Irrep_rph,Irrepx)

            Nrow_ph = Irpdpd(Irrep_lph,11)
            Ncol_ph = Irpdpd(Irrep_rph,11)
            Iend    = I000 + Ncol_ph*Nrow_ph
            If (Iend .Gt. Maxcor) 
     +          Call Insmem("Pccd_Pccd_symmetrize_phph",Iend,Maxcor)
            Call Transp(H_in(Ioff),Work(I000),Ncol_ph,Nrow_ph)
            Call Daxpy(Nrow_ph*Ncol_ph,One,Work(I000),1,H_in(Ioff),1)
            Ioff = Ioff + Nrow_ph*Ncol_ph
         Enddo 
        
      Endif 

      Return
      End
