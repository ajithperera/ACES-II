      Subroutine rcc_check_sym(Work,Maxcor,Ispin,List,Type,flag)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Character*6 Flag

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Character*4 type

      Write(6,"(a,a,i4,a,1x,6a)") " Antisymmetry/symmetry check",
     +                            " of list",List,
     +                            " and type",type
      Irrepx = 1
      Ioff   = 1

      If (Flag .Eq. "noread") Then

      Do Irrep_r = 1, Nirrep

         Irrep_l = Dirprd(Irrep_r, Irrepx)

         Nrow  = Irpdpd(Irrep_l,20+ispin)
         Ncol  = Irpdpd(Irrep_r,20+Ispin)

         Do icol = 1, Ncol  
         Write(6,*) 
         Write(6,"(a,i4)") " Printing <AB|IJ> for IJ pair: ", icol

         do irrep_sr =1, Nirrep
            Irrep_sl = Dirprd(Irrep_sr, Irrep_l)

            If (Type .eq. "pphh") Then
               Nrow_s = Vrt(irrep_sl,Ispin)
               Nrow_l = Vrt(irreP_sr,Ispin)

               call output(Work(ioff),1,Nrow_s,1,Nrow_l,
     +                     Nrow_s,Nrow_l,1)
           Endif
           Ioff = Ioff + Nrow_s*Nrow_l
         Enddo
         Enddo 
      Enddo

      Elseif (Flag .Eq. "doread") Then

      Do Irrep_r = 1, Nirrep

         Irrep_l = Dirprd(Irrep_r, Irrepx)

         Nrow  = Irpdpd(Irrep_l,Isytyp(1,List))
         Ncol  = Irpdpd(Irrep_r,Isytyp(2,List))
      
         Call Getlst(Work,1,Ncol,1,Irrep_r,List)

         Do irrep_sr =1, Nirrep
            Irrep_sl = Dirprd(Irrep_sr, Irrep_l)
 
            If (Type .eq. "pphh") Then
               Nrow_s = Vrt(irrep_sl,Ispin)
               Nrow_l = Vrt(irreP_sr,Ispin)
    
               call output(Work(ioff),1,Nrow_s,1,Nrow_l,
     +                     Nrow_s,Nrow_l,1)
           Endif
           Ioff = Ioff + Nrow_s*Nrow_l

         Enddo
      Enddo 

      Endif 

      Return
      End
