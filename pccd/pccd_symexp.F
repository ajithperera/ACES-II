      Subroutine Pccd_symexp(Doo,Dvv,Dvo,Work,Maxcor,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Logical Symmetry

      Dimension Doo(Nocc*Nocc)
      Dimension Dvv(Nvrt*Nvrt)
      Dimension Dvo(Nvrt*Nocc)
      Dimension Work(Maxcor)

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"

      Common /Symm/ Symmetry 
      
      Data Ione,Itwo,Inull/1,2,0/

      Ioffv   = Ione
      Ioffo   = Ione
      Ioffvo  = Ione
      Irrepx  = Ione

      If (Symmetry) Then
         Nocc2 = Irpdpd(Irrepx,21)
         Nvrt2 = Irpdpd(Irrepx,19)
         Nvovo = Irpdpd(Irrepx,9)
      Else
         Nvrt2   = Nvrt*Nvrt
         Nocc2   = Nocc*Nocc
         Nvovo   = Nocc*Nvrt
      Endif 

      Lenoo = Nocc*Nocc
      Lenvv = Nvrt*Nvrt
      Lenvo = Nvrt*Nocc
      Lenmx = Max(Lenoo,Lenvv,Lenvo)

      Do Ispin = Ione, Ione

         I000 = Ione 
         I010 = I000 + Lenmx 
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Doo(Ioffo),"OO",Work(I010),
     +                       Work(I020),Nocc,Irpdpd(Irrepx,20+Ispin),
     +                       Nocc2,Ispin)
         Call Dcopy(Lenoo,Work(I000),1,Doo,1)

         I000 = Ione
         I010 = I000 + Lenmx
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Dvv(Ioffv),"VV",Work(I010),
     +                       Work(I020),Nvrt,Irpdpd(Irrepx,18+Ispin),
     +                       Nvrt2,Ispin)
         Call Dcopy(Lenvv,Work(I000),1,Dvv,1)

         I000 = Ione
         I010 = I000 + Lenmx
         I020 = I010 + Lenmx
         Iend = I020 + Lenmx
         If (Iend .Gt. Maxcor) Call Insmem("pccd_symexp",Iend,Maxcor)
         Call Dzero(Work(I000),3*Lenmx)
         Call Pccd_blockcopy(Work(I000),Dvo(Ioffvo),"VO",Work(I010),
     +                       Work(I020),Nvrt,Irpdpd(Irrepx,8+Ispin),
     +                       Nvovo,Ispin)
         Call Dcopy(Lenvo,Work(I000),1,Dvo,1)

         Ioffo  = Ioffo  + Irpdpd(Irrepx,20+Ispin)
         Ioffv  = Ioffv  + Irpdpd(Irrepx,18+Ispin)
         Ioffvo = Ioffvo + Irpdpd(Irrepx,8+Ispin)
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Symmetry expanded OO,VV and VO matrices"
      Call output(Doo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Call output(Dvv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Call output(Dvo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
#endif
       Return
       End


