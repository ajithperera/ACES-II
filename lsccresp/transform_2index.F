      Subroutine Transform_2index(T1,Coo,Cvv,Work,Maxcor,List1,Iuhf,
     &                            T1ln,T1ln_aa,Cct1,Type)
      
      Implicit Double Precision(A-H,O-Z)
      Character*3 Type
      Integer T1ln,T1ln_aa,T1off
      Logical Cct1
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)
      Dimension T1(T1ln)

      Data Ione,Inull /1,0/

      Irrepx = Ione
      Do Ispin =1, 1+Iuhf
         
         T1off = Ione + (Ispin-1)*T1ln_aa

         If (Type .Eq. "C2N") Then
            Iget = Ispin
            Iput = Ispin + 2
         Elseif (Type .Eq. "N2C") Then
            Iget = Ispin + 2
            Iput = Ispin + 2
         Endif 

         Ioff   = Inull
         Joff   = Ione + Nfmi(1)*(Ispin-1)
         Koff   = Ione + Nfea(1)*(Ispin-1)
         Ndim   = Irpdpd(Irrepx,8+Ispin)

         I000 = Ione
         I010 = I000 + Ndim 
         I020 = I010 + Ndim 
         Iend = I020 + Ndim 

         If (Iend .Gt. Maxcor) Call Insmem("transform_2index",Iend,
     +                                      Maxcor)
         If (Cct1) Then
            Call Getlst(Work(I000),1,1,1,Iget,List1)
         Else
            Call Dcopy(Ndim,T1(T1off),1,Work(I000),1)
         Endif 

         do Irrep = 1, Nirrep
            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 
            
             Nrow = Vrt(Irrepl,Ispin)
             Ncol = Pop(Irrepr,Ispin)

             Call Trans(Work(I000+Ioff),Work(I010+Ioff),Coo(Joff),
     +                  Cvv(Koff),Nrow,Ncol,Type)
             
             Ioff = Ioff + Nrow*Ncol
             Joff = Joff + Ncol*Ncol
             Koff = Koff + Nrow*Nrow
         Enddo

#ifdef _DEBUG_LVL0
      call checksum("T1trn:",Work(I000),Ndim)
#endif 
         If (Cct1) Then
            Call Putlst(Work(I000),1,1,1,Iput,List1)
         Else
            Call Dcopy(Ndim,Work(I000),1,T1(T1off),1)
         Endif 

      Enddo 

      Return
      End
            
