      Subroutine Transform_2index_pt(T1,Coo,Cvv,Work,Maxcor,List1_h,
     &                               List1_t,T1ln,Iuhf,Irrepx,CCt1,
     &                               Type)
      
      Implicit Double Precision(A-H,O-Z)
      Character*3 Type
      Integer T1ln,T1off
      Logical CCt1
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "info.com"

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension T1(T1ln)
      Dimension Work(Maxcor)

      Data Ione,Inull /1,0/

      Do Ispin =1, 1+Iuhf

         Ioff   = Inull
         Ndim   = Irpdpd(Irrepx,8+Ispin)
         T1off  = (Ispin-1)*Irpdpd(Irrepx,9) + Ione

         I000 = Ione
         I010 = I000 + Ndim 
         Iend = I010 + Ndim 

         If (Iend .Gt. Maxcor) Call Insmem("transform_2index",Iend,
     +                                      Maxcor)
         If (CCt1) Then
             Call Getlst(Work(I000),1,1,1,Ispin,List1_h)
         ELse
             Call Dcopy(Ndim,T1(T1off),1,Work(I000),1)
         Endif 

         do Irrep = 1, Nirrep
            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 
            
             Nrow = Vrt(Irrepl,Ispin)
             Ncol = Pop(Irrepr,Ispin)

             Joff = Ioff_oo(Irrepr,Ispin)
             Koff = Ioff_vv(Irrepl,Ispin)

             Call Trans(Work(I000+Ioff),Work(I010+Ioff),Coo(Joff),
     +                  Cvv(Koff),Nrow,Ncol,Type)
             
             Ioff = Ioff + Nrow*Ncol
             Joff = Joff + Ncol*Ncol
             Koff = Koff + Nrow*Nrow
         Enddo

#ifdef _DEBUG_LVLM
      call checksum("T1trn:",Work(I000),Ndim)
#endif 
         If (CCt1) Then
            Call Putlst(Work(I000),1,1,1,Ispin,List1_t)
         Else
            Call Dcopy(Ndim,Work(I000),1,T1(T1off),1)
         Endif 

      Enddo 

      Return
      End
            
