      Subroutine Form_new_t1pamps(T1amps,Work,Iwork,Maxcor,Imaxcor,
     &                            List1,T1ln,T1ln_aa,Ipert,Irrepx,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)
      Character*3 Type
      Integer T1ln,T1off,T1ln_aa
      Integer Aend
      Logical CCt1
      
#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "active_space.com"


      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension T1amps(T1ln)
      Dimension Work(Maxcor)

      Data Ione,Inull /1,0/

      Do Ispin =1, 1+Iuhf

         Ioff   = Inull
         Ndim   = Irpdpd(Irrepx,8+Ispin)
         T1off  = (Ispin-1)*Ndim + Ione

         I000 = Ione
         I010 = I000 + Ndim 
         Iend = I010 + Ndim 

         If (Iend .Gt. Maxcor) Call Insmem("transform_2index",Iend,
     +                                      Maxcor)

C EOM-CCSD Tx on the disk while PEOM-CCSD is in memory.

         Call Dcopy(Ndim,T1amps(T1off),1,Work(I010),1)
         Call Getlst(Work(I000),1,1,1,Ispin,List1)

         do Irrep = 1, Nirrep
            Irrepr = Irrep
            Irrepl = Dirprd(Irrepr,Irrepx) 
            
             Nrow = Vrt(Irrepl,Ispin)
             Ncol = Pop(Irrepr,Ispin)

             Aend = Pactive_vv(Ipert,Irrepl,Ispin)
             Iend = Pactive_oo(Ipert,Irrepr,Ispin)

             Call Form(Work(I000+Ioff),Work(I010+Ioff),Iwork,
     +                      Imaxcor,Nrow,Ncol,Aend,Iend)

             Ioff = Ioff + Nrow*Ncol
             
         Enddo

         Call Dcopy(Ndim,Work(I010),1,T1amps,1)

#ifdef _DEBUG_LVLM
      call checksum("T1trn:",T1amps,Ndim)
#endif 
         Call Putlst(T1amps,1,1,1,Ispin,List1)
      Enddo 

      Return
      End
            
