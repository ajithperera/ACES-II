      Subroutine Form_t1(T1,Work,Maxcor,Pop,Vrt,Ispin,T1ln,T1ln_aa)

      Implicit Double Precision(A-H,O-Z)

      Integer Pop,Vrt 
      Integer I,A,T1ln,T1ln_aa

      Dimension Pop(8),Vrt(8)
      Dimension T1(T1ln)
      Dimension Work(Maxcor)

#include "syminf.com"
#include "machsp.com"

      Data Ione,Izero,Threshold /1,0,1.0D-09/
      Data Dzero /0.0D0/

      Ndim  = Izero 
      Ndimo = Izero 

      Do Irrep = 1, Nirrep 
         Ndim  = Ndim  + Pop(Irrep) + Vrt(Irrep)
         Ndimo = Ndimo + Pop(Irrep) 
      Enddo 
      
      I000 = Ione 
      If (Ispin .Eq. 1) Then
         Call Getrec(20,"JOBARC","SCFEVALA",Ndim*Iintfp,Work(I000))
      Else
         Call Getrec(20,"JOBARC","SCFEVALB",Ndim*Iintfp,Work(I000))
      Endif  
      Ind   = Izero 
      Indi  = Izero
      Inda0 = Izero 

      Do Irrep = 1, Nirrep
         Nrow = Vrt(Irrep)
         Ncol = Pop(Irrep)

         Do I = 1, Ncol
            Indi = Indi + Ione
            Inda = Inda0
            Do A = 1, Nrow 
               Inda = Inda + Ione
               Ind  = Ind + Ione 
               E = Work(Indi) - Work(Ndimo+Inda)
               If (Abs(E) .Gt. Threshold) Then
                  T1(Ind) = T1(Ind)/E
               Else 
                  T1(Ind) = Dzero 
               Endif 
            Enddo
         Enddo
         Inda0 = Inda0 + Nrow 
      Enddo 

      Return
      End


