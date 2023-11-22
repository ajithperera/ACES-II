      Subroutine Pccd_mkfull(Hpq_full,Hpq_pack,Work,Maxcor,Itrlen,Nbas,
     +                       Sympack)

      Implicit Double Precision(A-H,O-Z)
      Logical Sympack

#include "syminf.com"

      Dimension Hpq_full(Nbas,Nbas)
      Dimension Hpq_pack(Itrlen)
      Dimension Work(Maxcor)
      Dimension Itriln(9)
      Dimension Isqrln(9)
      Dimension Itriof(9)
      Dimension Isqrof(9)
      Dimension Irpoff(9)
      Dimension Nbfirr(8)
   
      Indx2(I,J,N) = I+(J-1)*N

      Call Getrec(20,"JOBARC",'NUMBASIR',Nirrep,Nbfirr)

      Irpoff(1) = 0
      Do Irrep = 1, Nirrep
         Itriln(Irrep)   = Nbfirr(Irrep)*(Nbfirr(Irrep)+1)/2
         Isqrln(Irrep)   = Nbfirr(Irrep)*(Nbfirr(Irrep))
         Irpoff(Irrep+1) = Irpoff(Irrep) + Nbfirr(Irrep)
      Enddo 

      Itriof(1) = 1
      Isqrof(1) = 1
      Irpoff(1) = 0
      Do Irrep = 1, Nirrep-1
         Itriof(Irrep+1) = Itriof(Irrep) + Itriln(Irrep)
         Isqrof(Irrep+1) = Isqrof(Irrep) + Isqrln(Irrep) 
         Irpoff(Irrep+1) = Irpoff(Irrep) + Nbfirr(Irrep)
      Enddo 

      If (Sympack) Then
         Do Irrep = 1, Nirrep
            Call Expnd2(Hpq_pack(Ioff),Work,Nbfirr(Irrep))
            N = Nbfirr(Irrep)
            Do I=1, Nbfirr(Irrep)
               Do J = 1, Nbfirr(Irrep)
                  Hpq_full(Irpoff(Irrep)+I,Irpoff(Irrep)+J) =
     +            Work(Indx2(I,J,Nbfirr(Irrep)))
              Enddo
           Enddo
         Enddo 
      Else
        Call Expnd2(Hpq_pack,Hpq_full,Nbas)
      Endif 

      Return
      End
