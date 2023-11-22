      Subroutine Write_brueckmos(Work,Maxcor,Iuhf,Nbas)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)

      Dimension Nbf4irrep(8), Nocc(8,2),Ireps(9)

#include "machsp.com"

      Ldim = Nbas*(Nbas+1)/2
      I000 = 1
      I010 = I000 + Nbas * Nbas 
      I020 = I010 + (Iuhf-1)*Nbas * Nbas
      I030 = I020 + Ldim 

      Call Getrec(20,"JOBARC","SCFEVCA0",Nbas*Nbas*IINTFP,
     &            Work(I000))
      If (Iuhf .NE. 0)
     &    Call Getrec(20,"JOBARC","SCFEVCB0",Nbas*Nbas*IINTFP,
     &                Work(I010))

      Call Getrec(20, "JOBARC", 'NIRREP  ',   1,      Nirrep)
      Call Getrec(20, "JOBARC", 'NUMBASIR',   Nirrep, Nbf4irrep)
      Call Getrec(20, "JOBARC", 'OCCUPYA0',   Nirrep, Nocc(1,1))
      If (Iuhf .NE. 0) Then
         Call Getrec(20, "JOBARC", 'OCCUPYB0',   Nirrep, Nocc(1,2)) 
      Else
         Call Icopy(8, Nocc(1,1),1,Nocc(1,2),1)   
      Endif 
C
      Ireps(1) = 1
      Do Irrep = 1, Nirrep
         Ireps(Irrep+1) = Ireps(Irrep)+Nbf4irrep(Irrep)
      Enddo

      Ioff = I020
      Do Ispin = 1, Iuhf + 1
         If (Ispin .Eq. 1) Joff = I000
         If (Ispin .Eq. 2) Joff = I010
      Do Irrep = 1, Nirrep
         Call Getblk(Work(I000),Work(Ioff),Nbf4irrep(Irrep),Nbas,
     &               Ireps(Irrep))
         Ioff = Ioff+ Nbf4irrep(Irrep)*Nbf4irrep(Irrep)
      Enddo
      Enddo 

      Call PutBrueckmos(Work(I020), Ldim, Iuhf, Nirrep, 
     &                     Nbf4irrep,Nocc)

      Return
      End 
 
