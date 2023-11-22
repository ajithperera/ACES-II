      Subroutine Get_mbpt2_nos(Doo,Dvv,Coo,Cvv,Work,Maxcor,Iuhf)
     &                         
      Implicit Double Precision(A-H,O-Z)

#include "sym.com"
#include "syminf.com"

      Dimension Doo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Dvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension Work(Maxcor)

      Data Ione,Done,Dnull,Inull,Mone /1,1.0D0,0.0D0,0,-1/
     
      Do Ispin = 1, (Iuhf+1)
         Ioff = (Ispin-1)*Nfmi(1) + Ione
         Joff = (Ispin-1)*Nfea(1) + Ione

         Do Irrep = 1, Nirrep
            Nd1 = pop(Irrep,Ispin)
            Nd2 = vrt(Irrep,Ispin)

            Call Eig(Doo(Ioff),Coo(Ioff),Ione,Nd1,Mone) 
            Call Eig(Dvv(Joff),Cvv(Joff),Ione,Nd2,Mone)

#ifdef _DEBUG_LVL0
      Write(6,*) 
      If (ispin .eq. 1) then
         Write(6,"(a)") " The occ-occ mbpt(2) alpha natural orbitals"
      Else
         Write(6,"(a)") " The occ-occ mbpt(2) beta natural orbitals"
      Endif 
      call checksum("COO  :",Coo(ioff),nd1*Nd1)
CSSS      Call output(Coo(ioff),1,nd1,1,nd1,nd1,nd1,1)
      If (ispin .eq. 1) then
      Write(6,"(a)") " The vrt-vrt mbpt(2) alpha natural orbitals"
      Else
      Write(6,"(a)") " The vrt-vrt mbpt(2) beta natural orbitals"
      Endif 
CSSS      Call output(Cvv(joff),1,nd2,1,nd2,nd2,nd2,1)
      call checksum("Cvv  :",Cvv(joff),nd2*Nd2)
#endif
            Ioff = Ioff + Nd1*Nd1
            Joff = Joff + Nd2*Nd2
         Enddo 
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " Check for unitary: O(t)O and V(t)V"
      I000 = Ione
      I010 = I000 + (Nfmi(1)+Iuhf*Nfmi(2))
      Iend = I010 + (Nfea(1)+Iuhf*Nfea(2))
      Ioff = Ione
      Joff = Ione
      do irrep =1,Nirrep
      Noo  = Pop(irrep,1)
      Nvv  = vrt(irrep,1)
      call Xgemm("T","N",Noo,Noo,Noo,Done,Coo(ioff),Noo,
     +            Coo(ioff),Noo,Dnull,Work(i000+ioff-1),Noo)
      call output(Work(i000+Ioff-1),1,noo,1,noo,noo,n00,1)
      Ioff = Ioff + noo*noo
      call Xgemm("T","N",Nvv,Nvv,Nvv,Done,Cvv(joff),Nvv,
     +            Cvv(joff),Nvv,Dnull,Work(i010+joff-1),Nvv)
      call output(Work(i010+joff-1),1,nvv,1,nvv,nvv,nvv,1)
      Joff = Joff + nvv*nvv
      enddo 
      Stop
#endif 

      Return
      End


