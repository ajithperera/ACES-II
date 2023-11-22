      Subroutine Get_perturb_nos(Doo,Dvv,Pcoo,Pcvv,Work,Maxcor,Iwork,
     +                           Imaxcor,Lenvv,Lenoo,Ioffset,Iuhf)

      Implicit Double Precision(A-H,O-Z)

#include "sym.com"
#include "syminf.com"

      Dimension Lenoo(2),Lenvv(2)
      Dimension Work(Maxcor)
      Dimension Dvv(Lenvv(1)+Iuhf*Lenvv(2))
      Dimension Doo(Lenoo(1)+Iuhf*Lenoo(2))
      Dimension Pcvv(Lenvv(1)+Iuhf*Lenvv(2))
      Dimension Pcoo(Lenoo(1)+Iuhf*Lenoo(2))

      Data Ione,Mone,Done,Dnull/1,-1,1.0D0,0.0D0/

C Obtain the transition natural orbitals for OO and VV density matrices.

      Do Ispin = 1, (Iuhf+1)
         Ioff = (Ispin-1)*Lenoo(1) + Ione
         Joff = (Ispin-1)*Lenvv(1) + Ione

         Do Irrep = 1, Nirrep
            Nd1 = pop(Irrep,Ispin)
            Nd2 = vrt(Irrep,Ispin)

            Call Eig(Doo(Ioff),Pcoo(Ioff),Ione,Nd1,Mone)
            Call Eig(Dvv(Joff),Pcvv(Joff),Ione,Nd2,Mone)

#ifdef _DEBUG_LVL0
      Write(6,*)
      If (ispin .eq. 1) then
         Write(6,"(2a)") " The occ-occ perturbed alpha",
     +                   " natural orbitals"
      Else
         Write(6,"(2a)") " The occ-occ perturbed beta",
     +                   " natural orbitals"
      Endif
      call checksum("PCOO :",Pcoo(ioff),nd1*Nd1)
CSSS      Call output(Pcoo(ioff),1,nd1,1,nd1,nd1,nd1,1)
      If (ispin .eq. 1) then
         Write(6,"(2a)") " The vrt-vrt perturbed alpha",
     +                   " natural orbitals"
      Else
         Write(6,"(2a)") " The vrt-vrt perturbed beta",
     +                   " natural orbitals"
      Endif
CSSS      Call output(Pcvv(joff),1,nd2,1,nd2,nd2,nd2,1)
      call checksum("Pcvv :",Pcvv(joff),nd2*Nd2)
#endif
            Ioff = Ioff + Nd1*Nd1
            Joff = Joff + Nd2*Nd2
         Enddo
      Enddo

      Call getrec(20,"JOBARC","ETA_4EXP",Ione,Eta_o)
      Eta_v = Eta_o

      Call Screen_peom_nos(Doo,Dvv,Pcoo,Pcvv,Iwork,Imaxcor,Lenvv,
     +                     Lenoo,Ioffset,Iuhf,Eta_o,Eta_v)
      Return
      ENd

