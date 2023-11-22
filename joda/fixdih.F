      Subroutine Fixdih(Iangle,Nx)

      Implicit Double Precision(A-H,O-Z)
      Integer Dep_ind 
      Logical Dep_dih  

#include "mxatms.par"
#include "coord.com"

      Data Done,Ione /1.0D0,1/

#ifdef _DEBUG_LVL0
      Write(6,"(a,6(1x,I3))") "The connectivity pattern: ",(Ncon(i),
     &                         i=1,Nx)
#endif 
      Icon1 = Ncon(Iangle)
      Icon2 = Ncon(Iangle-1)

#ifdef _DEBUG_LVL0
      Write(6,"(a,2i2)") " The atom making the angle: ", Icon1,Icon2
#endif 

      Ioff = 1
      Do Iatom = 2, (NX/3+1)
         Dep_dih = ((Ncon(Ioff+2) .Eq. Icon1 .And. Ncon(Ioff+3) .Eq.
     +               Icon2) .Or. (Ncon(Ioff+2) .Eq. Icon2 .And.
     +               Ncon(Ioff+3) .Eq.Icon1))
         If (Dep_dih) Dep_ind =  Ioff + 3
         Ioff = (Iatom - 1)*3 
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,"(a,2i2)") " The dependent diherdral: ", Dep_ind
#endif 
      R(Dep_ind) = Dacos(-Done) - R(Dep_ind)
      Call Putrec(20,"JOBARC","DIH_FIX ",Ione,Ione)

      Return
      End 
