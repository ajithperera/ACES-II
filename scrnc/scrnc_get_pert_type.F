      Subroutine Scrnc_get_pert_type(Work,Memleft,Nbfns,Label,Isympert,
     +                              Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Iaopop(8),Isympert(3)
      Character*8 Label(3)
      Dimension Iscr(Nbfns)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Parameter (Tol = 1.0D-05) 
      Label(1) = 'DIPOLE_X'
      Label(2) = 'DIPOLE_Y'
      Label(3) = 'DIPOLE_Z'

      Length = Nbfns * (Nbfns+1)/2

      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Iaopop)

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a,8(1x,i2))") "The Iaopop:",
     +                         (Iaopop(i),i=1,Nirrep)
#endif 
      Do Npert = 1, 3
         Call Getrec(20,"JOBARC",Label(Npert),Length*Iintfp,Work)
#ifdef _DEBUG_LVLM
      Write(6,*)
      write(6,"(a)") "The dipole integrals"
      write(6,"(6(1x,F12.6))") (work(i),i=1,Length)
      Write(6,*)
#endif 
         Ioff = 0
         Do Irrep = 1, Nirrep 
            Numirr = Iaopop(Irrep)
            Do I = 1, Numirr
               Ioff = Ioff + 1
               Iscr(Ioff) = Irrep
            Enddo
         Enddo 
             
         Ithru = 0
         Do Index1 = 1, Nbfns 
            Do Index2 = 1, Index1
               Ithru = Ithru + 1
               If (Dabs(Work(Ithru)) .GT. Tol) Then
                  Isym1 = Iscr(Index1)
                  Isym2 = Iscr(Index2)
                  Isym  = DIRPRD(Isym1,Isym2)
               Endif
            Enddo
         Enddo 

         Isympert(Npert) = Isym

      Enddo
#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a,3(1x,I2))") "The symmetry of dipole perturbations: ",
     +     (isympert(i),i=1,3)
#endif 

      Return
      End
