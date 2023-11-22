      Subroutine Process_uhf_rpa_roots(Zvec, Yvec, Eigval, Ndim, Irrep,
     &                                 Nsizea, Iuhf, Scr, Maxcor)

      Implicit Double Precision (A-H, O-Z)
      Integer Airrep, A 

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Dimension Nroot(8)
      Dimension Zvec(Ndim,Ndim),Yvec(Ndim,Ndim),Scr(Maxcor)
      Dimension Eigval(Ndim)
      Parameter(Thresh = 0.050D0)
      Character*7 Spin_state(40)
C
      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NRoot)
C
      If (Irrep .EQ. 1) then
      Write(6,*)
      Write(6,"(a,a)") " RPA excitation energies and RPA state vectors",
     &                 " (only coefs. which" 
      Write(6,"(a)")   " are greater than 0.05 are printed.)" 
      Write(6,"(a,a)")   " If the spin-state is labeled as ??????", 
     &                   " that means it is not assigned" 
      Endif 

      Call Assign_spin_state(Zvec,Nsizea,Spin_state,Nroot,Irrep)

      Write(6,"(a,a)")  " ------------------------------------",
     &                  "--------------------------------"
      Write(6,*)
      Write(6,"(a,1x,i1)") "The excited states of symmetry block:",
     &                      irrep

      If  (Nirrep .EQ. 1) Then

          Do Iroot = 1, Nroot(Irrep)
          Do Ispin = 1, 2
             Icount = (Ispin - 1) * Nsizea  + 1

             If (Ispin .Eq. 1) Then
             Write(6,*)
             Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                             Dsqrt(Eigval(Iroot))*27.2113961,
     &                             " eV"
             Write(6,"(a,7a)") " The spin-state        = ", 
     &                           Spin_state(Iroot)
             Write(6,*)
             Endif 
             If (Ispin .EQ. 1) Write(6,"(a)") " The AA coeffcients" 
             If (Ispin .EQ. 2) Write(6,"(a)") " The BB coeffcients" 
             Write(6,"(a,a)")  "------------------------------------",
     &                         "-------------------"
             Do i =1, Pop(Irrep, Ispin)
                Do a =1, Vrt(Irrep, Ispin)
 
                   If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                    
                      Write(6,1) I, Irrep, A, Irrep, 
     &                           Zvec(Icount,Iroot)
                    
CSSS                      Write(6,"(a,a)")  "------------------",
CSSS     &               "-------------------------------------"
                   Endif 
                   Icount = Icount + 1

                Enddo 
             Enddo 

          Enddo 
          Write(6,"(a,a)")  "------------------",
     &                      "-------------------------------------"
          Enddo 
      Else 
          Do Iroot = 1, Nroot(Irrep)
          Do Ispin = 1, 2
             Icount = (Ispin - 1) * Nsizea + 1
 
             If (Ispin .EQ. 1) then
             Write(6,*)
             Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                             Dsqrt(Eigval(Iroot))*27.2113961,
     &                             " eV"
             Write(6,"(a,7a)") " The spin-state        = ", 
     &                           Spin_state(Iroot)
             Endif 
             Write(6,*)
             If (Ispin .EQ. 1) Write(6,"(a)") " The AA coeffcients" 
             If (Ispin .EQ. 2) Write(6,"(a)") " The BB coeffcients" 
             Write(6,"(a,a)")  "------------------------------------",
     &                         "-------------------"
             Do Iirrep = 1, Nirrep
                Airrep = Dirprd(Iirrep, Irrep)

                Do i =1, Pop(Iirrep, Ispin)
                   Do a =1, Vrt(AIrrep, Ispin)
 
                      If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                       
                         Write(6,1) I, Iirrep, A, Airrep, 
     &                              Zvec(Icount,Iroot)
                    
CSSS                      Write(6,"(a,a)")  "------------------",
CSSS     &               "-------------------------------------"
                      Endif 
                      Icount = Icount + 1
                   Enddo
                Enddo 
             Enddo 

          Enddo 
          Write(6,"(a,a)")  "------------------",
     &                      "-------------------------------------"
          Enddo
      Endif 
C
C RPA densities and (TNOs) are made for the first root of each irrep
C (ie. lowest energy in a given irrep).

      Iroot = 1
      Call Rpa_dens_main(Irrep,Zvec(1,iroot),Ndim,Scr,Maxcor,Iuhf)

 1    Format(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     +       10X, F12.6)

      Return
      End
