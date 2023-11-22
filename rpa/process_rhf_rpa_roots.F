      Subroutine Process_rhf_rpa_roots(Zvec, Yvec, Eigval, Ndim, Irrep,
     &                             Ispin, Scr, Maxcor)
 
      Implicit Double Precision (A-H, O-Z)
      Integer Airrep, A 

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      integer Nroot(8)
      Dimension Zvec(Ndim,Ndim),Yvec(Ndim,Ndim),Scr(Maxcor) 
      Dimension Eigval(Ndim)
      Parameter(Thresh = 0.050D0)
C
      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NRoot)
C     Write(6,*) 'NIRREP:',NIRREP
C     Write(6,*) 'NRoot: ',NRoot

      Iroot = 1
      Iuhf  = 0
      If (Irrep .EQ. 1) then
      Write(6,*)
      Write(6,"(a,a)") " RPA excitation energies and RPA state vectors",
     &                 " (only coefs. which" 
      Write(6,"(a)")   " are greater than 0.05 are printed.)" 
      Write(6,"(a,a)")   " If the spin-state is labeled as ??????",
     &                   " that means it is not assigned"
      Endif 

      Write(6,"(a,a)")  "------------------------------------",
     &                  "--------------------------------"
      Write(6,*)
      Write(6,"(a,1x,i1)") "The excited states of symmetry block:",irrep

C     Test for the different RPA NTOs with EOM-CCSD 
C
C     IF (NRoot(Irrep).EQ.1) THEN
C      Write(6,*) "YCP, Change nroot into 3",Irrep
C      NRoot(Irrep)=3
C     ENDIF
     
      Do Iroot = 1, NRoot(Irrep)
         Write(6,*)
         Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                         Dsqrt(Eigval(Iroot))*27.2113961,
     &                         "eV"
         Write(6,"(a,a)") " The spin-state        = ", "Singlet"
         Write(6,"(a,a)")  "-----------------------",
     &                     "--------------------------------"
         Write(6,*) 
         Icount =  1
         Do Iirrep = 1, Nirrep
            Airrep = Dirprd(Iirrep, Irrep)

            Do i =1, Pop(Iirrep, Ispin)
               Do a =1, Vrt(AIrrep, Ispin)
 
                  If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                   
                     Write(6,1) I, Iirrep, A, Airrep, 
     &                          Zvec(Icount,Iroot)
                
CSSS                     Write(6,"(a,a)")  "-----------------------",
CSSS     &                           "--------------------------------"
                  Endif 
                  Icount = Icount + 1
               Enddo
            Enddo 
            Write(6,"(a,a)")  "------------------",
     &                        "-------------------------------------"
  
         Enddo 
C
C        RPA densities and (TNOs) are made for the last root of each irrep
C        (ie. highest energy in a given irrep).
C
C        If (Iroot .eq. 2) then
         If (Iroot .eq. Nroot(Irrep)) then
            Write(6,*) 
            Write(6,'(A,I3,A)') 'RPA_DENS_MAIN: NTOs for RPA state',
     &           Iroot,' will be saved.'
            Call Rpa_dens_main(Irrep,Zvec(1,iroot),Ndim,Scr,Maxcor,Iuhf)
         Endif
      Enddo 

 1    Format(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     +       10X, F12.6)

      Return
      End
