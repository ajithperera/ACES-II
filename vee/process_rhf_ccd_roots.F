      Subroutine Process_rhf_ccd_roots(Zvec,Yvec,Eigval,Dbls_corr,Ndim,
     &                                 Irrep, Ispin,Rpa,
     &                                 Eom_sf_rccd,Eom_sf_drccd,
     &                                 Eom_s_rccd,Eom_s_drccd,
     &                                 Imult,Jroot)
 
      Implicit Double Precision (A-H, O-Z)
      Integer Airrep, A,Tnroots
      Logical Rpa,Rccd,DrcCd,Eom_s_rccd,Eom_s_drccd
      Logical Eom_sf_rccd,Eom_sf_drccd,Flag
      Logical Syminfo_exist
      Dimension E0_keep(400,2),Ed_keep(400,2)
      Character*7 S_keep(400,2)
      Character*4 Irpchar(8)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      integer Nroot(8)
      Dimension Zvec(Ndim,Ndim),Yvec(Ndim,Ndim)
      Dimension Eigval(Ndim),Dbls_corr(Ndim,8) 
      Parameter(Thresh = 0.050D0)
      Data Au2ev /27.2113961D0/
C
      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NRoot)
      Tnroots = 0
      Flag    = .False.
      Do I = 1, Nirrep
         Tnroots = Tnroots + Nroot(I)
         If (Nroot(I) .Gt. 50) Flag = .True.
      Enddo
      If (Flag) Then
         Write(6,*)
         Write(6,"(3a,1x,3a)") " The number of roots requested for",
     +                         " at least one of the irreps is greater"
         Write(6,"(a,1x,3a)")  " than the maximum allowed","50."
         Write(6,*)
         Call Errex
     +
      Endif
      If (Tnroots .Gt. 400) Then
         Write(6,*)
         Write(6,"(2a,1x,i3,2a,1x,4a)") " The requested total number",
     +                  " of roots", Tnroots, " exceeded the maximum",
     +                  " allowed","800."
         Write(6,*)
         Call Errex 
      Endif

      Iroot = 1
      Iuhf  = 0
      If (Irrep .EQ. 1) then
      Write(6,*)
      If (Eom_sf_rccd) Then
      Write(6,"(a,a)") " EOM(SF)-RCCD (RPA) excitation energies",
     +                 " and state vectors" 
      Write(6,"(a,a)")   " (only coefs. which are greater than",
     +                   " 0.05 is printed.)"
      Elseif (Eom_sf_drccd) Then
      Write(6,"(a,a)") " EOM(SF)-DRCCD (DRPA) excitation energies",
     +                 " and  state vectors"
      Write(6,"(a,a)")   " (only coefs. which are greater than",
     +                   " 0.05 is printed.)"
      Elseif (Eom_s_rccd) Then
      Write(6,"(a,a)") " EOM(S)-RCCD excitation energies",
     +                 " and  state vectors" 
      Write(6,"(a,a)")   " (only coefs. which are greater than",
     +                   " 0.05 is printed.)"
      Elseif (Eom_s_drccd) Then
      Write(6,"(a,a)") " EOM(S)-DRCCD excitation energies",
     +                 " and  state vectors" 
      Write(6,"(a,a)")   " (only coefs. which are greater than",
     +                   " 0.05 is printed.)"
      Endif 
      Endif 

      Write(6,"(a,a)")  "---------------------------------------",
     &                  "-------------------"
      Write(6,*)
      Write(6,"(a,1x,i1)") "The excited states of symmetry block:",irrep

      Do Iroot = 1, NRoot(Irrep)
         Jroot = Jroot + 1
         Write(6,*)
         If (Eom_sf_rccd) Then 
         Write(6,"(a,F12.6,1x,a)") " The EOM(SF)-RCCD energy   = ",
     &                          Eigval(Iroot)*Au2ev,"eV"
         ElseIf (Eom_sf_drccd) Then 
         Write(6,"(a,F12.6,1x,a)") " The EOM(SF)-DRCCD energy  = ",
     &                          Eigval(Iroot)*Au2ev,"eV"
         ElseIf (Eom_s_rccd) Then 
         Write(6,"(a,F12.6,1x,a)") " The EOM(S)-RCCD energy    = ",
     &                          Eigval(Iroot)*Au2ev,"eV"
         ElseIf (Eom_s_drccd) Then 
         Write(6,"(a,F12.6,1x,a)") " The EOM(S)-DRCCD energy   = ",
     &                          Eigval(Iroot)*Au2ev,"eV"
         Endif 
         E0_keep(Jroot,Imult) = Eigval(Iroot)
         Ed_keep(Jroot,Imult) = Dbls_corr(Iroot,Irrep)

         Write(6,"(a,F12.6,1x,a)") " The doubles correction    = ",
     &                          Dbls_corr(Iroot,Irrep)*Au2ev,"eV"
         Write(6,"(a,F12.6,1x,a)") " The excitation energy     = ",
     &                          (Eigval(Iroot)+
     &                          Dbls_corr(Iroot,Irrep))*Au2ev,"eV"
         If (Imult .EQ. 1) Then
            Write(6,*)
            Write(6,"(a,a)") " The spin-state            = ", 
     &                       "Singlet"
            S_keep(Jroot,Imult)  = "Singlet"
         Elseif (Imult .EQ. 2) then 
            Write(6,*)
            Write(6,"(a,a)") " The spin-state            = ", 
     &                       "Triplet"
            S_keep(Jroot,Imult)  = "Triplet"
         Endif 

         Write(6,*)
         Icount =  1
         Do Iirrep = 1, Nirrep
            Airrep = Dirprd(Iirrep, Irrep)

            Do i =1, Pop(Iirrep, Ispin)
               Do a =1, Vrt(AIrrep, Ispin)
 
                  If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                   
                     Write(6,1) I, Iirrep, A, Airrep, 
     &                          Zvec(Icount,Iroot)
                
                  Endif 
                  Icount = Icount + 1
               Enddo
            Enddo 
         Enddo 
         Write(6,"(a,a)")  "---------------------------------------",
     &                     "-------------------"
      Enddo 

      If (Irrep .Eq. Nirrep) Then
         If (Eom_sf_rccd) Then
            Write(6,"(a)") " The summary of the EOM(SF)-RCCD results."
         ElseIf (Eom_sf_drccd) Then
            Write(6,"(a)") " The summary of the EOM(SF)-DRCCD results."
         ElseIf (Eom_s_rccd) Then
            Write(6,"(a)") " The summary of the EOM(S)-RCCD results."
         ElseIf (Eom_s_drccd) Then
            Write(6,"(a)") " The summary of the EOM(S)-DRCCD results."
         Endif

        Inquire(file='SYMINFO',Exist=Syminfo_exist)
        if (Syminfo_exist) then
        Open(unit=10, FILE='SYMINFO',STATUS='OLD')
        Do I = 1, Nirrep
           Read(10,599) J, Irpchar(i)
        Enddo
 599    format(i4,2x,4a)
        Close(unit=10, status='KEEP')
        Endif
        If (Imult .Eq. 2) Then

        Write(6,*)
        Write(6,"(2a)") "-----------------------------------------",
     +                   "----------------------"
        Write(6,"(2a,10x,a)") "Irrep","      Excitation Energy",
     +                         "       Multiplicity"
        Write(6,"(a)")          "                (eV)"
        Write(6,"(2a)") "-----------------------------------------",
     +                   "----------------------"
        Write(6,"(3a)") "        Base value","   Doubles",
     +                   "    Exc. Eng."
        Write(6,"(3a)") "                     correction"
        Write(6,"(2a)") "-----------------------------------------",
     +                   "----------------------"
        Do I =1, 2
        Kroot = 0
        Do Isym = 1, Nirrep
           Do Iroot = 1, Nroot(Isym)
            Kroot = Kroot + 1
            Write(6,"(a,1x,3(F12.6),7x,7a)") Irpchar(Isym),
     +                E0_keep(Kroot,I)*Au2ev,Ed_keep(Kroot,I)*Au2ev,
     +                (E0_keep(Kroot,I)+Ed_keep(Kroot,I))*Au2ev,
     +                S_keep(Kroot,I)

           Enddo
        Enddo
        Write(6,*)
        Enddo
        Write(6,"(2a)") "-----------------------------------------",
     +                  "----------------------"
      Endif
      Endif 
C
 1    Format(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     +       10X, F12.6)

      Return
      End
