#include "flags.h" 
      Subroutine Job_control(Occnums_file) 
      Logical Ca,Fc,Bca,Bfc,Occnums_file,Move
      Logical Nwo
      Logical Scf, Pscf
      Character*80 Fname 
      
#include "flags.com"

      Data Izz /100/

      Ca  = .False.
      Fc  = .False.
      BCa = .False.
      Bfc = .False.
      Nwo = .False.
      Fcr = .False.

      If (Iflags(h_IFLAGS_orbitals) .EQ. 2 .OR. 
     +    Iflags(h_IFLAGS_orbitals) .EQ. 0)     Then
          Ca = .True. 
      ElseIf (Iflags(h_IFLAGS_orbitals) .EQ. 3) Then
          Fc = .True.
      ElseIf (Iflags(h_IFLAGS_orbitals) .EQ. 4) Then
          Bca = .True.
      ElseIf (Iflags(h_IFLAGS_orbitals) .EQ. 5) Then
          Bfc = .True.
      ElseIf (Iflags(h_IFLAGS_orbitals) .EQ. 6) Then
          Nwo = .True.
      Endif 

      Scf  = .False.
      Pscf = .False.

      If (iflags(h_IFLAGS_calclevel) .Eq. 0) Scf  = .True.
      If (iflags(h_IFLAGS_calclevel) .Eq. 1) PScf = .True.

      If (Nwo) Then
         If (Scf) Then
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Call Runit("xvprops")
            Call Runit("xvscf")
            Call Runit("xget_acesinfo")
            Call Runit("python nwchem_mo.py nwchem.out dft")
            Call Runit("xget_acesmo")
            Call Runit("mv NEWMOS.work OLDMOS")
            Iflags(h_IFLAGS_props)      = 1
            Iflags(h_IFLAGS_scf_maxcyc) = 0
            Iflags(h_IFLAGS_guess)      = 4
            Call A2putrec(20, 'JOBARC', 'IFLAGS  ', Izz, IFLAGS)
            Call Runit("xvscf")
            Call Runit("xprops")
            Call Runit("xa2mix")
         Else 
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Call Runit("xvprops")
            Call Runit("xvscf")
            Call Runit("xvtran")
            Call Runit("xintprc")
            Call Runit("xvcc")
            Call Runit("xlambda")
            Call Runit("xdens")
            Call Runit("xprops")
            Call Runit("xa2mix")
         Endif
      Else 
         If (Iflags(h_IFLAGS_brueckner) .EQ. 1) Then
            Call Do_brueckner(BCa,Bfc,Occnums_file)

C After Brueckner is completed, we need to continue with the 
C Brueckner orbitals. Lets setup the flags for that. 
C Set GUESS=READ_SO_MOS and move NEWMOS to OLDMOS and remove the
C BRUECKNER flag. 

            Call Runit("mv NEWMOS OLDMOS")
            Iflags(h_IFLAGS_guess)      = 4
            Iflags(h_IFLAGS_brueckner)  = 0
            Iflags(h_IFLAGS_scf_maxcyc) = 0
            Call A2putrec(20, 'JOBARC', 'IFLAGS  ', Izz, IFLAGS)

            Call Gfname("OCCNUMS",Fname,Ilength)
            Inquire(File=Fname(1:7),Exist=Occnums_file)
            If (Occnums_file)  Call Runit("rm OCCNUMS")
            Move = (Occnums_file .and. Bfc .and. .not. Bca)
            If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xvscf")
            Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xccsd_light")
         Else 
            Call Runit("xvmol")
            Call Runit("xvmol2ja")
            Move = (Occnums_file .and. Fc .and. .not. Ca) 
            If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xvscf")
            Write(6,*)
            If (Occnums_file) Call Runit("cp OCCNUMS_DUP OCCNUMS")
            Call Runit("xccsd_light")
         Endif 
      Endif 

      return
      End
