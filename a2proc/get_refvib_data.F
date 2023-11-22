      SUBROUTINE GET_REFVIB_DATA(Vcoords, Norm_coords, Vomega, Omega,
     &                           AtmMass, AtmLabel, SymLabel, Coords,   
     &                           Btmp, A2grad, Grad, Hess, Vhess,  
     &                           Imap, Vib_Type, Nreals, Natoms, 
     &                           Nvibs, B2ang, Au2Invcm, Trns_state,
     &                           Mass_Weigh_nm, Mass_Weigh_gr,
     &                           Get_hess, Get_Grad, Ivib_level,
     &                           Igrad_calc,Icol)

      Implicit Double Precision (A-H, O-Z)
C
#include <machsp.com> 
#include <jodaflags.com> 
C
      Character*4 SymLabel(3*Natoms)
      Character*5 AtmLabel(Natoms)
      Character*11 Vib_Type(3*Natoms)
      Double Precision Norm_coords
C
      Logical Trns_state, Mass_weigh_nm, Get_Hess, Get_Grad, 
     &        Mass_Weigh_gr,Symmetry
      Logical HESSIANM_present
C
      Dimension VCoords(3*Natoms), Norm_coords(9*Nreals*Nreals), 
     &          Vomega(3*Nreals), Omega(3*Nreals), AtmMass(Natoms), 
     &          Coords(3*Nreals), Btmp(Natoms), Imap(Natoms),
     &          Hess(9*Natoms*Natoms), Vhess(9*Natoms*Natoms),
     &          Grad(3,Natoms), A2grad(3,Natoms)

      Data Thres /1.0D-09/

      Symmetry = (Iflags(h_IFLAGS_symmetry) .Gt. 0)
      Symmetry = .FALSE.
C The normal coordinates in NORMCORD are in M^(1/2)L dimensions.

      Call Getrec(20, "JOBARC", "REF_GEOM", Natoms*3*IINTFP,
     &            Vcoords)
      Call Getrec(20, "JOBARC", "NORMCORD", 3*Nreals*3*Nreals*
     &             IINTFP, Norm_coords)
#ifdef _DEBUG_LVL0
      Write(6,"(2a)") "The reference (with dummy atoms geometry)",
     &                " (in Angstrom)"
      Write(6, "(3F10.5)") (B2ang*VCoords(i),i=1,3*Natoms)
#endif 

C Remove the dummy atoms from the coordinates.

      Ioff = 1
      Joff = 1
      Ix   = 1
      Do Iatm = 1, Natoms
         Joff = Joff + (Iatm-1)*3
         If (AtmLabel(i)(1:5) .NE. 'X    ') Then
            Ioff = Ioff + (Ix-1)*3
            Call Dcopy(3,Vcoords(Joff),1,Vcoords(Ioff),1)
            Ix = Ix + 1
         Endif
      Enddo

#define _PROJECTED 
#ifdef _PROJECTED
C If projected frequencies are read, the VIB_TYPE record is 
C not in the correct order. Simply skip the first six frequencies.

      Call Getrec(20, "JOBARC", 'FORCECON', 3*Nreals*IINTFP,
     &            VOmega)
#else
      Call Getrec(20, "JOBARC", 'FREQUENC', 3*Nreals*IINTFP,
     &            VOmega)
#endif 
      Call Getcrec(20, "JOBARC", 'VIB_TYPE', 3*Nreals*11*IINTFP, 
     &             Vib_Type)

#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(1x,2a)") " Vibration types (correspond to unprojected",
     +                   " Hessian)"
      Write(6,"(6(1x,a))") (Vib_type(i),i=1,3*Nreals)
#endif 
C
#ifdef _PROJECTED 
      Ivibs = 0
#else 
      Ivibs = 0
#endif 
      Do Imode = 1, 3*Nreals
C
#ifdef _PROJECTED 
         If (Vomega(Imode) .lt. 0.0D0) Then
            Vomega(Imode) = - DSQRT(DABS(Vomega(Imode)))*Au2Invcm
         Else
            Vomega(Imode) = DSQRT(Vomega(Imode))*Au2Invcm
         Endif
         If (Imode .Gt. Icol) Then
             Ivibs = Ivibs + 1
             Omega(Ivibs) = Vomega(Imode)
         Endif 
#else
C
         If (Vib_Type(Imode) .EQ. ' VIBRATION ') Then
             Ivibs = Ivibs + 1
             Omega(Ivibs) = Vomega(Imode)
         Endif
#endif 
C
      Enddo
C
      If (Trns_state .AND. Omega(1) .GE. 0.0D0) Then
          Write(6, "(a)") "The IRC search must start from the TS."
CSSS          Call Errex
      Endif
C
      Call Getcrec(20, "JOBARC", "VIB_SYMS", 3*Nreals*4, SymLabel)
C
      Call Dcopy(Natoms, AtmMass, 1, Btmp, 1)
      Ireal = 0
      Icord = 1
      Jcord = 1
      Do Iatom = 1, Natoms
         Icord = Icord + (Iatom - 1)*3
         If (.Not. (Btmp(Iatom) .lt. 0.50D0)) Then
            Ireal = Ireal + 1
            Jcord = Jcord + (Ireal - 1)*3
            AtmMass(Ireal) = Btmp(Iatom)
            Call Dcopy(3, Vcoords(Icord), 1, Coords(Jcord), 1)
         Endif
         Icord = 1
         Jcord = 1
      Enddo
C
#ifdef _PROJECTED 
      Ivibs = 0
      Do Imode = 1, 3*Nreals
         If (Imode .Gt. Icol) Then
            Ivibs = Ivibs + 1
            SymLabel(Ivibs) = SymLabel(Imode)
         Endif
      Enddo
C
      Ivibs = 0
      Do Imodes = 1, 3*Nreals
         If (Imodes .Gt. Icol) Then
            Ivibs = Ivibs + 1
            Joff  = 3*Nreals*(Ivibs  - 1) + 1
            Ioff  = 3*Nreals*(Imodes - 1) + 1
            Call Dcopy(3*Nreals, Norm_coords(Ioff), 1,
     &                 Norm_coords(Joff), 1)
         Endif
      Enddo
#else
      Ivibs = 0
      Do Imode = 1, 3*Nreals
         If (Vib_Type(Imode) .EQ. ' VIBRATION ') Then
            Ivibs = Ivibs + 1
            SymLabel(Ivibs) = SymLabel(Imode)
         Endif
      Enddo
C
      Ivibs = 0
      Do Imodes = 1, 3*Nreals
C
         If (Vib_Type(Imodes) .EQ. ' VIBRATION ') Then 
            Ivibs = Ivibs + 1
            Joff  = 3*Nreals*(Ivibs  - 1) + 1
            Ioff  = 3*Nreals*(Imodes - 1) + 1
            Call Dcopy(3*Nreals, Norm_coords(Ioff), 1,
     &                 Norm_coords(Joff), 1)
         Endif
      Enddo
#endif 

#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,"(2a)") "@get_ref_vibdata; Cartesian/normal mode",
     &                   " transformation from JOBARC (dimension L)" 
   
         Call output(Norm_coords,1,3*nreals,1,3*nreals,3*nreals,
     &               3*nreals,1)
#endif 
C
      If (Mass_weigh_nm) Then 
          Call Get_massw_nrmlmodes(Norm_coords,AtmMass,Nvibs,Nreals)

#ifdef _DEBUG_LVLM
         Write(6,"(2a)") "@get_ref_vibdata; Cartesian/normal mode",
     &                   " transformation from JOBARC",
     &                     (dimension M^(1/2)L)" 
         Call output(Norm_coords,1,3*nreals,1,3*nreals,3*nreals,
     &               3*nreals,1)
#endif 
      Endif 
C The Hessian that is read is exactly what enters into vib1.F in 
C joda. This is ordered as in ZMAT order and the dummy atoms have
C been removed.

      If (Get_Hess) then
          Length = 9*Nreals*Nreals 
          Call Getrec(20,'JOBARC','CARTHESC',Length*IINTFP,Hess)

#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian from JOBARC"
         Call output(Hess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
#endif
      Endif 

#ifdef _NOSKIP
C These blocks of code are obsolete and some parts may even be not
C right. 

      If (Get_Hess) then
         If (Symmetry) Then
            Length = 9*Nreals*Nreals 
            Call Getrec(20,'JOBARC','CART_HES',Length*IINTFP,
     &                  VHess)
#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian from JARC"
         Call output(Vhess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
#endif
         Else
            Length = 9*Natoms*Natoms 
            Call Getrec(0,'JOBARC','HESSIANM',Length_int,0)
            HESSIANM_present = .False. 
            If (Length_int .Gt. 0) HESSIANM_present = .True.
            If (Length_int .Eq. Length) Then
               Call Getrec(20,'JOBARC','HESSIANM',Length*IINTFP,
     &                     VHess)
            Else
              Write(6,"(2a)") "Stored length is diffetent from"
     &                        " the requested length for the"
              Write(6,"(a)")  "HESSIANM record. "
              Call Errex
            Endif 
       
#ifdef _DEBUG_LVL0
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian from JARC"
         Call output(Vhess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
#endif
         Endif 
c
c-----Transform to the ZMAT order and remove dummy atoms.
C
         Call Zero(Hess, Length)
         Ioff = 1
         If (Ivib_level .Eq. 1) Then
         Do Ivmol = 1, Natoms
            Iatmzmat = Imap(Ivmol)
            Do Ixyz = 1, 3
               Icol = Ixyz + (Iatmzmat - 1)*3
               Do Jvmol = 1, Natoms
                  Jatmzmat = Imap(Jvmol)
                  If (Iatmzmat .NE. 0 .AND. Jatmzmat .NE. 0) Then
                      Irow = 1 + (Jatmzmat - 1)*3
                      Call Blkcpy(Vhess(ioff), 3, 1, Hess, 3*Natoms,
     &                            3*Natoms, Irow, Icol)
                      Ioff = Ioff + 3
                  Endif 
               Enddo
            Enddo
         Enddo
         Else 
         Call Dcopy(Length,Vhess,1,Hess,1)
         Endif 

#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,*) "@get_ref_vibdata; Hessian in ZMAT order"
         Call output(Hess,1,3*natoms,1,3*natoms,3*natoms,3*natoms,1)
#endif
C
C-----Eleminate the entry for dummy atoms
C
C This block is not needed
         Ix = 0
         Jx = 1
         Do Jatm = 1, 3*Natoms
            Do Iatm = 1, 3*Natoms
               Ink = 1 + (Iatm - 0.1)/3
               Jnk = 1 + (Jatm - 0.1)/3
               If ((AtmLabel(Ink)(1:5) .NE. 'X    ') .AND.
     &            ( AtmLabel(Jnk)(1:5) .NE. 'X    ')) Then
                   Ix = Ix + 1
                   Hess(Ix, Jx) = Hess(Iatm, Jatm)
                   If (Ix .EQ. 3*Natoms) Ix = 0
                   If (Ix .EQ. 0) Jx = Jx + 1
               Endif
           Enddo
         Enddo

      Endif
#endif 
C
#ifdef _DEBUG_LVLM
         Write(6,*)  
         Write(6,*) "@get_ref_vibdata; Hessian after removing dummies"
         Call output(hess,1,3*nreals,1,3*nreals,3*nreals,3*nreals,1)
#endif
C
      If (Get_Grad) Then

          Call Getrec(20,'JOBARC','GRADIENT',3*Nreals*IINTFP,
     &                A2Grad(1,1))
C
          Call Dzero(Grad, 3*Natoms)

C Remove dummy atoms 

          Do i=1, Natoms
             If (AtmLabel(i)(1:5) .NE. 'X    ') Then
                Do j=1,3
                   Grad(j,k)  = A2grad(j,i)
                Enddo
             Endif
          Enddo                  
C        
          Call Dcopy(3*Natoms, Grad, 1, A2grad, 1)

C Change to ZMAT order

          Do i=1, Natoms
             If (Imap(i) .NE. 0) Then
                 Call Dcopy(3, A2grad(1,Imap(i)), 1, Grad(1,i), 1)
             Endif
          Enddo
C
          If (Mass_Weigh_gr) Then 
             Do Iatom = 1, Nreals
                Do Ixyz = 1, 3
                   Grad(Ixyz, Iatom) = Grad(Ixyz, Iatom)/
     &                                 Dsqrt(AtmMass(Iatom))
                Enddo
             Enddo
          Endif
C
      Endif
C
      Return
      End

