      Subroutine Get_cfc(Cfc_Cart,Cfc,Omega,Nrmlmodes,Work,Maxcor,
     +                   Atmmass,Nreals,Nmodes,Iamlinear)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Nrmlmodes
      Dimension Cfc_cart(9*Nreals*Nreals,Nmodes)
      Dimension Cfc(Nmodes*Nmodes,Nmodes)
      Dimension Nrmlmodes(3*Nreals,3*Nreals) 
      Dimension Work(Maxcor)
      Dimension Omega(3*Nreals)
      Dimension Atmmass(Nreals)

      Data Done,Dnull/1.0D0,0.0D0/

      Data Done,Dnull,Ione,Thres,Three/1.0D0,0.0D0,1,1.0D-09,3.0D0/

      Ndim = 3*Nreals
      Onethird = Done/Three 

#ifdef _DEBUG_LVLM
      Write(6,*) 
      Do I = 1, Nmodes
         write(6,"(a,a,I3)") "The Cartesian cubic force constants"
     +                       " for mode ", I
         call output(Cfc_cart(1,I),1,Ndim,1,Ndim,Ndim,Ndim,1)
      Enddo
#endif

C Mass-weigh the Cubic force constant vectors

      Call Massw_cfc(Cfc_cart,Atmmass,Ndim,Nmodes,Nreals)

C Upoun completion of the following call the Cfc (Cubic fore constant) is a 
C ijk element where i,j,k represent the normal modes.

      Do I = 1, Nmodes
         Call Trans_Cart_2norml_cfc(Cfc_cart(1,I),Cfc(1,I),
     +                              Nrmlmodes,Work,Maxcor,Ndim,Nmodes,
     +                              Iamlinear)
      Enddo 

C Analyse the asymmetry the Cfc matrix and possibly symmetrize?

      Call Symmet_cfc(Cfc,Nmodes) 

C Convert to cm-1.
 
      Icol = 6
      if (iamlinear .Gt. 0) Icol = 5
      Call cfc_2cm(Cfc,Omega,Nmodes,Nreals,Icol)

#ifdef _DEBUG_LVL0
      write(6,*) 
      Write(6,"(a)") " The cubic force constant matrices (per mode)"
      do i=1,Nmodes
         call output(Cfc(1,i),1,Nmodes,1,Nmodes,Nmodes,Nmodes,1)
      enddo
#endif 

      Return
      End 
