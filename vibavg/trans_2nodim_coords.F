      Subroutine Trans_2nodim_coords(Prop_derv1,Prop_derv2,Omega,
     +                               Ndim,Nreals,Nmodes,Icol,
     +                               NMR_spincc,NMR_shifts) 


      Implicit Double Precision(A-H,O-Z)

      Dimension Prop_derv1(Ndim,Nmodes)
      Dimension Prop_derv2(Ndim,Nmodes)
      Dimension Omega(3*Nreals)
      Logical NMR_spincc
      Logical NMR_shifts

#include "const.h" 

      Data Done,Two /1.0D0,2.0D0/
     
      Pi = DAcos(-1.0D+00)
      C1 = Dsqrt(Two*Pi*CONST_c*1.0D02/CONST_hbar)
      C2 = CONST_a0
      C3 = Done/CONST_u
     
C Factor to convert first-derivative of propertieS to dimensionless
C normal coords. Their naked units are Hz Kg^(-1/2)m
    
      Factor_1der = (Done/(C1*C2))*Dsqrt(C3)

C Factor to convert second-derivative of couplings to dimensionless
C normal coords. Their naked units are Hz Kg^(-1)m^(2)

      Factor_2der = ((Done/(C1*C2))**2)*(C3)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a,ES18.9E2)") " The conversion factor for first",
     +                          " derivatives :",Factor_1der 
      Write(6,"(2a,ES18.9E2)") " The conversion factor for second",
     +                          " derivatives:", Factor_2der 
#endif
      Do I = 1, Nmodes
         D1 = Dsqrt(Omega(I+Icol))
         D2 = D1*D1
         If (NMR_spincc) Then
            Call Trans_derv_spncc(Prop_derv1(1,I),Factor_1der,D1,
     +                            Nreals)
            Call Trans_derv_spncc(Prop_derv2(1,I),Factor_2der,D2,
     +                            Nreals)
         Elseif (NMR_shifts) Then
            Call Trans_derv_shifts(Prop_derv1(1,I),Factor_1der,D1,
     +                             Nreals)
            Call Trans_derv_shifts(Prop_derv2(1,I),Factor_2der,D2,
     +                             Nreals)
         Endif 
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(2a)") " The first and second derivatives of ",
     +                " the property (in its respective units)"
      If (NMR_spincc) Then 
      do i=1,nmodes
      write(6,"(a,I2)") " The mode: ", i
      call output(Prop_derv1(1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      call output(Prop_derv2(1,i),1,Nreals,1,Nreals,Nreals,Nreals,1)
      enddo
      Elseif (NMR_shifts) Then 
      do i=1,nmodes
      write(6,"(a,I2)") " The mode: ", i
      write(6,"(5(1x,F15.8))") (Prop_derv1(j,i),j=1,Nreals)
      write(6,"(5(1x,F15.8))") (Prop_derv2(j,i),j=1,Nreals)
      enddo
      endif 
#endif

      Return 
      End  

