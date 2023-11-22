      Subroutine Get_fc(Fc,Nrmlmodes,P,Omega,WorK,maxcor,Atmmass,
     +                  Nreals,Nmodes,Iamlinear)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Max_offd,Max_diag,Max_asym 

      Dimension Fc(3*Nreals,Nmodes)
      Dimension P(3*Nreals,3*Nreals)
      Dimension Nrmlmodes(3*Nreals,3*Nreals) 
      Dimension Work(Maxcor)
      Dimension Atmmass(Nreals)
      Dimension Omega(3*Nreals)

      Data Done,Dnull,Ione,Thres,Au2cm1/1.0D0,0.0D0,1,1.0D-09,
     +                                  5.1404871D03/

      Ndim = 3*Nreals 

C Upoun finish the fc (fore constants) are in Cartesian representation,

      Do I = 1, Nmodes
         Do J = 1, Ndim
            Fact = Dsqrt(Atmmass(1+(J-1)/3))
            If (Fact .Lt. Thres) Then
               Fc(J,I) = Dnull
            Else
               Fc(J,I) = Fc(J,I)/Fact
            Endif
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,*) 
      write(6,"(a)") " The mass weighted Numerical Hessian"
      call output(Fc,1,Ndim,1,Nmodes,Ndim,Nmodes,1)
      write(6,"(a)") " The Normal to Cartesian transformation"
      call output(Nrmlmodes,1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Checking the projector P^2 = P"
      I000 = Ione 
      Call Dgemm("T","N",Ndim,Ndim,Ndim,Done,P,Ndim,P,Ndim,Dnull,
     +            Work(I000),Ndim)
      call Daxpy(Ndim*Ndim,-1.0D0,P,1,Work(I000),1) 
      call checkssum("P^2=P check: ", Work(I000),Ndim*Ndim,S)
#endif 

      Icol = 7
      If (Iamlinear.Gt.0) Icol = 6
      I000 = Ione 
      Iend = I000 + Ndim*Ndim
      If (Iend .Ge. Maxcor) Call Insmem("get_fc",Iend,Maxcor)

C Project the rotation and translations and then convert to normal 
C modes 

      Call Dgemm("N","N",Ndim,Nmodes,Ndim,Done,P,Ndim,Fc,Ndim,Dnull,
     +            Work(I000),Ndim)

      Call Dcopy(Ndim*Nmodes,Work(I000),1,Fc,1)

      Call Dgemm("T","N",Nmodes,Nmodes,Ndim,Done,Nrmlmodes(1,
     +            Icol),Ndim,Fc,Ndim,Dnull,Work(I000),Nmodes)

      Call Sqrt_diags(Work(I000),Nmodes,Max_asym)

#ifdef _DEBUG_LVL0
      Write(6,*)   
      write(6,"(2a)") " The sqrt of the Hessian in normal coord.",
     +                " Representation"
      call output(Work(I000),1,Nmodes,1,Nmodes,Nmodes,Nmodes,1)
      Write(6,*)
#endif
      Call Dscal(Ndim,Done/Au2cm1,Omega,1)

      Max_offd = Dnull
      Max_diag = Dnull 

      Do J = 1, Nmodes
         Do I = 1, Nmodes
            Ioff = (J-1)*Nmodes + (I-1)
            If (I .Ne. J) Then
               Max_offd = Max(Max_offd,Dabs(Work(Ioff+I000)))
            Elseif (I .Eq. J) Then
               Max_diag = Max(Max_diag,(Dabs(Work(Ioff+I000)-
     +                                  Omega(I+Icol-1))))
            Endif 
         Enddo
      Enddo 

      Write(6,"(2a,ES15.8E2)") " The abs. maximum deviation from the",
     +                         " exact eigenvalues: ",Max_diag
      Write(6,"(2a,ES15.8E2)") " The maximum abs. off-diagonal element",
     +                         "                : ",Max_offd
      Write(6,"(2a,ES15.8E2)") " The maximum assymmetry              ",
     +                         "                 : ",Max_asym

      Call Dscal(Ndim,Au2cm1,Omega,1)

      Return
      End 
