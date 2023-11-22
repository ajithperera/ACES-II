      Subroutine Get_massw_nrmlmodes(Norm_coords,AtmMass,Nvibs,
     +                               Nreals)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Norm_coords(3*Nreals,Nvibs) 
      Dimension AtmMass(Nreals)

      Data Done,Thres /1.0D0,1.0D-09/
     
      Ndims = 3*Nreals 

      Do I = 1, Nvibs
         Do J = 1, Ndims
            SqrtMass = Dsqrt(AtmMass(1+(J-1)/3))
            Norm_coords(J,I) = Norm_coords(J,I)*SqrtMass
         Enddo
      Enddo 

      Do I = 1, Nvibs
         Dnorm = Dsqrt(Ddot(Ndims,Norm_coords(1,I),1,
     +                      Norm_coords(1,I),1))
         If (Dnorm .Gt. Thres) Then 
             Dnorm = Done/Dnorm
             Call Dscal(Ndims,Dnorm,Norm_coords(1,I),1)
         Else
             Write(6,"(a)") "@-get_nomass_nrmlmodes, division",
     +                      " by zero error occurs" 
             Call Errex
         Endif
      Enddo

      Return
      End 


