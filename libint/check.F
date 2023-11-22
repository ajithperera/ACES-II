      Subroutine Check(C_f,I_f,Nbasis,K)

      Implicit Double Precision (A-H,O-Z)
      Double Precision I_f

      Dimension C_f(Nbasis,Nbasis,K)
      Dimension I_f(Nbasis,Nbasis,Nbasis,Nbasis)

      Data  DNull /0.0D0/

      Do I=1, Nbasis
      Do J=1, Nbasis 
      Do K=1, Nbasis 
      Do L=1, Nbasis 
         Sum = Dnull 
         Do M = 1, K
            Sum = Sum +C_f(K,L,M)*C_f(I,J,M) 
         Enddo 
         I_f(L,K,J,I)  = Sum
      Enddo
      Enddo
      Enddo
      Enddo

      Write(6,*) 
      Write(6,"(a)") " The final Cholesky derived integrals"

      Do I=1, Nbasis
      Do J=1, Nbasis 
      Write(6, "(5(1x,F15.10))") ((I_f(i,j,k,l),K=1,Nbasis),
     &                             L=1,Nbasis)
      Enddo
      Enddo 
      Return 
      End
   
