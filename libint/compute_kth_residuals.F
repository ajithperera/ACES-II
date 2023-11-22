      Subroutine Compute_kth_residuals(Ints_current,D_current,C_km1,
     &                                 Mu_max,Nu_max,Nbasis,Ncfps_a,
     &                                 Ncfps_b,Aa1,Aa2,Bb1,Bb2,Lmax,
     &                                 Kmax,K)

      Implicit Double Precision (A-H, O-Z)

      Double Precision Ints_current 
      Integer Aa1,Aa2,Bb1,Bb2

      Dimension Ints_current(Nbasis,Nbasis,Aa1:Aa2,Bb1:Bb2)
      Dimension D_current(Nbasis,Nbasis)

      Dimension C_km1(Nbasis,Nbasis,K)
      Dimension Mu_max(Nbasis)
      Dimension Nu_max(Nbasis)

      Data Dzero /0.0D0/

#ifdef _DEBUG_LVL0
      Print*
      Print*, "@-kth_residual C_k",k
      Do i = 1, k
         call output(C_km1(1,1,i),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Enddo 
      Print*
      Write(6,*) "The <mu nu|1/r12|lam 1> integral matrix"
      Write(6, "(5(1x,F15.10))") ((Ints_current(a,b,Lmax,Kmax),
     +                             a=1,Nbasis),b=1,Nbasis)
      Print*, "Mu_ and _Nu_max"
      Write(6,"(6(1x,i3))") (Mu_max(i),i=1,K)
      Write(6,"(6(1x,i3))") (Nu_max(i),i=1,K)
#endif 
      Do L = 1, Nbasis
         Do M = 1, Nbasis
            Sum = Dzero
             Do Iter = 1, K
                J = Nu_max(Iter)
                I = Mu_max(Iter)
                Sum = Sum +  C_km1(M,L,Iter)*C_km1(Lmax,Kmax,Iter)
              Enddo
              D_current(M,L) = Ints_current(M,L,Lmax,Kmax) - Sum
         Enddo 
      Enddo

#ifdef _DEBUG_LVL0
      Print*, "@-kth_residual D_k"
      call output(D_current,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
      Return
      End

