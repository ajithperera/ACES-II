      Subroutine Pccd_asymm(W,WA,Np,Nq)

      Implicit Double Precision(A-H,O-Z)
      Integer p,q

      Dimension W(Np,Nq),WA(Np,Nq)

      Do P = 1, Np
         Do Q = 1, Nq
            WA(P,Q)  = W(P,Q) - W(Q,P)
         Enddo
      Enddo 

      Return
      End
