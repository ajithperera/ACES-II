      Subroutine Set_diags2_one(Hess_vo,Hess_oo,Hess_vv,Lenoo,Lenvv,
     +                          Lenvo)

      Implicit Double Precision(A-H,O-Z)
      Integer A, B

      Dimension Hess_vo(Lenvo,Lenvo)
      Dimension Hess_oo(Lenoo,Lenoo)
      Dimension Hess_vv(Lenvv,Lenvv)

      Data Done /1.0D0/

      Call Dzero(Hess_vo,Lenvv*Lenoo)
      Call Dzero(Hess_oo,Lenoo*Lenoo)
      Call Dzero(Hess_vv,Lenvv*Lenvv)

      Do I = 1, Lenvo
          Hess_vo(I,I) = Done
      Enddo 

      Do I = 1, Lenoo
         Hess_oo(I,I) = Done
      Enddo 
 
      Do A = 1, Lenvv
         Hess_vv(A,A) = Done
      Enddo 

      Return
      End 

