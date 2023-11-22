      Subroutine compute_1st_residuals(Dk_m2,Dk_m1,Ck_m1,Nbasis,
     &                                 Nbasis2,Dmax,Index_a,Index_b)

      Implicit Double Precision (A-H, O-Z)

      Dimension Dk_m2(Nbasis*Nbasis)
      Dimension Dk_m1(Nbasis*Nbasis)
      Dimension Ck_m1(Nbasis*Nbasis,1:1)

      Data Dzero /0.0D0/


      Dmax_new = Dzero 
      Dmax_old = Dzero 

      Do I = 1, Nbasis2 
         Dk_m1(I) = Dk_m2(I) - Ck_m1(I,1)*Ck_m1(I,1)
         Dmax_new = Max(Dmax_old,Dk_m1(I))
         If (Dmax_new .Ne. Dmax_old) Then
             Imax = I
         Endif 
         Dmax_old = Dmax_new
      End do 

      I = Mod(Imax,Nbasis)
      J = (Imax-I)/Nbasis 

      Return
      End
