      Subroutine Pccd_sum_ppph_ph(Work,Hvo,Dpq,Ioo,Nv1,Nv2,Nv3,No1,
     +                            I_1,I_2,I_3,I_4,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Integer  A,B,C

      Dimension Hvo(Nvrt,Nocc)
      Dimension Dpq(Nbas,Nbas)
      Dimension Ioo(Nocc,Nocc)

      Dimension W(Nv1,Nv2,Nv3,No1)
      Data Dnull /0.0D0/

      Do I = 1, No1
      Do A = 1, Nv2
         I_i = I_i+I
         I_a = I_a+A
         Sum = Dnull
         Do C = 1, Nv3
         Do D = 1, Nv1
           I_c = I_c+C
           I_d = I_c+D
           Print*, Dpq(I_c,I_d),I_c,I_d
           Sum = Sum + Dpq(I_c,I_d)*Ioo(I_i,I_i)*
     +                 W(I_d,I_a,I_c,I_i)
         Enddo
         Enddo
        Hvo(A,I) = Hvo(A,I) + Sum
      Enddo 
      Enddo

      Return
      ENd 
  
     
