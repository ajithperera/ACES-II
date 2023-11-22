      Subroutine Pccd_sum_apaq(H_pq,Spq,Nn,Ni,Nj,P)

      Implicit Double Precision(A-H,O-Z)

      Dimension H_pq(Ni*Nj)
      Dimension Spq(Nn*Ni*Nn*Nj)
      Character*1 P

      call output(Spq,1,Nn*Ni,1,Nn*nj,Nn*Ni,Nn*nj,1)
      Call Pccd_sumpq(H_pq,Spq,Nn,Ni,Nj,P)
C      call output(h_pq,1,Ni,1,nj,Ni,nj,1)

      Return
      End

     
