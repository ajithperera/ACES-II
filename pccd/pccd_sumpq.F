      Subroutine Pccd_sumpq(H_pq,Spq,Nn,Ni,Nj,P)

      Implicit Double Precision(A-H,O-Z)

      Dimension H_pq(Ni,Nj)
      Dimension Spq(Nn,Ni,Ni,Nj)
      Character*1 P

      Do J = 1, Nj
         Do I = 1, Ni
            Do K = 1, Nn
               If (P .EQ. "N") Then
                   H_pq(I,J) = H_pq(I,J) + Spq(I,K,K,J)
               Else if (P .EQ. "T") Then
                   Print*,I,J,K,Spq(I,K,K,J)
                   H_pq(J,I) = H_pq(J,I) + Spq(I,K,K,J)
               Endif
            Enddo
         Enddo
      Enddo

      Return
      End

