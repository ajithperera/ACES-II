      Subroutine Pccd_insist_no_oovv(Htau_pq,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)
      Dimension Htau_pq(Nbas,Nbas)
      Integer P,Q,A,B,I,J
 
      Data Dnull /0.0D0/

      Do I = 1, Nocc
         Do J = 1, Nocc
            Htau_pq(I,J) = Dnull
            Htau_pq(J,I) = Dnull
         Enddo
      Enddo 

      Do A = 1+Nocc,Nbas
         Do B = 1+Nocc,Nbas 
            Htau_pq(A,B) = Dnull
            Htau_pq(B,A) = Dnull
         Enddo
      Enddo 

#ifdef _NOSKIP
      Do P = 1, Nbas
         Do Q = 1, Nbas
            Htau_pq(Q,P) = Htau_pq(Q,P) 
         Enddo 
      Enddo

      Do P = 1, Nbas
         Do Q = 1, P-1 
            Htau_pq(P,Q) = -Htau_pq(P,Q)
         Enddo
      Enddo 
#endif 

      Return
      End 
