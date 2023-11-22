      Subroutine Init_2abij(Resid_aa,Resid_bb,Resid_ab,W0_aa,
     +                      W0_bb,W0_ab,Nocc_a,Nocc_b,Nvrt_a,
     +                      Nvrt_b)

      Implicit Double Precision(A-H,O-Z)
      
      Dimension Resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension Resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension Resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W0_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension W0_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension W0_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Integer I,J,A,B

C AAAA block 

      DO J = 1, Nocc_a
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_a
      DO A = 1, Nvrt_a
         Resid_aa(A,B,I,J) = W0_aa(A,B,I,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C BBBB block 

      DO J = 1, Nocc_b
      DO I = 1, Nocc_b
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_b
       Resid_bb(a,b,i,j) = W0_bb(a,b,i,j)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C ABAB block 


      DO J = 1, Nocc_b
      DO I = 1, Nocc_a
      DO B = 1, Nvrt_b
      DO A = 1, Nvrt_a
         Resid_ab(A,b,I,j) = W0_ab(A,b,I,j) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

#ifdef _DEBUG_LVL0
      Write(6,*) 
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("Resid_aa  :",Resid_aa,L_aaaa)
      call checksum("Resid_bb  :",Resid_bb,L_bbbb)
      call checksum("Resid_ab  :",Resid_ab,L_abab)
#endif

      Return
      End
