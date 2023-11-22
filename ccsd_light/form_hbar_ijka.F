      Subroutine Form_hbar_ijka(Hbar_ijka1,Hbar_ijka2,Hbar_ijka3,
     +                          Hbar_ijka4,T1aa,T1bb,W_aa,W_bb,
     +                          W_ab,W5_aa,W5_bb,W5_ab,W5_ba,Nocc_a,
     +                          Nocc_b,Nvrt_a,Nvrt_b)
     
      Implicit Double Precision(A-H,O-Z)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

#include "ccsdlight_vars.com"

      Integer A,I,J,K,E

#ifdef _DEBUG_LVLM
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("W_aa      :",W_aa,L_aaaa)
      call checksum("W_bb      :",W_bb,L_bbbb)
      call checksum("W_ab      :",W_ab,L_abab)
      call checksum("T1aa      :",T1aa,Nvrt_a*Nocc_a)
      call checksum("T1bb      :",T1bb,Nvrt_b*Nocc_b)
#endif

C Hbar_ijka(J,K,I,A) = W(J,K,I,A) + T1(E,I)*<JK||EA>

      Do A = 1, Nvrt_a
      Do I = 1, Nocc_a 
      Do K = 1, Nocc_a
      Do J = 1, Nocc_a 
         C = Ocn_oa(I)
         Hbar_ijka1(J,K,I,A) = W5_aa(J,K,I,A)*C
      Do E = 1, Nvrt_a
         Hbar_ijka1(J,K,I,A) =  Hbar_ijka1(J,K,I,A) + 
     +                          T1aa(E,I)*W_aa(J,K,E,A)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_ijka(j,k,i,a) = W(j,k,i,a) + T1(e,i)*<jk||ea>

      Do A = 1, Nvrt_b
      Do I = 1, Nocc_b
      Do K = 1, Nocc_b
      Do J = 1, Nocc_b
         C = Ocn_ob(i)
         Hbar_ijka2(j,k,i,a) = W5_bb(j,k,i,a)*C
      Do E = 1, Nvrt_b
         Hbar_ijka2(j,k,i,a) =  Hbar_ijka2(j,k,i,a) +
     +                          T1bb(e,i)*W_bb(j,k,e,a)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_ijka(J,k,I,a) = W(J,k,I,a) + T1(E,I)*<Kj|Ea>

      Do A = 1, Nvrt_b
      Do I = 1, Nocc_a
      Do K = 1, Nocc_b
      Do J = 1, Nocc_a
         C = Ocn_oa(I)
         Hbar_ijka3(J,k,I,a) = W5_ab(J,k,I,a)*C
      Do E = 1, Nvrt_a
         Hbar_ijka3(J,k,I,a) =  Hbar_ijka3(J,k,I,a) +
     +                          T1aa(E,I)*W_ab(J,k,E,a)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C Hbar_ijka(j,K,i,A) = W(j,K,I,A) + T1(e,i)*<Jk|Ae> 

      Do A = 1, Nvrt_a
      Do K = 1, Nocc_a
      Do I = 1, Nocc_b
      Do J = 1, Nocc_b
         C = Ocn_ob(i)
         Hbar_ijka4(j,K,i,A) = W5_ba(K,j,A,i)*C 
      Do E = 1, Nvrt_b
         Hbar_ijka4(j,K,i,A) =  Hbar_ijka4(j,K,i,A) + 
     +                          T1bb(e,i)*W_ab(K,j,A,e)
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aaaa = Nocc_a*Nocc_a*Nocc_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nocc_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nocc_a*Nvrt_b
      L_baba = Nocc_b*Nocc_a*Nocc_b*Nvrt_a
      call checksum("Hbar_ijka1:",Hbar_ijka1,L_aaaa)
      call checksum("Hbar_ijka2:",Hbar_ijka2,L_bbbb)
      call checksum("Hbar_ijka3:",Hbar_ijka3,L_abab)
      call checksum("Hbar_ijka4:",Hbar_ijka4,L_baba)
#endif

      Return
      End

