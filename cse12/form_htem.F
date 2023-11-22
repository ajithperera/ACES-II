      Subroutine form_Htem(T1aa,T1bb,Htmi_aa,Htmi_bb,Htae_aa,Htae_bb,
     +                     Htme_aa,Htme_bb,Htem_aa,Htem_bb,W3_aaaa,
     +                     W3_bbbb,W3_abab,W3_baba,W3_abba,W3_baab,
     +                     Fockov_a,Fockov_b,Fockoo_a,Fockoo_b,
     +                     Fockvv_a,Fockvv_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension W3_aaaa(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W3_bbbb(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W3_abab(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W3_baba(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension W3_abba(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension W3_baab(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Dimension Htmi_aa(Nocc_a,Nocc_a),Htmi_bb(Nocc_b,Nocc_b)
      Dimension Htae_aa(Nvrt_a,Nvrt_a),Htae_bb(Nvrt_b,Nvrt_b)
      Dimension Htme_aa(Nocc_a,Nvrt_a),Htme_bb(Nocc_b,Nvrt_b) 
      Dimension Htem_aa(Nvrt_a,Nocc_a),Htem_bb(Nvrt_b,Nocc_b)

      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)
      Dimension Fockoo_a(Nocc_a,Nocc_a),Fockoo_b(Nocc_a,Nocc_a) 
      Dimension Fockvv_a(Nvrt_a,Nvrt_a),Fockvv_b(Nvrt_b,Nvrt_b)

      Integer A,I,B,E,J,F,K

#include "cse12_vars.com"
      
C Initilize to fock elements (assumed to be Hermitian).

C Htem_aa(A,I) <- Htmi_aa(M,I)*T(A,M) + fockov(A,I) 

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
         C = Ocn_oa(I)*(1.0D0-Ocn_va(A))
         Htem_aa(A,I) = Fockov_a(I,A)*C
      Do M = 1, Nocc_a
         T = T  + Htmi_aa(M,I)*T1aa(A,M)
      Enddo
         Htem_aa(A,I) = Htem_aa(A,I) - T
      Enddo
      Enddo

C Htem_bb(a,i) <- Htmi_bb(m,i)*T(a,m) + fockov(a,i)

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
         C = Ocn_ob(i)*(1.0D0-Ocn_vb(a))
         Htem_bb(a,i) = Fockov_b(i,a)*C
      Do M = 1, Nocc_b
         T = T + Htmi_bb(m,i)*T1bb(a,m)
      Enddo
         Htem_bb(a,i) = Htem_bb(a,i) - T
      Enddo
      Enddo

C Htem_aa(A,I) <- Htae_aa(A,E)*T(E,I)

      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do E = 1, Nvrt_a
         T = T + Htae_aa(A,E)*T1aa(E,I)
      Enddo
         Htem_aa(A,I) = Htem_aa(A,I) + T
      Enddo
      Enddo

C Htem_bb(a,i) <- Htae_bb(a,e)*T(e,i)

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do E = 1, Nvrt_b
         T = T + Htae_bb(a,e)*T1bb(e,i)
      Enddo
         Htem_bb(a,i) = Htem_bb(a,i) + T
      Enddo
      Enddo

C Htem_aa(A,I) <- +T1(F,N)*W(N,A,F,I)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))*Ocn_oa(I)
         T = T + T1aa(F,N)*W3_aaaa(N,A,F,I)*C
      Enddo
      Enddo
         Htem_aa(A,I) = Htem_aa(A,I) + T
      Enddo
      Enddo

C Htem_aa(A,I) <- +T1(f,n)*W(n,A,f,I)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         C = (1.0D0-Ocn_va(A))*Ocn_oa(I)
         T = T + T1bb(f,n)*W3_baba(n,A,f,I)*C
      Enddo
      Enddo
         Htem_aa(A,I) = Htem_aa(A,I) + T
      Enddo
      Enddo

C Htem_bb(a,i) <- +T1(f,n)*W(n,a,f,i)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do n = 1, Nocc_b
      Do f = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))*Ocn_ob(i)
         T = T + T1bb(f,n)*W3_bbbb(n,a,f,i)*C
      Enddo
      Enddo
         Htem_bb(a,i) = Htem_bb(a,i) + T
      Enddo
      Enddo

C Htem_bb(a,i) <- +T1(F,N)*W(N,a,F,i)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do n = 1, Nocc_a
      Do f = 1, Nvrt_a
         C = (1.0D0-Ocn_vb(a))*Ocn_ob(i)
         T = T + T1aa(F,N)*W3_abab(N,a,F,i)*C
      Enddo
      Enddo
         Htem_bb(a,i) = Htem_bb(a,i) + T
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("Htem_aa   :",Htem_aa,L_aa)
      call checksum("Htem_bb   :",Htem_bb,L_bb)
      Write(6,*) 
#endif

      Return
      End 
   


