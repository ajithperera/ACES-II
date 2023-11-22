      Subroutine T2_inrai(T2aa,T2bb,T2ab,W4_aa,W4_bb,W4_ab,W4_ba,
     +                    W5_aa,W5_bb,W5_ab,W5_ba,T1resid_aa,
     +                    T1resid_bb,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
     
      Dimension T1resid_aa(Nvrt_a,Nocc_a)
      Dimension T1resid_bb(Nvrt_b,Nocc_b)

      Dimension W4_aa(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension W4_bb(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension W4_ab(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension W4_ba(Nvrt_a,Nvrt_b,Nocc_a,Nvrt_b)

      Dimension W5_aa(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
      Dimension W5_bb(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
      Dimension W5_ab(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
      Dimension W5_ba(Nocc_a,Nocc_b,Nvrt_a,Nocc_b)

      Integer M,N,A,F,E

#include "cse12_vars.com"

C AA block 
C Rai_a(A,I)<- +1/2T2(E,F,I,M)*W(E,F,A,M)
      
      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do e = 1, Nvrt_a
      Do f = 1, Nvrt_a
         C = (1.0D0-Ocn_va(A))
         T = T + 0.50D0*T2aa(E,F,I,M)*W4_aa(E,F,A,M)*C
      Enddo
      Enddo
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo 
      Enddo 

C Rai_a(A,I)<- +1/2T2(E,f,I,m)*W(E,f,A,m)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do e = 1, Nvrt_a
      Do f = 1, Nvrt_b
         C = (1.0D0-Ocn_va(A))
         T = T + T2ab(E,f,I,m)*W4_ab(E,f,A,m)*C
      Enddo
      Enddo
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo
      Enddo

C Rai_a(A,I)<- +1/2T2(A,E,M,N)*W(N,M,I,E)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_a
      Do e = 1, Nvrt_a 
         C = Ocn_oa(I)
         T = T + 0.50D0*T2aa(A,E,M,N)*W5_aa(N,M,I,E)*C
      Enddo
      Enddo
      Enddo 
         T1resid_aa(A,I) = T1resid_aa(A,I) + T
      Enddo 
      Enddo 

C Rai_a(A,I) <- -1/2T2(A,e,N,m)*W(N,m,I,e)

      Do i = 1, Nocc_a
      Do a = 1, Nvrt_a
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_a
      Do e = 1, Nvrt_b
         C = Ocn_oa(I)
         T = T + T2ab(A,e,N,m)*W5_ab(N,m,I,e)*C
      Enddo
      Enddo
      Enddo
         T1resid_aa(A,I) = T1resid_aa(A,I) - T
      Enddo
      Enddo

C BB block
C Rai_a(a,i) <- +1/2T2(e,f,i,m)*W(e,f,a,m)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_b
         C = (1.0D0-Ocn_vb(a))
         T = T + 0.50D0*T2bb(e,f,i,m)*W4_bb(e,f,a,m)*C
      Enddo
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

C Rai(a,i) <- +1/2T2(e,F,i,M)*W(e,F,a,M)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do e = 1, Nvrt_b
      Do f = 1, Nvrt_a
         C = (1.0D0-Ocn_vb(a))
         T = T + T2ab(F,e,M,i)*W4_ba(F,e,M,a)*C 
      Enddo
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(A,I) + T
      Enddo
      Enddo

C Rai_a(a,i) <- +1/2T2(a,e,m,n)*W(n,m,i,e)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_b
      Do n = 1, Nocc_b
      Do e = 1, Nvrt_b
         C = Ocn_ob(i)
         T = T + 0.50D0*T2bb(a,e,m,n)*W5_bb(n,m,i,e)*C
      Enddo
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) + T
      Enddo
      Enddo

C Rai_a(a,i) <- -1/2T2(E,a,M,n)*W(M,n,E,i)

      Do i = 1, Nocc_b
      Do a = 1, Nvrt_b
         T = 0.0D0
      Do m = 1, Nocc_a
      Do n = 1, Nocc_b
      Do e = 1, Nvrt_a
         C = Ocn_ob(i)
         T = T + T2ab(E,a,M,n)*W5_ba(M,n,E,i)*C
      Enddo
      Enddo
      Enddo
         T1resid_bb(a,i) = T1resid_bb(a,i) - T
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
      L_aa = Nocc_a*Nvrt_a
      L_bb = Nocc_b*Nvrt_b
      call checksum("T1resid_aa:",T1resid_aa,L_aa)
      call checksum("T1resid_bb:",T1resid_bb,L_bb)
#endif

      Return 
      End 
