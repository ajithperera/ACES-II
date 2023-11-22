      subroutine Htmbej_inrabij(T2aa_old,T2bb_old,T2ab_old,
     +                          T2resid_aa,T2resid_bb,T2resid_ab,
     +                          Htmbej_1,Htmbej_2,Htmbej_3,Htmbej_4,
     +                          Htmbej_5,Htmbej_6,Nocc_a,Nocc_b,
     +                          Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa_old(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb_old(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab_old(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)
 
      Dimension T2resid_aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2resid_bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2resid_ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension Htmbej_1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Htmbej_2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Htmbej_3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Htmbej_4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
      Dimension Htmbej_5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
      Dimension Htmbej_6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

      Integer A,B,I,J,M,E

#include "cse12_vars.com"

C AAAA block 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a 
      Do E = 1, Nvrt_a 
         T = T + Htmbej_1(M,B,E,J)*T2aa_old(A,E,I,M)
      Enddo 
      Enddo 

      T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T
      T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) - T
      T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) - T
      T2resid_aa(B,A,J,I) = T2resid_aa(B,A,J,I) + T

      Enddo 
      Enddo 
      Enddo 
      Enddo 

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
         T = T + Htmbej_4(m,B,e,J)*T2ab_old(A,e,I,m)
      Enddo
      Enddo

      T2resid_aa(A,B,I,J) = T2resid_aa(A,B,I,J) + T 
      T2resid_aa(B,A,I,J) = T2resid_aa(B,A,I,J) - T 
      T2resid_aa(A,B,J,I) = T2resid_aa(A,B,J,I) - T 
      T2resid_aa(B,A,J,I) = T2resid_aa(B,A,J,I) + T 

      Enddo
      Enddo
      Enddo
      Enddo

C BBBB block 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
         T = T + Htmbej_2(m,b,e,j)*T2bb_old(a,e,i,m)
      Enddo
      Enddo

      T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T
      T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) - T
      T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) - T
      T2resid_bb(b,a,j,i) = T2resid_bb(b,a,j,i) + T

      Enddo
      Enddo
      Enddo
      Enddo 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         T  = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         T = T + Htmbej_3(M,b,E,j)*T2ab_old(E,a,M,i)
      Enddo
      Enddo

      T2resid_bb(a,b,i,j) = T2resid_bb(a,b,i,j) + T 
      T2resid_bb(b,a,i,j) = T2resid_bb(b,a,i,j) - T 
      T2resid_bb(a,b,j,i) = T2resid_bb(a,b,j,i) - T 
      T2resid_bb(b,a,j,i) = T2resid_bb(b,a,j,i) + T 

      Enddo
      Enddo
      Enddo
      Enddo


C ABAB block 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         T = T +  Htmbej_1(M,A,E,I)*T2ab_old(E,b,M,j)
      Enddo
      Enddo
      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 

      Enddo
      Enddo
      Enddo
      Enddo 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a 
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
        T = T +  Htmbej_2(m,b,e,j)*T2ab_old(A,e,I,m)
      Enddo
      Enddo

      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 

      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_a
         T = T +  Htmbej_3(M,b,E,j)*T2aa_old(A,E,I,M)
      Enddo
      Enddo

      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 

      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_b
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_b
        T = T + Htmbej_4(m,B,e,J)*T2bb_old(a,e,i,m)
      Enddo
      Enddo

      T2resid_ab(B,a,J,i) = T2resid_ab(B,a,J,i) + T 

      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_a
      Do E = 1, Nvrt_b
         T = T -  Htmbej_5(M,b,I,e)*T2ab_old(A,e,M,j)
      Enddo
      Enddo

      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T 

      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         T = 0.0D0
      Do M = 1, Nocc_b
      Do E = 1, Nvrt_a
         T = T - Htmbej_6(m,A,j,E)*T2ab_old(E,b,I,m)
      Enddo
      Enddo

      T2resid_ab(A,b,I,j) = T2resid_ab(A,b,I,j) + T

      Enddo
      Enddo
      Enddo
      Enddo

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Htmbej_inRabij"
      L_aaaa = Nocc_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb = Nocc_b*Nocc_b*Nvrt_b*Nvrt_b
      L_abab = Nocc_a*Nocc_b*Nvrt_a*Nvrt_b
      call checksum("T2resid_aa:",T2resid_aa,L_aaaa)
      call checksum("T2resid_bb:",T2resid_bb,L_bbbb)
      call checksum("T2resid_ab:",T2resid_ab,L_abab)
      Write(6,*) 
#endif

      Return
      End 
