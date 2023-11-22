      Subroutine Form_t2ijkl(T2aa,T2bb,T2ab,W_aa,W_bb,W_ab,
     +                       T2ijkl_1,T2ijkl_2,T2ijkl_3,Nocc_a,
     +                       Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension W_aa(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension W_bb(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
      Dimension W_ab(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Dimension T2ijkl_1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
      Dimension T2ijkl_2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
      Dimension T2ijkl_3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)

      Integer M,N,I,J,E,F
      Data Dzero /0.0D0/

#include "ccsdlight_vars.com" 

C T2ijkl_1(MN,IJ) = ) + 1/4T2(IJ,EF)*W(MN,EF) 
 
      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do N = 1, Nocc_a
      Do M = 1, Nocc_a
      T2ijkl_1(M,N,I,J) = Dzero 
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_a
      T2ijkl_1(M,N,I,J) = T2ijkl_1(M,N,I,J) + 
     +                    T2aa(E,F,I,J)*W_aa(M,N,E,F)*0.250D0
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2ijkl_2(mn,ij) =  1/4T2(ij,ef)*W2(mn,ef) 

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do N = 1, Nocc_b
      Do M = 1, Nocc_b
      T2ijkl_2(m,n,i,j) = Dzero 
      Do E = 1, Nvrt_b
      Do F = 1, Nvrt_b
      T2ijkl_2(m,n,i,j) = T2ijkl_2(m,n,i,j) + 
     +                    T2bb(e,f,i,j)*W_bb(m,n,e,f)*0.250D0
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

C T2ijkl_3(Mn,Ij) =  T2(Ij,Ef)*W2(Mn,Ef)

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do N = 1, Nocc_b
      Do M = 1, Nocc_a
      T2ijkl_3(M,n,I,j) = Dzero
      Do E = 1, Nvrt_a
      Do F = 1, Nvrt_b
      T2ijkl_3(M,n,I,j) = T2ijkl_3(M,n,I,j) + 
     +                    T2ab(E,f,I,j)*W_ab(M,n,E,f)*0.50D0
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Return
      End
