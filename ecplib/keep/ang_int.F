      Subroutine ang_int(Kx, Ky, Kz, Lamu, xhat, yhat, zhat, 
     &                   Sum1_xyz)
C
C Compute the angular integral Z_lam,muI_{lam,mu,l,m}^{nx,ny,nz} =
C 4pi(i-1)!!(j-1)!!(k-1)!!/(i+j+k+1)!! if all i,j,and k are even
C otherwise zero. 

      Implicit Double Precision(A-H, O-Z)

#include "ecp.par"
#include "ecp.com"

      Dimension Sum1_xyz(Maxang+1,3)

      Logical kxeven, kyeven, kzeven

      Call Dzero(Ang_sum,9*Lamu)

      Do N = 1, Kx

         If (Sum1_xyz(N, 1) .EQ. 0.0D0) Go To 10

         Do L = 1, Ky

            If (Sum1_xyz(N, 1) .EQ. 0.0D0) Go To 20
          
            Do M = 1, Kz

              If (Sum1_xyz(N, 1) .EQ. 0.0D0) Go To 30

                 NLM = N + L + M - 2
                 Lam_Lo = Mod(NLM - 1, 2) + 1
                 Lam_hi = Min(NLM, Lamu)
                 
                 Do Lam = Lam_Lo, Lam_hi, 2
 
                    L2 = 2*Lam - 1

                    Do Mu
      If (Sum1_xyz(
      kxeven = mod(kx,2) .EQ. 0 
      kyeven = mod(ky,2) .EQ. 0
      kzeven = mod(kz,2) .EQ. 0
C
      kxyz = kx + ky + kz

      If (kxeven .AND. kyeven .AND. kzeven) the

      angint = (dfact(kx-1)*dfact(ky-1)*dfact(kz-1))/dfact(kxyz+1)
      angint = 4*pi*angint

      Else
          angint = 0.0D0

      Endif

      Return
      End
