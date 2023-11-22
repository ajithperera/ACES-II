      subroutine Cmp_intm_I(Nxyz, Nx, Ny, Nz, CX, Cy, Cz, Mu) 
C
C This compute the intermediate I 
C Note that Nxyz = Nx + NY + NZ
C
      Implicit Double Precision (A-H, O-Z)

      Dimension CSxyz(0:Maxang,3)
C
c Compute CP_i^(n-k) where i=x,y,z; n=0,1,2;k=0,...n
C
      Call Make_xyz(Nxyz, Cx, CSxyz(1,1)) 
      Call Make_xyz(Nxyz, Cy, CSxyz(1,2)) 
      Call Make_xyz(Nxyz, Cz, CSxyz(1,3)) 

C Compute the  fact(n_i,k_i)CP_i(n_i-k_i) where i=x,y,z;n_1=1,2..
C k_i=0,1...n_i
C
      Call Sum_xyz(Nx, CSxyz(1,1), BCxyz(1,1)) 
      Call Sum_xyz(Ny, CSxyz(1,2), BCxyz(1,2)) 
      Call Sum_xyz(Nz, CSxyz(1,3), BCxyz(1.3)) 
C
C Compute the Sum_mu{-lam,+lam} 2*dsqrt(pi)Z_lammu(unit{Ca})Omega
C Where omega is the angualr integral. 
C
     Do Ix = 0, Nx
        Do Jy = 0, Ny
           Do Kz = 0, Nz
 
      Call Ang_int(
   

