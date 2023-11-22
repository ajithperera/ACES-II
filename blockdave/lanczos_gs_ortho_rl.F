      Subroutine Lanczos_gs_ortho_rl(Work,Memleft,Vecs,Nsize,Ncount,
     +                               Irrepx)
      
      Implicit Double Precision(A-H,O-Z)
       
      Dimension Work(Memleft), Vecs(Nsize)

C At this point Vecs corrspond to k^th Lanczos right-hand vector
C and there are (k-1) Left Lanczos vectors are on the disk. We
C need to maintain the Biorthogonality of the current vector to
C all the previous vectors on the disk. Here Ncount is the variable
C k. 
#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "----Entered Tdee_gs_ortho_rl----"
#endif 
CSS      Return
      I000 = 1
      Iend = I000 + Nsize
      call checksum("checksum of V(j)",Vecs,Nsize,s)
      If (Iend .GE. Memleft) Call Insmem("@-Tdee_gs_ortho_rl",
     +                                    Iend,Memleft)
CSS      RNorm=Dsqrt(Ddot(Nsize,Vecs,1,Vecs,1))
CSS      Sqrt_Norm = 1.0D0/Dsqrt(RNorm)
CSS      Call Dscal(Nsize,Sqrt_norm,Vecs,1) 

      Do K = 1, Ncount-1
         Call Getlst(Work,K,1,0,Irrepx,498)
         Rl_norm=-Ddot(Nsize,Vecs,1,Work,1)
         Call Daxpy(Nsize,Rl_norm,Work,1,Vecs,1)
      Enddo

      RNorm=Dabs(Ddot(Nsize,Vecs,1,Vecs,1))
      Write(6,*) "Rnorm-----", Rnorm 

      Sqrt_Norm = 1.0D0/Dsqrt(RNorm)
      Call Dscal(Nsize,Sqrt_norm,Vecs,1) 

      Return
      End 
      
