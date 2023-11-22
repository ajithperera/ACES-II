      Subroutine Lanczos_screen_eigs_light(Eig_vals,Leig_vecs,Reig_vecs,
     +                                     Residual_Norm,Rki_agg,Scr,
     +                                     Beta_k,Nrow,Ncol,Irrepx,
     +                                     Nsize,Iuhf,Memleft)

      Implicit Double Precision (A-H,O-Z)
      Double Precision Leig_vecs
      Double Precision Null

      Dimension Eig_vals(Nrow),Leig_vecs(Nrow,Ncol+1),Scr(Memleft),
     +          Reig_vecs(Nrow,Ncol+1),Residual_Norm(Nrow),
     +          Rki_agg(Nrow)

      I000 = 1
      I010 = I000 + Nsize
      Iend = I010

      If (Iend .GE. Memleft) Call Insmem("@-Tdee_screen_eigs",
     +                                   Iend,Memleft)

C Retrive the k+1 (k is the chain length) Right Lanczos vector
 
      Kp1 = Ncol+1
      Call Getlst(Scr(I000),Kp1,1,0,Irrepx,497)

      Q_kp1  = Ddot(Nsize,Scr(I000),1,Scr(I000),1)

#ifdef _DEBUG_LVL0
      Write(6,"(1x,a,F10.5)") "T(k,k+1)=B_k : ", Beta_k
      Write(6,"(1x,a,F10.5)") "Q_kp1        : ", Q_kp1
#endif
      Do I = 1, Ncol
         R_ki = Dabs(Reig_vecs(Nrow,I))
         Residual_norm(I) = Beta_k * R_ki * Q_kp1
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "The residual norms of Lanczos vectors,light"
      Write(6,*)
      Write(6,"(5(1x,F15.7))") (Residual_norm(I),I=1,Ncol)
#endif

      Return
      End 
