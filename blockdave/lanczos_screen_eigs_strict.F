      Subroutine Lanczos_screen_eigs_strict(Eig_vals,Leig_vecs,
     +                                      Reig_vecs,Residual_Norm,
     +                                      Scr,Nrow,Ncol,Irrepx,Nsize,
     +                                      Iuhf,Memleft)

      Implicit Double Precision (A-H,O-Z)
      Double Precision Leig_vecs
      Double Precision Null

      Dimension Eig_vals(Nrow),Leig_vecs(Nrow,Ncol),Scr(Memleft),
     +          Reig_vecs(Nrow,Ncol),Residual_Norm(Nrow)

      Data Null,One /0.0D0,1.0D0/

      I000 = 1 
      I010 = I000 + Nrow * Nsize
      I020 = I010 + Nrow * Nsize
      I030 = I020 + Nsize
      I040 = I030 + Nsize
      Iend = I040 

      If (Iend .GE. Memleft) Call Insmem("@-Tdee_screen_eigs",
     +                                   Iend,Memleft)

C Obtain the memory left for multiplication (the following Hbar
C multiplication use the memory left from the pointer Iend.

      Memleft = Memleft - Iend

C Built the subspapce Lanczos matrices that represent the full
C Hbar matirces.

      Ioff = I000
      Joff = I010
      Do Icount = 1, Nrow

         Call Getlst(Scr(Ioff),Icount,1,0,Irrepx,497)
         Call Getlst(Scr(Joff),Icount,1,0,Irrepx,498)

         Ioff = Ioff + Nsize 
         Joff = Joff + Nsize 
      Enddo 

#ifdef _DEBUG_LVLM
      Write(6,*)
      I030 = I020 + Nrow * Ncol 
      Iend = I030
      If (Iend .GE. Memleft) Call Insmem("@-Tdee_screen_eigs",
     +                                   Iend,Memleft)
      Call Xgemm("T","N",Nrow,Ncol,Nsize,One,Scr(I010),Nsize,
     +            Scr(I000),Nsize,Null,Scr(I020),Nrow)
      Call checksum("Trace of T=P^tQ matrix",Scr(I020),  
     +               Nrow*Ncol,S)
#endif 

c Loop over each eigenvector and compute (Hbar R - E R) = Res. If
c the eigenvalues are proper then the residue must fall below
c accepted threshold. This must be done with the left vector too. If
c one of them fails then that should be good enough to reject that 
c particular eigenvalue. The remaining proper ones must be checked 
c with the left multiplication. At the moment we do do not do this 
c to save time, but easier to extend. 

      Mubar_s_pq_t0 = 490
      Mubar_d_pq_t0 = 443
      Mubar_s_pq_tn = 493
      Mubar_d_pq_tn = 447 

      Do N = 1, Ncol

        Ioffr1 = 0
        Ioffr2 = 0
        Ioffsp = 0
     
         Call Lanczos_normaliz_init_vec(Reig_Vecs(1,N),Nrow)

         Call Xgemm("N","N",Nsize,1,Nrow,One,Scr(I000),Nsize,
     +               Reig_vecs(1,N),Nrow,Null,Scr(I020),Nsize)

#ifdef _DEBUG_LVL0
      call checksum("Q            :",Scr(I000),Nsize*Nrow,s)
      call checksum("Q*R_n        :",Scr(I020),Nsize,s)
#endif 
         Call Lanczos_dump_vec(Irrepx,Scr(I020),Nsize,
     +                         Mubar_s_pq_t0,Ioffr1,Ioffsp,
     +                         Mubar_d_pq_t0,Ioffr2,Iuhf,.False.)

         Call Hbarxc(Scr(Iend),Memleft,Iuhf,1,Irrepx)

         Call Lanczos_Load_vec(Irrepx,Scr(I030),Nsize,
     +                      Mubar_s_pq_tn,Ioffr1,Mubar_d_pq_tn,
     +                      Ioffr2,Iuhf,.False.)
         Call Lanczos_normaliz_init_vec(Scr(I030),Nsize)

#ifdef _DEBUG_LVLM
         call checksum("Hbar*Q*R_n:",Scr(I030),Nsize,s)
         write(6,*) 
#endif 

C Construct Hbar*Q*R_n - E * R_n 

         Call Daxpy(Nsize,-Eig_vals(N),Scr(I020),1,Scr(I030),1)

         Residual_norm(N) = Ddot(Nsize,Scr(I030),1,Scr(I030),1)

      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") "The residual norms of Lanczos vectors,Strict"
      Write(6,*)
      Write(6,"(5(1x,F15.7))") (Residual_norm(I),I=1,Ncol)
#endif 

      Return
      End 
     
