      Subroutine Tdee_lanczos_main(Vecs,Work,Memleft,Nsize,Irrepx,Iuhf,
     +                             Maxmem)
 
      Implicit Double Precision(A-H,O-Z)

#include "lanczos_vars.com"

C This  subroutine uses  the Lanczos  algorithm  with look-ahead to
C compute eigenvalue  estimates.  The algorithm was first described
C in  the  RIACS Technical  Report 90.45, `An Implementation of the
C Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices, Part I`,
C by R.W. Freund, M.H. Gutknecht and N.M. Nachtigal, November 1990.
C The routine will first run the look-ahead Lanczos algorithm for a
C number of steps.  It will then prompt the user  for the number of
C eigenvalues to compute.  In addition, the routine will  also find
C the common eigenvalues in two successive runs.

      Integer Revcom,Colx,Colb,Info(4)
      Dimension Vecs(Nsize,2*Maxvw+6)
      Dimension Work(Memleft)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "--- Entered Lanczos_main ---"
      Write(6,"(a)") " The input vars,memleft,Nsize,Irrepx,Iuhf"
      Write(6,"(1x,i8,1x,i8,1x,i2,1x,i2)")Memleft,Nsize,Irrepx,iuhf
#endif 

      If (Maxn .Gt. Maxn_int) Then
         Write(6,"(a,a)") " The number of requested Lanczos vectors",
     +                    " are more than the currently allowed"
         Write(6,"(a,I6)")" maximum of ", No_max_vecs
         Call Errex
      Endif

C To initilize read both left and right starting vectors 

      Mubar_s_pq_t0 = 390
      Mubar_d_pq_t0 = 313
      Mubar_s_pq_tn = 394
      Mubar_d_pq_tn = 333

C Load Mubar_R0

      Ioffr1 = 0
      Ioffr2 = 0
      Call Tdee_Load_vec(Irrepx,Vecs(1,1),Nsize,Mubar_s_pq_t0,Ioffr1,
     +                   Mubar_d_pq_t0,Ioffr2,Iuhf,.False.)

C Load Mubar_L0 

      Ioffr1 = 2
      Ioffr2 = 10
      Call Tdee_Load_vec(Irrepx,Vecs(1,2),Nsize,Mubar_s_pq_t0,Ioffr1,
     +                   Mubar_d_pq_t0,Ioffr2,Iuhf,.False.)


C Du-lal stands for Double Precision Lanczos, LAL refer to three
C term look-ahead Lanczos. 
#ifdef _DULAL
      Info(1) = 111314
#else
      Info(1) = 1
#endif 

      Info(2) = 0
      Ncount  = 0
      Icount  = 0
      Jcount  = 0

 10   Continue 

      Ncount = Ncount + 1

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "Input for Dulal or Lanczos_light"
      If (ncount .eq.1) then
         call checksum("Lanczos_main:",Vecs(1,1),nsize,s)
         call checksum("Lanczos_main:",Vecs(1,2),nsize,s)
      else 
         if (revcom.eq.1) call checksum("Lanczos_main:",
     +                     Vecs(1,colb),nsize,s)
         if (revcom.eq.2) call checksum("Lanczos_main:",
     +                     Vecs(1,colb),nsize,s)
      endif 
#endif 
#ifdef _DULAL
C Dulal is a QMRPACK routine (see comments within the file Dulal.F).
C This implement a "look-ahead" algorithm. Implemented in 1990s 
C and I have not been able to get this to working correctly.
      Call Dulal(Nsize,Nsize,Nlim,Maxn,Maxvw,M,Norm,Dwk,Idx,Iwk,Hwk,
     +           Vecs,Info)
#endif 
C Always do one more iteration than requested. This is useful to asses 
C the convegence of Lanczos chains toward convergence. 

      If (Ncount .gt. 2*Maxn+4) Go to  99
      Call Lanczos_light_v1(Nsize,Vecs,Info)

      Ierror = Info(1)
#ifdef _DEBUG_LVL0
      Write(6,"(a,i2)") "Monitoring variable ", Ierror 
      If (Ierror .EQ.0) Go to 99
#endif 
      
      Revcom = Info(2)
      Colx   = Info(3)
      Colb   = Info(4)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The Revcom,Colx,Colb vars returned"
      Write(6,"(5(1x,i2))") Info(2),Info(3),Info(4),Colx,Colb
      if (revcom.eq.1) call checksum("Lanczos_main:",
     +                                Vecs(1,colx),nsize,s)
      if (revcom.eq.2) call checksum("Lanczos_main:",
     +                                Vecs(1,colx),nsize,s)
#endif 
      If (Revcom .EQ. 1) Then

         If (Convrg .EQ. "Strict") Then 
            Icount = Icount + 1
#ifdef _DEBUG_LVL0
            Write(6,*) "Icount & Ncount", Icount, Ncount 
#endif 
C Do the external biorthogonalization of Lanczos vectors. This should
C prevent the deterioration of biorthogonality of Lanczos vectors. 
C First right vector with all the left vectors.

            If (Ncount .Gt. 1) Call Tdee_gs_ortho_rl(Work,Memleft,
     +                              Vecs(1,Colx),Nsize,Icount,
     +                              Irrepx)

            If (Ncount .eq.1) Call Updmoi(maxn+4,nsize,Irrepx,497,0,0)
            Call Putlst(Vecs(1,Colx),Icount,1,0,Irrepx,497)
  
         Else if (Convrg .EQ. "Light") Then
            Icount = Icount + 1
#ifdef _DEBUG_LVL0
            Write(6,*) "Icount & Ncount", Icount, Ncount 
#endif 
            If (Icount .EQ. Maxn+2) Then
               Call Updmoi(Icount-1,nsize,Irrepx,497,0,0)
               Call Putlst(Vecs(1,Colx),Icount-1,1,0,Irrepx,497)
            Endif 
 
         Endif 

         Ioffr1 = 0
         Ioffr2 = 0
         Ioffsp = 0

C Do the right hand multiplication and continue

         Call Tdee_dump_vec(Irrepx,Vecs(1,Colx),Nsize,Mubar_s_pq_t0,
     +                      Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                      Iuhf,.False.)

         Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,1) 

C This Returns Mubar_R(n+1), iside=1

         Call Tdee_Load_vec(Irrepx,Vecs(1,Colb),Nsize,Mubar_s_pq_tn,
     +                     Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)
         Go to 10 

      Elseif (Revcom .EQ. 2) Then

         If (Convrg .EQ. "Strict") Then 
            Jcount = Jcount + 1

C Do the external biorthogonalization of Lanczos vectors. This should
C prevent the deterioration of biorthogonality of Lanczos vectors.
C Now left vector with all the right vectors.

            If (Ncount .Gt. 1) Call Tdee_gs_ortho_lr(Work,Memleft,
     +                              Vecs(1,Colx),Vecs(1,3),Nsize,Jcount,
     +                              Irrepx)
            If (Ncount .eq.2) Call updmoi(maxn+4,nsize,Irrepx,498,0,0)
#ifdef _DEBUG_LVL0
         Write(6,*) "Jcount & Ncount",jcount, Ncount 
#endif 
           Call Putlst(Vecs(1,Colx),Jcount,1,0,Irrepx,498)

         Else if (Convrg .EQ. "Light") Then
            Jcount = Jcount + 1
#ifdef _DEBUG_LVL0
            Write(6,*) "Jcount & Ncount", Jcount, Ncount 
#endif 
            If (Jcount .EQ. Maxn+2) Then
               Call Updmoi(Jcount-1,nsize,Irrepx,498,0,0)
               Call putlst(Vecs(1,Colx),Jcount-1,1,0,Irrepx,498)
            Endif 
  
         Endif 

         Ioffr1 = 2
         Ioffr2 = 10
         Ioffsp = 0

         Call Tdee_dump_vec(Irrepx,Vecs(1,Colx),Nsize,Mubar_s_pq_t0,
     +                    Ioffr1,Ioffsp,Mubar_d_pq_t0,Ioffr2,
     +                    Iuhf,.False.)

         Call Tdee_hbar_mult(Work,Memleft,Iuhf,Irrepx,2) 

C This Returns Mubar_L(n+1), iside=2

         Call Tdee_Load_vec(Irrepx,Vecs(1,Colb),Nsize,Mubar_s_pq_tn,
     +                      Ioffr1,Mubar_d_pq_tn,Ioffr2,Iuhf,.False.)
         Go to 10 

      Endif 

  99  Continue

      Write(6,"(a)") "Lanczos iterations are completed!" 

C#ifdef _DEBUG_LVL0
      If (Convrg .EQ. "Strict") Then
          Call check_lanczos(Work,Maxmem,Nsize,irrepx,Maxn)
      Endif 
C#endif 
    
C At this point non-symmetric tridiagonal matirx is complete. Digonalize
C (using dgeev) and obtain left and right eigenvectors and eigenvalues.

      Call Tdee_diag_trd(Work,Maxmem,Ncount,Nsize,Iuhf,Irrepx)

      Return
      End
