      Subroutine Lanczos_diag_trd(Work,Memleft,Nrow,Nsize,Iuhf,Irrepx)

      Implicit Double Precision (A-H,O-Z)

#include "lanczos_vars.com"

      Dimension Work(Memleft)

C Nrow is the actual dimension of the tridiagonal matrix and lda is
C the lead dimension. The lead dimension is set to 10,000 (ie. That
C many lanczos vectors are alloeewed).

      Lda  = Maxn_int
      Ncol = (Nrow - 4)/2 
      Nrow = Ncol
      
      If (Ncol .NE. Maxn) Then
          Write(6,*)
          Write(6,"(a,a)") " The dimension of the tridigonal matrix",
     +                     " is different from the what is requested."
          Call Errex
      Endif

#ifdef _DEBUG_LVL0
      Write(6,*) 
      Write(6,"(a)") " The tridigonal Lanczos matrix"
      Write(6,"(a,a,3(1x,I5))") " The variabels governing tridiagonal",
     +                        " matrix:Nrow,ncol,lda:", Nrow,Ncol,Lda
      Call output(Trd,1,Nrow+1,1,Ncol+1,Lda,Lda,1) 
#endif 
      I000 = 1
      I010 = I000 + Nrow 
      I020 = I010 + Nrow 
      I030 = I020 + Nrow * Ncol
      I040 = I030 + Nrow * Ncol
      I050 = I040 + 4*Lda
      Iend = I050

      If (Iend .GE. Memleft) Call Insmem("@-Tdee_diag_trd",Iend,
     +                                      Memleft)

C Save the following element before diagonalization. Note that this
C element of the Trd matrix is not subjected to the diagonalization,
C but it is better to save it from here. 

      Beta_k = Trd(Nrow,Ncol+1) 

      Call Dgeev("V","V",Nrow,Trd,Lda,Work(I000),Work(I010),
     +            Work(I020),Nrow,Work(I030),Nrow,Work(I040),
     +            4*Lda,Ierror)

      If (Ierror .Ne. 0) then
         Write(6,*) 
         Write(6,"(a,a)") " The diagonalization of tridiagonal",
     +                    " subspace matirx failed."
         Call Errex
      Endif 

C The real and imaginary eigenvalues are stored at I000 and I010 locations.
C
C#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " The real eigenvalues of the Lanczos matrix"
      Write(6,"(5(1X,F15.7))") (Work(I000-1+i),i=1,Nrow) 
      Write(6,*)
      Write(6,"(a)") " The imaginary eigenvalues of Lanczos matrix"
      Write(6,"(5(1X,F15.7))") (Work(I010-1+i),i=1,Nrow) 
C#endif 

C The left and right eigenvectors are stored at I020 and I030.
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " The right eigenvectors of the Lanczos matrix"
      Call output(Work(I030),1,Nrow,1,Ncol,Nrow,Ncol,1)
      Write(6,*)
      Write(6,"(a)") " The left eigenvectors of the Lanczos matrix"
      Call output(Work(I020),1,Nrow,1,Ncol,Nrow,Ncol,1)
CSSS      call dgemm("T","N",Ncol,Ncol,Nrow,1.0D0,Work(I020),Nrow,
CSSS     +             Work(I030),Nrow,0.0D0,Work(I040),Nrow)
CSSS      Call output(Work(I040),1,Nrow,1,Ncol,Nrow,Ncol,1)
#endif 

      If (Convrg .EQ. "Strict") Then
         I050 = I040 + Nrow
         Iend = I050

         If (Iend .GE. Memleft) Call Insmem("@-Tdee_diag_trd",Iend,
     +                                        Memleft)
         Memleft = Memleft - Iend

C Save the righr eigenvectors to I050 before residual is built

         Call Lanczos_screen_eigs_strict(Work(I000),Work(I020),
     +                                   Work(I030),Work(I040),
     +                                   Work(Iend),Nrow,Ncol,
     +                                   Irrepx,Nsize,Iuhf,Memleft)

      Elseif (Convrg .EQ. "Light") Then 

C Access the element of the tridiagonal matrix that corresond to k+1
C where k is the Lanczos chain length. Lanczos loops are set up such
C that we always do a one more iteration than the requested.

         I050 = I040 + Nrow
         I060 = I050 + Nrow
         Iend = I060
         If (Iend .GE. Memleft) Call Insmem("@-Tdee_diag_trd",Iend,
     +                                        Memleft)
         Memleft = Memleft - Iend
         Call Lanczos_screen_eigs_light(Work(I000),Work(I020),
     +                                  Work(I030),Work(I040),
     +                                  Work(I050),
     +                                  Work(Iend),Beta_k,Nrow,Ncol,
     +                                  Irrepx,Nsize,Iuhf,Memleft)
       Endif
 
C Having established the convergence, built the dipole responce 
C functions.

      Call Lanczos_built_dp_respns(Work(I000),Work(I020),Work(I030),
     +                             Work(I040),Work(I050),Work(Iend),
     +                             Nrow,Ncol,Irrepx,Iuhf,Memleft)

      Return
      End


