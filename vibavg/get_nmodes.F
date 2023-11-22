      Subroutine Get_nmodes(Hess,Nrmlmodes,Vfrqs,Work,Maxcor,Atmmass,
     +                      Nreals)

      Implicit Double Precision(A-H,O-Z)

      Double Precision Nrmlmodes 

#include "const.h"
#include "machsp.com"

      Dimension Hess(3*Nreals,3*Nreals)
      Dimension Nrmlmodes(3*Nreals,3*Nreals)
      Dimension Vfrqs(3*Nreals)
      Dimension Atmmass(Nreals)
      Dimension Work(Maxcor)

C See vib1.F in jode for comments about the origin of the 
C Fact. This convert the sqrt of the eigenvalues of the 
C Hessian to cm ^(-1).

      Data Dnull,Done,Thres,Au2cm1/0.0D0,1.0D0,1.0D-6,5.1404871D03/

      Ndim = 3*Nreals 

      Do I = 1, Ndim
         Do J = 1, Ndim
            Fact = Dsqrt(Atmmass(1+(I-1)/3)*Atmmass(1+(J-1)/3))
            If (Fact .Lt. Thres) Then
               Hess(J,I) = Dnull
            Else 
               Hess(J,I) = Hess(J,I)/Fact
            Endif 
         Enddo
      Enddo

#ifdef _DEBUG_LVL0
      write(6,"(a)") " The mass-weighted Cartesian Hessian"
      call output(Hess,1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 

C Note that the digonalization leads to all the motions.
     
      Call Eig(Hess,Nrmlmodes,Ndim,Ndim,0)

#ifdef _DEBUG_LVL0
      Print*,Au2cm1
      write(6,"(a)") "The eigenvalus of the Hessain (in cm^-1)"
      Write(6,"(5(1x,F15.10))") (Dsqrt(Dabs(Hess(i,i)))*Au2cm1,
     +                          i=1,Ndim)
      write(6,"(a)") "The eigenvalus of the Hessain (in a.u.)"
      Write(6,"(5(1x,F15.10))") (Dsqrt(Dabs(Hess(i,i))),i=1,Ndim) 
#endif 

C The reference point Cartesian-Normal transformation from the
C above diagonalization works for the cases where the full point
C group is Abelian. For non-Abelian groups, we must resolve the
C degenerate modes to pure irreducible representations of the Abelian
C subgroup. This is a complicated procedure and hase already been
C done in joda (vib1.F). Instead of trying to repeat it, let's 
C read the modified transformation matrices. 

      Call Getrec(20,"JOBARC","SYM_NMDS",Ndim*Ndim*Iintfp,Nrmlmodes)

C The eigenvalues in (ML^2T^-2) L^-2 * M^-1 and the vectors are in 
C are dimensionless.

      Do I = 1, Ndim
         Vfrqs(I) = Dsqrt(Dabs(Hess(I,I)))*Au2cm1 
      Enddo
C 
#ifdef _DEBUG_LVL0
C Check the orthonormality of the eigenvectors.
      write(6,"(a)") " The normal mode transformation matrix"
      call output(Nrmlmodes,1,Ndim,1,Ndim,Ndim,Ndim,1)
      I000 = 1
      Iend = I000 + Ndim*Ndim 
      If (Iend .Ge. Maxcor) Call Insmem("get_nmodes",Iend,Maxcor)
      call Dgemm("T","N",Ndim,Ndim,Ndim,1.0D0,Nrmlmodes,Ndim,
     +            Nrmlmodes,Ndim,0.0D0,Work,Ndim)
      Write(6,*) 
      Write(6,"(a)") " Test for the orthonormality of the Normal modes"
      call output(Work(I000),1,Ndim,1,Ndim,Ndim,Ndim,1)
#endif 

#ifdef _NOSKIP
C This is something that I thought that we needed (from Stephan 
C Sauer's student's work, but it seems that it is wrong

C Undo the mass weighing in Normal modes and normalize them.

      Do I = 1, Ndim
         Do J = 1, Ndim
            Fact = Dsqrt(Atmmass(1+(J-1)/3))
            If (Fact .Gt. Thres) Then
               Nrmlmodes(J,I) = Nrmlmodes(J,I)/Fact
            ELse 
               Write(6,"(2a)") " A fault occured. The square root of",
     +                         " the atomic masses can not be zero."
               Call Errex
            Endif 
         Enddo
      Enddo

      Do J = 1, Ndim
         Dnorm = Dsqrt(Ddot(Ndim,Nrmlmodes(1,J),1,Nrmlmodes(1,J),1))
         If (Dnorm .Gt. Thres) Then
             Dnorm = Done/Dnorm
             Call Dscal(Ndim,Dnorm,Nrmlmodes(1,J),1)
         Else
             Write(6,"(a)") "@-get_nmodes, division by zero",
     +                      " error occurs." 
             Call Errex
         Endif
      Enddo 
#endif

      Return
      End
      


