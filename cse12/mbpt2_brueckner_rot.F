#include "flags.h"
      Subroutine Mbpt2_brueckner_rot(T1aa,T1bb,Work,Nocc_a,Nocc_b,
     +                               Nvrt_a,Nvrt_b,Maxcor)
 
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Character*1 SP(2)
      Character*8 String
      Double Precision One,Mone 

#include "cse12_vars.com"
#include "info.com"
#include "machsp.com"
#include "jodaflags.com"

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Data Ione,Izero,One,Mone,Zilch,Izz /1,0,1.0D0,-1.0D0,0.0D0,100/
      Data SP /"A","B"/

      Nbasis_a   = Nocco(1) + Nvrto(1)
      Nbasis_b   = Nocco(2) + Nvrto(2)
      Nbasis     = Max(Nbasis_a,Nbasis_b)
      
      Laa = Nvrt_a*Nocc_a 
      Lbb = Nvrt_b*Nocc_b

      X_aa = Fndlrgab(T1aa,Laa)
      X_bb = Fndlrgab(T1bb,Lbb)
      X    = Max(X_aa,X_bb)

      If (X .Lt. Brueck_tol) Then
          Call Putrec(20,"JOBARC","BRUKTEST",Ione,Ione)
          Write(6,"(2a)") " Frac. occ. MBPT(2) Brueckner iterations",
     +                    " have converged."
          Write(6,"(2a,F15.10,1x,E8.2E2)") " The largest t1 amplitude",
     +                                     " and the tolerance:", X,
     +                                       Brueck_tol
      Else
          Write(6,"(2a)") " Current T1 amplitudes are used to define",
     +                    " new set of orbitals." 

          Write(6,"(2a,F15.10,1x,E8.2E2)") " The largest t1 amplitude",
     +                                     " and the tolerance:", X,
     +                                       Brueck_tol
          Call Putrec(20,"JOBARC","BRUKTEST",Ione,Izero)
  
          Do Ispin = 1, 2 

             If (Ispin .Eq. 1) Then
                 Nvrt = Nvrt_a
                 Nocc = Nocc_a
             Elseif (Ispin .Eq. 2) Then
                 Nvrt = Nvrt_b
                 Nocc = Nocc_b
             Endif
             If (Nbasis .EQ. Nocc) Then
                Ioff = 0
             Else
                Ioff  = Nbasis*Nocc  
             Endif 

             Icmo  = 1
             Itcmo = Icmo  + Nbasis*Nbasis 
             Iovlp = Itcmo + Nbasis*Nbasis
             Iovlc = Iovlp + Nbasis*Nbasis 
             Itmp1 = Iovlc + Nbasis*Nbasis
             Itmp2 = Itmp1 + Nbasis*Nbasis
             Iend  = Itmp2 + Nbasis*Nbasis 

             String = "SCFEVC"//SP(Ispin)//"0"
             Call Getrec(20,"JOBARC",String,Nbasis*Nbasis*Iintfp,
     +                   Work(Icmo)) 

             Call Dcopy(Nbasis*Nbasis,Work(Icmo),1,Work(Itcmo),1)

             If (Ispin .EQ. 1) Then
C Generate rotated virtuals 

             Call Xgemm("N","T",Nbasis,Nvrt,Nocc,Mone,Work(Icmo),
     +                   Nbasis,T1aa,Nvrt,One,Work(Itcmo+Ioff),Nbasis)

C Generate rotated occupied 

             Call Xgemm("N","N",Nbasis,Nocc,Nvrt,One,Work(Icmo+Ioff),
     +                   Nbasis,T1aa,Nvrt,One,Work(Itcmo),Nbasis)
             Else 
C Generate rotated virtuals

             Call Xgemm("N","T",Nbasis,Nvrt,Nocc,MOne,Work(Icmo),
     +                   Nbasis,T1bb,Nvrt,One,Work(Itcmo+Ioff),Nbasis)

C Generate rotated occupied

             Call Xgemm("N","N",Nbasis,Nocc,Nvrt,One,Work(Icmo+Ioff),
     +                   Nbasis,T1bb,Nvrt,One,Work(Itcmo),Nbasis)
             Endif 

C Form S^-(1/2) in Mo basis 

             Call Getrec(20,"JOBARC","AOOVRLAP",Nbasis*Nbasis*Iintfp,
     +                   Work(Iovlp))

C Transform to MO basis using the rotated vectors. 

             Call Xgemm("N","N",Nbasis,Nbasis,Nbasis,One,Work(Iovlp),
     +                   Nbasis,Work(Itcmo),Nbasis,Zilch,Work(Itmp1),
     +                   Nbasis)

             Call Xgemm("T","N",Nbasis,Nbasis,Nbasis,One,Work(Itcmo),
     +                   Nbasis,Work(Itmp1),Nbasis,Zilch,Work(Itmp2),
     +                   Nbasis)

C Build the orthonormalization matrix.

             Call Eig(Work(Itmp2),Work(Iovlc),Nbasis,Nbasis,0)
             Call Invsqt(Work(Itmp2),Nbasis+1,Nbasis)
              
             Call Xgemm("N","N",Nbasis,Nbasis,Nbasis,One,Work(Iovlc),
     +                   Nbasis,Work(Itmp2),Nbasis,Zilch,Work(Itmp1),
     +                   Nbasis)
             Call Xgemm("N","T",Nbasis,Nbasis,Nbasis,One,Work(Itmp1),
     +                   Nbasis,Work(Iovlc),Nbasis,Zilch,Work(Iovlp),
     +                   Nbasis)

C Transform the rotated vectors to   

             Call Xgemm("N","N",Nbasis,Nbasis,Nbasis,One,Work(Itcmo),
     +                   Nbasis,Work(Iovlp),Nbasis,Zilch,Work(Icmo),
     +                   Nbasis)

#ifdef _DEBUG_LVLM
      Write(6,"(a)") " @-Mbpt2_brueckner_rot: The C_new(t)*S*C_new" 
      Call Getrec(20,"JOBARC","AOOVRLAP",Nbasis*Nbasis*Iintfp,
     +            Work(Iovlp))
      Call Xgemm("T","N",Nbasis,Nbasis,Nbasis,One,Work(Icmo),
     +            Nbasis,Work(Iovlp),Nbasis,Zilch,Work(Itmp1),
     +            Nbasis)
      Call Xgemm("N","N",Nbasis,Nbasis,Nbasis,One,Work(Itmp1),
     +            Nbasis,Work(Icmo),Nbasis,Zilch,Work(Itmp2),
     +            Nbasis)
      Call output(Work(Itmp2),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

             Call Putrec(20,"JOBARC",String,Nbasis*Nbasis*Iintfp,
     +                   Work(Icmo)) 
          Enddo 
      Endif 

C Since we have to build a fock matrix, it is better to transfer
C the control to standard SCF code. We need to set the flags 
C scf_maxcyc=0,non_hf=on. 
      
      Iflags(h_IFLAGS_scf_maxcyc) = 0
      Iflags(h_IFLAGS_nonhf)      = 1

      Call Putrec(20,"JOBARC","IFLAGS  ",Izz,Iflags) 

      Return
      End 
