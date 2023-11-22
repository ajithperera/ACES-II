#include "flags.h"
      Subroutine Pccd_rotgrd(Work,Maxcor,Nbas,Nocca,Noccb,Nvrta,Nvrtb,
     +                       Iuhf,Gradcnv,Tol,Y,X,Ncycle)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Tol,X,Y,Fndlrgab
      Double Precision Grad_stat
      Logical Symmetry
      Logical Diis 

      Dimension Grad_stata(6),Grad_statb(6)
      Dimension Work(Maxcor)
      Logical Gradcnv,B_like

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
      Common/Symm/Symmetry
      Common/Extrapo/Diis
      Data Ione, Izero /1,0/
#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,*) "-----------------pCCD_rotgrd-------------------"
#endif 

      Irrepx   = Ione
      Length   = Nbas*Nbas

      Ndim_ooa   = Nocca*Nocca
      Ndim_vva   = Nvrta*Nvrta
      Ndim_voa   = Nocca*Nvrta
      Ndim_oob   = Noccb*Noccb
      Ndim_vvb   = Nvrtb*Nvrtb
      Ndim_vob   = Noccb*Nvrtb

      Tol=10.0**-(Iflags((h_IFLAGS_bruck_conv)))

      Do Ispin = 1, Iuhf + 1

         If (Ispin .Eq. 1) Then 
            Nocc    = Nocca 
            Nvrt    = Nvrta
            Ndim_oo = Ndim_ooa
            Ndim_vv = Ndim_vva 
            Ndim_vo = Ndim_voa
         Elseif (Ispin .Eq. 2) Then
            Nocc    = Noccb
            Nvrt    = Nvrtb
            Ndim_oo = Ndim_oob
            Ndim_vv = Ndim_vvb
            Ndim_vo = Ndim_vob
         Endif 

         I000 = Ione
         I010 = I000 + Ndim_oo
         I020 = I010 + Ndim_vv
         I030 = I020 + Ndim_vo
         I040 = I030 + Ndim_vo
         Iend = I040 + Length
         Memleft = Maxcor - Iend 
         If (Iend .Ge. Maxcor) Call Insmem("pccd_rotg",Iend,Maxcor)

         If (Ispin .Eq. 1) Then
             Call Getrec(20,"JOBARC","OBRTGRDA",
     +                                Length*Iintfp,Work(I040))
             Call Pccd_vstat(Work(I040),Grad_stata,Nbas*Nbas)
             Call Putrec(20,"JOBARC","GRDSTATA",6*Iintfp,Grad_stata)
         Endif 
         If (Ispin .Eq. 2) Then 
            Call Getrec(20,"JOBARC","OBRTGRDB",
     +                               Length*Iintfp,Work(I040))
            Call Pccd_vstat(Work(I040),Grad_statb,Nbas*Nbas)
            Call Putrec(20,"JOBARC","GRDSTATB",6*Iintfp,Grad_statb)
         Endif 

#ifdef _DEBUG_LVL0
         write(6,*)
         If (Ispin .Eq. 1) write(6,"(2a)") "The alpha orbital rotation",
     +                                     " gradients"
         If (Ispin .Eq. 2) write(6,"(2a)") "The beta orbital rotation",
     +                                     " gradients"
         call output(Work(I040),1,Nbas,1,Nbas,Nbas,Nbas,1)
         Write(6,*)
         If (Ispin .Eq. 1) Write(6,"(2a)")  " The statistics of the ",
     +                     "alpha orbital rotation gradient matrix"
         If (Ispin .Eq. 2) Write(6,"(2a)")  " The statistics of the ",
     +                     "beta orbital rotation gradient matrix"
         Write(6,"(1x,2a)") "-----------------------------------------",
     +                      "----------"
         Write(6,"(5x,a,5xa,5xa)") "Minimum grad.", "Maximum grad.",
     +                              "RMS grad"
         If (Ispin .Eq. 1) Then 
             Write(6,*)
             Write(6,"(3(5x,E12.6))") Grad_stata(3), Grad_stata(4), 
     +                                Grad_stata(5)
         Elseif (Ispin .Eq. 2) Then
             Write(6,*)
             Write(6,"(3(5x,E12.6))") Grad_statb(3), Grad_statb(4), 
     +                                Grad_statb(5)
         Endif 
         Write(6,*)
         Write(6,"(1x,2a)") "-----------------------------------------",
     +                      "----------"
#endif
         Call Pccd_prep_grad(Work(I040),Work(I000),Work(I010),
     +                       Work(I020),Work(I030),Work(Iend),
     +                       Memleft,Ndim_oo,Ndim_vv,Ndim_vo,
     +                       Nocc,Nvrt,Nbas)

         Call Pccd_rotg(Work(I040),Work(I000),Work(I010),
     +                  Work(I020),Work(I030),Grad_stata,Grad_statb,
     +                  B_like,Work(Iend),Memleft,Lenoo,Lenvv,Lenvo,
     +                  Nocc,Nvrt,Nbas,Ispin,Iuhf,Ncycle,Tol)
      Enddo 

      If(Diis) Call Pccd_rotate_aux(Work(I000),Maxcor,Nbas,Iuhf,
     +                              Ncycle)

C B_like is true for Breuckner like optimiztion using orbital rotation
C gradients (i.e. fake).

      If (B_like) Return 

C Grad_stat(1) : Larget absolute 
C Grad_stat(2) : Smallest absolute 
C Grad_stat(3) : Largest value
C Grad_stat(4) : Smallest value 
C Grad_stat(5) : RMS gradient 
C Grad_stat(6) : Dynamic range (abs. min-abs. max)

#ifdef _DEBUG_LVLM
      write(6,"(a)") " The alpha gradient stats"
      Write(6,"(6(1x,F12.6))") (Grad_stata(i),i=1,6)
      write(6,"(a)") " The beta   gradient stats"
      Write(6,"(6(1x,F12.6))") (Grad_statb(i),i=1,6)
#endif 
      Write(6,*)
      Write(6,90)
90    FORMAT(T3,'Determination of pCCD/CCD/LCCD optimum orbitals.') 
      Write(6,300)
300   Format(T3,"Orbital rotation gradients are used to determine",
     +          " new orbitals.")
      Write(6,301)
301   Format(T3,"Newton-Raphson algorithm combined with DIIS is used",
     +          " to define the unitary rotation matrix.")
      Ya  = Grad_stata(1)
      Xa  = Grad_stata(5)
      Yb  = Ya
      Xb  = Xa
      If (Iuhf .Ne. 0) Then 
         Yb  = Grad_statb(1)
         Xb  = Grad_statb(5)
      Endif 
      Y = Max(Ya,Yb)
      X = Max(Xa,Xb)
      Call Getrec(20,"JOBARC","TOTENERG",Iintfp,E_new)
      If (Ncycle .Eq. 1 .or. Ncycle .Eq. 0) Then 
          Call Putrec(20,"JOBARC","SAVEENR ",Iintfp,E_new)
      Else 
          Call Getrec(20,"JOBARC","SAVEENR ",Iintfp,E_old)
      Endif 
      If (Ncycle .gt. 0) Ediff = Dabs(E_old-E_new)
      Write(6,"(a,F15.8,a,i2)")  "  |E(n-1)-E(n)| = ", Ediff, 
     +                           " where n = ", ncycle
      Write(6,100) Y,X,Tol
100   FORMAT(T3,"Largest and RMS gradients and tolerance are: ",
     +           E12.6,1x,E12.6,1x,E12.6,'.')
      Call Getrec(-20,'JOBARC','ORBOPITR',Ione,Ncycle)
      If (X .Lt. Tol) Then
         Call Putrec(20,'JOBARC','GRADTEST',Ione,Ione)
         Call Putrec("20","JOBARC","MICRO_CV",Iintfp,Ione)
         Write(6,200) Ncycle 
200      Format(T3,"Pccd/CCD/LCCD gradient rotation converged in ", i3, 
     +             " optimization cycles.")
         GradCNV = .TRUE.
CThis is a dumb idea. The NON-HF nature of the fock matrix has to be recognized
C during the orbital rotation gradient evaluation.
CSSSS         Call Pccd_t2tot1(Work(I040),Memleft,0)
      Else
         Call Putrec(20,'JOBARC',"GRADTEST",IONE,0)
         Gradcnv = .FALSE.
      Endif 

      Call Pccd_orbopt_summary(Iuhf,Ncycle)

      RETURN
      END
